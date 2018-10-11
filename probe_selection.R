#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(logging))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(Seurat))

# Set up logging
# Logging level choices
C_LEVEL_CHOICES <- names(loglevels)
#initialize to info setting.
logging::basicConfig(level = "INFO")


check_arguments <- function(arguments){
# Check arguments supplied by user meet requirements.
  # Require the name of an output pdf file
  if (( ! ("input_dir" %in% names(arguments) )) ||
       (arguments$input_dir == "") ||
       (is.na(arguments$input_dir))) {
       logging::logerror(paste(":: --input_dir: Please enter a file path to ",
                               "the input directory.", sep = ""))
        stop("error, no --input_dir")
  }
  # Make sure the input directory exists
  if ( ! file.exists(arguments$input_dir)){
    error_message <- paste0("Provided input directory, ",
    arguments$input_dir, ",  does not exist.",
    "Please check your input and try again")
    stop(error_message)
  }
  # Parse inputs where a range may be supplied by user
  for ( i in c("pcs_print", "jackstraw_pcs", "dims_use")) {
    if ( ( ! is.numeric(arguments[[i]])) ){
      arguments[[i]] <- eval(parse(text = arguments[[i]]))
    }
    if ( ( ! is.numeric(arguments[[i]])) || (arguments[[i]] == "") ||
         (is.na(arguments[[i]])) || (!( i %in% names(arguments))) ) {
         logging::logerror(paste(" --", i,
                                 ": Please provide a vector of integers",
                                 "for the dimensions to be used"))
          stop("error, ", i, " invalid")
    }
  }
  return(arguments)
}


classify_data <- function(seuratObj, seed){
# Generate random forest classifier
  classifier_input_x <- t(seuratObj@data)
  classifier_input_y <- as.factor(as.numeric(seuratObj@ident))
  set.seed(seed)
  trained_classifier <- randomForest(x = as.matrix(classifier_input_x),
                                    y = factor(classifier_input_y),
                                    importance = TRUE )
  print(trained_classifier)
  return(trained_classifier)
}


output_gene_names <- function(trained_classifier, nGenes, filename) {
# Save gene names for top 100 factors from the classification model
# in order of importance as tsv file
  ranked_list <- trained_classifier$importance[order(
    -trained_classifier$importance[, "MeanDecreaseAccuracy"]), ]

  dput(rownames(ranked_list)[1:nGenes])
  top_markers <- rownames(ranked_list)[1:nGenes]

  write.table(top_markers, file = filename, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = TRUE )
}

report_runtime_variables <- function(arguments){
  print(arguments)
}

# All functions below and their descriptive comments derive from
# Seurat tutorial https://satijalab.org/seurat/pbmc3k_tutorial.html


filter_cells <- function(seuratObj, filt_genes_low, filt_genes_high,
                         filt_mito_low, filt_mito_high){
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratObj@data),
                     value = TRUE)
  percent.mito <- Matrix::colSums(
    seuratObj@raw.data[mito.genes, ]) / Matrix::colSums(seuratObj@raw.data)
  seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito,
                           col.name = "percent.mito")
  vln <- VlnPlot(object = seuratObj, features.plot =
                   c("nGene", "nUMI", "percent.mito"), nCol = 3)
  print(vln)
  par(mfrow = c(1, 2))
  GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "nGene")

# Filter out cells that have unique gene counts over 2,500 or less than
# 200 (defaults) Note that low.thresholds and high.thresholds are used
# to define a 'gate'.  -Inf and Inf should be used if you don't want a
# lower or upper threshold.
  seuratObj <- FilterCells(object = seuratObj, subset.names =
                           c("nGene", "percent.mito"),
                           low.thresholds = c(filt_genes_low, filt_mito_low),
                           high.thresholds = c(filt_genes_high, filt_mito_high))
  seuratObj <- NormalizeData(object = seuratObj,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)
  return(seuratObj)
}

find_variable_genes <- function(seuratObj, x_low_cutoff,
                                x_high_cutoff, y_low_cutoff){

  seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean,
                                 dispersion.function = LogVMR,
                                 x.low.cutoff =  x_low_cutoff,
                                 x.high.cutoff = x_high_cutoff,
                                 y.cutoff = y_low_cutoff)
  seuratObj <- ScaleData(object = seuratObj, vars.to.regress =
                         c("nUMI", "percent.mito"))
  return(seuratObj)
}

run_PCA <- function(seuratObj, pcs_print, genes_print){
  seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes,
                      do.print = TRUE, pcs.print = pcs_print,
                      genes.print = genes_print)
  PrintPCA(object = seuratObj, pcs.print = pcs_print,
           genes.print = genes_print,
           use.full = FALSE)
  viz <- VizPCA(object = seuratObj, pcs.use = pcs_print )
  print(viz)
  for (i in pcs_print) {
    if ( i == 1 ) {
      next
    } else {
    PCAPlot(object = seuratObj, dim.1 = (i - 1), dim.2 = i )
    }
  }
# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
  seuratObj <- ProjectPCA(object = seuratObj, do.print = FALSE)
  return(seuratObj)
}

run_JackStraw <- function(seuratObj, num_replicate, pcs){
  seuratObj <- JackStraw(object = seuratObj, num.replicate = num_replicate,
                         display.progress = FALSE)
  jackstraw <- JackStrawPlot(object = seuratObj, PCs = pcs)
  print(jackstraw)
  return(seuratObj)
}

run_elbow <- function(seuratObj){
  elbow <- PCElbowPlot(object = seuratObj)
  print(elbow)
}

find_clusters <- function(seuratObj, dims_use, resolution){
  seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca",
                            dims.use = dims_use, resolution = resolution,
                            print.output = 0, save.SNN = TRUE)
  PrintFindClustersParams(object = seuratObj)
  seuratObj <- RunTSNE(object = seuratObj, dims.use = dims_use, do.fast = TRUE)
  TSNEPlot(object = seuratObj)
  return(seuratObj)
}

find_all_markers <- function(seuratObj, thresh_use){
# find markers for every cluster compared to all remaining cells, report
# only the positive ones
  seuratObj.markers <- FindAllMarkers(object = seuratObj, only.pos = TRUE,
                                      min.pct = 0.25,
                                      thresh.use = thresh_use)
  seuratObj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
  return(seuratObj)
}

# Command line arguments
pargs <- optparse::OptionParser(usage = paste("%prog [options]",
                                            "--input_dir directory"))


# Path to Cell Ranger output directory to be used as analysis input
pargs <- optparse::add_option(pargs, c("--input_dir"),
            type = "character",
            action = "store",
            dest = "input_dir",
            metavar = "directory",
            help = paste("Input data for analysis.",
                         "Expected format is 10x Genomics CellRanger pipeline",
                         "filtered gene-barcode matrices MEX output directory",
                         "[Default %default][REQUIRED]"))

pargs <- optparse::add_option(pargs, c("--plotfile"),
            type = "character",
            default = "ps_plots.pdf",
            action = "store",
            dest = "plotfile",
            help = paste("Filename for plotted output.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--stdoutfile"),
            type = "character",
            default = "ps_stdout.txt",
            action = "store",
            dest = "stdoutfile",
            help = paste("Filename for messages captured from stdout.",
                         "[Default %default]"))
                         
pargs <- optparse::add_option(pargs, c("--outputfile"),
            type = "character",
            default = "FINAL_LIST.tsv",
            action = "store",
            dest = "outputfile",
            help = paste("Filename for top markers output.",
                         "[Default %default]"))
                         
pargs <- optparse::add_option(pargs, c("--seurat_object"),
            type = "character",
            default = NULL,
            action = "store",
            dest = "save_obj",
            help = paste("Filename for resulting Seurat object",
                        "[Default Seurat object not saved]"))

pargs <- optparse::add_option(pargs, c("--run_jackstraw"),
            type = "logical",
            default = FALSE,
            action = "store",
            dest = "run_js",
            help = paste("Control jackstraw plotting.",
                         "Please note enabling this option increases runtime.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--run_elbow_plot"),
            type = "logical",
            default = TRUE,
            action = "store",
            dest = "run_elbow_plot",
            help = paste("This argument controls elbow plotting.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--num_genes"),
            type = "integer",
            default = 100,
            action = "store",
            dest = "num_genes",
            help = paste("Number of top factor genes to print from",
                         "classification model.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--seed"),
            type = "integer",
            default = 1,
            action = "store",
            dest = "seed",
            help = paste("Set seed for Random Forest Classifier generation",
                         "to enable reproducible behavior",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--min_cells"),
            type = "integer",
            default = 3,
            action = "store",
            dest = "min_cells",
            help = paste("Seurat.CreateSeuratObject, when filtering,",
                         "include genes with detected expression",
                         "in at least this many cells.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--min_genes"),
            type = "integer",
            default = 200,
            action = "store",
            dest = "min_genes",
            help = paste("Seurat.CreateSeuratObject: when filtering,",
                         "include genes where at least",
                         "this many genes are detected.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--filter_genes_low_threshold"),
            type = "double",
            default = 200,
            action = "store",
            dest = "filt_genes_low",
            help = paste("Seurat.FilterCells: low cutoff for number of genes.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--filter_genes_high_threshold"),
            type = "double",
            default = 2500,
            action = "store",
            dest = "filt_genes_high",
            help = paste("Seurat.FilterCells: high cutoff for number of genes.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--filter_mito_low_threshold"),
            type = "double",
            default = -Inf,
            action = "store",
            dest = "filt_mito_low",
            help = paste("Seurat.FilterCells: low cutoff for",
                         "percent mitochondrial genes.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--filter_mito_high_threshold"),
            type = "double",
            default = 0.05,
            action = "store",
            dest = "filt_mito_high",
            help = paste("Seurat.FilterCells: high cutoff for",
                         "percent mitochondrial genes.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--x_low_cutoff"),
            type = "double",
            default = 0.0125,
            action = "store",
            dest = "x_low_cutoff",
            help = paste("Seurat.FindVariableGenes: Bottom cutoff on x-axis",
                         "for identifying variable genes.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--x_high_cutoff"),
            type = "double",
            default = 3,
            action = "store",
            dest = "x_high_cutoff",
            help = paste("Seurat.FindVariableGenes: Top cutoff on x-axis",
                         "for identifying variable genes.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--y_low_cutoff"),
            type = "double",
            default = 0.5,
            action = "store",
            dest = "y_low_cutoff",
            help = paste("Seurat.FindVariableGenes: Bottom cutoff on y-axis",
                         "for identifying variable genes.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--pcs_print"),
            type = "character",
            default = "1:5",
            action = "store",
            dest = "pcs_print",
            help = paste("Seurat.PrintPCA: Principle components to",
                         "print genes for.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--genes_print"),
            type = "integer",
            default = 5,
            action = "store",
            dest = "genes_print",
            help = paste("Seurat.PrintPCA: Number of genes to print for",
                         "each principle component.",
                         "[Default %default]"))
                         
                         
pargs <- optparse::add_option(pargs, c("--jackstraw_replicates"),
            type = "integer",
            default = 100,
            action = "store",
            dest = "jackstraw_repl",
            help = paste("Seurat.jackstraw: Number of replicate samplings",
                         "to perform for JackStraw.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--jackstraw_pcs"),
            type = "character",
            default = "1:12",
            action = "store",
            dest = "jackstraw_pcs",
            help = paste("Seurat.jackstraw: Number of PCs to compute",
                         "significance for. [Default %default]"))

pargs <- optparse::add_option(pargs, c("--dims_use"),
            type = "character",
            default = "1:10",
            action = "store",
            dest = "dims_use",
            help = paste("Seurat.FindClusters: A vector of dimensions",
                         "to use for the SNN graph. [Default %default]"))


pargs <- optparse::add_option(pargs, c("--resolution"),
            type = "double",
            default = 0.6,
            action = "store",
            dest = "resolution",
            help = paste("Seurat.FindClusters: resolution parameter,",
                         "use a value above(below) 1.0 if you want",
                         "to obtain a larger(smaller) number of communities.",
                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--thresh_use"),
            type = "double",
            default = 0.25,
            action = "store",
            dest = "thresh_use",
            help = paste("Seurat.FindMarkers: Limit testing to genes",
                         "with average of at least X-fold difference.",
                         "Increasing thresh_use speeds up the function,",
                         "but can miss weaker signals. [Default %default]"))

#process user-submitted inputs and set defaults
args_parsed <- optparse::parse_args(pargs, positional_arguments = TRUE)
args <- args_parsed$options
logging::loginfo(paste("Checking input arguments", sep = "" ))
args <- check_arguments(args)

# setup to capture outputs
pdf(args$plotfile)
sink(args$stdoutfile, append = FALSE, split = TRUE)

logging::loginfo(paste("Reading in 10X data", sep = ""))
tenX.data <- Read10X( data.dir = args$input_dir)

logging::loginfo(paste("Creating Seurat data object", sep = ""))
data_obj <- CreateSeuratObject(raw.data = tenX.data,
                               min.cells = args$min_cells,
                               min.genes = args$min_genes,
                               project = "10X_data")

logging::loginfo(paste("Filtering cells", sep = ""))
data_obj <- filter_cells(data_obj, args$filt_genes_low,
                             args$filt_genes_high,
                             args$filt_mito_low,
                             args$filt_mito_high)

logging::loginfo(paste("Finding variable genes", sep = ""))
data_obj <- find_variable_genes(data_obj,
                                 args$x_low_cutoff,
                                 args$x_high_cutoff,
                                 args$y_low_cutoff)

logging::loginfo(paste("Run PCA", sep = ""))
data_obj <- run_PCA(data_obj, args$pcs_print, args$genes_print)

if (args$run_js) {
  logging::loginfo(paste("Running and plotting jackstraw analysis (optional)"))
  data_obj <- run_JackStraw(data_obj, args$jackstraw_repl, args$jackstraw_pcs)
}

if (args$run_elbow_plot) {
  logging::loginfo(paste("Running elbow plot", sep = ""))
  run_elbow(data_obj)
}

logging::loginfo(paste("Finding clusters, generating tSNE", sep = ""))
data_obj <- find_clusters(data_obj, args$dims_use, args$resolution)

logging::loginfo(paste("Report all positive markers for each cluster"))
#reset output redirection to avoid capturing progress meter output
sink()
data_obj <- find_all_markers(data_obj, args$thresh_use)

#resume output capture
sink(args$stdoutfile, append = TRUE, split = TRUE)

# classify data
logging::loginfo(paste("Generate random forest classifier", sep = ""))
classifier <- classify_data(data_obj, args$seed)

logging::loginfo(paste("Output gene names for top markers",  sep = ""))
output_gene_names(classifier, args$num_genes, args$outputfile)

if (!is.null(args$save_obj)) {
  logging::loginfo(paste("Saving final Seurat object as", args$save_obj))
  saveRDS(data_obj, args$save_obj)
}

logging::loginfo(paste("R Session info, including Seurat version"))
sessionInfo()


logging::loginfo(paste("List of runtime variables (user inputs and defaults)"))
report_runtime_variables(args)

if (!is.null(args$save_obj)) {
  logging::loginfo(paste("Results in", args$outputfile, args$save_obj,
                       args$plotfile, "and", args$stdoutfile))
} else {
  logging::loginfo(paste("Results in", args$outputfile, args$plotfile, "and",
                        args$stdoutfile))
}

#reset output redirection
sink()
dev.off()