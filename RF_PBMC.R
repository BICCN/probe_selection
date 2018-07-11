#!/usr/bin/env Rscript


# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(logging))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(randomForest))

# # Logging level choices
# C_LEVEL_CHOICES <- names(loglevels)



# logging::basicConfig(level='INFO') #initialize to info setting.  

process_data <- function(input.data){
	dense.size <- object.size(x = as.matrix(x = input.data))
	print(dense.size)
# is ARGUMENT needed for min.cells, min.genes, project name
	seuratObj <- CreateSeuratObject(raw.data = input.data, min.cells = 3, min.genes = 200, project = "10X_data")
	print(seuratObj)
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratObj@data), value = TRUE)
	percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
	seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")
	vln <- VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
	print(vln)
	par(mfrow = c(1, 2))
	GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "percent.mito")
	GenePlot(object = seuratObj, gene1 = "nUMI", gene2 = "nGene")
	seuratObj <- FilterCells(object = seuratObj, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
	seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
	seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
	length(x = seuratObj@var.genes)
	seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("nUMI", "percent.mito"))
	seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
#    PrintPCA(object = seuratObj, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
#    viz <- VizPCA(object = seuratObj, pcs.use = 1:2)
#    print(viz)
#	pca <- PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2)
#	print(pca)
	seuratObj <- ProjectPCA(object = seuratObj, do.print = FALSE)
	seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, display.progress = FALSE)
	JackStrawPlot(object = seuratObj, PCs = 1:12)
#	jackstraw <- JackStrawPlot(object = seuratObj, PCs = 1:12)
#	print(jackstraw)
#	elbow <- PCElbowPlot(object = seuratObj)
#	print(elbow)
	seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#	PrintFindClustersParams(object = seuratObj)


	seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:10, do.fast = TRUE)
	TSNEPlot(object = seuratObj)
	
# # find all markers of cluster 1
# cluster1.markers <- FindMarkers(object = seuratObj, ident.1 = 1, min.pct = 0.25)
# print(x = head(x = cluster1.markers, n = 5))
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(object = seuratObj, ident.1 = 5, ident.2 = c(0, 3), 
#     min.pct = 0.25)
# print(x = head(x = cluster5.markers, n = 5))
# # find markers for every cluster compared to all remaining cells, report
# # only the positive ones
# seuratObj.markers <- FindAllMarkers(object = seuratObj, only.pos = TRUE, min.pct = 0.25, 
#     thresh.use = 0.25)
# seuratObj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
# cluster1.markers <- FindMarkers(object = seuratObj, ident.1 = 0, thresh.use = 0.25, 
#     test.use = "roc", only.pos = TRUE)
# top10 <- seuratObj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# 
# classifier_input_x = t(seuratObj@data)
# classifier_input_x
# classifier_input_y = as.factor(as.numeric(seuratObj@ident))
# trained_classifier = randomForest(x = as.matrix(classifier_input_x) , y=factor(classifier_input_y) , importance = TRUE )
# trained_classifier 
# ranked_list = trained_classifier$importance[order(-trained_classifier$importance[,'MeanDecreaseAccuracy']),]
# dput(rownames(ranked_list)[1:300]) 
# 
# ### end of cell 1

}

# #' From https://github.com/broadinstitute/inferCNV/blob/master/scripts/inferCNV.R
# #'
# #' Check arguments and make sure the user input meet certain 
# #' additional requirements.
# #'
# #' Args:
# #'    @param arguments: Parsed arguments from user.
# #'
# #' Returns:
# #'    @return: Updated arguments
# check_arguments <- function(arguments){
# 
#     logging::loginfo(paste("::check_arguments:Start", sep=""))
#     # Require the name of a output pdf file
#     if ( (!( "output_dir" %in% names(arguments))) || (arguments$output_dir == "") || (is.na(arguments$output_dir)) ) {
#         logging::logerror(paste(":: --output_dir: Please enter a file path to ",
#                                 "save the heatmap.",
#                                  sep=""))
# 
#         stop("error, no --output_dir")
#     }
# 
#     # Require the cut off to be above 0
#     if (arguments$cutoff < 0){
#         logging::logerror(paste(":: --cutoff: Please enter a value",
#                                 "greater or equal to zero for the cut off.",
#                                 sep=""))
# 
#         stop("error, no --cutoff")
#     }
# 
#     # Require the logging level to be one handled by logging
#     if (!(arguments$log_level %in% C_LEVEL_CHOICES)){
#         logging::logerror(paste(":: --log_level: Please use a logging level ",
#                                 "given here: ", C_LEVEL_CHOICES,
#                                 collapse=",", sep=""))
#         stop("error, not recognizing log level")
#     }
#     
#     # Require the visualization outlier detection to be a correct choice.
#     if (!(arguments$bound_method_vis %in% C_VIS_OUTLIER_CHOICES)){
#         logging::logerror(paste(":: --vis_bound_method: Please use a method ",
#                                 "given here: ", C_VIS_OUTLIER_CHOICES,
#                                 collapse=",", sep=""))
#         stop("error, must specify acceptable --vis_bound_method")
#     }
# 
#     if (! (arguments$ref_subtract_method %in% C_REF_SUBTRACT_METHODS) ) {
#         logging::logerror(paste(":: --ref_subtract_method: acceptable values are: ",
#                                 paste(C_REF_SUBTRACT_METHODS, collapse=","), sep="") )
#         stop("error, must specify acceptable --ref_subtract_method")
#     }
# 
# 
#     if (! (arguments$hclust_method %in% C_HCLUST_METHODS) ) {
#         logging::logerror(paste(":: --hclust_method: acceptable values are: ",
#                                 paste(C_HCLUST_METHODS, collapse=","), sep="") )
#         stop("error, must specify acceptable --hclust_method")
#     }
#         
#     # Warn that an average of the samples is used in the absence of
#     # normal / reference samples
#     if (is.null(arguments$reference_observations)){
#         logging::logwarn(paste(":: --reference_observations: No reference ",
#                       "samples were given, the average of the samples ",
#                       "will be used.",
#                       sep=""))
#     }
# 
#     # Make sure the threshold is centered.
#     arguments$max_centered_expression <- abs(arguments$max_centered_expression)
#     arguments$magnitude_filter <- abs(arguments$magnitude_filter)
# 
#     # Require the contig tail to be above 0
#     if (is.na(arguments$contig_tail)){
#         arguments$contig_tail <- (arguments$window_length - 1) / 2
#     }
# 
#     if (arguments$contig_tail < 0){
#         logging::logerror(paste(":: --tail: Please enter a value",
#                                 "greater or equal to zero for the tail.",
#                                 sep=""))
# 
#         stop(980)
#     }
# 
#     if (! is.na(suppressWarnings(as.integer(arguments$num_groups)))){
#         arguments$num_groups <- list(as.integer(arguments$num_groups))
#     } else {
#         # Warn references must be given.
#         if (is.null(arguments$reference_observations)){
#             logging::logerror(paste(":: --ref_groups to use this function ",
#                                     "references must be given. "))
#         }
# 
#         # TODO need to check and make sure all reference indices are given.
#         num_str <- unlist(strsplit(arguments$num_groups,","))
#         if (length(num_str) == 1){
#             logging::logerror(paste(":: --ref_groups. If explicitly giving ",
#                                     "indices, make sure to give atleast ",
#                                     "two groups", sep =""))
#             stop(990)
#         }
# 
#         num_groups <- list()
#         for (num_token in num_str){
#             token_numbers <- unlist(strsplit(num_token, ":"))
#             number_count <- length(token_numbers)
#             if (number_count == 1){
#                 singleton <- as.integer(number_count)
#                 num_groups[[length(num_groups) + 1]] <- singleton
#             } else if (number_count == 2){
#                 from <- as.integer(token_numbers[1])
#                 to <- as.integer(token_numbers[2])
#                 num_groups[[length(num_groups) + 1]] <- seq(from, to)
#             } else {
#                 logging::logerror(paste(":: --ref_groups is expecting either ",
#                                         "one number or a comma delimited list ",
#                                         "of numbers or spans using ':'. ",
#                                         "Examples include: --ref_groups 3 or ",
#                                         " --ref_groups 1,3,5,6,3 or ",
#                                         " --ref_groups 1:5,6:20 or ",
#                                         " --ref_groups 1,2:5,6,7:10 .", sep=""))
#                 stop(999)
#             }
#         }
#         arguments$num_groups <- num_groups
#     }
#     return(arguments)
# }


# Command line arguments
pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                            "--input_dir directory"))

# pargs <- optparse::add_option(pargs, c("--color_safe"),
#                               type="logical",
#                               default=FALSE,
#                               action="store_true",
#                               dest="use_color_safe",
#                               metavar="Color_Safe",
#                               help=paste("To support the needs of those who see ",
#                                          "colors differently, use this option to",
#                                          "change the colors to a palette visibly ",
#                                          "distinct to all color blindness. ",
#                                          " [Default %default]"))
# 
# pargs <- optparse::add_option(pargs, c("--contig_lab_size"),
#                               type="integer",
#                               action="store",
#                               default=1,
#                               dest="contig_label_size",
#                               metavar="Contig_Label_Size",
#                               help=paste("Used to increase or decrease the text labels",
#                                          "for the X axis (contig names).",
#                                          "[Default %default]"))
# 
# pargs <- optparse::add_option(pargs, c("--cutoff"),
#                               type="numeric",
#                               default=0,
#                               action="store",
#                               dest="cutoff",
#                               metavar="Cutoff",
#                               help=paste("A number >= 0 is expected. A cut off for",
#                                          "the average expression of genes to be used",
#                                          "for CNV inference (use the value before log2 transformation). [Default %default]"))
# 
# 
# pargs <- optparse::add_option(pargs, c("--log_file"),
#                               type="character",
#                               action="store",
#                               default=NA,
#                               dest="log_file",
#                               metavar="Log",
#                               help=paste("File for logging. If not given,",
#                                          "logging will occur to console.",
#                                          "[Default %default]"))
# 
# pargs <- optparse::add_option(pargs, c("--delim"),
#                               type="character",
#                               action="store",
#                               default="\t",
#                               dest="delim",
#                               metavar="Delimiter",
#                               help=paste("Delimiter for reading expression matrix",
#                                         " and writing matrices output.",
#                                          "[Default %default]"))
# 
# pargs <- optparse::add_option(pargs, c("--log_level"),
#                               type="character",
#                               action="store",
#                               default="INFO",
#                               dest="log_level",
#                               metavar="LogLevel",
#                               help=paste("Logging level. Valid choices are",
#                                          paste(C_LEVEL_CHOICES,collapse=", "),
#                                          "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--input_dir"),
                              type="character",
                              action="store",
                              dest="input_dir",
                              metavar="Input_Directory",
                              help=paste("Input directory for analysis.",
                                         "[Default %default][REQUIRED]"))


# pargs <- optparse::add_option(pargs, c("--save"),
#                               type="logical",
#                               action="store_true",
#                               default=FALSE,
#                               dest="save",
#                               metavar="save",
#                               help="Save workspace as infercnv.Rdata")


args_parsed <- optparse::parse_args(pargs)

# , positional_arguments=2)
args <- args_parsed$options
print(args_parsed)
print(args_parsed["input_dir"])
# args["input_matrix"] <- args_parsed$args[1]
# args["gene_order"] <- args_parsed$args[2]

# # Check arguments
# args <- check_arguments(args)


# Make sure the input directory exists
if(!file.exists(args_parsed$input_dir[1])){
    error_message <- paste0("Provided input directory, ",
							args_parsed$input_dir[1], ",  does not exist.",
    						"Please check your input and try again")
    stop(error_message)
}

# # Parse bounds
# bounds_viz <- c(NA,NA)
# if (!is.na(args$bound_threshold_vis)){
#     bounds_viz <- as.numeric(unlist(strsplit(args$bound_threshold_vis,",")))
# }
# if (length(bounds_viz) != 2){
#     error_message <- paste("Please use the correct format for the argument",
#                            "--vis_bound_threshold . Two numbers seperated",
#                            "by a comma is expected (lowerbound,upperbound)",
#                            ". As an example, to indicate that outliers are",
#                            "outside of -1 and 1 give the following.",
#                            "--vis_bound_threshold -1,1")
#     stop(error_message)
# }

# # Set up logging file
# logging::basicConfig(level=args$log_level)
# if (!is.na(args$log_file)){
#     logging::addHandler(logging::writeToFile,
#                         file=args$log_file,
#                         level=args$log_level)
# }
# 
# # Log the input parameters
# logging::loginfo(paste("::Input arguments. Start.")) 
# for (arg_name in names(args)){
#     logging::loginfo(paste(":Input_Argument:",arg_name,"=",args[[arg_name]],
#                            sep="")) 
# }
# logging::loginfo(paste("::Input arguments. End.")) 

# # Manage inputs
# logging::loginfo(paste("::Reading data matrix.", sep=""))
# # Row = Genes/Features, Col = Cells/Observations
tenX.data <- Read10X( data.dir = args_parsed$input_dir[1])
# logging::loginfo(paste("Original matrix dimensions (r,c)=",
#                  paste(dim(expression_data), collapse=",")))

# process data
pdf("output.pdf")
procd_data <- process_data(tenX.data)
dev.off()

# # Read in the gen_pos file
# input_gene_order <- seq(1, nrow(expression_data), 1)
# if (args$gene_order != ""){
#     input_gene_order <- read.table(args$gene_order, header=F, row.names=1, sep="\t")
#     names(input_gene_order) <- c(CHR, START, STOP)
# }
# logging::loginfo(paste("::Reading gene order.", sep=""))
# logging::logdebug(paste(head(args$gene_order[1]), collapse=","))
# 
# # Default the reference samples to all
# input_reference_samples <- colnames(expression_data)
# if (!is.null(args$reference_observations)){
#     # This argument can be either a list of column labels
#     # which is a comma delimited list of column labels
#     # holding a comma delimited list of column labels
#     refs <- args$reference_observations
#     if (file.exists(args$reference_observations)){
#         refs <- scan(args$reference_observations,
#                      what="character",
#                      quiet=TRUE)
#         refs <- paste(refs, collapse=",")
#     }
#     # Split on comma
#     refs <- unique(unlist(strsplit(refs, ",", fixed=FALSE)))
#     # Remove multiple spaces to single spaces
#     refs <- unique(unlist(strsplit(refs, " ", fixed=FALSE)))
#     refs <- refs[refs != ""]
#     # Normalize names with make.names so they are treated
#     # as the matrix column names
#     refs <- make.names(refs)
#     if (length(refs) > 0){
#         input_reference_samples <- refs
#     }
#     logging::logdebug(paste("::Reference observations set to: ", input_reference_samples, collapse="\n"))
# }

# # Make sure the given reference samples are in the matrix.
# if (length(input_reference_samples) !=
#     length(intersect(input_reference_samples, colnames(expression_data)))){
#     missing_reference_sample <- setdiff(input_reference_samples,
#                                         colnames(expression_data))
#     error_message <- paste("Please make sure that all the reference sample",
#                            "names match a sample in your data matrix.",
#                            "Attention to: ",
#                            paste(missing_reference_sample, collapse=","))
#     logging::logdebug(paste("::colnames(expression_data): ", colnames(expression_data), collapse="\n"))
#     logging::logerror(error_message)
#     stop(error_message)
# }

# # Order and reduce the expression to the genomic file.
# order_ret <- infercnv::order_reduce(data=expression_data,
#                                     genomic_position=input_gene_order)
# expression_data <- order_ret$expr
# input_gene_order <- order_ret$order
# if(is.null(expression_data)){
#     error_message <- paste("None of the genes in the expression data",
#                            "matched the genes in the reference genomic",
#                            "position file. Analysis Stopped.")
#     stop(error_message)
# }


# if (args$save) {
#     logging::loginfo("Saving workspace")
#     save.image("infercnv.Rdata")
# }

# # Run CNV inference
# ret_list = infercnv::infer_cnv(data=expression_data,
#                                gene_order=input_gene_order,
#                                cutoff=args$cutoff,
#                                reference_obs=input_reference_samples,
#                                transform_data=args$log_transform,
#                                window_length=args$window_length,
#                                max_centered_threshold=args$max_centered_expression,
#                                noise_threshold=args$magnitude_filter,
#                                num_ref_groups=args$num_groups,
#                                out_path=args$output_dir,
#                                k_obs_groups=args$num_obs,
#                                plot_steps=args$plot_steps,
#                                contig_tail=args$contig_tail,
#                                method_bound_vis=args$bound_method_vis,
#                                lower_bound_vis=bounds_viz[1],
#                                upper_bound_vis=bounds_viz[2],
#                                ref_subtract_method=args$ref_subtract_method,
#                                hclust_method=args$hclust_method)

# # Log output
# logging::loginfo(paste("::infer_cnv:Writing final data to ",
#                        file.path(args$output_dir,
#                        "expression_pre_vis_transform.txt"), sep="_"))
# # Output data before viz outlier
# write.table(ret_list["PREVIZ"], sep=args$delim,
#             file=file.path(args$output_dir,
#                        "expression_pre_vis_transform.txt"))
# # Output data after viz outlier
# write.table(ret_list["VIZ"], sep=args$delim,
#             file=file.path(args$output_dir,
#                        "expression_post_viz_transform.txt"))
# logging::loginfo(paste("::infer_cnv:Current data dimensions (r,c)=",
#                        paste(dim(ret_list[["VIZ"]]), collapse=","), sep=""))
# 
# logging::loginfo(paste("::infer_cnv:Drawing plots to file:",
#                            args$output_dir, sep=""))


# if (args$save) {
#     logging::loginfo("Saving workspace")
#     save.image("infercnv.Rdata")
# }

