#!/usr/bin/env Rscript


# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(logging))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(Seurat))

 
# process_data function follows Seurat tutorial for processing PBMC data set made available by 10x Genomics
# dataset https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
# Seurat tutorial https://satijalab.org/seurat/pbmc3k_tutorial.html
process_data <- function(input.data){
	dense.size <- object.size(x = as.matrix(x = input.data))
	print(dense.size)
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
	PrintPCA(object = seuratObj, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
    viz <- VizPCA(object = seuratObj, pcs.use = 1:2)
    print(viz)
	pca <- PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2)
	seuratObj <- ProjectPCA(object = seuratObj, do.print = FALSE)
	seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, display.progress = FALSE)
	JackStrawPlot(object = seuratObj, PCs = 1:12)
	jackstraw <- JackStrawPlot(object = seuratObj, PCs = 1:12)
	print(jackstraw)
	elbow <- PCElbowPlot(object = seuratObj)
	print(elbow)
	seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
	PrintFindClustersParams(object = seuratObj)
	seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:10, do.fast = TRUE)
	TSNEPlot(object = seuratObj)
	# find all markers of cluster 1
 	cluster1.markers <- FindMarkers(object = seuratObj, ident.1 = 1, min.pct = 0.25)
	print(x = head(x = cluster1.markers, n = 5))
	# find all markers distinguishing cluster 5 from clusters 0 and 3
 	cluster5.markers <- FindMarkers(object = seuratObj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
 	print(x = head(x = cluster5.markers, n = 5))
	# find markers for every cluster compared to all remaining cells, report
	# only the positive ones
	seuratObj.markers <- FindAllMarkers(object = seuratObj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
	seuratObj.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
	cluster1.markers <- FindMarkers(object = seuratObj, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
	top10 <- seuratObj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
	return(seuratObj)
}

classify_data <- function(seuratObj){
	classifier_input_x = t(seuratObj@data)
	print(classifier_input_x)
	classifier_input_y = as.factor(as.numeric(seuratObj@ident))
	trained_classifier = randomForest(x = as.matrix(classifier_input_x) , y=factor(classifier_input_y) , importance = TRUE )
	print(trained_classifier) 
	ranked_list = trained_classifier$importance[order(-trained_classifier$importance[,'MeanDecreaseAccuracy']),]
	dput(rownames(ranked_list)[1:300]) 
	FINAL_LIST = rownames(ranked_list)[1:100]

	write.table(FINAL_LIST, file='FINAL_LIST.tsv', quote=FALSE, sep='\t', row.names = FALSE , col.names = TRUE )
}



# Command line arguments
pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                            "--input_dir directory"))



pargs <- optparse::add_option(pargs, c("--input_dir"),
                              type="character",
                              action="store",
                              dest="input_dir",
                              metavar="Input_Directory",
                              help=paste("Input directory for analysis.",
                                         "[Default %default][REQUIRED]"))



args_parsed <- optparse::parse_args(pargs)



# Make sure the input directory exists
if(!file.exists(args_parsed$input_dir[1])){
    error_message <- paste0("Provided input directory, ",
							args_parsed$input_dir[1], ",  does not exist.",
    						"Please check your input and try again")
    stop(error_message)
}


# # Manage inputs
tenX.data <- Read10X( data.dir = args_parsed$input_dir[1])

# process data
pdf("output.pdf")
procd_data <- process_data(tenX.data)
dev.off()

# classify data
classify_data(procd_data)