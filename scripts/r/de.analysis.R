#load the required packages and functions
source("prepare.counts.for.de.R")
source("prepare.metadata.for.de.R")
source("get.gene.annotations.R")
library(DESeq2) #v1.28.1 
library(ggplot2) #v3.3.3
library("vsn") #3.56.0
library("pheatmap") #1.0.12
library("RColorBrewer")
setwd(paste0("../../results/",
             "exp01/5_count/"))

#---------------------------------------------------------------#

#prepare count matrix and metadata
count.mat <- prep_count(paste0("../../results/exp01/5_de/",
                                  "/data/counts.txt"))
meta.sample <- prep_metadata(count.mat,"dn")


#get gene annotations
gene.names <- get_gene_label(paste0("../../results/",
                                    "exp01/5_count/data/counts.txt"))

outfile <- "data/read.counts.gene.deseq2.Rdata"


#---------------------------------------------------------------#

#create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count.mat,
                              colData = meta.sample,
                              design= formula(~ batch + condition)
                              )
#save the count matrix and DESeqDataSet
save(count.mat, meta.sample, dds, file=outfile)


#---------------------------------------------------------------#


#iterate over the conditions to produce all combinations of DE lists
for(i in unique(dds$condition)){
  #relevel the conditions
  dds$condition <- relevel(dds$condition, ref = i)
  #run the DESeq function
  dds <- DESeq(dds)

  #Get DE lists for different comparisons across conditions
  for(j in resultsNames(dds)[-3:-1]){
    out.file <- substr(j,11,nchar(j))
    
    res <- results(dds, name=j, alpha=0.05)
    jpeg(paste0("plots/",out.file,"_MAplot_before",".jpg"), 
         width = 600, height = 600, quality = 100)
    plotMA(res, ylim=c(-5,5), cex.lab=1.5)
    dev.off()
    
    #perform LFC shrinkage using apeglm 
    #as suggested in the DESeq2 tutorial
    res <- lfcShrink(dds, coef=j, type="apeglm")
    jpeg(paste0("plots/",out.file,"_MAplot_after",".jpg"), 
         width = 600, height = 600, quality = 100)
    plotMA(res, ylim=c(-5,5), cex.lab=1.5)
    dev.off()
    
    res <- data.frame(ensembl_gene_id=rownames(res), res)
    resanno <- merge(res, gene.names[,c("gene_id", "hgnc_symbol")], 
                     by.x= "ensembl_gene_id",
                     by.y="gene_id", 
                     all.x = T)
    resanno <- resanno[order(resanno$padj),]
    out.file <- paste0("data/",out.file, "_de.txt")
    write.table(resanno, file=out.file, sep="\t", col.names=T, row.names=F, quote=F)
  }
}


#---------------------------------------------------------------#

#Data transformations and visualization

rld <- rlog(dds, blind=FALSE)

#effects of transformation on variance
jpeg(paste0("plots/","mean_sd_plot",".jpg"), 
     width = 600, height = 600, quality = 100)
meanSdPlot(assay(rld))
dev.off()

#heatmap of the count matrix
select_rows <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","batch")])
jpeg(paste0("plots/","heatmap_count_matrix",".jpg"), 
     width = 800, height = 800, quality = 100)
pheatmap(assay(rld)[select_rows,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

#heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$batch, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
jpeg(paste0("plots/","heatmap_sample_distances",".jpg"), 
     width = 800, height = 800, quality = 100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#PCA plot of samples
jpeg(paste0("plots/","pca_plot_samples",".jpg"), 
     width = 800, height = 800, quality = 100)
pcaData <- plotPCA(rld, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  stat_ellipse(type = "t", size=1) +
  geom_point(size=4.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  theme(legend.title = element_text(size = 18), legend.key.size = unit(1, "cm"), 
        legend.text = element_text(size = 18), legend.key = element_rect(size = 6),
        axis.title=element_text(size=18), axis.text=element_text(size=15))
dev.off()

#sessionInfo
writeLines(capture.output(sessionInfo()), "data/sessionInfo.txt")

#------------------------ END ----------------------------------#
#---------------------------------------------------------------#
