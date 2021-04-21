#load required packages
library(DESeq2)
library(ggplot2)
require(data.table)
require(yaml)

#set the owrking directory
setwd("//hydrogen/chembio_datasets/csdev/AA/betacellproliferation/rnaseq_fucci/results/exp01/5_count")


## FUNCTIONS ------------------------------------------------------------------

prep.counts.deseq2 <- function(count.matrix, meta.sample, design.formula, outfile){
  
  dds <- DESeqDataSetFromMatrix(count.matrix, meta.sample, design = formula(design.formula))
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  
  colData(dds)$name <- colData(dds)$sample
  
  normLogCounts <- log2(fpm(dds)+0.25)
  colnames(normLogCounts) <- colData(dds)$name
  
  save(count.matrix, meta.sample, normLogCounts, dds, file=outfile)
}


## SCRIPT ---------------------------------------------------------------------

count.file <- "counts.txt"
model.factor <- " batch + condition"
ref.level <- "dn"
outfile <- "read.counts.gene.deseq2.Rdata"


## Read data
count.matrix <- as.matrix(read.table(count.file, header=TRUE, row.names=1, sep="\t"))
colnames(count.matrix)[1:4] <- c("DoubleNeg2", "DoubleNeg3", "DoublePos2", "DoublePos3")

## prepare meta data
meta.sample = data.frame(
  row.names = colnames( count.matrix ),
  sample_id = colnames(count.matrix),
  condition = c(rep("dn",2),rep("dp",2),rep("gr",3),rep("or",3)),
  batch = substr(colnames(count.matrix), nchar(colnames(count.matrix)), nchar(colnames(count.matrix))),
  libType = rep("paired-end", 10) )

# make experimental variable a factor and relevel 
meta.sample$condition <- as.factor(meta.sample$condition)
meta.sample$condition <- relevel(meta.sample$condition, "dn")
meta.sample$batch <- as.factor(meta.sample$batch)
meta.sample$batch <- relevel(meta.sample$batch, "3")

# remove "unmapped" counts etc. all genes that were never detected
count.matrix <- subset(count.matrix, !grepl("N_", rownames(count.matrix)))
count.matrix <- count.matrix[apply(count.matrix > 0, 1, any),]
count.matrix <- count.matrix[,order(colnames(count.matrix))]

# DESeq2
stopifnot(colnames(count.matrix) == meta.sample$sample_id)
design.formula <- paste0("~", model.factor)
prep.counts.deseq2(count.matrix, meta.sample, design.formula, outfile)

## view sample counts ----------------------------------------------------------
pseudoCount = log2(count.matrix + 1)
df <- reshape2::melt(pseudoCount)
colnames(df)[2] <- "Samples"
df = data.frame(df, Condition = substr(df$Samples, 1, nchar(as.character(df$Samples))-1))
ggplot(df, aes(x = Samples, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1))) 

pseudoCount = log2(counts(dds, normalized = TRUE) + 1)
df <- reshape2::melt(pseudoCount)
colnames(df)[2] <- "Samples"
df = data.frame(df, Condition = substr(df$Samples, 1, nchar(as.character(df$Samples))-1))
ggplot(df, aes(x = Samples, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1)))

## Prepare DE lists-------------------------------------------------------------
count.file <- "read.counts.gene.deseq2.Rdata"

load(count.file)
## calculate differential expression for all samples
nl <- nlevels(colData(dds)$condition) - 2
fl <- c("gr","or")#levels(colData(dds)$condition)
for(j in 1:nl){
  for (i in c(1:nl)[-j]){
    res <- results(dds, contrast = c("condition", fl[i], fl[j]))
    res <- data.frame(gene_id=rownames(res), res)
    res <- res[order(res$pvalue),]
    out.file <- paste0(fl[i],"_vs_",fl[j], ".de.txt")
    write.table(res, file=out.file, sep="\t", col.names=T, row.names=F, quote=F)
  }
}


## Get gene annotations --------------------------------------------------------
library(biomaRt)
library(plyr)
library(data.table)

count.file <- "counts.txt"
outfile.raw <- "gene.annotation.raw.txt"
outfile.final <- "gene.annotation.txt"

# get gene ids
gene_id <- fread(count.file, select = "gene_id")
gene_id[, ensembl_gene_id:=sub("\\..*", "", gene_id)]

# define biomart object
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")

# query biomart (split in to queries because of 3-external-attribute limit)
gene.names.mrna <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "refseq_mrna", "refseq_mrna_predicted"),
                         filters = "ensembl_gene_id", values = gene_id$ensembl_gene_id,
                         mart = mart)

gene.names.nc <- getBM(attributes = c("ensembl_gene_id", "refseq_ncrna", "refseq_ncrna_predicted"),
                       filters = "ensembl_gene_id", values = gene_id$ensembl_gene_id,
                       mart = mart)

# merge mRNA and ncRNA results
gene.names <- merge(gene_id, gene.names.mrna, by="ensembl_gene_id", all.x=T, allow.cartesian=TRUE)
gene.names <- merge(gene.names, gene.names.nc, by="ensembl_gene_id", all.x=T, allow.cartesian=TRUE)
gene.names[is.na(gene.names)] <- " "

#sort on gene id
gene.names <- gene.names[order(gene.names$gene_id),]

# write all annotations to file
write.table(gene.names, outfile.raw, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)



# prepare single refseq id
gene.names$refseq_evidence <- 0
gene.names$refseq_evidence[gene.names$refseq_ncrna_predicted != ""] <- 1
gene.names$refseq_evidence[gene.names$refseq_mrna_predicted != ""] <- 2
gene.names$refseq_evidence[gene.names$refseq_ncrna != ""] <- 3
gene.names$refseq_evidence[gene.names$refseq_mrna != ""] <- 4

gene.names <- ddply(gene.names, .(ensembl_gene_id), function(x) x[which.max(x$refseq_evidence),])

gene.names$refseq <- gene.names$refseq_mrna
gene.names$refseq[gene.names$refseq == ""] <- gene.names$refseq_ncrna[gene.names$refseq == ""]
gene.names$refseq[gene.names$refseq == ""] <- gene.names$refseq_mrna_predicted[gene.names$refseq == ""]
gene.names$refseq[gene.names$refseq == ""] <- gene.names$refseq_ncrna_predicted[gene.names$refseq == ""]

#sort on gene id
gene.names <- gene.names[order(gene.names$gene_id),]

# write cleaned annotations to file
write.table(gene.names[,c("gene_id", "ensembl_gene_id", "hgnc_symbol", "refseq", "refseq_evidence")], 
            outfile.final, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


## annotate de list ---------------------------------------------------------------------
require(data.table)

infiles <- c("gr_vs_or.de.txt","or_vs_gr.de.txt")#list.files(pattern = ".de.txt")

for(i in 1:length(infiles)){
  outfile <- paste0(substr(infiles[i],1,nchar(infiles[i])-4),".annotated.txt")
  res <- fread(infiles[i])
  resanno <- merge(res, gene.names[,c("gene_id", "hgnc_symbol")], by="gene_id", all.x = T)
  resanno <- resanno[order(resanno$padj),]
  write.table(resanno, outfile, row.names=F, quote=F, sep="\t")
  
  #plot the volcano plot
  #volcano plots
  jpeg(paste0(substr(infiles[i],1,nchar(infiles[i])-4),".jpg"), 
       width = 1800, height = 1200)
  par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
  topT <- as.data.frame(resanno)
  #Adjusted P values (FDR Q values)
  with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=2, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
  with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=2))
  with(subset(topT, padj<0.05 & abs(log2FoldChange)>2)[1:10,], text(log2FoldChange, -log10(padj), 
                                                              labels=subset(topT$hgnc_symbol, topT$padj<0.05 & abs(topT$log2FoldChange)>2)[1:10], cex=1.5, pos=3))
  #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
  abline(v=0, col="black", lty=3, lwd=1.0)
  abline(v=-2, col="black", lty=4, lwd=2.0)
  abline(v=2, col="black", lty=4, lwd=2.0)
  abline(h=-log10(0.05), col="black", lty=4, lwd=2.0)
  dev.off()
}

## plot counts of proliferation genes ----------------------------------------------------

for(i in proli.genes){
  geneid <- gene.names[gene.names$hgnc_symbol == i,'gene_id']
  jpeg(paste0("normcounts_",i,".jpg"),
       width = 1200, height = 1500)
  plotCounts(dds, gene=i, intgroup="condition")
  dev.off()
}

#pca plots
rld <- rlog(dds, blind=FALSE)
plotPCA( rld, intgroup = "condition") + 
  #stat_ellipse(type = "norm", linetype = 2, size =1) + 
  stat_ellipse(type = "t", size=1) + 
  geom_point(size=4.5) + 
  theme(legend.title = element_text(size = 18), legend.key.size = unit(1, "cm"), 
        legend.text = element_text(size = 18), legend.key = element_rect(size = 6),
        axis.title=element_text(size=18), axis.text=element_text(size=15))

