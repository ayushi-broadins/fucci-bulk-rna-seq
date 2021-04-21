## Get gene annotations --------------------------------------------------------
library(biomaRt)
library(plyr)
library(data.table)

get_gene_label<- function(count.file){
  outfile.raw <- "data/gene.annotation.raw.txt"
  outfile.final <- "data/gene.annotation.txt"
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
  return(gene.names)
}