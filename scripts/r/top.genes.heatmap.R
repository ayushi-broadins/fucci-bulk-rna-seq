library("gplots")
library(RColorBrewer)

source(paste0("//hydrogen/chembio_datasets/csdev/",
              "AA/betacellproliferation/rnaseq_fucci/",
              "scripts/r/prepare.counts.for.de.new.R"))

setwd("//hydrogen/chembio_datasets/proj/bcproli/rnaseq/bulk/results/exp01/5_de/data")

#read in the raw counts file
raw_counts <- prep_count("counts.txt")
out.path <- "//hydrogen/chembio_datasets/proj/bcproli/rnaseq/bulk/results/exp01/5_de/sandbox/"

#get the names of all the DE gene lists
lst <- list.files(pattern='_de.txt')
#genes <- NULL

for(i in lst){
  x <- read.csv(i, sep='\t', stringsAsFactors = F)
  #get only abs(loffoldchange) >= 2
  x <- x[abs(x$log2FoldChange) >= 2,]
  #order by padj
  x <- x[order(x$padj),]
  #get the top 10 genes
  x.genes <- x[1:10,]
  print(x.genes)
  dat <- as.data.frame(raw_counts[unique(x.genes$ensembl_gene_id),])
  dat <- dat[,c("DoubleNeg2",
                "DoubleNeg3",
                "Orange1",
                "Orange2",
                "Orange3",
                "DoublePos2",
                "DoublePos3",
                "Green1",
                "Green2",
                "Green3",
                "DoubleNeg2",
                "DoubleNeg3")]
  dat$ensembl_gene_id <- rownames(dat)
  dat <- merge(dat,x[,c('ensembl_gene_id','hgnc_symbol')], by='ensembl_gene_id')
  rownames(dat) <- dat$hgnc_symbol
  dat <- as.matrix(dat[,c(2:13)])
  colnames(dat) <- c('DN2','DN3','Or1','Or2','Or3','DP2','DP3','Gr1','Gr2','Gr3','DN2','DN3')
  x <- scale(t(as.matrix(dat)))
  file.name <- substr(i,1,8)
  jpeg(paste0(out.path,file.name,".jpeg"),
       width = 700, height = 500,
       pointsize = 15)
  heatmap.2(t(x),
            main = paste0("Top 10 DE Genes (",file.name,")"),
            scale     = "none",
            trace     = "none",
            #margins =c(12,9), 
            col       = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
            distfun   = function(x) as.dist(1-cor(t(x))), 
            hclustfun = function(x) hclust(x, method="ave"),
            Colv = FALSE,
            dendrogram = "row",
            ColSideColors = c(rep('grey',2),
                              rep('orange',3),
                              rep('yellow',2),
                              rep('green',3),
                              rep('grey',2))
            
  )
  
  legend("topright",    
         xpd=TRUE,
         inset=c(0,-0.07),
         legend = c("G0 phase",
                    "G1 phase (mature FUCCI cells)",
                    "G1/S phase transition",
                    "S/G2/M (proliferating FUCCI)"),
         col = c("grey",
                 "orange",
                 "yellow",
                 "green"), 
         lty= 1,             
         lwd = 5,           
         cex=0.8,
         bty="o"

  )
  
dev.off()

}  



#############################################
#GeLiNea results
#############################################

library(ggplot2)
setwd("//hydrogen/chembio_datasets/proj/bcproli/rnaseq/bulk/results/exp01/5_de/sandbox/AA/gelinea_results")
lst <- list.files()
x1 <- read.csv(lst[1])
x2 <- read.csv(lst[2])
x3 <- read.csv(lst[3])
x <- rbind.data.frame(x1,x2,x3)
x <- x[x$GeLiNEA_p_adj <= 0.05,]
x <- x[order(x$GeLiNEA_p_adj),]

ggplot(data=x[1:20,], aes(x=reorder(id,gene_list_overlap), y=gene_list_overlap, fill = GeLiNEA_p_adj)) +
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw()+
  xlab("gene set")+
  ylab("gene list overlap") +
  ggtitle("GeLiNEA pathway enrichment (Green vs Orange)")+
  scale_fill_continuous(type = "viridis")
