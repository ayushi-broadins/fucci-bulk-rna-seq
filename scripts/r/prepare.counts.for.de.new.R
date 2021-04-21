#preparation of count data
prep_count <- function(path.count.file){
  ## Read data
  count.matrix <- as.matrix(read.table(path.count.file, header=TRUE, row.names=1, sep="\t"))
  colnames(count.matrix)[1:4] <- c("DoubleNeg2", "DoubleNeg3", "DoublePos2", "DoublePos3")
  # remove "unmapped" counts etc. all genes that were never detected
  count.matrix <- subset(count.matrix, !grepl("N_", rownames(count.matrix)))
  count.matrix <- count.matrix[apply(count.matrix > 0, 1, any),]
  count.matrix <- count.matrix[,order(colnames(count.matrix))]
  return(count.matrix)
}
