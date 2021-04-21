#preparation of metadata
prep_metadata <- function(count.matrix, ref.condition="dn"){
  meta.sample <- data.frame(
                            row.names = colnames( count.matrix ), 
                            sample_id = colnames( count.matrix ), 
                            condition = c(rep("dn",2),
                                          rep("dp",2),
                                          rep("gr",3),
                                          rep("or",3)),
                            batch = substr(colnames(count.matrix), 
                                           nchar(colnames(count.matrix)), 
                                           nchar(colnames(count.matrix))),
                            libType = rep("paired-end", 10) 
                            )
  
  # make experimental variables a factor and relevel 
  meta.sample$condition <- factor(meta.sample$condition)
  meta.sample$condition <- relevel(meta.sample$condition, ref = ref.condition)
  meta.sample$batch <- factor(meta.sample$batch)
  
  return(meta.sample)
}



