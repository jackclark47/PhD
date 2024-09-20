#Script to find how mnay times each isolate is min or max expressed and the average expression of each isolate over all loci

library(openxlsx)

#Load rna data
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx', sheet = 8)[-1,]
rna$`1717.genes` <- str_trim(rna$`1717.genes`)
rna$X13 <- as.numeric(rna$X13)
rna$`log2(EXP/Min)` <- as.numeric(rna$`log2(EXP/Min)`)

#find min and max
mincols <- c()
maxcols <- c()
for(i in 1:nrow(rna)){
  counts <- unlist(rna[i, 15:22])
  mincount <- min(counts)
  min_index <- which(counts == mincount)[1]
  mincols <- c(as.numeric(min_index), mincols)
  
  maxcount <- max(counts)
  max_index <- which(counts == maxcount)[1]
  maxcols <- c(as.numeric(max_index), maxcols)
}

minfreq <- as.data.frame(table(mincols))
maxfreq <- as.data.frame(table(maxcols))

minfreq[,1] <- c('27509', '27553', '28269', '28262', '53930', '28287', '53948', '53951')
maxfreq[,1] <- c('27509', '27553', '28269', '28262', '53930', '28287', '53948', '53951')
sum(minfreq[,2])
sum(maxfreq[,2])

min_percent <- minfreq
max_percent <- maxfreq

min_percent[,2] <- 100* (min_percent[,2]/2271)
max_percent[,2] <- 100* (max_percent[,2]/2271)

rna[,15]
#calculate average expression
mean_exp <- function(colnum){
  return(mean(as.numeric(rna[,colnum])))
}

isolates <- c('27509', '27553', '28269', '28262', '53930', '28287', '53948', '53951')
mean_df <- as.data.frame(matrix(NA, nrow = 8, ncol = 2))
colnames(mean_df) <- c('isolate', 'mean_exp')
for(i in 15:22){
  avg <- mean_exp(i)
  mean_df[i-14,2] <- avg
}

mean_df[,1] <- isolates
