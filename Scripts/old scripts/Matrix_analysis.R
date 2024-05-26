###Script to make matrices for each locus for both the igr and rna seq datasets and then use ANOVA tests
###to compare their structure. 
library(vegan)
library(openxlsx)
library(ggplot2)
#56,   32, 39
#55,   32, 39
#54    32, 39
#53    32, 39
#52.   42, 40
#51.   34, 34
#50    45, 40
#49.   40  42
#48    

setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data")
#Load igr data of only igrs with variation
up_igrs <- read.xlsx("igr_new_12102023.xlsx", sheet = 3)[, c(1,3, 12:19, 21)]
down_igrs <- read.xlsx("igr_new_12102023.xlsx", sheet = 3)[, c(1, 2, 4:11, 20)]
#Load expression data
rna <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx", sheet = 8)[, c(8:10, 101:156)]

#Change structure of rna object so that the isolate order in that is the same as in the igr objects
#Swap N188.1 (col17:23) and N222.1 (col24:30) as well as N222.2 (col31:37) and N445.1 (col38:44)
rna <- rna[,c(1:17, 25:31, 18:24, 39:45, 32:38, 46:ncol(rna))]
#this changes the order of missing isolates from 1, 10, 19, 28, 37, 46, 55, 64, to 1, 10, 20, 27, 38, 47, 55, 64
#So swap the columns within each pairwise comparison as well
rna <- rna[,c(1:4, 6,5,8,7,9:11, 13,12,15,14, 16:20, 22,21, 23:27,29,28,30:33,35,34,
              36:40, 42,41, 43:47, 49,48, 51,50, 52:54, 56,55, 58,57, 59)]

#Filter rna dataset to remove rows which are all 1s for their q values
removables <- c()
for(i in 2:nrow(rna)){
  if(sum(as.numeric(rna[i,4:ncol(rna)])) < 56){
    removables <- c(removables, rna$`11-7_gene`[i])
  }
}

rna <- rna[which(rna$`11-7_gene` %in% removables),]
rownames(rna) <- 1:nrow(rna)
# rna <- rna[-c(removables),]
# ncol(rna)


#define functions for the loop
igr_matrices <- function(dataset, j){
  skip = FALSE
  igr <- dataset[j, ]
  locus <- as.character(igr[1])
  igr <- igr[-c(1, 2, 11)]
  if (!all(is.na(igr))) {
    igr_matrix <- matrix(data = NA,
                         nrow = 8,
                         ncol = 8)
    for (k in 1:8) {
      for (l in 1:8) {
        if (is.na(igr[k]) | is.na(igr[l])) {
          next
        }
        if (igr[k] == igr[l]) {
          igr_matrix[k, l] = 1
        }
        if (igr[k] != igr[l]) {
          igr_matrix[k, l] = 0
        }
      }
    }
    colnames(igr_matrix) <- mnames
    rownames(igr_matrix) <- mnames
  }
  if (all(is.na(igr)) | sum(!is.na(igr)) == 1) {
    skip = TRUE
  }
  return(list(igr_matrix, skip, locus))
}
rna_matrices <- function(dataset, j){
  rna_locus <- rna[j, ][-c(1:3)]
  rna_locus <- as.numeric(rna_locus)
  rna_matrix <- matrix(data = NA,
                       nrow = 8,
                       ncol = 8)
  #need the names of the test matrix in pubmlst id form
  colnames(rna_matrix) <- mnames
  rownames(rna_matrix) <- mnames
  
  row = 1
  col = 1
  counter = 1
  for (k in 1:64) {
    row = ceiling(k / 8)
    col = k %% 8
    if (k %% 8 == 0) {
      col = 8
    }
    #Handle the cases where the pairwise comparison is between the same isolate
    if (k == seq(1, 64, 9)[counter]) {
      counter = counter + 1
      rna_matrix[row, col] = 1
      next
    }
    if (rna_locus[k - counter + 1] < 0.05 |
        rna_locus[k - counter + 1] > 1) {
      rna_matrix[row, col] = 0
    }
    if (rna_locus[k - counter + 1] >= 0.05) {
      rna_matrix[row, col] = 1
    }
  }
  return(rna_matrix)
}
matrix_stats <- function(dataset){
  
  rna_matrix2 <-
    vegdist(rna_matrix, method = 'jaccard', binary = TRUE)
  
  igr_matrix2 <-
    vegdist(
      up_igr_matrix,
      method = 'jaccard',
      binary = TRUE,
      na.rm = TRUE
    )
  
  if (skip_up == FALSE) {
    p <- mantel(rna_matrix2, igr_matrix2, na.rm = TRUE)[4]
    output <- as.numeric(p)
    locus <- locus
    return(list(output, locus))
  }
}
overlap_matrices <- function(dataset, threshold){
  overlap <- matrix(data = NA,
                    nrow = 8,
                    ncol = 8)
  for (i in 1:8) {
    for (j in 1:8) {
      if (i == j) {
        overlap[i, j] = NA
        next
      }
      if (is.na(rna_matrix[i, j]) | is.na(dataset[i, j])) {
        next
      }
      if (rna_matrix[i, j] == dataset[i, j]) {
        overlap[i, j] = 1
      }
      if (rna_matrix[i, j] != dataset[i, j]) {
        overlap[i, j] = 0
      }
    }
  }
  if (sum(overlap, na.rm = TRUE) >= threshold) {
    hit <- locus
    return(hit)
  }
}

#First make the igr matrix for a gene

test_igr <- down_igrs[926,]

test_locus <- as.character(test_igr[1])
test_locus
test_igr <- test_igr[-c(1,2,11)]
#make matrix of the test_igr
test_matrix <- matrix(data =NA, nrow = 8, ncol = 8)
for(i in 1:8){
  for(j in 1:8){
    if(test_igr[i] == test_igr[j]){
      test_matrix[i,j] = 1
    }
    if(test_igr[i] != test_igr[j]){
      test_matrix[i,j] = 0
    }
  }
}
colnames(test_matrix) <- c('27509', "27553", "28262", "28269", "28287", "53930", "53948", "53951")
rownames(test_matrix) <- c('27509', "27553", "28262", "28269", "28287", "53930", "53948", "53951")

data_melt <- reshape::melt(test_matrix)
data_melt$X1 <- as.character(data_melt$X1)
data_melt$X2 <- as.character(data_melt$X2)
ggplot(data_melt, aes(X1, X2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "purple3", high ="yellow2") +
  labs(x = "Isolate", y = "Isolate", title = test_locus)
igr_matrix = test_matrix

test <- matrix(data = c(1,0), nrow = 8, ncol = 8)
test
#Now make an expression matrix for the same gene

test_rna <- rna[38,] #NEIS2430
test_rna$`1717.genes`
test_rna <- test_rna[-c(1:3)]
test_rna <- as.numeric(test_rna)
test_matrix <- matrix(data=NA, nrow = 8, ncol = 8)
#need the names of the test matrix in pubmlst id form
colnames(test_matrix) <- c('27509', "27553", "28262", "28269", "28287", "53930", "53948", "53951")
rownames(test_matrix) <- c('27509', "27553", "28262", "28269", "28287", "53930", "53948", "53951")

row = 1
col = 1
counter = 1
for(i in 1:64){
  print(i)
  print(paste('index is', i-counter+1))
  row = ceiling(i/8)
  col = i%%8
  if(i %% 8 == 0){
    col = 8
  }
  #Handle the cases where the pairwise comparison is between the same isolate
  if(i == seq(1,64, 9)[counter]){
    counter = counter + 1
    test_matrix[row,col] = 1
    next
  }
  if(test_rna[i-counter+1] < 0.05 | test_rna[i-counter+1] > 1){
    test_matrix[row, col] = 0
  }
  if(test_rna[i-counter+1] >= 0.05){
    test_matrix[row, col] = 1
  }
}

rna_matrix = test_matrix
test_rna <- as.numeric(test_rna)
data_melt <- reshape::melt(test_matrix)
data_melt$X1 <- as.character(data_melt$X1)
data_melt$X2 <- as.character(data_melt$X2)
ggplot(data_melt, aes(X1, X2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "purple3", high ="yellow2") +
  labs(x = "Isolate", y = "Isolate", title = test_locus)

#Run stats on the structure of the matrices
rna_matrix2 <- vegdist(rna_matrix, method = 'jaccard', binary = TRUE)
igr_matrix2 <- vegdist(igr_matrix, method = 'jaccard', binary = TRUE)
x <- mantel(rna_matrix2, igr_matrix2)[4]
x

#Method to superimpose two matrices and convert to a single matrix of TRUE/FALSE depending on if they have identical values
overlap <- matrix(data = NA, nrow = 8, ncol = 8)
for(i in 1:8){
  for(j in 1:8){
    if(i == j){
      overlap[i,j] = NA
      next
    }
    if(rna_matrix[i,j] == igr_matrix[i,j]){
      overlap[i,j] = 1
    }
    if(rna_matrix[i,j] != igr_matrix[i,j]){
      overlap[i,j] = 0
    }
  }
}

data_melt <- reshape::melt(overlap)
data_melt$X1 <- as.character(data_melt$X1)
data_melt$X2 <- as.character(data_melt$X2)
ggplot(data_melt, aes(X1, X2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "red", high ="green") +
  labs(x = "Isolate", y = "Isolate", title = test_locus)

sum(overlap, na.rm = T)

#turn into a loop
out_up <- c()
out_down <- c()
loci_up <- c()
loci_down <- c()
no_match <- c()
hits_up <- c()
hits_down <- c()
mnames <- c('27509', "27553", "28262", "28269", "28287", "53930", "53948", "53951") #replace the below loop colnames = x entries with mnames
for(i in 2:length(rna$`1717.genes`)) {
  match = FALSE
  if (is.na(rna$`1717.genes`[i])) {
    next
  }
  for (j in 1:length(up_igrs$Locus)) {
    #if same locus, make matrices and perform stats, saving in table
    if (rna$`1717.genes`[i] == up_igrs$Locus[j]) {
      match = TRUE
      print('match')
      
      ###IGR_up matrix first
      out <- igr_matrices(up_igrs, j)
      up_igr_matrix <- out[[1]]
      skip_up <- out[[2]]
      locus <- out[[3]]
      
      ###IGR_down matrix next
      out <- igr_matrices(down_igrs, j)
      down_igr_matrix <- out[[1]]
      skip_down <- out[[2]]
      locus <- out[[3]]
      
      ###RNA Matrix finally
      rna_matrix <- rna_matrices(rna, j)
      
      ###Now perform the stats and save to dataframe
      out <- matrix_stats(up_igr_matrix)
      out_up <- c(out_up, out[[1]])
      loci_up <- c(loci_up, out[[2]])
      
      out <- matrix_stats(down_igr_matrix)
      out_down <- c(out_down, out[[1]])
      loci_down <- c(loci_down, out[[2]])

      
      #Now compute the overlap matrix
      hits_down <- c(hits_down, overlap_matrices(down_igr_matrix, 48))
      hits_up <- c(hits_up, overlap_matrices(up_igr_matrix, 48))

    }
    
  }
  if (j == length(up_igrs$Locus) && match == FALSE) {
    no_match <- c(no_match, rna$`1717.genes`[i])
  }
}

length(no_match)
length(no_match) + length(unique(c(loci_up, loci_down)))
length(rna$`1717.genes`)
#soem of the no)match entries are in both datasets. Not sure why theyre failing

test <- as.data.frame(cbind(loci_up, out_up))
test_down <- as.data.frame(cbind(loci_down, out_down))

p1
