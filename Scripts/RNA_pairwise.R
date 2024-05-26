#Script to find number of up and down expressed genes per every pairwise comparison 

library(openxlsx)
library(pheatmap)
library(ggplot2)
library(grid)
#Read data
rna <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx", sheet = 8)

#Only columns needed are the gene names in NMB and NEIS, their products, expression values, and pairwise qvalues
rna <- rna[, c(8:10,12, 15:22, 101:156) ]

colnames(rna) <- c("CDS_code", "NEIS_id", "NMB_id", "Product", 
                   "exp_20578", "exp_22689", "exp_N188.1", "exp_N222.1", "exp_N222.2", "exp_N445.1", "exp_N459.3", "exp_N459.6",
                   "q20578_22689", "q20578_N188.1", "q20578_N222.1", "q20578_N222.2", "q20578_N445.1", "q20578_N459.3", "q20578_N459.6",
                   "q22689_20578", "q22689_N188.1", "q22689_N222.1", "q22689_N222.2", "q22689_N445.1", "q22689_N459.3", "q22689_N459.6",
                   "qN188.1_20578", "qN188.1_22689", "qN188.1_N222.1", "qN188.1_N222.2", "qN188.1_N445.1", "qN188.1_N459.3", "qN188.1_N459.6",
                   "qN222.1_20578", "qN222.1_22689", "qN222.1_N188.1", "qN222.1_N222.2", "qN222.1_N445.1", "qN222.1_N459.3", "qN222.1_N459.6",
                   "qN222.2_20578", "qN222.2_22689", "qN222.2_N188.1", "qN222.2_N222.1", "qN222.2_N445.1", "qN222.2_N459.3", "qN222.2_N459.6",
                   "qN445.1_20578", "qN445.1_22689", "qN445.1_N188.1", "qN445.1_N222.1", "qN445.1_N222.2", "qN445.1_N459.3", "qN445.1_N459.6",
                   "qN459.3_20578", "qN459.3_22689", "qN459.3_N188.1", "qN459.3_N222.1", "qN459.3_N222.2", "qN459.3_N445.1", "qN459.3_N459.6",
                   "qN459.6_20578", "qN459.6_22689", "qN459.6_N188.1", "qN459.6_N222.1", "qN459.6_N222.2", "qN459.6_N445.1", "qN459.6_N459.3"
                   )
#Remove first row
rna = rna[-c(1),]
rna[,5:68] <- as.numeric(rna[,5:68])

cols <- colnames(rna)[5:68]
rna[cols] <- sapply(rna[cols], as.numeric)


###Main part
#This is where we check every pairwise combination for a signficant logfold change and tally it

#Initialise matrix to store counts
res <- as.data.frame(matrix(NA, nrow = 56, ncol = 3))
colnames(res) <- c("Comparison", "upreg", 'downreg')
res$Comparison <- substr(colnames(rna)[13:68], 2, nchar(colnames(rna)[13:68]))

#Principle of the analysis:
#For a pairwise comparison:
#If there is a qvalue <0.05 at a locus
#And logfold change is >= 1
#Add 1 to the count for upreg or downreg, depending on direction of change
#All changes are relative to the first isolate in the rowname

q_check <- function(isolate1, isolate2, locus){
  colname = paste("q", isolate1, '_', isolate2, sep='')
  loc <- rna[[colname]][locus]
  if(loc < 0.05){
    return(TRUE)
  }
  return(FALSE)
}

exp_check <- function(isolate1, isolate2, locus){
  col1 = paste("exp_", isolate1, sep = '')
  col2 = paste("exp_", isolate2, sep = '')
  exp1 = rna[[col1]][locus]
  exp2 = rna[[col2]][locus]
  change = log2(exp1/exp2)
  if(abs(change) >= 1){
    if(change < 0){
      return("downreg")
    }
    return("upreg")
  }
  return(FALSE)
}

locus_loop <- function(isolate1, isolate2){
  upreg = 0
  downreg = 0
  for(i in 1:nrow(rna)){
    if(q_check(isolate1, isolate2, i) == TRUE && exp_check(isolate1, isolate2, i) == 'downreg'){
      downreg = downreg + 1
    }
    if(q_check(isolate1, isolate2, i) == TRUE && exp_check(isolate1, isolate2, i) == 'upreg'){
      upreg = upreg + 1
    }
    
  }
  
  return(c(upreg, downreg))
}

main <- function(isolates, res){
  r = 0
  for(i in isolates){
    print(r)
    for(j in isolates){
      if(i == j){
        next
      }
      r = r + 1
      out <- locus_loop(i,j)
      res[r,2] = out[1]
      res[r,3] = out[2]
      
      
    }
  }
  return(res)
}


isolates <- c("20578", "22689", "N188.1", "N222.1", "N222.2", "N445.1", "N459.3", "N459.6")
res <- main(isolates, res)
rownames(res) <- res[,1]
pheatmap(res[,2:3], scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = T, color= hcl.colors(50, "Blues 3"),
        angle_col = 45, fontsize_col = 13, fontsize_row = 11, display_numbers = T, fontsize_number = 9, number_color = 'black')
grid.text("Comparisons", x=-0.07, rot=90, gp=gpar(fontsize=16))
#visualise again removing all duplicated comparisons
res2 <- res
res2 <- res2[c(1:7, 9:14, 17:22, 25:28, 33:35, 41,42,49),]

pheatmap(res2[,2:3], scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = T, color= hcl.colors(50, "Blues 3"),
         angle_col = 45, fontsize_col = 13, fontsize_row = 11, display_numbers = T, fontsize_number = 9, number_color = 'black')
grid.text("Comparisons", x=-0.07, rot=90, gp=gpar(fontsize=16))