###Script to quantify nucleotide diversity between the 8 isolates for each igr for each locus
###APPROACH
#Use RNA data to find the reference isolate for a given locus
#this is chosen based on which isolate had the lowest expression at that locus
#Then using this isolate as a reference, measure the number of nucleotide differences in the up and downstream igrs at that locus,
#relative to that isolate. 

#Load libraries
library(seqinr)
library(pegas)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(stringr)

#load data
rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 6)[, c(8,9,15:22, 25:32)]
up_igrs <- read.FASTA("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_up.xmfa")
down_igrs <- read.FASTA("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_down.xmfa")

#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

#Order of isolates in both datasets:
#RNA: 20578 22698 N188.1  N222.1, N222.2, N445.1, N459.3, N459.6
#IGR: 20578 22698 N222.1  N188.1  N445.1, N222.2  N459.3, N459.6

#So rearrange columns in the rna dataset so theyre the same order
rna <- rna[,c(1:4, 6,5, 8,7, 9:12, 14,13,16,15, 17, 18)]

#Create empty matrix to store pi values
pi_values <- as.data.frame(matrix(data = NA, nrow = (length(up_igrs)/8), ncol = 5))
colnames(pi_values) <- c('Locus', 'pi_up', 'variance_up', 'pi_down', 'variance_down')

count = 0
for(i in seq(1,length(up_igrs), 8)){
  count = count+1
  #store the locus' name as an object
  locusname <- names(up_igrs[i])
  locusname <- substr(locusname, start = 1, stop = nchar(locusname)-6)
  
  #extract the igrs from each isolate for the current locus
  testlocus_up <- up_igrs[i:(i+7)]
  testlocus_down <- down_igrs[i:(i+7)]
  
  #if the upstream igrs are identical in length, calculate the nucleotide diversity for them
  if(length(unique(sapply(testlocus_up, length))) == 1){
    out <- nuc.div(testlocus_up, variance = TRUE)
    pi_values$Locus[count] <- locusname
      pi_values$pi_up[count] <- out[1]
    pi_values$variance_up[count] <- out[2]
  }
  
  #if the up_igrs are not all the same length, record only the locus name
  if(length(unique(sapply(testlocus_up, length))) > 1){
    pi_values$Locus[count] <- locusname
  }
  
  #if the down_igrs are all the same length, calculate the nucletodie diversity 
  if(length(unique(sapply(testlocus_down, length))) == 1){
    out <- nuc.div(testlocus_down, variance = TRUE)
    pi_values$pi_down[count] <- out[1]
    pi_values$variance_down[count] <- out[2]
  }
  #if not then record only the locus name
  if(length(unique(sapply(testlocus_down, length))) > 1){
    pi_values$Locus[count] <- locusname
  }
  #Note, nuc.div can only be calculated for sequences of the same length, hence the above
}
#save this table to a file
write.xlsx(pi_values, "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/pi_values.xlsx")

#Select a locus and identify the reference isolate
testlocus_up <- up_igrs[8585:8592]
nuc.div(testlocus_up, variance = TRUE)
names(up_igrs[1])

hist(pi_values$pi_down)
hist(pi_values$pi_up[pi_values$pi_up != 0])
hist(pi_values$pi_down[pi_values$pi_down != 0])
#Using that reference isolate, calculate nucleotide diversity of the other 7 isolates against it. 
#Store this in a table with format:

#Locus    Isolate1    Isolate2      Isolate3  
#NEIS0001     0           23            23      #isolate 1 was the reference here, isolates 2 and 3 have the same number of diffs, so may have the same allele
#NEIS0002     0           0             0       #all same allele
#This table would be a way to get a quantification of pairwise nucleotide differences for each igr for each locus. A regression can be done on this table.


# I need to identify the most variable loci and see if they are significant in the linear regression
#Load significant loci from the regression analysis
sigtable <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/linearregression_igrRNA_pvalues.xlsx")
up_sig <- sigtable[which(sigtable$p_up < 0.05),]
down_sig <- sigtable[which(sigtable$p_down < 0.05),]

pi_values_upsig <- pi_values[which(pi_values$Locus %in% up_sig$Locus),]
pi_values_downsig <- pi_values[which(pi_values$Locus %in% down_sig$Locus),]

#average diversity of the hits
mean_pi_up_sig_up <- mean(pi_values_upsig$pi_up) #0.0444
mean_pi_down_sig_up <- mean(pi_values_upsig$pi_down) #0.0488
mean_pi_down_sig_down <- mean(pi_values_downsig$pi_down) #0.0527
mean_pi_up_sig_down <- mean(pi_values_downsig$pi_up) #0.0492

mean_sig_up <- (mean_pi_up_sig_up + mean_pi_down_sig_up)/2 #0.0466
mean_sig_down <- (mean_pi_up_sig_down + mean_pi_down_sig_down)/2 #0.509

#average diversity of the full dataset
mean_total_up <- mean(pi_values$pi_up, na.rm = TRUE) #0.0285
mean_total_down <- mean(pi_values$pi_down, na.rm = TRUE) #0.0198
mean_total <- (mean_total_up+mean_total_down)/2 #0.0242

avgs <- data.frame(samples = c('up_significant', 'down_significant', 'all_loci'), averages = c(mean_sig_up, mean_sig_down, mean_total))


ggplot(data=avgs, aes(x=samples, y= averages)) +
  geom_bar(stat='identity')



#prep data for a boxplot. Need to add a column defining the group
pi_values$Group="nonsignificant"
pi_values$Group[which(pi_values$Locus %in% up_sig$Locus)] = "up_significant"
pi_values_upsig <- pi_values
pi_values$Group = "nonsignificant"
pi_values$Group[which(pi_values$Locus %in% down_sig$Locus)] = "down_significant"
pi_values_downsig <- pi_values

#comparinh up pi in up_igrs between loci signficiant in the up_igrg in the linear regression against any that were nonsig
my_comparisons <- list(c("nonsignificant", "up_significant"))
ggplot(data=subset(pi_values_upsig, !is.na(pi)), aes(x=Group, y=pi_up, na.rm = TRUE)) +
  geom_boxplot(outlier.size= 0.2, outlier.color = 'red') + 
  stat_summary(fun.y=mean, geom="point", shape = 23, size = 4) +
  stat_compare_means(method = 't.test', comparisons = my_comparisons)

#The same but comparing pi in down_igrs sig in their down_igr in the linear regression against any loci that were nonsig
my_comparisons <- list(c("nonsignificant", "down_significant"))
ggplot(data=subset(pi_values_downsig, !is.na(pi)), aes(x=Group, y=pi_up, na.rm = TRUE)) +
  geom_boxplot(outlier.size= 0.2, outlier.color = 'red') + 
  stat_summary(fun.y=mean, geom="point", shape = 23, size = 4) +
  stat_compare_means(method = 't.test', comparisons = my_comparisons)

#looking at the pi_up difference betwen loci that were significant in either igr 
pi_values$Group="nonsignificant"
pi_values$Group[which(pi_values$Locus %in% up_sig$Locus)] = "significant"
pi_values$Group[which(pi_values$Locus %in% down_sig$Locus)] = "significant"
my_comparisons <- list(c("nonsignificant", "significant"))
ggplot(data=subset(pi_values, !is.na(pi)), aes(x=Group, y=pi_up, na.rm = TRUE)) +
  geom_boxplot(outlier.size= 0.2, outlier.color = 'red') + 
  stat_summary(fun.y=mean, geom="point", shape = 23, size = 4) +
  stat_compare_means(method = 't.test', comparisons = my_comparisons)

#looking at the pi_down difference between loci that were significant in either igr
ggplot(data=subset(pi_values, !is.na(pi)), aes(x=Group, y=pi_down, na.rm = TRUE)) +
  geom_boxplot(outlier.size= 0.2, outlier.color = 'red') + 
  stat_summary(fun.y=mean, geom="point", shape = 23, size = 4) +
  stat_compare_means(method = 't.test', comparisons = my_comparisons)



#average nucleotide diversity of whole genome"
wg_pi <- function(isolate){
  seq <- read.FASTA(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/", isolate, ".fasta", sep = ""))
  
  
  
  return(nuc.div(seq))
}


isolate <- "27509"
seq <- read.FASTA(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/", isolate, ".fasta", sep = ""))

wg_join <- function(sequence){
  isolate<- read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/", sequence, ".fasta",sep=""), as.string = TRUE)
  seq <- ""
  for(i in 1:length(isolate)){
    seq <- paste(seq, isolate[i], sep = '')
  }
  return(seq)
}

wg_join("27509")
isolate


i_27509 <- wg_join("27509")
i_27553 <- wg_join("27553")
i_28262 <- wg_join("28262")
i_28269 <- wg_join("28269")
i_28287 <- wg_join("28287")
i_53930 <- wg_join("53930")
i_53948 <- wg_join("53948")
i_53951 <- wg_join("53951")
samples <- list(i_27509, i_27553, i_28262, i_28269, i_28287, i_53930, i_53948, i_53951)
write.fasta(samples, names = c("27509", "27553", "28262", "28269", "28287", "53930", "53948", "53951"), "~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/each_wg_concatenated.fasta")

i_27509 <- read.FASTA("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/each_wg_concatenated.fasta")
#test <- nuc.div(i_27509)
#I need pi graphs for pi_up vs pi_up_sig and pi_down vs pi_down_sig for sig linear regression loci
#Also pi_up vs pi_up_sig_exp and pi_down vs pi_down_sig_exp for loci with significant expression
#and maybe the same for loci with log2fold change above 1. 

#Also need a pi value across the hwole genome to find some genome-wide mean for each isolate and then get a mean for the entire dataset.



#Per isolate pi values:



#To find loci with nucleotide diversity substantially above the average we'll use the boxplot outlier equation:
#(Q1- 1.5*IGQ) or (Q3 +1.5*IQR).   Note - values of 0 pi are excluded from the caluclation of quartiles

#IQR of the loci:
#up_igrs:
iqrange_up <- IQR(pi_values$pi_up[which(pi_values$pi_up != 0)], na.rm = T) #0.1245
q3_up <- quantile(pi_values$pi_up[which(pi_values$pi_up != 0)], probs = c(0.75), na.rm = T) #0.1255
#down_igrs:
iqrange_down <- IQR(pi_values$pi_down[which(pi_values$pi_down != 0)], na.rm = T) #0.0832858
q3_down <- quantile(pi_values$pi_down[which(pi_values$pi_down != 0)], probs = c(0.75), na.rm = T) #0.0847858

#label loci as high variation if both their up and downstream igrs have pi values high enough to be considered outliers in a boxplot. 
pi_values$Var <- "low"
for(i in 1:length(pi_values$Locus)){
  if(is.na(pi_values$pi_up[i]) | is.na(pi_values$pi_down[i])){
    next
  }
  if(pi_values$pi_up[i] > (q3_up + (1.5*iqrange_up)) && pi_values$pi_down[i] > (q3_down + (1.5*iqrange_down))){
    pi_values$Var[i] <- 'high'
  }
}

#or try it so that the definition of high variation is those where at least one of the igr pair would be an outlier, and the other at least has as much variation as the mean 
for(i in 1:length(pi_values$Locus)){
  if(is.na(pi_values$pi_up[i]) | is.na(pi_values$pi_down[i])){
    next
  }
  if(pi_values$pi_up[i] > (q3_up + (1.5*iqrange_up)) && pi_values$pi_down[i] > (median(pi_values$pi_down[which(pi_values$pi_down != 0)]))){
    pi_values$Var[i] <- 'high'
  }
  if(pi_values$pi_down[i] > (q3_down + (1.5*iqrange_down)) && pi_values$pi_up[i] > (median(pi_values$pi_up[which(pi_values$pi_up != 0)]))){
    pi_values$Var[i] <- 'high'
  }
}


#or if i just take highVar to be twice the mean:
mean_total_up #0.0285
mean_total_down #0.0198
mean_total #0.0242
pi_values$Var <- "low"
for(i in 1:length(pi_values$Locus)){
  if(is.na(pi_values$pi_up[i]) | is.na(pi_values$pi_down[i])){
    next
  }
  if(pi_values$pi_up[i] > 2*mean_total && pi_values$pi_down[i] > mean_total){
    pi_values$Var[i] <- 'high'
  }
  if(pi_values$pi_down[i] > 2*mean_total && pi_values$pi_up[i] > mean_total){
    pi_values$Var[i] <- 'high'
  }
}
table(pi_values$Var)
length(pi_values$Locus[which(pi_values$pi_down > 0 & pi_values$pi_up > 0)])

seqs <- list("c", "g", "a")
as.alignment(seqs)
x <- as.SeqFastadna(seqs)
seqs <- as.DNAbin(seqs)
seqs <-list(c('g', 'c', 'g', 'c','a','t','a','t','c','g','a','t','c'), 
     c('g', 'c', 'g', 'c','t','t','a','t','c','t','a','t','c'),
     c('g', 'c', 'g', 'c','a','t','a','t','c','g','a','t','c'),
     c('g', 'c', 'g', 'c','t','t','a','t','c','g','a','t','c'))
align
nuc.div(seqs)


example <- DNAStringSet(c("GCGCATATCGATC", 
                          "GCGCTTATCTATC", 
                          "GCGCATATCGATC", 
                          "GCGCTTATCGATC"))
example <- AlignSeqs(example)
example
nuc.div(as.DNAbin(example))

6 comparisons each of length 13
3 differences
3 diffs
1 diff

example <- DNAStringSet(c("GGGC", 
                          "GGGG", 
                          "GGGG", 
                          "GGGG"))
example <- AlignSeqs(example)
example
nuc.div(as.DNAbin(example))

3 diffs
6 of 4*2

