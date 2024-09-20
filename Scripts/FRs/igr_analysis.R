#Script to combine the pubmlst and blast igrs together into a single file.
#Check theres no overlap between those two datasets
#then assign alleles to them
#Then conduct nucleotide diversity analysis
#Then split into groups based on nucleotide diversity values
#Then conduct some basic analysis with that data. 

#####Contents (CMD-F to go to each section)
###1 - COMBINING IGRS
###2 - ALLELE ASSIGNMENT
###3 - NUCLEOTIDE DIVERSITY
###4 - CORRELATING EXPRESSION
###5 - LINEAR REGRESSION

library(seqinr)
library(pegas)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(stringr)

###1 - COMBINING IGRS

#Load in the blast igrs
igrs_up_BLAST <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/up_igr_BLAST.xmfa', as.string = TRUE)
igrs_down_BLAST <-read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/down_igr_BLAST.xmfa', as.string = TRUE)
#load in the pubmlst igrs
igrs_up_MLST <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/Old/PubMLST_upstream_igr.xmfa', as.string = TRUE, whole.header = TRUE)
igrs_down_MLST <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/Old/PubMLST_downstream_igr.xmfa', as.string = TRUE, whole.header= TRUE)

#Change names of igrs_up_MLST and igrs_down_MLST to contain just the isolate and locus name
for(i in 1:length(igrs_up_MLST)){
  locus <- tail(str_split(names(igrs_up_MLST[i]), ' ')[[1]], n=1)
  isolate <- substr(names(igrs_up_MLST[i]), 1, 5) 
  
  names(igrs_up_MLST)[i] <- paste(locus, '_', isolate, sep ='')
  names(igrs_down_MLST)[i] <- paste(locus, '_', isolate, sep ='')
  
}

#remove .fasta from the names of igrs_up_BLAST and igrs_down_BLAST
for(i in 1:length(igrs_up_BLAST)){
  names(igrs_up_BLAST)[i] <- substr(names(igrs_up_BLAST[i]), 1, (nchar(names(igrs_up_BLAST[i]))-6 ) )
  names(igrs_down_BLAST)[i] <- substr(names(igrs_down_BLAST[i]), 1, (nchar(names(igrs_down_BLAST[i]))-6 ) )
}

#Check for duplicate loci
for(i in 1:length(igrs_down_BLAST)){
  if(names(igrs_down_BLAST[i]) %in% names(igrs_down_MLST)){
    print(names(igrs_down_BLAST[i]))
  }
} #no duplicates found

#Combine the two datasets and write to files
igrs_up <- c(igrs_up_BLAST, igrs_up_MLST)
igrs_down <- c(igrs_down_BLAST, igrs_down_MLST)

#write.fasta(igrs_up, names(igrs_up), '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_up.xmfa')
#write.fasta(igrs_down, names(igrs_down), '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_down.xmfa')

#Convert igrs with nchar values < 100 to NA
for(i in 1:length(igrs_down)){
  if(nchar(igrs_down[i]) < 100){
    igrs_down[i] <- NA
  }
  if(nchar(igrs_up[i]) < 100){
    igrs_up[i] <- NA
  }
}





#2 - ALLELE ASSIGNMENT
###First create the dataframe to hold the data
out <- matrix(nrow = length(igrs_up)/8, ncol = 21)
colnames(out) <- c("Locus", "Length_down", "Length_up", "27509_down", "27553_down", "28262_down", "28269_down", "28287_down",
                   "53930_down", "53948_down", '53951_down', "27509_up", "27553_up", "28262_up", "28269_up",
                   "28287_up", "53930_up", '53948_up', "53951_up", "Alleles-Down", "Alleles-Up")

out <- as.data.frame(out)

locus_names <- substr(names(igrs_up)[seq(1,length(igrs_up),8)], 1, (nchar(names(igrs_up)[seq(1,length(igrs_up),8)])-6))
row = 0
for(i in seq(from = 1, to = length(igrs_up), by = 8)){

  row = row + 1
  #add locus info to the dataframe
  out$Locus[row] = locus_names[row] 
  
  #make vector of the 8 sequences at that locus
  all_sequences <- c(igrs_up[[i]], igrs_up[[i+1]], igrs_up[[i+2]], igrs_up[[i+3]],
                     igrs_up[[i+4]], igrs_up[[i+5]], igrs_up[[i+6]], igrs_up[[i+7]])
  
  #Check how many unique sequences there are
  sequencesold <- unique(all_sequences)
  #remove nas from sequences
  sequences <- c()
  for(j in 1:length(sequencesold)){
    if(sequencesold[j] == 'na' | is.na(sequencesold[j])){
      next
    }
    if(sequencesold[j] != 'na' && !is.na(sequencesold[j])){
      sequences <- c(sequences, sequencesold[j])
    }
  }
  if(length(sequences) == 0){
    out[row, (12:19)] <- NA
  }
  out$Length_up[row] <- nchar(sequences)[1]
  #Add the number of unique sequences to the Alleles-up column
  out$`Alleles-Up`[row] <- length(sequences)
  #Track the unique sequences for each locus so that they can be mapped back to the individual isolates and numbered individually
  for(j in 1:length(all_sequences)){
    for(k in 1:length(sequences)){ 
      if(all_sequences[j] == 'na' | is.na(all_sequences[j])){
        out[row,(j+11)] <- NA
        next
      }
      if(all_sequences[j] == sequences[k]){
        out[row,(j+11)] <- k
      }

    }
  }
  
  
  
  #repeat for downstream
  #make vector of the 8 sequences at that locus
  all_sequences <- c(igrs_down[[i]], igrs_down[[i+1]], igrs_down[[i+2]], igrs_down[[i+3]],
                     igrs_down[[i+4]], igrs_down[[i+5]], igrs_down[[i+6]], igrs_down[[i+7]])
  #Check how many unique sequences there are
  sequencesold <- unique(all_sequences)
  #remove nas from sequences
  sequences <- c()
  for(j in 1:length(sequencesold)){
    if(sequencesold[j] == 'na' | is.na(sequencesold[j])){
      next
    }
    if(sequencesold[j] != 'na' && !is.na(sequencesold[j])){
      sequences <- c(sequences, sequencesold[j])
    }
  }
  if(length(sequences) == 0){
    out[row, (4:11)] <- NA
  }
  #Add the number of unique sequences to the Alleles-down column
  out$`Alleles-Down`[row] <- length(sequences)
  out$Length_down[row] <- nchar(sequences)[1]
  #Track the unique sequences for each locus so that they can be mapped back to the individual isolates and numbered individually
  for(j in 1:length(all_sequences)){
    for(k in 1:length(sequences)){
      if(all_sequences[j] == 'na' | is.na(all_sequences[j])){
        out[row,(j+3)] <- NA
        next
      }
      if(all_sequences[j] == sequences[k]){
        out[row,(j+3)] <- k
      }
    }
  }
}
allele_data <- out
#save spreadsheet
write.xlsx(out, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_29112023.xlsx', keepNA = TRUE)


#load old spreadsheet and compare allele counts and variation
old_igr <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/NEIS_IGR_Alleles-cdb-21092022.xlsx", sheet = 4)
table(out$`Alleles-Down`, useNA = 'always')
table(old_igr$`Alleles-Down`, useNA = 'always')
table(out$`Alleles-Up`, useNA = 'always')
table(old_igr$`Alleles-Up`, useNA = 'always')


###3 - NUCLEOTIDE DIVERSITY
#Create empty matrix to store pi values
pi_values <- as.data.frame(matrix(data = NA, nrow = (length(up_igrs)/8), ncol = 5))
colnames(pi_values) <- c('Locus', 'pi_up', 'variance_up', 'pi_down', 'variance_down')

up_igrs <- readDNAStringSet("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_up.xmfa")
down_igrs <- readDNAStringSet("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_down.xmfa")

count = 0
for(i in seq(1,length(up_igrs), 8)){
  print(i)
  count = count+1
  #store the locus' name as an object
  locusname <- names(up_igrs[i])
  locusname <- substr(locusname, start = 1, stop = nchar(locusname)-6)
  
  #extract the igrs from each isolate for the current locus
  testlocus_up <- up_igrs[i:(i+7)]
  testlocus_down <- down_igrs[i:(i+7)]
  
  ######################Remove entries that contain NAs
  up_removables <- c()
  down_removables <- c()
  for(j in 1:length(testlocus_up)){
    if(is.na(testlocus_up[j]) | nchar(testlocus_up[j]) < 100){
      up_removables <- c(up_removables, j)
    }
    if(is.na(testlocus_down[j]) | nchar(testlocus_down[j]) < 100){
      down_removables <- c(down_removables, j)
    }
  }
  
  if(!(is.null(up_removables))){
    testlocus_up <- testlocus_up[-c(up_removables)]
  }
  if(!(is.null(down_removables))){
    testlocus_down <- testlocus_down[-c(down_removables)]
  }
  
  
  #if the upstream igrs are identical in length, calculate the nucleotide diversity for them
  
  if(length(unique(sapply(testlocus_up, length))) == 1 && length(testlocus_up) > 1){
    testlocus_up <- AlignSeqs(testlocus_up)
    out <- my.nuc.div(testlocus_up)
    pi_values$Locus[count] <- locusname
    pi_values$pi_up[count] <- out
    #pi_values$variance_up[count] <- out[2]
  }
  
  #if the up_igrs are not all the same length, record only the locus name
  if(length(unique(sapply(testlocus_up, length))) > 1){
    print(locusname)
    pi_values$Locus[count] <- locusname
  }
  
  #if the down_igrs are all the same length, calculate the nucletodie diversity 
  if(length(unique(sapply(testlocus_down, length))) == 1 && length(testlocus_down) > 1){
    testlocus_down <- AlignSeqs(testlocus_down)
    out <- my.nuc.div(testlocus_down)
    pi_values$pi_down[count] <- out
    #pi_values$variance_down[count] <- out[2]
  }
  #if not then record only the locus name
  if(length(unique(sapply(testlocus_down, length))) > 1){
    pi_values$Locus[count] <- locusname
    print(locusname)
  }
  #Note, nuc.div can only be calculated for sequences of the same length, hence the above
}
#save this table to a file
# write.xlsx(pi_values, "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/pi_values.xlsx")

#pi_values = read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/pi_valuesalign.xlsx")
#see distribution of pi values
hist(pi_values$pi_down)
hist(pi_values$pi_up[pi_values$pi_up != 0])
hist(pi_values$pi_down[pi_values$pi_down != 0])


#Grouping of loci using pi values
#Definitions of the 4 groups to use:
#>2x mean pi in both igrs
#>2x mean pi in up igr, < 2* mean in down
#>2x mean pi in down igr, <2* mean in up
#<2x mean pi in both igrs

#Identify mean pi values
up_pi_mean <- mean(pi_values$pi_up, na.rm = T) #0.0189
down_pi_mean <- mean(pi_values$pi_down, na.rm = T) #0.0.0211
total_pi_mean <- mean(c(up_pi_mean, down_pi_mean)) #0.0200

HighVar <- c()
UpVar <- c()
DownVar <- c()
LowVar <- c()
NoVar <- c()

###ERROR HERE need to index allele data after the loop by name, not by index, as loci in allele data and pi_values 
#are in different indices. 


# up_igrs["NEIS0025_28269"]
# 
# nuc.div(muscle5(as.DNAbin(up_igrs[14657:14661]), exec = "/Users/u5501917/opt/anaconda3/pkgs/muscle-5.1-h7d747a2_3/bin/muscle"))
# nuc.div(as.DNAbin(AlignSeqs()))
# NucleotideDiversity(DNA.dist(up_igrs[14657:14664]))
# up_igrs[91:96]
# up_igrs[14:14664]
# nuc.div()
# s1 <- up_igrs[92]
# s2 <- up_igrs[93]
# 
# count = 0
# for(char in 1:nchar(s1)){
#   if(substr(s1, char,char) != substr(s2, char,char)){
#     count = count + 1
#   }
# }
# 
# pi <- count/nchar(s1)/2
# pi
# nuc.div(as.DNAbin(AlignSeqs(up_igrs[92:96])), pairwise.deletion = T)
# nuc.div(as.DNAbin(list(c('g','c','t', 'g'), c('g', 'c', 't', 't'))))
# NucleotideDiversity(DNA.dist(up_igrs[92:96]))
# 
# write.fasta(up_igrs[92:96], file.out = 'nucdivtest.fasta', names = c('1', '2', '3', '4', '5'))
# up_igrs[95]

# AlignSeqs(up_igrs[92:96])[5] == AlignSeqs(up_igrs[92:96])[2]
# my.nuc.div(AlignSeqs(up_igrs[92:96]))
# DNA.dist(up_igrs[92:96])
# dist.dna(as.DNAbin(up_igrs[92:96]))
#now split loci into groups and write each group to its own spreadsheet
for(i in 1:length(pi_values$Locus)){
  
  #Dealing with loci with NA pi values in either igr. NA pi values are 0 
  if(is.na(pi_values$pi_up[i])){
    pi_values$pi_up[i] = 0
  }
  
  if(is.na(pi_values$pi_down[i])){
    pi_values$pi_down[i] = 0
  }
  
  #NoVar group
  if(pi_values$pi_up[i] == 0 && pi_values$pi_down[i] == 0){
    NoVar <- c(NoVar, pi_values$Locus[i])
    next
  }
  
  #HighVar group
  if(pi_values$pi_up[i] >= 2* total_pi_mean && pi_values$pi_down[i] >= 2*total_pi_mean){
    HighVar <- c(HighVar, pi_values$Locus[i])
    next
  }
    
  #UpVar group
  if(pi_values$pi_up[i] > 0 && pi_values$pi_down[i] == 0){
    UpVar <- c(UpVar, pi_values$Locus[i])
    next
  }
    
  #DownVar group
  if(pi_values$pi_down[i] > 0 && pi_values$pi_up[i] == 0){
    DownVar <- c(DownVar, pi_values$Locus[i])
    next
  }
    
  #LowVar group
  if(pi_values$pi_up[i] < 2* total_pi_mean | pi_values$pi_down[i] < 2*total_pi_mean){
    LowVar <- c(LowVar, pi_values$Locus[i])
    next
  }
}

length(LowVar) + length(HighVar) + length(UpVar) + length(DownVar)
length(Var)
length(NoVar)
allele_data <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx')
HighVardata <- allele_data[which(allele_data$Locus %in% HighVar),]
UpVardata <- allele_data[which(allele_data$Locus %in% UpVar),]
DownVardata <- allele_data[which(allele_data$Locus %in% DownVar),]
LowVardata <- allele_data[which(allele_data$Locus %in% LowVar),]
NoVardata <- allele_data[which(allele_data$Locus %in% NoVar),]
Vardata <- allele_data[which(!allele_data$Locus %in% NoVar),]

#Write each to a file
write.xlsx(HighVardata, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/HighVar.xlsx')
write.xlsx(UpVardata, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/UpVar.xlsx')
write.xlsx(DownVardata, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/DownVar.xlsx')
write.xlsx(LowVardata, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/LowVar.xlsx')
write.xlsx(NoVardata, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/NoVar.xlsx')
write.xlsx(Vardata, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/Var.xlsx')

#load rna data
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx', sheet = 6)[, c(8,9,15:22, 25:32)]

#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

#Order of isolates in both datasets:
#RNA: 20578 22698 N188.1  N222.1, N222.2, N445.1, N459.3, N459.6
#IGR: 20578 22698 N222.1  N188.1  N445.1, N222.2  N459.3, N459.6

#So rearrange columns in the rna dataset so theyre the same order
rna <- rna[,c(1:4, 6,5, 8,7, 9:12, 14,13,16,15, 17, 18)]

###graphing expression of each group 
#Refer to script RNA_IGR_correlations.R

#Look at log2fold change between all NonVar and Var
#Refer to script RNA_IGR_correlations.R

#Look at number of loci with log2fold change > 1 which are Var and how many are NonVar
#Refer to script RNA_IGR_correlations.R

#Look at log2fold change between all 5 groups
#Look


#graphing allele counts of each group


#graphing expression differences between variable and nonvariable loci



###4 - CORRELATING EXPRESSION
#Correlating expression by isolate
#Refer to RNA_correlations.R script



###5 - LINEAR REGRESSION
#linear regression of loci using rna and igr data
#Refer to script Linregression_igrvsrna.R


#graphing pi 
pi_values <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/pi_values.xlsx")
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

