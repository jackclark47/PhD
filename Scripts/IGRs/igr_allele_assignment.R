#Script to align igrs for each locus, calculate their nucleotide diversity, and call allele numbers. Outputs a spreadsheet of allele numbers and classes igrs by up,down, low, high var

library(openxlsx)
library(magrittr)
library(stringr)
library(seqinr)
library(Biostrings)
library(msa)
library(pegas)
#Read in each isolates igr datafile and get a list of unique loci across all of them
isolates <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
get_loci <- function(isolates){
  
  loci <- c()
  for(isolate in isolates){
    ids <- read.xlsx(paste('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/isolate_igr_initial/igrs_', isolate, '.xlsx', sep = ''))[,1]
    loci <- c(loci, ids)
  }
  loci <- unique(loci)
  return(loci)
}

loci <- get_loci()

#extract the igr sequences for each locus for each isolate and save as fasta files. 
write_igrs <- function(isolates){
  for(isolate in isolates){
    igr_data <- read.xlsx(paste('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/isolate_igr_initial/igrs_', isolate, '.xlsx', sep = ''))
    for(j in 1:nrow(igr_data)){
      igr_down <- igr_data$down_sequence[j]
      igr_up <- igr_data$up_sequence[j]
      id <- igr_data$pubmlst_id[j]
      
      up_name <-paste(id, 'up', sep='_')
      down_name <- paste(id, 'down', sep='_')
      
      print(paste('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igr_seqs/', id, '_up_', isolate, '.fasta', sep = ''))
      
      write.fasta(igr_up, up_name, paste('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igr_seqs/individual/', id, '_up_', isolate, '.fasta', sep = ''))
      write.fasta(igr_down, down_name, paste('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igr_seqs/individual', id, '_down_', isolate, '.fasta', sep = ''))
    }
  }
}

write_igrs(isolates)

#Combine individual fasta files into megafastas
dir <- "~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igr_seqs/individual/"
combine_igr_seqs <- function(dir){
  
  setwd(dir)
  files <- list.files()
  loci <- unlist(str_extract_all(files, 'NEIS[:alnum:]{4,}(?=_)'))
  
  for(locus in unique(loci)){
    print(locus)
    
    for(i in c('_up', '_down')){
      query <- paste(locus, i, '_[:digit:]{5}.fasta', sep = '')
      
      file_indices <- str_which(files, query)
      filenames <- files[file_indices]
      
      fastas <- sapply(filenames, read.fasta, as.string = TRUE)
      names(fastas) %<>% sapply(str_extract, pattern = '[:graph:]+(?=\\.fasta)')
      
      write.fasta(fastas, names = names(fastas), file.out = paste('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igr_seqs/multifastas/', locus, i, '.fasta', sep = ''))
    }
  }
  
}

combine_igr_seqs(dir)

###align igrs, call allele numbers, calculate nucleotide diversity
#align igrs - read multifasta of given locus, align the up or the down stream. return the alignment
# align_igrs <- function(file){
#   sequences <- read.FASTA(file)
#   sequences <- read.fasta(file, as.string = T)
#   nuc.div(sequences, variance = F, pairwise.deletion = FALSE)
#   DNAbin
#   alignment <- msa(sequences)
#   return(alignment)
# }

# b <- hclust(stringDist(a))
# c <- stringDist(a)
# 
# plot(b)

#assign allele numbers based on the alignment. return an updated dataframe with the information
assign_alleles <- function(df, file, locus, strand, count){
  
  count = ceiling(count / 2)
  
  df$Locus[count] <- locus
  
  sequences <- read.fasta(file, as.string = TRUE)
  queries <- unique(unlist(sequences))
  
  print(locus)
  if(strand == 'down'){
    
    seq_length <- sort(table(sapply(sequences, nchar)), decreasing = TRUE)[1]
    seq_length %<>% names() %>% as.numeric()
    
    df$Length_down[count] <- seq_length
    
    for(i in 1:length(sequences)){
      sequence <- sequences[i]
      name <- names(sequences)[i]
      isolate <- str_extract(name, "(?<=_)[:digit:]{5}")
      columnname <- paste(isolate, 'down', sep = '_')
      print(columnname)
      df[count, columnname] <- which(queries == sequence)
    }
  }
  
  if(strand == 'up'){
    seq_length <- sort(table(sapply(sequences, nchar)), decreasing = TRUE)[1]
    seq_length %<>% names() %>% as.numeric()
    
    df$Length_up[count] <- seq_length
    
    for(i in 1:length(sequences)){
      sequence <- sequences[i]
      name <- names(sequences)[i]
      isolate <- str_extract(name, "(?<=_)[:digit:]{5}")
      columnname <- paste(isolate, 'up', sep = '_')
      print(columnname)
      df[count, columnname] <- which(queries == sequence)
    }
  }
  
  return(df)
}

sequences <- unlist(sequences)
length(unique(unlist(sequences)))
unique(unlist(sequences))[1]
sequences[1] == unique(unlist(sequences))[1]

allele_summaries <- function(df){
  for(i in 1:nrow(df)){
    df$`Alleles-Down`[i] <- max(as.numeric(na.omit(as.numeric(df[i, 4:11]))))
    df$`Alleles-Up`[i] <- max(as.numeric(na.omit(as.numeric(df[i,12:19]))))
    
    df$Combined[i] <- df$`Alleles-Down`[i] + df$`Alleles-Up`[i]
  }
  return(df)
}

df <- allele_summaries(df_alleles)

group_by_var <- function(df){
  Var <- c()
  NoVar <- c()
  upVar <- c()
  downVar <- c()
  lowVar <- c()
  highVar <- c()
  
  #full vargroups
  for(i in 1:nrow(df)){
    locus <- df$Locus[i]
    print(locus)
    if(df$Combined[i] > 2){
      Var <- c(Var, locus)
    } else if(df$Combined[i] == 2){
      NoVar <- c(NoVar, locus)
    }
    
    if(df$`Alleles-Down`[i] == 1 && df$`Alleles-Up`[i] > 1){
      upVar <- c(upVar, locus)
    } else if(df$`Alleles-Down`[i] > 1 && df$`Alleles-Up`[i] == 1){
      downVar <- c(downVar, locus)
    } else if(df$Combined[i] >= 10){
      highVar <- c(highVar, locus)
    } else if(df$Combined[i] < 10 && df$Combined[i] > 3){
      lowVar <- c(lowVar, locus)
    }
  }
  
  #subset df
  upVar_igrs <- df[which(df$Locus %in% upVar),]
  downVar_igrs <- df[which(df$Locus %in% downVar),]
  Var_igrs <- df[which(df$Locus %in% Var),]
  NoVar_igrs <- df[which(df$Locus %in% NoVar),]
  lowVar_igrs <- df[which(df$Locus %in% lowVar),]
  highVar_igrs <- df[which(df$Locus %in% highVar),]
  #write spreadsheets
  write.xlsx(upVar_igrs, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/upvar_igrs.xlsx')
  write.xlsx(downVar_igrs, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/downvar_igrs.xlsx')
  write.xlsx(Var_igrs, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/var_igrs.xlsx')
  write.xlsx(NoVar_igrs, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/novar_igrs.xlsx')
  write.xlsx(lowVar_igrs, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/lowvar_igrs.xlsx')
  write.xlsx(highVar_igrs, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/highvar_igrs.xlsx')
}
group_by_var(df)


df$Combined[i] <- df$`Alleles-Down`[i] + df$`Alleles-Up`[i]


#which alleles are significantlt expressed?




#plot vargroups by exp
varplot <- function(){
  
}


#calculate nucleotide diversity based on the alingment. return a different, updated dataframe with the information
nucleotide_div <- function(df, alignment, locus, strand, count){
  
}

#Bring it all together. For each file, loop through the functions.
main <- function(dir){
  setwd(dir)
  files <- list.files()
  
  df_alleles <- as.data.frame(matrix(data = NA, nrow = length(files), ncol = 22))
  #df_nuc_div <- as.data.frame(matrix(data = NA, nrow = length(files), ncol = 5))
  
  colnames(df_alleles) <- c('Locus', 'Length_down', 'Length_up', '27509_down', '27553_down', 
                            '28262_down', '28269_down', '28287_down', '53930_down', '53948_down', 
                            '53951_down', '27509_up', '27553_up', '28262_up', '28269_up', '28287_up', 
                            '53930_up', '53948_up', '53951_up', 'Alleles-Down', 'Alleles-Up', 'Combined')
  #colnames(df_nuc_div) <- c('Locus', 'pi_up', 'variance_up', 'pi_down', 'variance_down')
  count = 1
  for(file in files){
    locus <- str_extract(string = file, pattern ="[:alnum:]+(?=_)")
    strand <- str_extract(string = file, pattern = "[:alnum:]+(?=\\.)")
      
    #alignment <- align_igrs(file)
    
    df_alleles <- assign_alleles(df_alleles, file, locus, strand, count)
    df_alleles <- allele_summaries(df_alleles)
    #df_nuc_div <- nucleotide_div(df_nuc_div, alignment, locus, strand, count)
    count = count + 1
  }
  
  write.xlsx(x = df_alleles, file = '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs.xlsx')
  write.xlsx(x = df_nuc_div, file = '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct')
}


dir <- '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igr_seqs/multifastas'




