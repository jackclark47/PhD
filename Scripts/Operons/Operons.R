#This script uses the N222.2 hybrid reference to identify operons based on their coordinates
#First, a variety of igr thresholds are tested
#30 nt is selected as a cutoff IGR length and operons extracted
#Then genes in each operon have their PubMLST IDs identified by running BLAST searches against my local copy of 
#the PubMLST database, which contains a single allele of every N.meningitidis locus in the database
#Once annotations are made, operons are compared to the expression and igr datasets
#Each operon is also run as a query in a BLAST search against each isolate to check its presence/absence and integrity. 500 bases up and downstream of the operon 
#in each isolate are taken and the entire sequence aligned for all 8 isolates to identify variation.

#to do:
#add direction of each operon, forwards or reverse
#add function of each operon. 


#libraries
library(gggenomes)
library(openxlsx)
library(rBLAST)
library(seqinr)

#Function that runs a BLAST search of every gene against a blast database to annotate with pubmlst ids
annotate_genes <- function(reference){
  reference$pubmlst_ID <- NA
  blast_db <- blast(db="~/Documents/PhD/PhD/N_meningitidis_loci/translations/BLAST_db/coding.fsa", type = 'blastp')
  genome_seq <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222.2.fasta', as.string = TRUE)
  for(i in 1:nrow(reference)){
    print(i)
    gene_start <- reference$start[i]
    gene_end <- reference$end[i]
    
    gene_sequence <- DNAStringSet(substr(genome_seq, gene_start, gene_end))
    prot_sequence1 <- Biostrings::translate(gene_sequence)
    prot_sequence2 <- Biostrings::translate(reverseComplement(gene_sequence))
    
    res1 <- predict(blast_db, prot_sequence1, type = 'blastp')[1,]
    res2 <- predict(blast_db, prot_sequence2, type = 'blastp')[1,]
    res <- rbind(res1, res2)
    res <- res[order(res$evalue),]
    res <- res[1,]
    #Incorporate the reverse complement processing and choose the most significant hit of the two. 
    
    #If the percent identity and evalue are low in the top result, discard
    if(is.na(res$pident) | is.na(res$evalue)){
      next
    }
    #Otherwise, assign the sequence the same name as the top blast hit
    if(res$pident > 80 && res$evalue < 0.005){
      reference$pubmlst_ID[i] <- res$sseqid
    }

  }
  print('annotation complete')
  return(reference)
}


#Function that checks if any operon number has fewer entries than the operon threshold
operon_filter <- function(data, operon_threshold){
  removables <- c()
  operon_nums <- data$operon
  for(i in 1:length(unique(operon_nums))){
    if(sum(operon_nums == unique(operon_nums)[i]) < operon_threshold) {
      removables <- c(removables, unique(operon_nums)[i])
    }
  }
  data <- data[which(!(data$operon %in% removables)),]
  return(data)
}

#Function that identifies all operons in the dataset
operons_predict <- function(data, igr_threshold, operon_threshold, sheetnum, file){
  operon_num <- 1
  data$operon[1] <- operon_num
  for(feature in 2:length(reference$type)){
    
    igr_length = abs(reference$start[feature] - reference$end[feature-1])
    strand_current = reference$strand[feature]
    strand_previous = reference$strand[feature-1]
    #Accounting for cases where the start of end of the current CDS is contained within the coords of the previous CDS, 
    if(strand_current == strand_previous){
      if( (reference$start[feature-1] < reference$start[feature] && reference$start[feature] < reference$end[feature-1]) | (reference$start[feature-1] < reference$end[feature] && reference$end[feature] < reference$end[feature-1]) ){
        data$operon[feature] <- operon_num
        next
      }
      
      if(igr_length <= igr_threshold){
        data$operon[feature] <- operon_num
        next
      }
      if(igr_length > igr_threshold){
        operon_num = operon_num + 1
        data$operon[feature] <- operon_num
        next
      }
    }
    else if(strand_current != strand_previous){
      operon_num = operon_num + 1
      data$operon[feature] <- operon_num
    }

  }
  #Remove operons with fewer genes than the operon threshold
  data <- operon_filter(data, operon_threshold)

  #write this to an excel sheet
  wb <- loadWorkbook(file)
  addWorksheet(wb, sheetnum)
  writeData(wb, sheet = sheetnum, data)
  saveWorkbook(wb, file = file, overwrite = TRUE)
  print(paste('run ', sheetnum, ' complete', sep = ''))
  return()
}

#First load in the genbank file
reference <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222.2.gbk")[-1,]

#Attempt to annotate every gene in the reference
reference <- annotate_genes(reference)
reference <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/MC58.gbff")[-1,]
reference <- reference[which(reference$type == 'CDS'),]

main <- function(reference, igr_thresholds, operon_threshold, file){
  
  write.xlsx(x = NA, file = file)
  sheetnums <- seq(1, length(igr_thresholds))
  for(i in 1:length(igr_thresholds)){
    data <- as.data.frame(matrix(data = NA, nrow = nrow(reference), ncol = 5))
    colnames(data) <- c('gene_id', 'operon', 'start', 'stop', 'strand')
    data$gene_id <- reference$feat_id
    data$start <- reference$start 
    data$stop <- reference$end
    data$strand <- reference$strand
    #data$pubmlst_ID <- reference$pubmlst_ID

    operons_predict(data, igr_thresholds[i], operon_threshold, sheetnum = sheetnums[i], file = file)
  }
  return()
}

#this code will miss genes whose start is contained within the end of the previous gene!


res <- main(reference, c(10, 20, 30, 40, 50), operon_threshold = 2, file = '~/Documents/PhD/PhD/RNA_IGR/MC58_operon_test.xlsx')

#tabulate results and save in an xlsx

