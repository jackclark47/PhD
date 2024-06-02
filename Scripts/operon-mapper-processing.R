#This script corrects the formatting of operon_mapper outputs, before annotating all hits with PUBMLST IDs

#Function to correct formatting and remove single gene operons
formatter <- function(file){
  removables <- c()
  #Add operon numbers to each row and identify rows containing only operon numbers
  for(i in 1:nrow(file)){
    if(!(is.na(file$Operon[i]))){
      current_operon <- file$Operon[i]
    }
    if(is.na(file$Operon[i])){
      file$Operon[i] <- current_operon
    }
    if(length(na.omit(unlist(file[i,]))) == 1){
      removables <- c(removables, i)
    }
  }
  
  #Remove rows containing only operon numbers
  out <- file[-removables,]
  
  #Remove operons containing only a single gene
  operon_counts <- count(out, Operon)
  removable_ops <- c()
  for(i in 1:nrow(operon_counts)){
    if(operon_counts$n[i] == 1){
      removable_ops <- c(removable_ops, operon_counts$Operon[i])
    }
  }
  
  out <- out[which(!(out$Operon %in% removable_ops)),]
  
  print(paste('Formatting complete', sep = ''))
  n_operons <- length(unique(out$Operon))
  n_genes <- nrow(out$Operon)
  print(paste(n_operons, 'operons identified'))
  print(paste(n_genes, 'genes identified'))
  print(paste(n_genes/n_operons, 'average genes per operon'))
  return(out)
}

#file <- read.xlsx('~/Downloads/N222.2_operons/list_of_operons_3477439.xlsx')
#test <- formatter(file)

#Function to annotate an operon-mapper output file with pubmlst ids instead of randomly assigned prokka hits.
gene_anno <- function(file, isolate){
  
  #Create output dataframe
  out <- file
  out$pubmlst_id <- NA
  
  genome_seq <- read.fasta(paste('~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/', isolate, '.fasta', sep = ''), as.string = TRUE)
  blast_db <- blast(db="~/Documents/PhD/PhD/N_meningitidis_loci/translations/BLAST_db/coding.fsa", type = 'blastp')
  for(i in 1:nrow(file)){
    
    #Identify the contig containing the operon
    contig_number <- contig_finder(file, file$IdGene[i], isolate)
    gene_start <- file$PosLeft[i]
    gene_end <- file$postRight[i]
    gene_sequence <- DNAStringSet(substr(genome_seq[contig_number], gene_start, gene_end))
    if(isolate %in% c('53930', '53948', '53951')){
      gene_sequence <- DNAStringSet(substr(genome_seq[which(names(genome_seq) %in% contig_number)], gene_start, gene_end))
    }
    prot_sequence1 <- Biostrings::translate(gene_sequence)
    prot_sequence2 <- Biostrings::translate(reverseComplement(gene_sequence))
    
    res1 <- predict(blast_db, prot_sequence1, type = 'blastp')[1,]
    res2 <- predict(blast_db, prot_sequence2, type = 'blastp')[1,]
    res <- rbind(res1, res2)
    res <- res[order(res$evalue),]
    res <- res[1,]
    
    #If the percent identity and evalue are low in the top result, discard
    if(is.na(res$pident) | is.na(res$evalue)){
      next
    }
    #Otherwise, assign the sequence the same name as the top blast hit
    if(res$pident > 80 && res$evalue < 0.005){
      out$pubmlst_id[i] <- res$sseqid
    }
    #print(paste('Gene number ', i, 'assigned id: ', out$pubmlst_id[i]))
  }
  
  #Write output to an xlsx with the name of the isolate
  wb <- loadWorkbook('~/Documents/PhD/PhD/operon_mapper_res/operon_mapper_annos.xlsx')
  addWorksheet(wb, isolate)
  writeData(wb, sheet = isolate, out)
  saveWorkbook(wb, file = '~/Documents/PhD/PhD/operon_mapper_res/operon_mapper_annos.xlsx', overwrite = TRUE)
  print(paste('Writing operons for isolate ', isolate, 'to operon_mapper_annos.xlsx', sep = ''))
  print('iteration finished')
  print('===============')
  print('===============')
  return()
}

#gene_anno(formatted_file, '53930')
#Function to identify which contig of an open assembly contains the operon for the current iteration
contig_finder <- function(file, gene_ID, isolate){
  gene_ID_number <- as.numeric(substr(gene_ID, start = nchar(gene_ID)-4, stop = nchar(gene_ID)))
  gene_coords <- read.delim(paste('~/Documents/PhD/PhD/operon_mapper_res/', isolate, '_operons/ORFs_coordinates_', isolate, sep = ''), header = FALSE)
  
  #Obtain a list of contigs with the gene numbers contained in each
  for(i in 1:nrow(gene_coords)){
    number <- gene_coords$V9[i]
    number <- as.numeric(substr(number, start= 13, stop = 17))
    if(number == gene_ID_number){

      contig_number <- gene_coords$V1[i]
      if(isolate %in% c('53930', '53948', '53951')){
        contigs <- unique(gene_coords$V1)
        contig_name <- contigs[which(contigs %in% contig_number)]
        return(contig_name)
      }
      contig_number <- as.numeric(substr(contig_number, start = nchar(contig_number)-4, stop =nchar(contig_number)))
      return(contig_number)
    }
  }
}

#contig_finder(formatted_file, formatted_file$IdGene[i], '53948')

#test <- read.delim('~/Documents/PhD/PhD/operon_mapper_res/27509_operons/ORFs_coordinates_1658703', header = FALSE)

#isolate <- '53930'

main <- function(isolate_list){
  write.xlsx(x = NA, file = '~/Documents/PhD/PhD/operon_mapper_res/operon_mapper_annos.xlsx')
  for(isolate in isolate_list){
    print(paste('Beginning iteration for isolate', isolate))
    file <- read.xlsx(paste('~/Documents/PhD/PhD/operon_mapper_res/', isolate, '_operons/list_of_operons_', isolate, '.xlsx', sep = ''))
    formatted_file <- formatter(file)
    gene_anno(formatted_file, isolate)
  }
}

isolate_list <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
main(isolate_list)