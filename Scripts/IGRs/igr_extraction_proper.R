#Script to read isolate assemblies, convert the gbk gene names to NEIS codes via BLAST search against a shortened pubmlst database copy
#then save as new gbks, and extract igrs based on that information

library(gggenomes)
library(seqinr)
library(rBLAST)
library(Biostrings)
library(stringr)
library(openxlsx)

#Load in an annotated genome and its fasta file
test_gb <- read_gbk('~/Documents/PhD/RNA_IGR/Isolate_fastas/Annotations/53930/53930.gbff')
test_fa <- read.fasta('~/Documents/PhD/RNA_IGR/Isolate_fastas/Annotations/53930/53930.fna', as.string = T, seqonly = F)

#Loop through each CDS in the genbank file and BLAST it against the shortened pubmlst database
get_seq <- function(feature, fa){

  direction = feature$strand
  contig <- feature$seq_id
  
  start = feature$start
  stop = feature$end
  
  
  contig_seq <- fa[contig]
  CDS_seq <- DNAStringSet(substr(contig_seq, start = start, stop = stop))
  
  if(direction == '-'){
    CDS_seq <- reverseComplement(CDS_seq)
  }
  
  #get protein sequence for the blast search (blast db is a protein sequence db)
  aa_seq <- translate(CDS_seq)
  
  return(aa_seq)
}


testseq <- get_seq(test_gb[104,],test_fa)
testseq
blast_seq(testseq, db)

blast_seq <- function(aa_seq, db){
  
  res <- predict(object = db, newdata = aa_seq)[1,]
  
  print(res)
  if(is.na(res$evalue) | is.na(res$pident)){
    return(NA)
  }
  
  coverage <- (res$length/nchar(aa_seq)) * 100
  print(paste('coverage is:', coverage))
  if(res$evalue < 0.05 && res$pident >= 80 && coverage >= 80){
    pubmlst_id <- res$sseqid
    return(pubmlst_id)
  } 
  
  return(NA)
}
aa_seq <- AAStringSet('MKLEASKQKFKKIIYYKSIFFYSLYLSAFGC')


update_gb <- function(feat_id, pubmlst_id, gb){
  gb$feat_id <- str_replace_all(gb$feat_id, feat_id, pubmlst_id)
  gb$parent_ids <- str_replace_all(gb$parent_ids, feat_id, pubmlst_id)
  gb$locus_tag <- str_replace_all(gb$locus_tag, feat_id, pubmlst_id)
  gb$protein_id <- str_replace_all(gb$protein_id, feat_id, pubmlst_id)
  gb$geom_id <- str_replace_all(gb$geom_id, feat_id, pubmlst_id)
  
  return(gb)
}


#####CRITICAL
###Havent taken into account strandedness when extracting the igr sequences. 
###Havent taken into account cases where the igr length is negative. substr will return nothing if the substr() stop coord is lower than the start so this may work as is - TEST!
extract_igrs <- function(gb, fa, pubmlst_id, igrs, count){
  index <- which(gb$feat_id == paste('gene-', pubmlst_id,sep=''))
  
  #to account for cases where Blast identifies the same pubmlst id for a locus. always take the latest row.
  index <- index[length(index)]
  
  contig <- gb$seq_id[index]
  print(index)
  igrs$pubmlst_id[count] <- pubmlst_id
  igrs$strand[count] <- gb$strand[index]
  #print(gb$strand[index])
  #take into account the strand of the gene. + strand is the default forwards sequence
  if(gb$strand[index] == '+'){
    #print('up strand')
    #up_igr stop and down_igr start is just one less or more than the start/stop coords of the current locus
    igrs$up_stop[count] <- as.numeric(gb$start[index]) -1
    igrs$down_start[count] <- as.numeric(gb$end[index]) + 1
    
    #up_igr start is the stop coordinate of the previous entry, plus one. Unless the previous entry is on a different contig, meaning our gene is the first on that contig
    #so in that case the start coordinate is just 1
    if(gb$seq_id[index-2] != contig){
      igrs$up_start[count] <- 1
      igrs$up_neighbour[count] <- NA
    } else{
      igrs$up_start[count] <- as.numeric(gb$end[index-2])+1
      igrs$up_neighbour[count] <- substr(gb$feat_id[index-2], 6, nchar(gb$feat_id[index-2]))
      
    }
    
    #down_igr stop is the start coordinate of the next 'gene' type entry, minus one. Unless the next entry is on a different contig, meaning our gene is the last on the contig
    #so the stop coordinate would just be the contig length.
    if(gb$seq_id[index+2] != contig){
      igrs$down_stop[count] <- gb$end[which(gb$type == 'region' & gb$seq_id == contig)]
      igrs$down_neighbour[count] <- NA
    } else{
      igrs$down_stop[count] <- as.numeric(gb$start[index+2])-1
      igrs$down_neighbour[count] <- substr(gb$feat_id[index+2], 6, nchar(gb$feat_id[index+2]))
    }
    
    #Using the coordinates, extract the sequences of each igr
    #print('extracting up')
    igrs$up_sequence[count] <- substr(fa[contig], igrs$up_start[count], igrs$up_stop[count])
    igrs$down_sequence[count] <- substr(fa[contig], igrs$down_start[count], igrs$down_stop[count])
    
  } else{
    #the reverse strand case, need to reverse complement the sequence at the end. 
    igrs$up_stop[count] <- as.numeric(gb$end[index]) +1
    igrs$down_start[count] <- as.numeric(gb$start[index]) - 1
    
    if(gb$seq_id[index-2] != contig){
      igrs$down_stop[count] <- 1
      igrs$down_neighbour[count] <- NA
    } else{
      igrs$down_stop[count] <- as.numeric(gb$end[index-2])+1
      igrs$down_neighbour[count] <- substr(gb$feat_id[index-2], 6, nchar(gb$feat_id[index-2]))
    }
    
    
    if(gb$seq_id[index+2] != contig){
      igrs$up_start[count] <- gb$end[which(gb$type == 'region' & gb$seq_id == contig)]
      igrs$up_neighbour[count] <- NA
    } else{
      igrs$up_start[count] <- as.numeric(gb$start[index+2])-1
      igrs$up_neighbour[count] <- substr(gb$feat_id[index+2], 6, nchar(gb$feat_id[index+2]))
    }

    #Using the coordinates, extract the sequences of each igr
    up_seq <- substr(fa[contig], igrs$up_stop[count], igrs$up_start[count])
    down_seq <- substr(fa[contig], igrs$down_stop[count], igrs$down_start[count])
    
    #take reverse complement because these are on the reverse strand
    #print(up_seq)
    if(nchar(up_seq) > 0){
      print('100')
      print(reverseComplement(DNAString(up_seq)))
      igrs$up_sequence[count] <- as.character(reverseComplement(DNAString(up_seq)))
    }

    #print(down_seq)
    if(nchar(down_seq) > 0){
      igrs$down_sequence[count] <- as.character(reverseComplement(DNAString(down_seq)))
    }

    
  }
  
  igrs$up_length[count] <- nchar(igrs$up_sequence[count])
  igrs$down_length[count] <- nchar(igrs$down_sequence[count])

  return(igrs)
}


main <- function(isolate){
  
  gb <- read_gbk(paste('~/Documents/PhD/RNA_IGR/Isolate_fastas/Annotations/', isolate, '/', isolate, '.gbff', sep = ''))
  fa <- read.fasta(paste('~/Documents/PhD/RNA_IGR/Isolate_fastas/Annotations/', isolate, '/', isolate, '.fna', sep = ''), as.string = T, seqonly = F)
  
  new_gb <- gb
  count = 0
  #Initialise a df to store igrs
  igrs <- as.data.frame(matrix(data = NA, nrow = length(which(gb$type == 'CDS')), ncol = 12))
  colnames(igrs) <- c('pubmlst_id', 'strand', 'up_start', 'up_stop', 'up_sequence', 'up_length', 'up_neighbour', 'down_start', 'down_stop', 'down_sequence', 'down_length', 'down_neighbour')
  #load blast database
  db <- blast(db="~/Documents/PhD/PhD/N_meningitidis_loci/translations/BLAST_db/coding.fsa", type = 'blastp')
  
  for(row in 4:nrow(gb)){
    feature <- gb[row,]
    if(feature$type != 'CDS'){
      next
    }
    count = count + 1
    feat_id <- substr(feature$feat_id, start = 5, stop = nchar(feature$feat_id))
    
    #print('feature id extracted')
    
    #get protein sequnece of the feature
    aa_seq <- get_seq(feature, fa)
    #print('aa sequence obtained')
    
    #Use BLAST to identify the pubmlst id of the feature
    #print(aa_seq)
    pubmlst_id <- blast_seq(aa_seq, db)
    if(is.na(pubmlst_id)){
      next
    }
    #print(paste('pubmlst id is:', pubmlst_id))
    
    #Overwrite the gbk object with the pubmlst id in the feat_id, parent_ids, locus_tag, protein_id and geom_id
    #look for every instance of the name and replace it with the pubmlst id
    new_gb <- update_gb(feat_id, pubmlst_id, new_gb)
    #print(new_gb[row,])
    #print('new_gb updated')
    #Extract up and downstream igrs for the current feature and write to a dataframe
    igrs <- extract_igrs(new_gb, fa, pubmlst_id, igrs, count)
    #print('igrs extracted')
  }
  
  return(igrs)
  
}


testdf <- main('27509')
write.xlsx(testdf, "~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_test_27509.xlsx")

isolates <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
for(isolate in isolates){
  df <- main(isolate)
  write.xlsx(df, file = paste("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_", isolate, '.xlsx', sep = ''))
}


df <- main('53930')
write.xlsx(df, "~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_53930")
# for(feature in 1:nrow(test_gb)){
#   if(test_gb$type[feature] != 'CDS'){
#     next
#   }
#   direction = test_gb$strand[feature]
#   contig <- test_gb$seq_id[feature]
#   
#   start = test_gb$start[feature]
#   stop = test_gb$start[feature]
# 
#   
#   contig_seq <- test_fa[contig]
#   CDS_seq <- substr(contig_seq, start = start, stop = stop)
#   
#   if(direction == '-'){
#     CDS_seq <- reverseComplement(CDS_seq)
#   }
#   
# }
