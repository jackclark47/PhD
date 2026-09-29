#Script to test running IGR annotation on the command line
library(stringr)
library(seqinr)
library(gggenomes)
library(DECIPHER)
library(openxlsx)
library(parallelly)
library(doParallel)
library(optparse)

###IDEAS:
#Make blasted igrbins and igrflankbins be written to an intermediate directory so future runs can skip blast if a set of sequences have already been blasted


#############################

#Define arguments for command line usage
option_list = list(
  make_option(c('-a', '--refanno'), type='character', default=NULL, 
              help='Full path to gff3 file of a reference genome', metavar='character'),
  make_option(c('-o', '--out'), type='character', default=NULL,
              help='Full path to an output directory', metavar='character'),
  make_option(c('-r', '--refseqs'), type='character', default=NULL,
              help='Full path to a fasta-formatted reference genome sequence', metavar='character'),
  make_option(c('-n', '--ncores'), type='integer', default=1,
              help='Number of cores to run. [DEFAULT = 1]', metavar='integer'),
  make_option(c('-m', '--metadat'), type='character', default=NULL,
              help='Full path to a metadata file with isolate ids and disease information', metavar='character'),
  make_option(c('-g', '--genomeseqs'), type='character', default=NULL,
              help='Full path to a directory containing fasta-formatted query genome sequences for IGR annotation', metavar='character'),
  make_option(c('-i', '--id'), type='integer', default=90,
              help='Sequence identity threshold for BLAST searches. [DEFAULT = 97]', metavar='integer'),
  make_option(c('-c', '--cov'), type='integer', default=70,
              help='Sequence coverage threshold for BLAST searches. [DEFAULT = 70]', metavar='integer'),
  make_option(c('-f', '--flanklen'), type='integer', default=1000,
              help='Length (in base pairs) of sequence to take on both flanks of extracted IGRs for synteny calculations. [DEFAULT = 500]', metavar='integer'),
  make_option(c('-x', '--skipgff'), type='character', default=FALSE,
              help='Logical. If the user wants to provide a GFF3 file already produced by this tool they can skip the initial GFF3 creation step. [DEFAULT = FALSE]', metavar='character'),
  make_option(c('-s', '--syn'), default=0.4,
              help='Synteny score threshold for binning similar IGRs together. Values range from 0 to 1', metavar='character')
)

#Create parser and read the arguments 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

###################################


make_igr_gff <- function(ref_anno, ref_seqs, out){
  
  igr_gff <- as.data.frame(matrix(NA, nrow = 0, ncol = 13))
  colnames(igr_gff) <- c('seq_id', 'start', 'end', 'strand', 'type', 'source', 'locus_tag', 'upstream_gene', 'upstream_gene_strand', 'downstream_gene', 'downstream_gene_strand', 'introns', 'geom_id')
  
  contigs <- unique(ref_anno$seq_id)
  index <- 0
  for(i in 1:length(contigs)){
    contig <- ref_anno[which(ref_anno$seq_id == contigs[i]),]
    if(nrow(contig)==1) next
    #Do not get igrs at the ends of contigs
    for(j in 2:nrow(contig)){
      
      index <- index + 1
      igr_gff[index,] <- NA
      
      igr_gff$locus_tag[index] <- paste('IGR_', index, sep='')
      igr_gff$geom_id[index] <- paste('IGR_', index, sep='')
      igr_gff$seq_id[index] <- contig$seq_id[j]
      
      #get start and stop coords of the igr before the current entry in the gff file
      igr_gff$start[index] <- contig$end[j-1]
      igr_gff$end[index] <- contig$start[j]
      
      igr_gff$upstream_gene[index] <- contig$locus_tag[j-1]
      igr_gff$upstream_gene_strand[index] <- contig$strand[j-1]
      igr_gff$downstream_gene[index] <- contig$locus_tag[j]
      igr_gff$downstream_gene_strand[index] <- contig$strand[j]
      
    }
  }
  
  igr_gff$strand <- '+'
  igr_gff$type <- 'intergenic_region'
  igr_gff$source <- 'Neisseria isolates'
  
  write_gff3(igr_gff, file = paste(out, 'igr_annos.gff3', sep=''), id_var = 'locus_tag')
  return(igr_gff)
}

#Extract IGR sequences and IGR sequences + flanks
extract_igrs <- function(igr_gff, ref_seqs, out, flanklen){
  
  outdir <- paste(out, 'ref_IGRs', sep='')
  dir.create(outdir)
  
  igrs <- list()
  igr_and_flanks <- list()
  for(i in 1:nrow(igr_gff)){
    contigseq <- ref_seqs[[which(names(ref_seqs) == igr_gff$seq_id[i])]]
    
    startpos <- igr_gff$start[i]
    endpos <- igr_gff$end[i]
    
    if(startpos > endpos | endpos-startpos < 30) next #dont take igrs with length < 30
    
    igr <- substr(contigseq, startpos, endpos)
    names(igr) <- paste('ref_', igr_gff$locus_tag[i], sep='')
    igrs <- append(igrs, igr)
    
    startpos <- startpos-flanklen
    endpos <- endpos+flanklen
    
    if(startpos < 1) startpos <- 1
    if(endpos > nchar(contigseq)) endpos <- nchar(contigseq)
    
    igrflank <- substr(contigseq, startpos, endpos)
    names(igrflank) <- paste('ref_', igr_gff$locus_tag[i], sep='')
    igr_and_flanks <- append(igr_and_flanks, igrflank)
    
  }
  
  #write each to a file
  write.fasta(sequences = igrs, names = names(igrs), file.out = paste(outdir, '/IGRs.fasta', sep =''))
  write.fasta(sequences = igr_and_flanks, names = names(igr_and_flanks), file.out = paste(outdir, '/IGRs_and_flanks.fasta', sep =''))
  
  #Make files to store all BLAST hits from each IGR across all genomes
  dir.create(paste(out, 'BLAST_out', sep=''))
  dir.create(paste(out, 'BLAST_out/igrs', sep=''))
  dir.create(paste(out, 'BLAST_out/igrs_and_flanks', sep=''))
  
  #make one fasta file for each IGR
  for(i in 1:length(igrs)){
    write.fasta(sequences = igrs[i], names = names(igrs)[i], file.out = paste(out, 'BLAST_out/igrs/', substr(names(igrs)[i], 5, nchar(names(igrs[i]))), '.fasta', sep=''))
    write.fasta(sequences = igr_and_flanks[i], names = names(igr_and_flanks)[i], file.out = paste(out, 'BLAST_out/igrs_and_flanks/', substr(names(igrs)[i], 5, nchar(names(igrs[i]))), '.fasta', sep=''))
  }
}


#blast igrs against query genomes
blast_igrs <- function(out, genome_seqs, n_cores, identity = 90, coverage = 70, flanklen=1000) {
  
  #Create temporary directory to store blast databases
  dir.create(paste(out, 'db', sep =''))
  
  #dir location to store any fastas that need their names changed
  newdir <- paste(out, 'updated_fastas/', sep='')
  
  igrs <- read.fasta(paste(out, 'ref_IGRs/IGRs.fasta', sep=''), as.string = T)
  igr_and_flanks <- read.fasta(paste(out, 'ref_IGRs/IGRs_and_flanks.fasta', sep=''), as.string = T)
  
  iter <- length(genome_seqs)/10
  iters <- round(seq(iter, length(genome_seqs), by = iter))
  print('0% of query genomes processed')
  
  #load one query genome at a time and blast every igr against it
  for (i in 1:length(genome_seqs)) {
    
    filename <- genome_seqs[i]
    
    seq <- read.fasta(filename, as.string = T)
    id <- str_extract(filename, '(?<=/)[:digit:]+')
    
    #Check for special characters in the contig names
    
    if(any(stringi::stri_enc_mark(names(seq)) != 'ASCII' | any(stringi::stri_detect_regex(names(seq), '#')))){
      
      if(!dir.exists(newdir)) dir.create(newdir)
      
      names(seq) <- stringi::stri_trans_general(names(seq), 'latin-ascii')
      names(seq) <- stringi::stri_replace_all(names(seq), regex = '#', replacement = '_')
      filename <- paste(newdir, id, '.fasta', sep='')
      write.fasta(seq, names(seq), filename)
    }
    
    dbfile <- paste(out, 'db/', id, sep = '')
    
    command <- paste('makeblastdb -in ', filename,
                     ' -dbtype nucl ',
                     '-out ', dbfile, 
                     ' -title ', id,
                     sep='')
    system(command, ignore.stdout = T)
    
    if(i %in% iters){
      percent <- which(iters == i)*10
      print(paste(percent, '% of query genomes processed', sep=''))
    }
    
    command <- paste('blastn -query ', out, 'ref_IGRs/IGRs.fasta', 
                     ' -db ', dbfile,
                     ' -out ', dbfile, 'res.csv',
                     ' -max_target_seqs 1',
                     ' -max_hsps 1',
                     ' -perc_identity ', identity,
                     ' -qcov_hsp_perc ', cov,
                     ' -num_threads ', n_cores,
                     ' -outfmt "6 qseqid sseqid sstart send sstrand pident length evalue"',
                     sep = '')
    system(command, ignore.stdout = T, ignore.stderr = T)
    
    columnnames <- c('qseqid', 'sseqid', 'sstart', 'send', 'sstrand', 'pident', 'length', 'evalue')
    results <- read.table(paste(dbfile, 'res.csv', sep=''), col.names = columnnames)
    
    if(nrow(results) == 0) next
    
    
    #Process blast hits
    for (j in 1:nrow(results)) {
      res <- results[j,]
      
      sseq <- seq[which(names(seq) == res$sseqid)]
      if(length(sseq) > 1){
        print(paste('Multiple sequences obtained when querying ', res$qseqid, ' against isolate ', filename, ', skipping. The duplicated sequence id is ', res$sseqid, sep=''))
        next
      }
      coords <- c(res$sstart, res$send)
      startpos <- min(coords)
      endpos <- max(coords)
      
      #determine strandedness
      rev <- FALSE
      if (res$sstrand == 'minus') {
        rev <- TRUE
      }
      
      #get just the igr seqs
      igrseq <- substr(sseq, startpos, endpos)
      if (rev) {
        igrseq <- as.character(reverseComplement(DNAString(igrseq)))
      }
      
      igrname <- substr(res$qseqid, 5, nchar(res$qseqid))
      
      names(igrseq) <- paste(id, igrname, sep = '_')
      outfile <- paste(out, 'BLAST_out/igrs/', igrname, '.fasta', sep='')
      write.fasta(sequences = igrseq, names = names(igrseq), file.out = outfile, open = 'a')
      
      #get the igrs + flanks
      startpos <- startpos - flanklen
      endpos <- endpos + flanklen
      if (startpos < 1)
        startpos <- 1
      if (endpos > nchar(sseq))
        endpos <- nchar(sseq) - 1
      
      sseq <- substr(sseq, startpos, endpos)
      if (rev){
        sseq <- as.character(reverseComplement(DNAString(sseq)))
      }
      
      names(sseq) <- paste(id, igrname, sep = '_')
      outfile <- paste(out, 'BLAST_out/igrs_and_flanks/', igrname, '.fasta', sep='')
      write.fasta(sequences = sseq, names = names(sseq), file.out = outfile, open = 'a')
      
    }
    #delete the blast database for the isolate
    unlink(paste(dbfile, '*', sep=''), expand=TRUE)
  }
  
  #delete the folder of blast databases
  unlink(paste(out, 'db/',sep=''))
}


get_syntenies <- function(out, n_cores){
  
  bins <- list.files(paste(out, 'BLAST_out/igrs_and_flanks/', sep=''), full.names=T)

  cluster <- makeCluster(n_cores, outfile='')
  registerDoParallel(cluster)
  
  clusterCall(cluster, function()
    suppressPackageStartupMessages(library(DECIPHER, quietly=T)))
  clusterCall(cluster, function()
    suppressPackageStartupMessages(library(stringr, quietly=T)))
  clusterCall(cluster, function()
    suppressPackageStartupMessages(library(seqinr, quietly=T)))
  
  all_syntenies <- foreach(i = 1:length(bins), .export = c("process_synteny")) %dopar% {
    
    iter <- length(bins)/10
    iters <- round(seq(iter, length(bins), by = iter))
    
    if(i %in% iters){
      percent <- which(iters == i)*10
      print(paste('Syntenies calculated for ', percent, '% of IGRs', sep=''))
    }
    
    db <- paste(out, 'temp_synteny_db', i, sep ='')
    
    res_tab <- data.frame("isolate1" = character(),
                          "isolate2" = character(),
                          "length1" = integer(),
                          "length2" = integer(),
                          "overlap" = integer(),
                          "blocks" = integer(),
                          "synteny_score" = integer())
    
    seqs <- read.fasta(bins[i], as.string = T)
    seqs <- DNAStringSet(unlist(seqs))
    if(length(seqs) > 1){
      for(j in 2:length(seqs)){
        
        Seqs2DB(seqs[c(1,j)], type='XStringSet', identifier = names(seqs[c(1, j)]), dbFile = db, verbose=F, replaceTbl=T)
        syn_res <- FindSynteny(db,
                               maxSep=15, #maximum gap size between hits within a block
                               maxGap=15, #maximum number of gaps between hits within a block
                               useFrames = F,
                               processors = 1,
                               verbose = F
        )
        if(nrow(syn_res[2,1][[1]]) == 0) next
        
        res_tab <- process_synteny(syn_res, res_tab)
        
      }
      
      unlink(db)
    }
    
    res_tab
  }
  
  stopCluster(cl = cluster)
  
  nams <- str_extract(bins, '(?<=/)([:alnum:]|_)+(?=.fasta)')
  names(all_syntenies) <- nams
  return(all_syntenies)
}

process_synteny <- function(syn_res, res_tab){
  for(j in 2:ncol(syn_res)){
    for(k in 1:(j-1)){
      isolate1 <- as.character(str_split(colnames(syn_res)[k], "_")[[1]][1])
      isolate2 <- as.character(str_split(colnames(syn_res)[j], "_")[[1]][1])
      
      #only want synteny vs the ref isolate
      if(isolate1 != 'ref' & isolate2 != 'ref') next
      length1 <- syn_res[k,k][[1]]
      length2 <- syn_res[j,j][[1]]
      overlap <- sum(syn_res[k,j][[1]][,4])
      blocks <- nrow(syn_res[j,k][[1]])
      
      lengths <- c(length1, length2)
      syn_score <- 1+log10((overlap/min(lengths))/blocks)
      
      temprow <- data.frame(isolate1, isolate2, length1, length2, overlap, blocks, syn_score)
      res_tab <- rbind(res_tab, temprow)
    }
  }
  return(res_tab)
}

#if synteny is high enough, the sequences are the same IGR.
find_syntenic <- function(all_syntenies, out, syn_threshold=0.9){
  
  igrs <- read.fasta(paste(out, 'ref_IGRs/IGRs.fasta', sep=''), as.string = T)
  dir.create(paste(out, 'syntenic_IGR_seqs', sep=''))
  
  for(i in 1:length(all_syntenies)){
    if(is.null(all_syntenies[[i]])) next
    igrbin <- all_syntenies[[i]]
    igrbin <- igrbin[which(igrbin$syn_score > syn_threshold),]
    if(nrow(igrbin)==0) next
    
    igr_id <- names(all_syntenies)[[i]]
    
    #write a fasta file of each sequence of the same IGR above a synteny score threshold
    seqs <- list()
    seqs[1] <- igrs[igr_id]
    names(seqs) <- paste('ref_', igr_id, sep='')
    igrseqs <- read.fasta(paste(out, 'BLAST_out/igrs/', igr_id, '.fasta', sep=''), as.string = T)
    
    for(j in 1:nrow(igrbin)){
      if(is.na(igrbin$syn_score[j])) next
      isolate <- igrbin$isolate2[j]
      seq_id <- paste(isolate, igr_id, sep='_')
      igrseq <- igrseqs[which(names(igrseqs) == seq_id)]
      seqs <- append(seqs, igrseq)
      names(seqs)[length(seqs)] <- seq_id
    }
    
    out_file <- paste(out, 'syntenic_IGR_seqs/', igr_id, '.fasta', sep='')
    write.fasta(seqs, names(seqs), file.out = out_file)
    
  }
}

#annotate alleles of each IGR
assign_alleles <- function(filelist, metadat){
  
  #for every fasta file in the directory, read the igr sequences, align them, and assign allele numbers
  igr_ids <- str_extract(filelist, '(?<=/)([:alnum:]|_)+(?=.fasta)')
  
  #initialise annotation table
  
  igr_annos <- as.data.frame(matrix(NA, nrow = nrow(metadat), ncol = 2+length(filelist)))
  colnames(igr_annos) <- c('id', 'disease', igr_ids)
  igr_annos$id <- metadat$id
  igr_annos$disease <- metadat$disease
  
  #populate annotation table
  for(i in 1:length(filelist)){
    seqs <- read.fasta(filelist[i], as.string = T)
    unique_seqs <- unique(as.character(seqs))
    #write alleles and their sequence to a file
    
    igr_id <- igr_ids[i]
    
    alleles <- match(seqs, unique_seqs)
    ids <- str_extract(names(seqs), '^[:alnum:]+')
    
    for(j in 1:length(ids)){
      if(ids[j] %in% igr_annos$id){
        igr_annos[which(igr_annos$id == ids[j]),igr_id] <- alleles[j]
      } else{
        next
      }
    }
    
  }
  return(igr_annos)
}

#full pipeline:
find_synteny <- function(ref_anno, ref_seqs, genome_seqs, n_cores, skipgff, out, metadat, syn_threshold = 0.5, flanklen, id, cov){
  
  ref_anno <- ref_anno[which(ref_anno$type == 'CDS'),]
  
  if(skipgff != 'TRUE'){
    cat('\nCreating GFF file of IGR features within the reference isolate\n')
    ref_anno <- make_igr_gff(ref_anno, ref_seqs, out)
  }
  
  cat('\nExtracting IGR sequences from reference isolate\n')
  extract_igrs(ref_anno, ref_seqs, out, flanklen)
  cat('\nReference IGR sequences written to: ', out, 'ref_IGRs/IGRs.fasta\n' , sep='')
  cat('\nReference IGRs and up to ',  flanklen, ' bp flanking sequences written to: ', out, 'ref_IGRs/IGRs_and_flanks.fasta\n', sep='')
  
  #turn all query genomes into blast databases and run BLAST searches
  cat('\nBLASTing reference IGR sequences against query genomes\n')
  start_time <- Sys.time()
  blast_igrs(out, genome_seqs, n_cores, identity=id, coverage=cov, flanklen = flanklen)
  end_time <- Sys.time()
  tt <- round(end_time - start_time)
  cat('\nCompleted in', tt, units(tt),'\n')
  
  
  cat('\nCalculating synteny scores\n')
  start_time <- Sys.time()
  all_syntenies <- get_syntenies(out, n_cores)
  end_time <- Sys.time()
  tt <- round(end_time - start_time)
  cat('\nCompleted in', tt, units(tt),'\n')
  
  
  
  cat('\nIdentifying syntenic IGRs\n')
  find_syntenic(all_syntenies, out, syn_threshold)
  cat('\nSequences of syntenic IGRs were written to: ', out, 'syntenic_IGR_seqs/\n', sep='')
  
  cat('\nAnnotating alleles of syntenic IGRs\n')
  start_time <- Sys.time()
  outdir <- paste(out, 'syntenic_IGR_seqs/', sep='')
  filelist <- list.files(outdir, full.names = T)
  igr_annos <- assign_alleles(filelist, metadat)
  end_time <- Sys.time()
  tt <- round(end_time - start_time)
  cat('\nCompleted in', tt, units(tt),'\n')
  
  filename <- paste(out, 'igr_annos_synthres-', syn_threshold, '_flank-', flanklen, '_id-', id, '_cov-', cov, '.xlsx', sep='' )
  
  write.xlsx(igr_annos, filename)
  cat('\nIGR annotations written to:', filename, '\n')
}

########################
#Assign arguments passed from the command line
# ref_anno <- read_gff3(opt$refanno)
# out <- opt$out
# ref_seqs <- read.fasta(opt$refseqs, as.string = T)
# n_cores <- opt$ncores
# metadat <- read.xlsx(opt$metadat)
# genome_seqs <- list.files(path = opt$genomeseqs, full.names = T)
# id <- opt$id
# cov <- opt$cov
# flanklen <- opt$flanklen
# skipgff <- opt$skipgff
# syn_threshold <- opt$syn


########################
ref_anno <- read_gff3('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/Annotation/references/240.gff3')
out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/Annotation/out/'
ref_seqs <- read.fasta('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/B/cc41_44/IGRs/Annotation/references/240.fasta', as.string = T)
n_cores <- 7
metadat <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/Annotation/metadata.xlsx')
genome_seqs <- list.files(path = '/Users/u5501917/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/Annotation/test_fastas', full.names = T)
id <- 90
cov <- 70
flanklen <- 500
skipgff <- FALSE
syn_threshold <- 0.4
########################
#/Users/u5501917/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/Annotation/test_fastas

dir.create(out)
cat('Number of cores available is', detectCores(), '\n')

#Create log file
# seq_ids <- str_extract(genome_seqs, '(?<=/)([:digit:]+)(?=.fasta)')
# genome_seqs <- genome_seqs[which(seq_ids %in% data_init$id)]

#annotate igrs
find_synteny(ref_anno = ref_anno, ref_seqs = ref_seqs, genome_seqs = genome_seqs, 
             n_cores = n_cores, skipgff = skipgff, out = out, metadat = metadat,
             syn_threshold = syn_threshold, flanklen = flanklen, id = id, cov = cov)

