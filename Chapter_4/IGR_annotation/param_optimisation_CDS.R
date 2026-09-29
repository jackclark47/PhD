#Script to test running IGR annotation on the command line
library(stringr)
library(seqinr)
library(gggenomes)
library(rBLAST)
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
  make_option(c('-i', '--id'), type='integer', default=97,
              help='Sequence identity threshold for BLAST searches. [DEFAULT = 97]', metavar='integer'),
  make_option(c('-c', '--cov'), type='integer', default=70,
              help='Sequence coverage threshold for BLAST searches. [DEFAULT = 70]', metavar='integer'),
  make_option(c('-f', '--flanklen'), type='integer', default=500,
              help='Length (in base pairs) of sequence to take on both flanks of extracted IGRs for synteny calculations. [DEFAULT = 500]', metavar='integer'),
  make_option(c('-x', '--skipgff'), type='character', default=FALSE,
              help='Logical. If the user wants to provide a GFF3 file already produced by this tool they can skip the initial GFF3 creation step. [DEFAULT = FALSE]', metavar='character')
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
  
  dir.create(paste(out, 'IGR_sequences', sep=''))
  dir.create(paste(out, 'IGR_sequences_and_flanks', sep=''))
  
  igrs <- list()
  igr_and_flanks <- list()
  for(i in 1:nrow(igr_gff)){
    contigseq <- ref_seqs[[which(names(ref_seqs) == igr_gff$seq_id[i])]]
    
    startpos <- igr_gff$start[i]
    endpos <- igr_gff$end[i]
    
    if(startpos > endpos | endpos-startpos < 30) next #dont take igrs with length < 30
    
    igr <- substr(contigseq, startpos, endpos)
    names(igr) <- igr_gff$locus_tag[i]
    igrs <- append(igrs, igr)
    
    startpos <- startpos-flanklen
    endpos <- endpos+flanklen
    
    if(startpos < 1) startpos <- 1
    if(endpos > nchar(contigseq)) endpos <- nchar(contigseq)
    
    igrflank <- substr(contigseq, startpos, endpos)
    names(igrflank) <- igr_gff$locus_tag[i]
    igr_and_flanks <- append(igr_and_flanks, igrflank)
    
  }
  
  #write each to a file
  for(i in 1:length(igrs)){
    
    igrfile <- paste(out, 'IGR_sequences/', names(igrs)[i],'.fasta', sep='')
    write.fasta(sequences = igrs[i], names = names(igrs)[i], file.out = igrfile)
    
    flankfile <- paste(out, 'IGR_sequences_and_flanks/', names(igrs)[i],'.fasta', sep='')
    write.fasta(sequences = igr_and_flanks[i], names = names(igr_and_flanks)[i], file.out = flankfile)
  }
  returnable <- list(igrs, igr_and_flanks)
  return(returnable)
}

#blast igrs against query genomes
blast_igrs <- function(out, igrs, igr_and_flanks, genome_seqs, n_cores, identity = 97, coverage = 70, flanklen=500) {
  
  #Create temporary directory to store blast databases
  dir.create(paste(out, 'db', sep =''))
  
  #initialise list containing igrs + flanks
  bins <- replicate(length(igrs), list())
  #the first entry in every list should be the reference igr
  for (i in 1:length(bins)) {
    bins[[i]] <- igr_and_flanks[i]
    names(bins)[i] <- names(igrs)[i]
    nam <- paste('ref_', names(igrs)[i], sep = '')
    names(bins[[i]]) <- nam
  }
  
  #initialise list containing igrs
  igr_bins <- replicate(length(igrs), list())
  #the first entry in every list should be the reference igr
  for (i in 1:length(igr_bins)) {
    igr_bins[[i]] <- igrs[i]
    names(igr_bins)[i] <- names(igrs)[i]
    nam <- paste('ref_', names(igrs)[i], sep = '')
    names(igr_bins[[i]]) <- nam
  }
  
  iter <- length(genome_seqs)/10
  iters <- round(seq(iter, length(genome_seqs), by = iter))
  print('0% of query genomes processed')
  #load one query genome at a time and blast every igr against it
  for (i in 1:length(genome_seqs)) {
    
	print(i)
    seq <- read.fasta(genome_seqs[i], as.string = T)
    id <- str_extract(genome_seqs[i], '(?<=/)[:digit:]+')
    
    dbfile <- paste(out, 'db/', id, sep = '')
    makeblastdb(genome_seqs[i], db_name = dbfile, dbtype='nucl', verbose=F)
    blastdb <- blast(db = dbfile, type = 'blastn')
    
    if(i %in% iters){
      percent <- which(iters == i)*10
      print(paste(percent, '% of query genomes processed', sep=''))
    }
    
    cluster <- makeCluster(n_cores, outfile='')
    registerDoParallel(cluster)
    clusterCall(cluster, function()
      library(rBLAST, quietly=TRUE))
    
    results <- foreach(j = 1:length(igrs), .errorhandling = 'pass') %dopar% {
	  print(j)
      res <- predict(object = blastdb, newdata = DNAStringSet(igrs[[j]]))[1, ]
      res
    }
    
    stopCluster(cl = cluster)
    
    names(results) <- names(igrs)
    
    #Process blast hits
    for (j in 1:length(results)) {
      res <- results[[j]]
      
      if (all(is.na(res))) next
      
      rescov <- 100 * (res$length / nchar(igrs[[j]]))
      if (res$pident > identity & rescov > coverage) {
        sseq <- seq[which(names(seq) == res$sseqid)]
        
        coords <- c(res$sstart, res$send)
        startpos <- min(coords)
        endpos <- max(coords)
        
        if (coords[1] > coords[2]) {
          rev <- TRUE
        } else{
          rev <- FALSE
        }
        
        #get just the igr seqs
        igrseq <- substr(sseq, startpos, endpos)
        if (rev)
          igrseq <- as.character(reverseComplement(DNAString(igrseq)))
        
        names(igrseq) <- paste(id, names(igrs[j]), sep = '_')
        igr_bins[[j]] <- append(igr_bins[[j]], igrseq)
        
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
        
        names(sseq) <- paste(id, names(igrs[j]), sep = '_')
        bins[[j]] <- append(bins[[j]], sseq)
        
      }
    }
    #delete the blast database for the isolate
    unlink(paste(dbfile, '*', sep=''), expand=TRUE)
  }
  
  #delete the folder of blast databases
  unlink(paste(out, 'db/',sep=''))
  
  returnable <- list(bins, igr_bins)
  return(returnable)
}


get_syntenies <- function(out, bins, n_cores){
  
  # all_syntenies <- list()
  
  cluster <- makeCluster(n_cores, outfile='')
  registerDoParallel(cluster)
  clusterCall(cluster, function()
    library(DECIPHER, quietly=T))
  clusterCall(cluster, function()
    library(stringr, quietly=T))
  
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
    
    if(length(bins[[i]]) > 1){
      for(j in 2:length(bins[[i]])){
        
        Seqs2DB(DNAStringSet(unlist(bins[[i]])[c(1,j)]), type='XStringSet', identifier = names(bins[[i]][c(1,j)]), dbFile = db, verbose = F, replaceTbl = T)
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
  
  names(all_syntenies) <- names(bins)
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
find_syntenic <- function(all_syntenies, igrs, igr_bins, out, syn_threshold=0.9){
  dir.create(paste(out, 'syntenic_IGR_seqs', sep=''))
  for(i in 1:length(all_syntenies)){
    if(is.null(all_syntenies[[i]])) next
    igrbin <- all_syntenies[[i]]
    igr_id <- names(all_syntenies)[i]
    if(nrow(igrbin)==0) next

    #write a fasta file of each sequence of the same IGR above a synteny score threshold
    seqs <- list()
    seqs[1] <- igrs[igr_id]
    names(seqs)[1] <- paste('ref', igr_id, sep = '_')
    for(j in 1:nrow(igrbin)){
      if(is.na(igrbin$syn_score[j])) next
      if(igrbin$syn_score[j] > syn_threshold){
        isolate <- igrbin$isolate2[j]
        seq_id <- paste(isolate, igr_id, sep='_')
        igrseq <- igr_bins[[igr_id]][[seq_id]]
        seqs <- append(seqs, igrseq)
        names(seqs)[length(seqs)] <- seq_id
      }
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
find_synteny <- function(out, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.9, flanklen, id, cov){
  
  cat('\nIdentifying syntenic IGRs\n')
  find_syntenic(all_syntenies, igrs, igr_bins, out, syn_threshold)
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
ref_anno <- read_gff3(opt$refanno)
out <- opt$out
ref_seqs <- read.fasta(opt$refseqs, as.string = T)
n_cores <- opt$ncores
metadat <- read.xlsx(opt$metadat)
genome_seqs <- list.files(path = opt$genomeseqs, full.names = T)
id <- opt$id
cov <- opt$cov
flanklen <- opt$flanklen
skipgff <- opt$skipgff

print(class(skipgff))
print(skipgff)


########################

########################

dir.create(out)

cat('Number of cores available is', detectCores(), '\n')
#Create log file



#Run initial steps up to and including BLAST
ref_anno <- ref_anno[which(ref_anno$type == 'CDS'),]

if(skipgff == 'FALSE\r'){
  cat('\nCreating GFF file of IGR features within the reference isolate\n')
  ref_anno <- make_igr_gff(ref_anno, ref_seqs, out)
}

cat('\nExtracting IGR sequences from reference isolate\n')
res <- extract_igrs(ref_anno, ref_seqs, out, flanklen)
igrs <- res[[1]]
igr_and_flanks <- res[[2]]
cat('\nReference IGR sequences written to:', out, 'IGR_sequences/\n')
cat('\nReference IGRs and up to ',  flanklen, ' bp flanking sequences written to: ', out, 'IGR_sequences_and_flanks/\n', sep='')

#turn all query genomes into blast databases and run BLAST searches
cat('\nBLASTing reference IGR sequences against query genomes\n')
start_time <- Sys.time()
res <- blast_igrs(out, igrs, igr_and_flanks, genome_seqs, n_cores, identity=id, coverage=cov)
bins <- res[[1]]
igr_bins <- res[[2]]
end_time <- Sys.time()
tt <- round(end_time - start_time)
cat('\nCompleted in', tt, units(tt),'\n')

cat('\nCalculating synteny scores\n')
start_time <- Sys.time()
all_syntenies <- get_syntenies(out, bins, n_cores)
end_time <- Sys.time()
tt <- round(end_time - start_time)
cat('\nCompleted in', tt, units(tt),'\n')



#########################################
#test each synteny threshold
root <- paste(out, 'synteny_results/',sep='')
dir.create(root)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_10/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.1, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_20/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.2, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_30/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.3, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_40/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.4, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_50/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.5, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_60/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.6, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_70/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.7, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_80/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.8, flanklen, id, cov)

outmain <- paste(root, 'vs_CDS_feats_threshold_0_90/',sep='')
dir.create(outmain)
find_synteny(outmain, all_syntenies, igrs, igr_bins, metadat, syn_threshold = 0.9, flanklen, id, cov)
