#' Merges ids of synonymous alleles
#'
#' @param df A dataframe object containing allele ids for a set of genes akin to Genome Comparator output from PuBMLST. Should already be filtered with filter_data()
#' @param allele_dir A directory containing a set of fasta files, each of which contain all the allele sequences for a given gene. Allele ids in df will be matched to sequences in these files.
#' @param out The path and filename to write a .xlsx file of the output to.
#'
#' @export
#'
merge_synonymous <- function(df, allele_dir, out){

  skipcols <- c('id', 'isolate', 'disease')
  all_data <- list()
  for(i in 1:ncol(df)){

    gene <- colnames(df)[i]
    print(gene)
    if(gene %in% skipcols){
      next
    }

    alleles <- stats::na.omit(unique(df[,gene]))
    alleles <- alleles[which(alleles != 0)]

    #get each allele sequence
    gene_file <- seqinr::read.fasta(paste(allele_dir, gene, '.fasta', sep = ''))
    allele_data <- list()
    for(allele in alleles){

      gene_file_ids <- as.numeric(stringr::str_extract(names(gene_file), '(?<=_)[:digit:]+$'))

      allele_seq <- gene_file[which(gene_file_ids == allele)]

      if(length(allele_seq) == 0) next

      #translate sequence
      allele_seq <- seqinr::c2s(seqinr::translate(allele_seq[[1]]))

      names(allele_seq) <- allele
      allele_data <- append(allele_data, allele_seq)

    }

    if(length(allele_data) == 0) next

    all_data[[length(all_data)+1]] <- allele_data
    names(all_data)[length(all_data)] <- gene
  }

  #now align all allele seqs to see if they match. If so, assign same number
  for(i in 1:length(all_data)){
    gene_seqs <- all_data[[i]]
    gene <- names(all_data)[i]
    print(gene)
    unique_seqs <- unique(gene_seqs)

    #if all seqs are unique then no need to edit the genome comparator output
    if(length(gene) == length(unique_seqs)){
      next
    }

    #get a copy of the gc column
    temp_col <- df[,c('id', gene)]
    temp_col[,2] <- NA

    for(j in 1:length(gene_seqs)){
      allele_id <- names(gene_seqs)[j]
      index <- which(unique_seqs == gene_seqs[[j]])

      temp_col[which(df[,gene] == allele_id),2] <- index
    }

    df[,gene] <- temp_col[,2]

  }

  openxlsx::write.xlsx(df, out)

  return(df)
}
