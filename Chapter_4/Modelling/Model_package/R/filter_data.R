#' @importFrom magrittr %<>%

filter_genes <- function(df, gene_freq_cutoff){

  removables <- c()
  threshold <- nrow(df) * gene_freq_cutoff
  for(i in 2:ncol(df)){
    if(sum(is.na(df[,i])) > threshold){
      removables <- c(removables, i)
    } else if(length(unique(df[,i])) == 1){
      removables <- c(removables, i)
    }
  }

  df <- df[,-removables]
  return(df)
}

code_disease <- function(data, i){
  if(data$disease[i] == 'carrier'){
    data$disease[i] <- 0
  } else{
    data$disease[i] <- 1
  }
  return(data)
}

main_filtering <- function(data,missing_gene_cutoff=0.8){

  #set cutoff for removing isolates with NA or 0 values at more than 80% of genes in the dataset -  removes MLST only isolates etc
  n_features <- (ncol(data) - 1)
  cutoff <- missing_gene_cutoff*n_features

  removables <- c()
  for(i in 1:nrow(data)){
    count <- 0
    for(j in 1:ncol(data)){
      cell <- data[i,j]
      #convert disease and carrier to 1 and 0 values
      if(colnames(data)[j] == 'disease'){
        data <- code_disease(data, i)
        next
      }
      if(is.na(cell)){
        count = count + 1
        data[i,j] <- 0 #convert NA values to 0
        next
      }
      #convert cells with multiple alleles to 0 values #maybe this should be NAs if a model can handle NA values
      if(stringr::str_detect(cell, ';')){
        data[i,j] <- 0
      }
    }
    #if an isolate has NA values at more loci than the cutoff, remove the isolate
    if(count >= cutoff){
      removables <- c(removables, i)
    }
  }
  return(data)
}

correct_classes <- function(data){
  for(i in 2:ncol(data)){
    #CHANGE BACK TO NUMERIC IF ISSUES
    data[,i] %<>% as.character()
  }
  data$disease %<>% as.factor()
  return(data)
}

edit_names <- function(data){
  for(i in 1:ncol(data)){
    query <- colnames(data)[i]
    newname <- stringr::str_replace(query, "'|\\(|\\)", '')
    colnames(data)[i] <- newname

    query <- colnames(data)[i]
    newname <- stringr::str_replace(query, '\\)|-', '')
    colnames(data)[i] <- newname

    query <- colnames(data)[i]
    newname <- stringr::str_replace(query, '_', '')
    colnames(data)[i] <- newname
  }

  return(data)
}

#' Filter Genome Comparator data
#'
#' Removes rubbish data...
#'
#' @param data A data frame of Genome Comparator output with disease as a column
#'
#' @param gene_freq_cutoff Genes in less than this proportion of isolates in the dataset will be removed
#' @param missing_gene_cutoff Genes missing in more than this many isolates will be removed
#' @param out Filepath to a directory for writing output files
#'
#' @export
filter_data <- function(data, gene_freq_cutoff = 0.2, missing_gene_cutoff= 0.8, out=NA){

  l0 <- print(paste('The dataset contains ', nrow(data), ' isolates and ', ncol(data), ' features, including the disease class.', sep = ''))

  outfile <- as.data.frame(matrix(nrow=4))
  #remove isolates for which disease data is missing
  start_n <- nrow(data)
  keepers <- c('carrier', 'invasive (unspecified/other)', 'meningitis', 'meningitis and septicaemia', 'septicaemia', 'conjunctivitis')
  data <- data[which(data$disease %in% keepers),]
  removed <- start_n - nrow(data)
  l1 <- print(paste(removed, 'isolates lack disease data and have been removed from the dataset.'))

  #remove genes missing in more than 20% of isolates and genes with only one allele
  start_genes <- ncol(data)
  data <- filter_genes(data, gene_freq_cutoff)
  removed <- start_genes - ncol(data)
  percent <- gene_freq_cutoff *100
  l2 <- print(paste(removed, ' genes are missing in more than ', percent, '% of isolates and have been removed from the dataset.', sep = ''))


  #remove isolates with NAs in > 80% genes. Convert disease, carrier to 1, 0. Convert NAs to 0. Convert multiple alleles to 0
  #Target encode the predictors
  #####Should these be NA's instead???? I think so
  print('Binary encoding disease states. Converting NA gene values to 0. Converting multiple allele values to 0')
  start_n <- nrow(data)
  data <- main_filtering(data, missing_gene_cutoff)
  removed <- start_n - nrow(data)
  percent <- missing_gene_cutoff*100
  l3 <- print(paste(removed, ' isolates are missing data for more than ', percent,'% of genes in the dataset and have been removed.', sep = ''))

  #correct classes of all columns
  print('Setting predictors to numeric. Setting disease column to factor')
  data <- correct_classes(data)

  #edit column names to remove certain problematic characters
  print('Removing problematic characters from column names')
  data <- edit_names(data)

  if(is.character(out)){
    outfile[,1] <- c(l0, l1, l2, l3)
    utils::write.csv(as.vector(outfile$V1), file = paste(out, '/filter_log.csv', sep = ''), quote = F, row.names = F)
  }

  return(data)
}
