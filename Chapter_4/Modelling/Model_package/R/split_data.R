#' @importFrom graphics text
#' @import ggplot2
#' @import ROSE
data.selection <- function(df, train_proportion = 0.7){
  nrows <- nrow(df)
  indices <- sample(1:nrows, train_proportion*nrows)

  train <- df[indices,]
  test <- df[-indices,]
  returnable <- list(train, test)

  return(returnable)
}

downsample <- function(dataset, data_name, fix_imbalance = T, out){

  freqs <- as.data.frame(table(dataset["disease"]))
  freqs$disease <- c('Carriage', 'Disease')

  outcsv <- as.data.frame(matrix(nrow = 4))
  l1 = print(paste('Number of carriage isolates is ', freqs$Freq[1], sep=''))
  l2 = print(paste('Number of invasive isolates is ', freqs$Freq[2], sep=''))
  outcsv[c(1,2),1] <- c(l1, l2)
  grDevices::pdf(file = paste(out, '.pdf', sep =''),
      width = 4,
      height = 4)
  p <- ggplot2::ggplot(data=freqs, aes(x=disease, y=Freq)) +
    geom_bar(stat = 'identity') +
    theme_classic() +
    ggokabeito::scale_color_okabe_ito() +
    xlab('Class') +
    ylab('Frequency') +
    ggtitle(paste('Class frequencies for ', data_name, 'ing data.', sep =''))
  print(p)
  grDevices::dev.off()

  n_maj_old <- sort(freqs$Freq)[2]

  if(fix_imbalance){
    #N <- 2*sort(freqs$Freq)[1]

    p_minority <- sort(freqs$Freq)[1]/(sum(freqs$Freq))
    if(p_minority < 0.3) p=0.3 else p=(p_minority+0.05)


    data <- dataset
    data <- ovun.sample(formula = disease~., data=data,
                        method='under', p=p)$data

    freqs <- as.data.frame(table(data["disease"]))
    freqs$disease <- c('Carriage', 'Disease')

    l3 = print(paste('Downsampled number of carriage isolates is ', freqs$Freq[1], sep=''))
    l4 = print(paste('Downsampled number of invasive isolates is ', freqs$Freq[2], sep=''))
    outcsv[c(3,4),1] <- c(l3,l4)

    grDevices::pdf(file = paste(out, '.pdf', sep =''),
                   width = 4,
                   height = 4)
    p <- ggplot2::ggplot(data=freqs, aes(x=disease, y=Freq)) +
      geom_bar(stat = 'identity') +
      theme_classic() +
      ggokabeito::scale_color_okabe_ito() +
      xlab('Class') +
      ylab('Frequency') +
      ggtitle(paste('Corrected class frequencies for ', data_name, 'ing data.', sep =''))
    print(p)
    grDevices::dev.off()
  }

  n_maj_new <- sort(freqs$Freq)[2]
  downsample_factor <- n_maj_old/n_maj_new


  utils::write.csv(as.vector(outcsv$V1), file = paste(out, '.csv', sep=''), quote = F, row.names = F)

  return(list(data, downsample_factor))
}

#' Split data into train and test datasets
#'
#' @param df Input data frame of filtered and feature selected Genome Comparator output
#'
#' @param train_proportion Proportion of the dataset to use for training. The remainder will be used as test data
#' @param balance_test_data Should the proportion of carriage and disease isolates in the test dataset be equalised?
#' @param out Filepath to write files to
#'
#' @export
split_data <- function(df, train_proportion = 0.7, balance_test_data = FALSE, out){
  sets <- data.selection(df, train_proportion)

  #check class imbalances. Fix imbalances with downsampling, as typically better when lots of observations available
  trainout <- paste(out, '/train_class_freqs', sep = '')
  testout <- paste(out, '/test_class_freqs', sep ='')

  train <- downsample(dataset=sets[[1]], data_name='train', fix_imbalance = T, out = trainout)
  ds_fac_train <- train[[2]]
  train <- train[[1]]
  test <- downsample(dataset=sets[[2]], data_name='test', fix_imbalance = balance_test_data, out = testout)[[1]]

  returnable <- list(train, test, ds_fac_train)
  return(returnable)
}
