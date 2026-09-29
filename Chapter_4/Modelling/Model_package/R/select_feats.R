#' Reduce feature space with Boruta
#'
#' @param df Input data frame of filtered Genome Comparator output.
#'
#' @param maxRuns Maximum number of feature selection iterations for Boruta to perform
#' @param out File path to a directory to save files to
#'
#' @import Boruta
#'
#' @export
select_feats_Bor <- function(df, maxRuns = 1000, out){
  selection <- Boruta::Boruta(disease ~ ., data = df, pValue = 0.05, maxRuns = maxRuns, doTrace = 2, getImp = Boruta::getImpRfZ, respect.unordered.factors = T)
  utils::write.csv(selection$finalDecision, paste(out, '/boruta_out.csv', sep=''),
            quote = F)

  grDevices::pdf(file = paste(out, '/boruta_plot.pdf', sep =''),
      width = 7,
      height = 5)
  p1 <- plot(selection, cex.axis=0.7, las=2)
  print(p1)
  grDevices::dev.off()

  keep <- names(selection$finalDecision[which(selection$finalDecision == 'Confirmed')])
  keep <- c(keep, 'disease')
  df <- df[, keep]

  utils::write.csv(x = df, file = paste(out, '/selected_features.csv', sep = ''),quote = F, row.names = F)

  return(df)
}


fit_anova <- function(x, y) {
  anova_res <- apply(x, 2, function(f) {caret::anovaScores(f, y)})
  return(anova_res)
}

#' Reduce feature space with univariate anova selection
#' Calculates ANOVA-derived p values for each feature against the predictor and keeps only the most significant features
#'
#' @param data - data frame of filtered Genome Comparator output used as training data.
#'
#' @param size Number of the most significant features to keep
#' @param out File path to a directory to save files to
#'
#' @export
select_feats_aov <- function(data, size, out){
  aov_res <- fit_anova(x = data[,-1],
                       y = data[,1])

  tophits <- c('disease', names(aov_res)[order(aov_res)[1:size]])
  return(tophits)
}


