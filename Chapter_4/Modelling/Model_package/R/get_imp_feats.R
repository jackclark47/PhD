#' Identify and plot important predictors
#' @param model A model object
#'
#' @param annotate Logical. Annotate PubMLST ids or plot the ids themselves
#' @param annotations A dataframe consisting of PuBMLST ids and the annotations to convert them to
#' @param topimportant Number of features to plot
#' @param out File path to write output files to
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @export
get_imp_feats <- function(model, annotate = F, annotations, topimportant = 20, out){

  metric <- model$importance.mode

  important_feats <- model$variable.importance %>%
    data.frame()
  colnames(important_feats) <- 'imp'
  important_feats$Locus <- names(model$variable.importance)

  important_feats <- important_feats[order(important_feats$imp, decreasing = TRUE),]
  top_important <- important_feats[1:topimportant,]

  for(i in 1:nrow(top_important)){
    top_important$Anno[i] <- annotations$Annotation[which(annotations$Locus == top_important$Locus[i])][1]
  }
  openxlsx::write.xlsx(important_feats, paste(out,'/important_feats.xlsx', sep=''))

  if(annotate) filename='annos' else filename='ids'
  if(annotate) colnames(top_important)[3]='Plot' else colnames(top_important)[2]='Plot'

  grDevices::pdf(file = paste(out, '/top_feats_', filename, '.pdf', sep =''),
                 width = 5.5,
                 height = 4)

  p <- ggplot(top_important, aes(x = stats::reorder(Plot, imp),
                                     y = imp)) +
    geom_bar(stat='identity', fill = '#E69F00') +
    coord_flip() +
    theme_classic() +
    labs(
      x = 'Feature',
      y = paste('Importance - ', metric, sep='')
    )
  print(p)

  grDevices::dev.off()

}


