#' @importFrom graphics text
#' @importFrom stats predict
#' @import ggplot2
gen_truths <- function(data, pred, out, filename){

  truths <- as.data.frame(matrix(nrow=nrow(data), ncol=2))
  colnames(truths) <- c('Truth', 'Prediction')
  truths$Truth <- data$disease
  truths$Prediction <- apply(pred$predictions, 1, which.max)
  carriage <- which(colnames(pred$predictions) == '0')
  disease <- which(colnames(pred$predictions) == '1')

  truths$Prediction[which(truths$Prediction == carriage)] <- 0
  truths$Prediction[which(truths$Prediction == disease)] <- 1

  mean(truths$Prediction == data$disease)
  yardstick::pr_auc(truths, Truth, Prediction)
  ######## testing probabilities to get more granular plots instead of discrete values Prediction values
  # truths$Prediction <- pred$predictions[,1]
  # yardstick::pr_auc(truths, Truth, Prediction)
  # ggplot2::autoplot(yardstick::pr_curve(truths, Truth, Prediction))
  ########

  grDevices::pdf(file = paste(out, filename, sep =''),
      width = 7,
      height = 5)
  truths$Truth <- factor(truths$Truth, levels = c('1', '0'))
  p <- ggplot2::autoplot(yardstick::pr_curve(truths, Truth, Prediction))
  print(p)
  grDevices::dev.off()

  return(truths)
}

#calculate metrics
model_perf <- function(data, truths, out, prefix){

  #calculate auc
  rf_perf <- pROC::auc(data$disease, truths$Prediction)
  roc_val <- pROC::roc(data$disease, truths$Prediction)

  grDevices::pdf(file = paste(out, '/', prefix, '_roc.pdf', sep =''),
      width = 7,
      height = 5)
  plot(roc_val, lwd = 1.2, ylim = c(0,1), xlim = c(1,0))
  label = paste('AUC: ', round(rf_perf, digits = 2), sep='')
  text(x = 0.80, y = 0.850, labels = label)
  grDevices::dev.off()

  y <- caret::confusionMatrix(factor(truths$Prediction, levels = c('1', '0')), data$disease, positive = '1')
  y
  utils::write.csv(y$byClass, file = paste(out, '/', prefix, '_metrics.csv', sep = ''))
  utils::write.csv(y$table, file = paste(out, '/', prefix, '_confusionMatrix.csv', sep = ''))
  return(y)
}

#' Predict disease
#' @param train Training data
#'
#' @param test Test data
#' @param params Optimised parameters for the model to use
#' @param importance Importance metric - either 'impurity' or 'permutation'
#' @param out filepath for output files to be written to
#' @param ds_fac_train downsampling factor for the training data, by which the class weights will be modified.
#'
#' @export
run_model <- function(train, test, params, importance = importance, out, ds_fac_train){

  model <- ranger::ranger(disease ~ .,
                          data = train,
                          num.trees = params$num.trees,
                          mtry = params$mtry,
                          importance = importance,
                          write.forest = T,
                          min.node.size = params$min.node.size,
                          classification = T,
                          probability = T,
                          num.threads = 0,
                          class.weights = c(1, ds_fac_train),
                          respect.unordered.factors = T
  )


  #predict for train data
  pred <- predict(model, train, type = 'response')

  truths <- gen_truths(train, pred, out, filename = '/pr_curve_train.pdf')
  model_perf(train, truths, out, prefix = 'train')

  #test on test data
  pred <- predict(model, test, type = 'response')

  truths <- gen_truths(test, pred, out, filename = '/pr_curve_test.pdf')
  model_perf(test, truths, out, prefix = 'test')

  return(model)
}
