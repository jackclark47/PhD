#' Define the parameter space for model optimisation
#'
#' @param n_trees_min Min number of trees
#'
#' @param n_trees_max Max number of trees
#' @param n_trees_inc Increment for trees to test
#' @param mtry_min Min number of mtry
#' @param mtry_max Max number of mtry
#' @param mtry_inc Increment for mtry values to test
#' @param nodesize_min Min node size
#' @param nodesize_max Max node size
#' @param nodesize_inc Increment for node size values to test
#' @param searchsize Number of searches in the grid to perform.
#'
#' @importFrom stats predict
#' @export
set_params <- function(n_trees_min = 500, n_trees_max = 5000, n_trees_inc = 500,
                       mtry_min = 5, mtry_max = 30, mtry_inc = 1,
                       nodesize_min = 1, nodesize_max = 50, nodesize_inc = 5,
                       searchsize = 100){

  params <- expand.grid(
    num_trees <- seq(n_trees_min, n_trees_max, by=n_trees_inc),
    mtries <- seq(mtry_min, mtry_max, by=mtry_inc),
    nodesizes <- seq(nodesize_min, nodesize_max, by=nodesize_inc)
  )

  names(params) <- c('num.trees', 'mtry', 'min.node.size')

  #Randomly sample conditions
  params <- params[sample(x=1:nrow(params), size=searchsize),]

  params$auc <- NA
  params$f1 <- NA
  params$accuracy <- NA
  params$ppv <- NA
  params$mcc <- NA
  params$pr_auc <- NA

  return(params)
}

#' Optimise a random forest model
#' @param data input training dataset
#'
#' @param params params to test
#' @param ds_fac_train Downsampling factor for training data by which to modify class weights for model training
#'
#' @export
opt_model <- function(data, params, ds_fac_train){

  for(i in 1:nrow(params)){
    print(paste('running grid search, iteration #', i, sep=''))

    #train model with kfold cv
    rows <- sample(nrow(data), size = nrow(data))

    #initalise metrics for each iteration across the k-fold cvs
    metrics <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
    colnames(metrics) <- c('accs', 'ppvs', 'aucs', 'pr_aucs', 'f1s', 'mccs')

    count = 0
    for(j in seq(0,0.9, by=0.1)){
      count = count + 1
      total_rows <- length(rows)
      start = total_rows*j
      stop = total_rows *(j+0.1)
      validation <- data[rows[start:stop],]
      train_cv <- data[-rows[start:stop],]

      if(length(unique(validation$disease)) != 2) next

      model <- ranger::ranger(disease ~ .,
                      data = train_cv,
                      num.trees = params$num.trees[i],
                      mtry = params$mtry[i],
                      importance = 'permutation',
                      write.forest = T,
                      min.node.size = params$min.node.size[i],
                      classification = T,
                      probability = T,
                      num.threads = 0,
                      class.weights = c(1, ds_fac_train),
                      respect.unordered.factors = T
      )

      #get predicted probabilites
      pred <- predict(model, validation)

      #find auc
      error <- pROC::auc(validation$disease, pred$predictions[,2])
      metrics$aucs[count] <- error

      #convert probs to 1s or 0s depending on if the prob of them being 1 is > 0.5
      rf_classes <- as.factor(dplyr::if_else(pred$predictions[,2] < 0.5, 0, 1))
      #####testing
      rf_classes <- factor(rf_classes, levels = c('1', '0'))
      ######
      y <- caret::confusionMatrix(rf_classes, validation$disease, positive = '1')
      #calculate metrics
      metrics$accs[count] <- y$overall["Accuracy"]
      metrics$ppvs[count] <- y$byClass["Pos Pred Value"]
      metrics$f1s[count] <- y$byClass["F1"]
      metrics$pr_aucs[count] <- yardstick::pr_auc_vec(truth = validation$disease, estimate = pred$predictions[,2], event_level = 'second')
      metrics$mccs[count] <- mltools::mcc(rf_classes, validation$disease)

    }

    #average metrics for the iteration
    params$auc[i] <- mean(metrics$aucs)
    params$f1[i] <- mean(metrics$f1s)
    params$accuracy[i] <- mean(metrics$accs)
    params$ppv[i] <- mean(metrics$ppvs)
    params$mcc[i] <- mean(metrics$mccs)
    params$pr_auc[i] <- mean(metrics$pr_aucs)

  }

  params <- params[order(params$f1, decreasing = TRUE),]
  return(params)
}
