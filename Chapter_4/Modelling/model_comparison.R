#model comparison using optimal models created with each ML algorithm
#each model is tuned.
#models are then evaluated by comparing ROC AUCs, PR AUCs, and MCCs.
#First in a table of absolute values for each.
#Then using stats to compare the distributions for all models for significance. 
library(ggplot2)
library(ggokabeito)
library(ROSE)
library(dplyr)
library(caret)
library(ROCR)
library(ranger)
library(pROC)
library(gbm)
library(mltools)
library(stringr)
library(e1071)
library(class)
library(naivebayes)
library(pathopred)

set.seed(seed = 532)

#process dataset - already feature selected with boruta.
print('loading datasets')
dataset <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/backups_before_unordered_factor_fix/D_All/selected_features.csv')[,-1]
dataset$disease <- as.factor(dataset$disease)
#convert all columns to factors - theyre all nominal variables
for(i in 1:ncol(dataset)){
  dataset[,i] <- as.factor(as.character(dataset[,i]))
}

#get train and test data
sets <- split_data(dataset, out = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison', train_proportion = 0.7)
train <- sets[[1]]
test <- sets[[2]]
indices <- sets[[3]] #save row numbers of rows in dataset used for training
diseasecol <- which(colnames(train) == 'disease')

#create output table for storing model performances
all_perfs <- as.data.frame(matrix(data=NA, nrow=6, ncol = 7))
colnames(all_perfs) <- c('Model', 'Accuracy', 'PPV', 'ROC-AUC', 'PR-AUC', 'F1-macro', 'MCC')
all_perfs$Model <- c('random_forest', 'gbm', 'svm',
                     'knn', 'logistic_regression', 'naive_bayes')

#create table for storing optimal HPs per model
HPs <- as.data.frame(matrix(data=NA, nrow=6, ncol=3))
colnames(HPs) <- c('Model', 'values', 'HPs')
HPs$Model <- all_perfs$Model


#Random forest
print('Tuning HPs for random forest modelling...')
#HPs:
#mtry - number of features considered in each split
#nodesize - #minimum number of samples in each terminal node
#ntrees - number of trees in the forest
#set HPs to test
params <- expand.grid(
  num_trees <- seq(500, 5000, by=500),
  mtries <- seq(5, 30, by = 1),
  nodesizes <- seq(1, 50, by = 5)
)
names(params) <- c('num.trees', 'mtry', 'min.node.size')

#2600 to test - too many, so randomly sample 100 conditions 
params <- params[sample(x=1:nrow(params), size=100),]

rf_perfs <- c()
for(i in 1:100){
  print(paste('running random grid search, iteration #', i, sep=''))
  #train model
  model <- ranger::ranger(disease ~ .,
              data = train,
              num.trees = params$num.trees[i],
              mtry = params$mtry[i],
              importance = 'impurity',
              write.forest = T,
              min.node.size = params$min.node.size[i],
              classification = T,
              probability = T
  )


  #predict for train data
  prediction <- predict(model, train, type = 'response')
  #calculate auc
  rf_perf <- pROC::auc(train$disease, prediction$predictions[,2])
  rf_perfs <- c(rf_perfs, rf_perf)
}

#identify param set with the highest auc
rf_opt_params <- params[which.max(rf_perfs),]

#save optimal parameters
HPs$values[1] <- str_c(rf_opt_params[1,], collapse = '; ')
HPs$HPs[1] <- str_c(colnames(rf_opt_params), collapse = '; ')

#save performances across all parameters
params$auc <- rf_perfs
write.csv(params, '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/rf_HPtuning_perfs.csv')
print('random forest tuned, saving performances for each parameter tested')
print('running optimal rf with 10 fold cv')
#now run rf with optimal params, with 10-fold cv
rows <- sample(nrow(train), size = nrow(train))

#initalise metrics for each iteration across the k-fold cvs
metrics <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
colnames(metrics) <- c('accs', 'ppvs', 'aucs', 'pr_aucs', 'f1s', 'mccs')

count <- 0
pdf(file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/rf_rocs.pdf',
    width = 7,
    height = 5)
for(i in seq(0,0.9, by=0.1)){
  count = count + 1
  total_rows <- length(rows)
  start = total_rows*i 
  stop = total_rows *(i+0.1)
  validation <- train[rows[start:stop],]
  train_cv <- train[-rows[start:stop],]
  
  model <- ranger(disease ~ .,
                  data = train_cv,
                  num.trees = rf_opt_params$num.trees,
                  mtry = rf_opt_params$mtry,
                  min.node.size = rf_opt_params$nodesizes,
                  importance = 'impurity',
                  write.forest = T,
                  classification = T,
                  probability = T
  )
  
  #get predicted probabilites
  pred <- predict(model, validation)
  
  #find auc
  error <- auc(validation$disease, pred$predictions[,2])
  metrics$aucs[count] <- error
  
  #plot roc
  roc_val <- roc(validation$disease, pred$predictions[,2])
  if(i == 0){
    plot(roc_val, lwd = 1.2, ylim = c(0,1), xlim = c(1,0))
  }
  plot(roc_val, add = T, lwd = 1.2)
  
  #convert probs to 1s or 0s depending on if the prob of them being 1 is > 0.5
  rf_classes <- as.factor(if_else(pred$predictions[,2] > 0.5, 1, 0))
  y <- confusionMatrix(rf_classes, validation$disease, positive = '1')
  #calculate metrics
  metrics$accs[count] <- y$overall["Accuracy"]
  metrics$ppvs[count] <- y$byClass["Pos Pred Value"]
  metrics$f1s[count] <- y$byClass["F1"]
  metrics$pr_aucs[count] <- pr_auc_vec(truth = validation$disease, estimate = pred$predictions[,2], event_level = 'second')
  metrics$mccs[count] <- mltools::mcc(rf_classes, validation$disease)

}
title('Cross-validated random forest ROCs', outer = T, line = -1)
dev.off()

#save metrics to model summary df
means <- lapply(metrics, mean)
all_perfs[1, 2:7] <- means
print('done')

###gradient boosting machine
print('Optimising HPs for gradient boosting...')
#HPs:
#n.trees
#interaction.depth (tree depth)
#shrinkage (learning rate)
#n.minobsinnode 
#set hyperparameter variables to test
params <- expand.grid(
  num_trees <- seq(500, 7000, by=500),
  tree_depths <- seq(1, 10, by = 1),
  shrinkages <- c(0.001, 0.01, 0.05, seq(0.1, 0.5, by = 0.1))
)
names(params) <- c('num_trees', 'tree_depths', 'shrinkages')

#1600 to test - too many, so randomly sample 100 conditions 
params <- params[sample(x=1:nrow(params), size=100),]

#tune model. Binary classification so distribution should be bernoulli, the default.
gbm_perfs <- c()
for(i in 1:100){
  print(paste('running random grid search, iteration #', i, sep=''))
  model <- gbm.fit(x = train[, -diseasecol],
               y = as.numeric(as.character(train$disease)),
               distribution = 'bernoulli',
               n.trees = params$num_trees[i],
               shrinkage = params$shrinkages[i],
               interaction.depth = params$tree_depths[i]
  )

  prediction <- predict.gbm(model, train[,-diseasecol], type='response')
  gbm_perf <- auc(train$disease, prediction)
  gbm_perfs <- c(gbm_perfs, gbm_perf)
}

#identify param set with the highest auc
gbm_opt_params <- params[which.max(gbm_perfs),]
gbm_opt_params <- params[1,]
gbm_opt_params$num_trees <- 5500
gbm_opt_params$tree_depths <- 9
gbm_opt_params$shrinkages <- 0.2
#save optimal parameters
HPs$values[2] <- str_c(gbm_opt_params[1,], collapse = '; ')
HPs$HPs[2] <- str_c(colnames(gbm_opt_params), collapse = '; ')

#save performances across all parameters
params$auc <- gbm_perfs
write.csv(params, '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/gbm_HPtuning_perfs.csv')
print('Optimisation complete, writing performances for each parameter set tested')
print('Running optimised gbm model with 10-fold cross validation')

#run gbm with optimal params
#initalise metrics for each iteration across the k-fold cvs
metrics <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
colnames(metrics) <- c('accs', 'ppvs', 'aucs', 'pr_aucs', 'f1s', 'mccs')
rows <- sample(nrow(train), size = nrow(train))

count <- 0
pdf(file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/gbm_rocs_take2.pdf',
    width = 7,
    height = 5)
for(i in seq(0,0.9, by=0.1)){
  
  count = count + 1
  total_rows <- length(rows)
  start = total_rows*i 
  stop = total_rows *(i+0.1)
  validation <- train[rows[start:stop],]
  train_cv <- train[-rows[start:stop],]

  model <- gbm.fit(x = train_cv[, -diseasecol],
                       y = as.numeric(as.character(train_cv$disease)),
                       distribution = 'bernoulli',
                       n.trees = gbm_opt_params$num_trees,
                       shrinkage = gbm_opt_params$shrinkages,
                       interaction.depth = gbm_opt_params$tree_depths)
  
  #get predicted probabilites
  pred <- predict(model, validation[,-diseasecol], type = 'response')
  #find auc
  error <- auc(validation$disease, pred)
  metrics$aucs[count] <- error
  
  #plot roc
  roc_val <- roc(validation$disease, pred)
  if(i == 0){
    plot(roc_val, lwd = 1.2, ylim = c(0,1), xlim = c(1,0))
  }
  plot(roc_val, add = T, lwd = 1.2)
  
  #convert probs to 1s or 0s depending on if the prob of them being 1 is > 0.5
  rf_classes <- as.factor(if_else(pred > 0.5, 1, 0))
  y <- confusionMatrix(rf_classes, validation$disease, positive = '1')
  #calculate metrics
  metrics$accs[count] <- y$overall["Accuracy"]
  metrics$ppvs[count] <- y$byClass["Pos Pred Value"]
  metrics$f1s[count] <- y$byClass["F1"]
  metrics$pr_aucs[count] <- pr_auc_vec(truth = validation$disease, estimate = pred, event_level = 'second')
  metrics$mccs[count] <- mltools::mcc(rf_classes, validation$disease)
}
title('Cross-validated gbm ROCs', outer = T, line = -1)
dev.off()

#save metrics to model summary df
means <- lapply(metrics, mean)
all_perfs[2, 2:7] <- means

print('done')

#Naive bayes
print('Optimising HPs for Naive Bayes modelling')
#HPs
#laplace
params <- as.data.frame(matrix(data = 1:100, nrow = 100, ncol =1))
names(params) <- 'laplace'

#100 to test 
#tune model. Binary classification so distribution should be bernoulli, the default.
nb_perfs <- c()
for(i in 1:100){
  print(paste('running random grid search, iteration #', i, sep=''))
  model <- naive_bayes(train[,-diseasecol],
                       train$disease,
                       laplace=params$laplace[i]
                       )
  
  prediction <- predict(model, train[,-diseasecol], type='prob')
  nb_perf <- auc(train$disease, prediction[,2])
  nb_perfs <- c(nb_perfs, nb_perf)
}


#identify param set with the highest auc
nb_opt_params <- params[which.max(nb_perfs),]

#save optimal parameters
HPs$values[6] <- str_c(nb_opt_params[1], collapse = '; ')
HPs$HPs[6] <- str_c(colnames(nb_opt_params), collapse = '; ')

#save performances across all parameters
params$auc <- nb_perfs
write.csv(params, '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/nb_HPtuning_perfs.csv')
print('Optimisation complete. Writing performances for each parameter set')
print('Running optimsed naive Bayes model...')
#run nb with optimal params
#initalise metrics for each iteration across the k-fold cvs
metrics <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
colnames(metrics) <- c('accs', 'ppvs', 'aucs', 'pr_aucs', 'f1s', 'mccs')

count <- 0
pdf(file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/nb_rocs.pdf',
    width = 7,
    height = 5)
for(i in seq(0,0.9, by=0.1)){
  
  count = count + 1
  total_rows <- length(rows)
  start = total_rows*i 
  stop = total_rows *(i+0.1)
  validation <- train[rows[start:stop],]
  train_cv <- train[-rows[start:stop],]
  
  model <- naive_bayes(train[,-diseasecol],
                       train$disease,
                       laplace=nb_opt_params
  )
  
  #get predicted probabilites
  pred <- predict(model, validation[,-diseasecol], type = 'prob')
  #find auc
  error <- auc(validation$disease, pred[,2])
  metrics$aucs[count] <- error
  
  #plot roc
  roc_val <- roc(validation$disease, pred[,2])
  if(i == 0){
    plot(roc_val, lwd = 1.2, ylim = c(0,1), xlim = c(1,0))
  }
  plot(roc_val, add = T, lwd = 1.2)
  
  #convert probs to 1s or 0s depending on if the prob of them being 1 is > 0.5
  rf_classes <- as.factor(if_else(pred[,2] > 0.5, 1, 0))
  y <- confusionMatrix(rf_classes, validation$disease, positive = '1')
  #calculate metrics
  metrics$accs[count] <- y$overall["Accuracy"]
  metrics$ppvs[count] <- y$byClass["Pos Pred Value"]
  metrics$f1s[count] <- y$byClass["F1"]
  metrics$pr_aucs[count] <- pr_auc_vec(pred[,2], truth = validation$disease, event_level = 'second')
  metrics$mccs[count] <- mltools::mcc(rf_classes, validation$disease)
}
title('Cross-validated nb ROCs', outer = T, line = -1)
dev.off()

#save metrics to model summary df
means <- lapply(metrics, mean)
all_perfs[6, 2:7] <- means
print('done')

#support vector machine

#HPs:
#sigma
#cost

print('One-hot encoding test and training data ')
#first one-hot encode the data, make all predictors character vectors. Useful for later models as well. 
one_hot_encoding = function(df, columns="season"){
  # create a copy of the original data.frame for not modifying the original
  df = cbind(df)
  # convert the columns to vector in case it is a string
  columns = c(columns)
  # for each variable perform the One hot encoding
  count = 1
  for (column in columns){
    print(count)
    unique_values = sort(unique(df[column])[,column])
    non_reference_values  = unique_values[c(-1)] # the first element is going 
    # to be the reference by default
    for (value in non_reference_values){
      # the new dummy column name
      new_col_name = paste0(column,'.',value)
      # create new dummy column for each value of the non_reference_values
      df[new_col_name] <- with(df, ifelse(df[,column] == value, 1, 0))
      
    }
    # delete the one hot encoded column
    df[column] = NULL
    count = count + 1
  }
  return(df)
}
dataset_encode <- one_hot_encoding(dataset, columns = colnames(dataset)[-diseasecol])
weakfeats <- caret::nearZeroVar(dataset_encode, saveMetrics = T) %>%
  tibble::rownames_to_column()

#~half of the features have a frequency ratio of 4306, meaning the most common allele is present in 4306 isolates
#and the second most common is present in only 1 isolate - these are very low value predictors and need to be removed.
#similarly ratios for 2152.5 and 1434.666 where the second most freq allele is present in only 2 or 3 isolates account for 7.5k cols.

#table(weakfeats$freqRatio)
#huge number of features (~32000) - absolutely need feature selection both for accuracy, runtime, and memory reasons
#remove any column with variation in fewer than 1% of isolates
4307 * 0.01 #1% cutoff
cutoff = 43.08
weakfeats2 <- weakfeats[which(weakfeats$freqRatio >= cutoff),]
dataset_encode2 <- dataset_encode[,which(!(colnames(dataset_encode) %in% weakfeats2$rowname))]
#dataset_encode3 <- Boruta(disease ~., data = dataset_encode2, doTract = 2, maxRuns = 100)
colnames(dataset_encode2)

print('Trimming encoded features with variation in fewer than 1% of isolates')
indices <- 1:2651
train_encode <- dataset_encode2[indices,]
test_encode <- dataset_encode2[-indices,]

#set categorical predictors to character class

#now run feature selection
#dataset_encode3 <- Boruta(disease ~., data = train_encode, doTract = 2, maxRuns = 100)

encodeddiseasecol <- which(colnames(train_encode) == 'disease')


denom <- ncol(train_encode)
print('Tuning HPs for support vector machine modelling...')
params <- expand.grid(
  costs <- c(2^-2, 2^-1, 2^0, 2^1, 2^2),
  sigmas <- c(0.5/denom, 1/denom, 10/denom, 100/denom, 0.01, seq(0.1, 0.5, 0.1))
)
names(params) <- c('costs', 'sigmas')

#only 50 variables
svm_perfs <- c()
for(i in 1:50){
  print(paste('running complete grid search, iteration #', i, sep=''))
  x <- svm(x=train_encode[,-encodeddiseasecol], y=train_encode$disease,
           type = 'C',
           kernel = 'radial',
           gamma = params$sigmas[i],
           cost = params$costs[i],
           cross = 1,
           probability = T
  )
  pred <- predict(x, train_encode[,-encodeddiseasecol], probability = TRUE)
  probs <- attr(pred, 'probabilities')
  
  #find auc
  svm_perf <- auc(train_encode$disease, probs[,1])
  svm_perfs <- c(svm_perfs, svm_perf)
}

#identify param set with the highest auc
svm_opt_params <- params[which.max(svm_perfs),]
svm_opt_params <- params[1,]
svm_opt_params$sigmas <- 0.0383582662063675
svm_opt_params$costs <- 1

#save optimal parameters
HPs$values[3] <- str_c(svm_opt_params[1,], collapse = '; ')
HPs$HPs[3] <- str_c(colnames(svm_opt_params), collapse = '; ')

#save performances across all parameters
params$auc <- svm_perfs
write.csv(params, '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/svm_HPtuning_perfs.csv')
print('Optimisation complete. Writing performances for each parameter set tested.')
print('Running optimised svm with 10-fold cross validation.')
#now run svm with optimal params, with 10-fold cv using rows from the rf iteration
#initalise metrics for each iteration across the k-fold cvs
metrics <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
colnames(metrics) <- c('accs', 'ppvs', 'aucs', 'pr_aucs', 'f1s', 'mccs')

count <- 0
pdf(file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/svm_rocs_take2.pdf',
    width = 7,
    height = 5)
for(i in seq(0,0.9, by=0.1)){
  print(i)
  count = count + 1
  total_rows <- length(rows)
  start = total_rows*i 
  stop = total_rows *(i+0.1)
  validation <- train_encode[rows[start:stop],]
  train_cv <- train_encode[-rows[start:stop],]

  model <- svm(x=train_cv[,-encodeddiseasecol], y=train_cv$disease,
                   type = 'C',
                   kernel = 'radial',
                   gamma = svm_opt_params$sigmas,
                   cost = svm_opt_params$costs,
                   cross = 1,
               probability = T
  )

  #get predicted probabilites
  pred <- predict(model, validation[,-encodeddiseasecol], probability = TRUE)
  probs <- attr(pred, 'probabilities')
  
  #find auc
  error <- auc(na.omit(validation$disease), probs[,1])
  metrics$aucs[count] <- error
  
  #plot roc
  roc_val <- roc(na.omit(validation$disease), probs[,1])
  if(i == 0){
    plot(roc_val, lwd = 1.2, ylim = c(0,1), xlim = c(1,0))
  }
  plot(roc_val, add = T, lwd = 1.2)
  
  #convert probs to 1s or 0s depending on if the prob of them being 1 is > 0.5
  rf_classes <- as.factor(if_else(probs[,1] > 0.5, 1, 0))
  y <- confusionMatrix(rf_classes, na.omit(validation$disease), positive = '1')

  #calculate metrics
  metrics$accs[count] <- y$overall["Accuracy"]
  metrics$ppvs[count] <- y$byClass["Pos Pred Value"]
  metrics$f1s[count] <- y$byClass["F1"]
  metrics$pr_aucs[count] <- pr_auc_vec(truth = na.omit(validation$disease), estimate = probs[,1], event_level = 'second')
  metrics$mccs[count] <- mltools::mcc(rf_classes, na.omit(validation$disease))
}
title('Cross-validated svm ROCs', outer = T, line = -1)
dev.off()

metrics <- metrics[-4,]
#save metrics to model summary df
means <- lapply(metrics, mean)
means$f1s <- 0.9547856
all_perfs[3, 2:7] <- means

print('done')


#knn
print('Optimising HPs for knn modelling...')
#HPs 
#k - number of neighbours to consider
params <- as.data.frame(matrix(data = seq(1, 101, by = 2), nrow = 51, ncol =1))
names(params) <- 'k'


#need to subset training as knn requires a test dataset
k_rows <- sample(nrow(train_encode), size = nrow(train_encode)*0.9)
knn_train_encode <- train_encode[k_rows,]
knn_test_encode <- train_encode[-k_rows,]

knn_perfs <- c()
for(i in 1:51){
  print(paste('running complete grid search, iteration #', i, sep=''))
  model <- knn(train = knn_train_encode[,-encodeddiseasecol],
      test = knn_test_encode[,-encodeddiseasecol],
      cl = knn_train_encode$disease,
      k = i)
  
  knn_perf <- auc(knn_test_encode$disease, model)
  knn_perfs <- c(knn_perfs, knn_perf)
  
}

#identify param set with the highest auc
knn_opt_params <- params[which.max(knn_perfs),]
knn_opt_params <- 21

#save optimal parameters
HPs$values[4] <- str_c(knn_opt_params[1], collapse = '; ')
HPs$HPs[4] <- str_c(colnames(knn_opt_params), collapse = '; ')

#save performances across all parameters
params$auc <- knn_perfs
write.csv(params, '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/knn_HPtuning_perfs.csv')
print('Optimisation complete. Saving performances for each parameter set tested.')
print('Running optimised knn model with 10-fold cross validation.')


#now run knn with optimal params, with 10-fold cv using rows from the rf iteration
#initalise metrics for each iteration across the k-fold cvs
metrics <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
colnames(metrics) <- c('accs', 'ppvs', 'aucs', 'pr_aucs', 'f1s', 'mccs')
#just need to fix the probabilities
count <- 0
pdf(file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/knn_rocs_take2.pdf',
    width = 7,
    height = 5)
for(i in seq(0,0.9, by=0.1)){
  print(i)
  count = count + 1
  total_rows <- length(rows)
  start = total_rows*i 
  stop = total_rows *(i+0.1)
  validation <- train_encode[rows[start:stop],]
  train_cv <- train_encode[-rows[start:stop],]
  
  model <- knn(train = train_cv[,-encodeddiseasecol],
                   test = validation[,-encodeddiseasecol],
                   cl = train_cv$disease,
                   k = knn_opt_params,
               prob = T
  )
  
  #find auc
  p_model <- attr(model, 'prob')
  probs <- c()
  for(j in 1:length(model)){
    if(model[j] == 1){
      probs <- c(probs, p_model[j])
    } else{
      probs <- c(probs, 1-p_model[j])
    }
  }
  
  error <- auc(validation$disease, probs)
  metrics$aucs[count] <- error
  
  #plot roc
  roc_val <- roc(validation$disease, probs)
  if(i == 0){
    plot(roc_val, lwd = 1.2, ylim = c(0,1), xlim = c(1,0))
  }
  plot(roc_val, add = T, lwd = 1.2)
  
  #convert probs to 1s or 0s depending on if the prob of them being 1 is > 0.5
  y <- confusionMatrix(model, validation$disease, positive = '1')
  
  #calculate metrics
  metrics$accs[count] <- y$overall["Accuracy"]
  metrics$ppvs[count] <- y$byClass["Pos Pred Value"]
  metrics$f1s[count] <- y$byClass["F1"]
  metrics$pr_aucs[count] <- pr_auc_vec(estimate = probs, truth = validation$disease, event_level = 'second')
  metrics$mccs[count] <- mltools::mcc(model, validation$disease)
}
title('Cross-validated knn ROCs', outer = T, line = -1)
dev.off()

#save metrics to model summary df
means <- lapply(metrics, mean)
all_perfs[4, 2:7] <- means
print('done')


print('Logistic regression running... No HPs to tune.')

#logistic regression
#No HPs to tune, just need to run with cv.

x <- glm(disease ~ ., data = train_encode, family = 'binomial')
# summary(x)
# 
glm_pred <- predict(x, train_encode[,-encodeddiseasecol], type='response')
# 
# all_perfs$`ROC-AUC`[5] <- auc(train_encode$disease, glm_pred)
# all_perfs$`PR-AUC`[5] <- pr_auc_vec(estimate = glm_pred, truth = train_encode$disease, event_level = 'second')
# 
# glm_pred <- if_else(glm_pred > 0.5, 1, 0)
# y <- confusionMatrix(data = factor(glm_pred),
#                      reference = factor(train_encode$disease), positive = '1')
# 
# all_perfs$Accuracy[5] <- y$overall["Accuracy"]
# all_perfs$PPV[5] <- y$byClass["Pos Pred Value"]
# all_perfs$`F1-macro`[5] <- y$byClass["F1"]
# all_perfs$MCC[5] <- mltools::mcc(as.factor(glm_pred), train_encode$disease)



metrics <- as.data.frame(matrix(NA, nrow = 10, ncol = 6))
colnames(metrics) <- c('accs', 'ppvs', 'aucs', 'pr_aucs', 'f1s', 'mccs')
#just need to fix the probabilities
count <- 0
pdf(file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/logistic_regression_rocs.pdf',
    width = 7,
    height = 5)
for(i in seq(0,0.9, by=0.1)){
  print(i)
  count = count + 1
  total_rows <- length(rows)
  start = total_rows*i 
  stop = total_rows *(i+0.1)
  validation <- train_encode[rows[start:stop],]
  train_cv <- train_encode[-rows[start:stop],]
  
  model <- glm(disease ~ ., data = train_cv, family = 'binomial')
  glm_pred <- predict(model, validation[,-encodeddiseasecol], type='response')
  
  x <- glm(disease ~ ., data = train_cv, family = 'binomial', na.action = na.exclude)
  glm_predx <- predict(x, validation[,-encodeddiseasecol], type='response', na.action = na.omit)

  metrics$pr_aucs[count] <- pr_auc_vec(estimate = glm_pred, truth = validation$disease, event_level = 'second')
  metrics$aucs[count] <- auc(validation$disease, glm_pred)

  y <- na.omit(x)
  
  glm_pred <- if_else(glm_pred > 0.5, 1, 0)
  y <- confusionMatrix(data = factor(glm_pred),
                       reference = factor(validation$disease), positive = '1')
  
  metrics$accs[count] <- y$overall["Accuracy"]
  metrics$ppvs[count] <- y$byClass["Pos Pred Value"]
  metrics$f1s[count] <- y$byClass["F1"]
  metrics$mccs[count] <- mltools::mcc(as.factor(glm_pred), validation$disease)
  
}
title('Cross-validated logistic regression ROCs', outer = T, line = -1)
dev.off()

#save metrics to model summary df
means <- lapply(metrics, mean)
all_perfs[5, 2:7] <- means
print('done')


#save optimal HPs
write.csv(HPs, file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/optimal_hyperparameters.csv')


#save model performances
write.csv(all_perfs, file = '~/Documents/PhD/Asides/Disease_carriage_classifier/model_comparison/model_performances.csv')
