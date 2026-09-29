library(pathopred)
library(openxlsx)
library(magrittr)
library(gggenomes)

#devtools::install('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/Bioinformatics/pathopred')

main <- function(data, outdir){
  
  # data <- filter_data(data, out = outdir)
  # 
  # data <- select_feats(data, maxRuns = 200, out = outdir)
  
  #1 - remove isolates lacking disease data
  keepers <- c('carrier', 'invasive (unspecified/other)', 'meningitis', 'meningitis and septicaemia', 'septicaemia', 'conjunctivitis')
  data <- data[which(data$disease %in% keepers),]
  
  #2 - Code disease as 0, 1
  data$disease[which(data$disease != 'carrier')] <- 1
  data$disease[which(data$disease == 'carrier')] <- 0
  
  #3 - Code cells with multiple alleles or no entry as 0
  removables <- c()
  for(i in 2:ncol(data)){
    column <- data[,i]
    column[which(is.na(column))] <- 0
    column[which(stringr::str_detect(column, ';'))] <- 0
    data[,i] <- column
    
    #4 - Remove genes absent in more than 20% of isolates
    if(sum(column == '0') > (0.2*length(column))){
      removables <- c(removables, i)
    }
  }
  
  if(length(removables) > 0){
    data <- data[,-removables]
  }

  #5 Remove isolates with no data in more than 20% of genes
  removables <- c()
  for(i in 1:nrow(data)){
    row <- data[i,-1]
    if(sum(row == '0') > 0.5*length(row)){
      removables <- c(removables, i)
    }
  }
  
  if(length(removables) > 0){
    data <- data[,-removables]
  }
  
  #6 - remove genes where the 2 most common alleles account for more than 80% of all isolates
  removables <- c()
  for(i in 2:ncol(data)){
    column <- data[,i]
    tab <- as.data.frame(table(column))
    tab <- tab[order(tab$Freq, decreasing = T),]
    
    if(nrow(tab) < 2){
      removables <- c(removables, i)
    }
    
    else if(tab$Freq[1] + tab$Freq[2] > (length(column)*0.8)){
      removables <- c(removables, i)
    } 
    #7 - convert rare alleles to the code 999.
    else{
      minors <- tab$column[which(tab$Freq < (length(column)*0.01))]
      column[which(as.numeric(column) %in% minors)] <- '999'
      data[,i] <- column
      
      #8 - Remove genes with many rare alleles
      if(sum(column == '999') > 0.2*length(column)){
        removables <- c(removables, i)
      }
    }
  }
  
  if(length(removables) > 0){
    data <- data[,-removables]
  }
  
  data <- split_data(df = data, balance_test_data = T, out = outdir)
  train <- data[[1]]
  test <- data[[2]]
  ds_fac_train <- data[[3]]
  
  #univariate feature selection
  fit_chi <- function(x, y){
    chi_res <- apply(x, 2, function(f) {chisq.test(f, y, simulate.p.value = F, rescale = T)})
  }
  
  chi_res <- fit_chi(x = train[,-1],
                     y = train[,1])
  
  
  df <- as.data.frame(matrix(data = NA, nrow = length(chi_res), ncol = 2))
  colnames(df) <- c('Locus', 'p')
  df$Locus <- names(chi_res)
  for(i in 1:length(chi_res)){
    truep <- pchisq(q = chi_res[[i]]$statistic,
                    df = chi_res[[i]]$parameter,
                    lower.tail = F)
    df$p[i] <- truep
  }
  
  df <- df[order(df$p),]
  tophits <- c('disease', df$Locus[1:100])
  # 
  # fit_anova <- function(x, y) {
  #   anova_res <- apply(x, 2, function(f) {caret::anovaScores(f, y)})
  #   return(anova_res)
  # }
  # aov_res <- fit_anova(x = train[,-1],
  #                      y = train[,1])
  # tophits <- c('disease', names(aov_res)[order(aov_res)[1:100]])

  train <- train[, which(colnames(train) %in% tophits)]
  test <- test[, which(colnames(test) %in% tophits)]
  
  train$disease <- as.factor(train$disease)
  test$disease <- as.factor(test$disease)
  
  ps <- ncol(train)-1
  
  n_trees_min <- 10*ps
  n_trees_max <- 30*ps
  n_trees_inc <- round((n_trees_max - n_trees_min)/10)
  mtry_min <- 1
  mtry_max <- round(ps/2)
  mtry_inc <- round((mtry_max - mtry_min)/10)
  nodesize_min <- 1
  nodesize_max <- 11
  nodesize_inc <- 2
  searchsize <- 50
  
  n_trees_min <- 1000
  n_trees_max <- 1100
  n_trees_inc <- 10
  
  params <- set_params(n_trees_min = n_trees_min, n_trees_max = n_trees_max, n_trees_inc = n_trees_inc,
                       mtry_min = mtry_min, mtry_max = mtry_max, mtry_inc = mtry_inc,
                       nodesize_min = nodesize_min, nodesize_max = nodesize_max, nodesize_inc = nodesize_inc,
                       searchsize = searchsize)
  
  print('optimising model')
  params <- opt_model(train, params, ds_fac_train=ds_fac_train)
  params <- params[order(params$f1, decreasing = T),]
  print(params$f1[1])
  opt_params <- params[1,]
  
  model <- run_model(train, test, opt_params, importance = 'permutation', out=outdir, ds_fac_train=ds_fac_train)
  
  annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
  feats <- get_imp_feats(model = model, annotate = F, annotations = annotations, out = outdir)
  
  return(model)
}

set.seed(202)
####################################################################################
#load data
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data_init <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/dataset_isolates_under100_contigs_all_loci.xlsx')[,-c(1,2)]
data <- data_init
data <- data[,which(!(colnames(data) %in% removables$Locus))]
data <- data[,which(colnames(data) != 'capsule_group')]

out <- '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/D_All/chisqselection/testing_times/'

#Run model on full dataset
start_time <- Sys.time()
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)
end_time <- Sys.time()
tt <- round(end_time - start_time)
cat('\nCompleted in', tt, units(tt),'\n')

#Run on same dataset with shuffled disease information
shuffled_df <- data
states <- shuffled_df$disease
shuffled <- sample(states, length(states))
shuffled_df$disease <- shuffled
#shuffled_df$disease <- as.factor(shuffled_df$disease)

out <- '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/shuffled_disease'
model <- main(shuffled_df, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


#Run on IGR data
data_init <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/IGR_modelling/igr_annos_synthres-0.4_flank-1000_id-90_cov-70.xlsx')[,-c(1)]
out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/IGR_modelling/unordered_fix/chisq'
model <- main(data_init, out)

treeInfo(model)
model$forest$covariate.levels

igrs <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/IGR_modelling/unordered_fix/chisq/important_feats.xlsx')
annos <- read_gff3('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/IGR_modelling/igr_annos.gff3')

igrs$Up_gene <- NA
igrs$Down_gene <- NA
for(i in 1:nrow(igrs)){
  locusname <- substr(igrs$Locus[i], 4, nchar(igrs$Locus[i]))
  locusname <- paste('IGR', locusname, sep ='')
  igrs$Up_gene[i] <- annos$upstream_gene[which(annos$feat_id == locusname)]
  igrs$Down_gene[i] <- annos$downstream_gene[which(annos$feat_id == locusname)]
}

write.xlsx(igrs, '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/IGR_modelling/unordered_fix/chisq/important_feats_annotated.xlsx')




#How many missing datapoints?
igrannotations <- igrs[,c(2,4)]
colnames(igrannotations) <- c('Locus', 'Annotation')
for(i in 1:nrow(igrannotations)){
  print(i)
  #if(i == 30) break
  
  try(igrannotations$Annotation[i] <- annotations$Annotation[which(annotations$Locus == igrannotations$Annotation[i])])
}

get_imp_feats(model=model, annotate=T, annotations = igrannotations, out = out)

################################################
#run on GWAS datasets
####c41_44 replication study####
#load data
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/GWAS_validationset2.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]
data <- data[,which(colnames(data) != 'capsule_group')]

out <- '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


######Czechia######
data <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Czechia/GWAS_validationset1.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]
data <- data[,which(colnames(data) != 'capsule_group')]

out <- '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Czechia/unordered_fix'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


#####Sweden######
data <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Sweden/GWAS_Sweden_isolates.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]
data <- data[,which(colnames(data) != 'capsule_group')]

out <- '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Sweden/unordered_fix'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)




#####Farzand_MenW######
data <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Farzand_MenW/GWAS_Farzand_MenW_isolates.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]
data <- data[,which(colnames(data) != 'capsule_group')]

out <- '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Farzand_MenW/unordered_fix'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)





######Farzand_MenY#####
data <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Farzand_MenY/GWAS_Farzand_MenY_isolates.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]
data <- data[,which(colnames(data) != 'capsule_group')]

out <- '~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/Farzand_MenY/unordered_fix'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)



################################################
#Run on specific lineages
#MenB - cc41/44
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/B/cc41_44/MenB_cc41_44_dataset.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]
out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/B/cc41_44'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)

##########
##########
#MenW:cc11
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/W/cc11/MenW_cc11_dataset.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]

out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/W/cc11'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)

##########
##########
#MenY:cc23
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/Y/cc23/MenY_cc23_dataset.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]

out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/Y/cc23'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


##########
##########
#MenA:cc5
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/A/cc5/MenA_cc5_dataset.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]

out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/A/cc5'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


##########
##########
#MenC:cc11
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/C/cc11/MenC_cc11_dataset.xlsx')[,-c(1,2)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]

out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/C/cc11'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


##########
##########
#MenY:cc23 UK vs USA
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/Y/cc23_UK_USA/MenY_cc23_UK_USA.xlsx')[,-c(1,2,4)]
data <- data[,which(!(colnames(data) %in% removables$Locus))]

#convert country to disease type, where UK is invasive (1) and USA is carriage (0)
data$country[which(data$country == 'UK')] <- 'invasive (unspecified/other)'
data$country[which(data$country == 'USA')] <- 'carrier'
colnames(data)[1] <- 'disease'

out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/Y/cc23_UK_USA'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


##########
##########
#Strep pneumo 19A carriage vs disease
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/strep_pneumo/imp_feats_annotations.xlsx')
#removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/strep_pneumo/strep_19A.xlsx')[,-c(1,2)]
#data <- data[,which(!(colnames(data) %in% removables$Locus))]

#convert diagnosis to disease
colnames(data)[1] <- 'disease'
diseases <- c('bacteraemia', 'meningitis', 'pneumonia')
data$disease[which(data$disease %in% diseases)] <- 'invasive (unspecified/other)'
data$disease[which(data$disease == 'carriage')] <- 'carrier'


out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/strep_pneumo/19A'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)


#Strep pneumo major STs carriage vs disease
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/strep_pneumo/imp_feats_annotations.xlsx')
#removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')

data <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/strep_pneumo/major_STs.xlsx')[,-c(1,2)]
#data <- data[,which(!(colnames(data) %in% removables$Locus))]

#convert diagnosis to disease
colnames(data)[1] <- 'disease'
diseases <- c('bacteraemia', 'meningitis', 'pneumonia')
data$disease[which(data$disease %in% diseases)] <- 'invasive (unspecified/other)'
data$disease[which(data$disease == 'carriage')] <- 'carrier'


out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/strep_pneumo/major_STs'

#Run model on full dataset
model <- main(data, out)
get_imp_feats(model=model, annotate=T, annotations = annotations, out = out)



################
################
#Do model on cc41/44 GWAS IGRs
#load data
annotations <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/imp_feats_annotations.xlsx')
removables <- read.csv('~/Documents/PhD/Asides/Disease_carriage_classifier/columnnames.csv')


#Run on IGR data
data_init <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/Annotation/out/igr_annos_synthres-0.4_flank-500_id-90_cov-70.xlsx')[,-c(1)]
out <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs'
model <- main(data_init, out)

treeInfo(model)
model$forest$covariate.levels

igrs <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/important_feats.xlsx')
annos <- read_gff3('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/Annotation/out/igr_annos.gff3')

igrs$Up_gene <- NA
igrs$Down_gene <- NA
for(i in 1:nrow(igrs)){
  locusname <- substr(igrs$Locus[i], 4, nchar(igrs$Locus[i]))
  locusname <- paste('IGR', locusname, sep ='')
  igrs$Up_gene[i] <- annos$upstream_gene[which(annos$feat_id == locusname)]
  igrs$Down_gene[i] <- annos$downstream_gene[which(annos$feat_id == locusname)]
}

write.xlsx(igrs, '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/GWAS_validation/c41_44_replication/unordered_fix/IGRs/important_feats_annotated.xlsx')


#How many missing datapoints?
igrannotations <- igrs[,c(2,4)]
colnames(igrannotations) <- c('Locus', 'Annotation')
for(i in 1:nrow(igrannotations)){
  print(i)
  if(i == 30) break
  
  try(igrannotations$Annotation[i] <- annotations$Annotation[which(annotations$Locus == igrannotations$Annotation[i])])
}

get_imp_feats(model=model, annotate=T, annotations = igrannotations, out = out)



