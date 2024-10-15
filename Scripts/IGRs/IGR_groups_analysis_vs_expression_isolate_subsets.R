igrs <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs.xlsx')[, -c(4,5,7, 8, 12, 13, 15, 16)]

upvar <- c()
highvar <- c()
downvar <- c()
lowvar <- c()
var <- c()
novar <- c()

for(i in 1:nrow(igrs)){
  print(i)
  igr_up <- as.vector(igrs[i,4:7])
  igr_down <- as.vector(igrs[i, 8:11])
  
  igrs$`Alleles-Down`[i] <- length(unique(as.numeric(na.omit(as.numeric(igrs[i, 4:7])))))
  igrs$`Alleles-Up`[i] <- length(unique(as.numeric(na.omit(as.numeric(igrs[i, 8:11])))))
  igrs$Combined[i] <- igrs$`Alleles-Down`[i] + igrs$`Alleles-Up`[i]
  
  if(length(unique(igr_up)) == 1 & length(unique(igr_down)) == 1){
    print('yes1')
    novar <- c(novar, i)
  } else if(length(unique(igr_up)) == 1 & length(unique(igr_down)) > 1){
    print('yes2')
    downvar <- c(downvar, i)
    var <- c(var, i)
  } else if(length(unique(igr_up)) > 1 & length(unique(igr_down)) == 1){
    print('yes3')
    upvar <- c(upvar, i)
    var <- c(var, i)
  } else if(igrs$Combined[i] < 5){
    print('yes4')
    lowvar <- c(lowvar, i)
    var <- c(var, i)
  } else if(igrs$Combined[i] >= 5){
    print('yes5')
    highvar <- c(highvar, i)
    var <- c(var, i)
  }
}


upvar <- igrs[upvar,]
highvar <- igrs[highvar,]
downvar <- igrs[downvar,]
lowvar <- igrs[lowvar,]
var <- igrs[var,]
novar <- igrs[novar,]

nrow(upvar) + nrow(highvar) + nrow(downvar) + nrow(lowvar) == nrow(var)
nrow(upvar) + nrow(highvar) + nrow(downvar) + nrow(lowvar) + nrow(novar) == nrow(igrs)

          #excl 28269     #samehost only
nrow(upvar) #323.     #269
nrow(highvar) #20     #39
nrow(downvar) #291    #245
nrow(lowvar) #222    #120
nrow(var) #856.      #673
nrow(novar) #928    #1111


#load rna data
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx', sheet = 8)[-1, c(8,9,13:22, 25:32)]

#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

#which loci in rna are up,downvar etc?
vargroup_plotdata <- as.data.frame(matrix(data = NA, nrow = nrow(rna), ncol = 3))
var_plotdata <- as.data.frame(matrix(data = NA, nrow = nrow(rna), ncol = 3))
colnames(vargroup_plotdata) <- c('Locus', 'log2fold.change', 'VarGroup')
colnames(var_plotdata) <- c('Locus', 'log2fold.change', 'Var')


#Recaclculate LFC for just the within host isolates.
#28262, 53930, 53948, 53951

rna_samehost<- rna[,c(2,3 ,8, 9, 11,12)]
colnames(rna_samehost) <- c('1717.genes', 'X13', '28262', '53930','53948', '53951')
rna_samehost$`log2(EXP/Min)` <- NA
get_LFC <- function(df){
  for(row in 1:nrow(df)){
    min_exp <- min(as.numeric(df[row,3:6]))
    max_exp <- max(as.numeric(df[row,3:6]))
    change <- max_exp/min_exp
    log_change <- log2(change)
    df$`log2(EXP/Min)`[row] <- log_change
    
    if(min_exp == 0){
      df$`log2(EXP/Min)`[row] <- NA
    }
  }
  return(df)
}

rna <- get_LFC(rna_samehost)

for(i in 1:nrow(rna)){
  locus <- rna$`1717.genes`[i]
  if(locus %in% upvar$Locus){
    print('1')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'UpVar'
    
  }else if(locus %in% downvar$Locus){
    print('2')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'DownVar'
  } else if(locus %in% lowvar$Locus){
    print('3')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'LowVar'
  } else if(locus %in% highvar$Locus){
    print('4')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'HighVar'
  } 
  
  if(locus %in% var$Locus){
    print('5')
    var_plotdata$Locus[i] <- locus
    var_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    var_plotdata$Var[i] <- 'Var'
  }else if(locus %in% novar$Locus){
    print('6')
    var_plotdata$Locus[i] <- locus
    var_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    var_plotdata$Var[i] <- 'NoVar'
    
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'NoVar'
  }
}

#remove rows with NA values for lfc
removable <- c()
for(i in 1:nrow(vargroup_plotdata)){
  if(is.na(vargroup_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}

vargroup_plotdata <- vargroup_plotdata[-removable,]


removable <- c()
for(i in 1:nrow(var_plotdata)){
  if(is.na(var_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}
var_plotdata <- var_plotdata[-removable,]

#plot by vargroup
var_plotdata$log2fold.change %<>% as.numeric()
vargroup_plotdata$log2fold.change %<>% as.numeric()
wilcox.test((var_plotdata$log2fold.change[var_plotdata$Var == 'Var']), var_plotdata$log2fold.change[var_plotdata$Var == 'NoVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'UpVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'LowVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'HighVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'UpVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'UpVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'LowVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'UpVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'HighVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'LowVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'HighVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'LowVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'HighVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])


mycomparisons <- list(c('DownVar', 'UpVar'), c('DownVar', 'NoVar'), c('DownVar', 'LowVar'), c('DownVar', 'HighVar'), c('UpVar', 'NoVar'), c('UpVar', 'LowVar'), c('UpVar', 'HighVar'), c('LowVar', 'HighVar'), c('LowVar', 'NoVar'), c('HighVar', 'NoVar'))
#Make plot
ggplot(na.omit(vargroup_plotdata), aes(x = VarGroup, y = log2fold.change, color = VarGroup)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom='point', shape = 0, size = 1,  color = 'black') +
  geom_jitter(shape = 16, size = 1.2, position=position_jitter(0.2)) +
  geom_signif(comparisons = mycomparisons, map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05), step_increase = 0.1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size = 15),
        legend.text = element_text(size=10), legend.title = element_text(size = 15)) +
  labs(y = 'log2-fold change')


#plot by var
mycomparisons <- list(c('NoVar', 'Var'))

ggplot(var_plotdata[!is.na(var_plotdata$Var),], aes(x = Var, y = log2fold.change, color = Var, na.rm = TRUE)) +
  geom_boxplot(na.rm = T) +
  stat_summary(fun.y=mean, geom='point', shape = 5, size = 2,  color = 'black') +
  geom_jitter(shape = 16, position=position_jitter(0.2)) +
  geom_signif(comparisons = list(c("Var", "NoVar")), map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size = 20),
        legend.text = element_text(size=15), legend.title = element_text(size = 15))


#find sig loci
sigloci <- c()
for(i in 1:nrow(rna)){
  
  if(is.na(rna$`log2(EXP/Min)`[i]) | is.na(rna$X13[i])){
    next
  }
  
  if(as.numeric(rna$`log2(EXP/Min)`[i]) >= log2(1.5/1) && as.numeric(rna$X13[i]) < 0.05){
    sigloci <- c(sigloci, rna$`1717.genes`[i])
  }
}


sigvar_plotdata <- var_plotdata[which(var_plotdata$Locus %in% sigloci), ]
sigvargroup_plotdata <- vargroup_plotdata[which(vargroup_plotdata$Locus %in% sigloci), ]

sigvar_plotdata <- sigvar_plotdata[which(sigvar_plotdata$log2fold.change >= log2(1.5/1)), ]
sigvargroup_plotdata <- sigvargroup_plotdata[which(sigvargroup_plotdata$log2fold.change >= log2(1.5/1)), ]

removable <- c()
for(i in 1:nrow(sigvargroup_plotdata)){
  if(is.na(sigvargroup_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}

sigvargroup_plotdata <- sigvargroup_plotdata[-removable,]


removable <- c()
for(i in 1:nrow(sigvar_plotdata)){
  if(is.na(sigvar_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}
sigvar_plotdata <- sigvar_plotdata[-removable,]





#sigplot data
wilcox.test((sigvar_plotdata$log2fold.change[sigvar_plotdata$Var == 'Var']), sigvar_plotdata$log2fold.change[sigvar_plotdata$Var == 'NoVar'])
wilcox.test((sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$Var == 'DownVar']), sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$VarGroup == 'UpVar'])
wilcox.test((sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$Var == 'DownVar']), sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$VarGroup == 'NoVar'])
wilcox.test((sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$Var == 'UpVar']), sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$VarGroup == 'NoVar'])




mycomparisons <- list(c('DownVar', 'UpVar'), c('DownVar', 'NoVar'), c('DownVar', 'LowVar'), c('DownVar', 'HighVar'), c('UpVar', 'NoVar'), c('UpVar', 'LowVar'), c('UpVar', 'HighVar'), c('LowVar', 'HighVar'), c('LowVar', 'NoVar'), c('HighVar', 'NoVar'))
#Make plot
ggplot(na.omit(sigvargroup_plotdata), aes(x = VarGroup, y = log2fold.change, color = VarGroup)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom='point', shape = 0, size = 1,  color = 'black') +
  geom_jitter(shape = 16, size = 1.2, position=position_jitter(0.2)) +
  geom_signif(comparisons = mycomparisons, map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05), step_increase = 0.1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size = 15),
        legend.text = element_text(size=10), legend.title = element_text(size = 15)) +
  labs(y = 'log2-fold change')


#plot by var
mycomparisons <- list(c('NoVar', 'Var'))

ggplot(sigvar_plotdata[!is.na(sigvar_plotdata$Var),], aes(x = Var, y = log2fold.change, color = Var, na.rm = TRUE)) +
  geom_boxplot(na.rm = T) +
  stat_summary(fun.y=mean, geom='point', shape = 5, size = 2,  color = 'black') +
  geom_jitter(shape = 16, position=position_jitter(0.2)) +
  geom_signif(comparisons = list(c("Var", "NoVar")), map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size = 20),
        legend.text = element_text(size=15), legend.title = element_text(size = 15))


#plot sig loci only
mean(sigvar_plotdata$log2fold.change[which(sigvar_plotdata$Var == 'Var')])
mean(sigvar_plotdata$log2fold.change[which(sigvar_plotdata$Var == 'NoVar')])
mean(vargroup_plotdata$log2fold.change[which(sigvar_plotdata$Var == 'NoVar')])





#====================
#====================
#====================
###Repeated, excluding 28269

igrs <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs.xlsx')[,-c(7, 15)]

upvar <- c()
highvar <- c()
downvar <- c()
lowvar <- c()
var <- c()
novar <- c()

for(i in 1:nrow(igrs)){
  print(i)
  igr_up <- as.vector(igrs[i,4:10])
  igr_down <- as.vector(igrs[i, 11:17])
  
  igrs$`Alleles-Down`[i] <- length(unique(as.numeric(na.omit(as.numeric(igrs[i, 4:10])))))
  igrs$`Alleles-Up`[i] <- length(unique(as.numeric(na.omit(as.numeric(igrs[i, 11:17])))))
  igrs$Combined[i] <- igrs$`Alleles-Down`[i] + igrs$`Alleles-Up`[i]
  
  if(length(unique(igr_up)) == 1 & length(unique(igr_down)) == 1){
    print('yes1')
    novar <- c(novar, i)
  } else if(length(unique(igr_up)) == 1 & length(unique(igr_down)) > 1){
    print('yes2')
    downvar <- c(downvar, i)
    var <- c(var, i)
  } else if(length(unique(igr_up)) > 1 & length(unique(igr_down)) == 1){
    print('yes3')
    upvar <- c(upvar, i)
    var <- c(var, i)
  } else if(igrs$Combined[i] < 9){
    print('yes4')
    lowvar <- c(lowvar, i)
    var <- c(var, i)
  } else if(igrs$Combined[i] >= 9){
    print('yes5')
    highvar <- c(highvar, i)
    var <- c(var, i)
  }
}


upvar <- igrs[upvar,]
highvar <- igrs[highvar,]
downvar <- igrs[downvar,]
lowvar <- igrs[lowvar,]
var <- igrs[var,]
novar <- igrs[novar,]

nrow(upvar) + nrow(highvar) + nrow(downvar) + nrow(lowvar) == nrow(var)
nrow(upvar) + nrow(highvar) + nrow(downvar) + nrow(lowvar) + nrow(novar) == nrow(igrs)

#excl 28269     #samehost only
nrow(upvar) #323.     #269
nrow(highvar) #20    #39
nrow(downvar) #291    #245
nrow(lowvar) #222    #120
nrow(var) #856.      #673
nrow(novar) #928    #1111


#load rna data
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx', sheet = 8)[-1, c(8,9,13:22, 25:32)]

#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

#which loci in rna are up,downvar etc?
vargroup_plotdata <- as.data.frame(matrix(data = NA, nrow = nrow(rna), ncol = 3))
var_plotdata <- as.data.frame(matrix(data = NA, nrow = nrow(rna), ncol = 3))
colnames(vargroup_plotdata) <- c('Locus', 'log2fold.change', 'VarGroup')
colnames(var_plotdata) <- c('Locus', 'log2fold.change', 'Var')


#recalculate LFC(MAX/MIN) to remove 28269.
#find minimum expressed without 28269.
#Find maximum expressed without 28269
#take LFC base 2 of max/min
rna_7 <- rna[,c(2,3 , 5,6, 8:12)]
colnames(rna_7) <- c('1717.genes', 'X13','27509','27553', '28262', '53930', '28287', '53948', '53951')
rna_7
rna_7$`log2(EXP/Min)` <- NA
get_LFC <- function(df){
  for(row in 1:nrow(df)){
    
    min_exp <- min(as.numeric(df[row,3:9]))
    max_exp <- max(as.numeric(df[row,3:9]))
    change <- max_exp/min_exp
    log_change <- log2(change)
    df$`log2(EXP/Min)`[row] <- log_change
    
    if(min_exp == 0){
      df$`log2(EXP/Min)`[row] <- NA
    }
  }
  return(df)
}

rna <- get_LFC(rna_7)

for(i in 1:nrow(rna)){
  locus <- rna$`1717.genes`[i]
  if(locus %in% upvar$Locus){
    print('1')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'UpVar'
    
  }else if(locus %in% downvar$Locus){
    print('2')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'DownVar'
  } else if(locus %in% lowvar$Locus){
    print('3')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'LowVar'
  } else if(locus %in% highvar$Locus){
    print('4')
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'HighVar'
  } 
  
  if(locus %in% var$Locus){
    print('5')
    var_plotdata$Locus[i] <- locus
    var_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    var_plotdata$Var[i] <- 'Var'
  }else if(locus %in% novar$Locus){
    print('6')
    var_plotdata$Locus[i] <- locus
    var_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    var_plotdata$Var[i] <- 'NoVar'
    
    vargroup_plotdata$Locus[i] <- locus
    vargroup_plotdata$log2fold.change[i] <- rna$`log2(EXP/Min)`[i]
    vargroup_plotdata$VarGroup[i] <- 'NoVar'
  }
}

#remove rows with NA values for lfc
removable <- c()
for(i in 1:nrow(vargroup_plotdata)){
  if(is.na(vargroup_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}

vargroup_plotdata <- vargroup_plotdata[-removable,]


removable <- c()
for(i in 1:nrow(var_plotdata)){
  if(is.na(var_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}
var_plotdata <- var_plotdata[-removable,]

#plot by vargroup
var_plotdata$log2fold.change %<>% as.numeric()
vargroup_plotdata$log2fold.change %<>% as.numeric()
wilcox.test((var_plotdata$log2fold.change[var_plotdata$Var == 'Var']), var_plotdata$log2fold.change[var_plotdata$Var == 'NoVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'UpVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'LowVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'DownVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'HighVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'UpVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'UpVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'LowVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'UpVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'HighVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'LowVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'HighVar'])
wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'LowVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])

wilcox.test((vargroup_plotdata$log2fold.change[vargroup_plotdata$Var == 'HighVar']), vargroup_plotdata$log2fold.change[vargroup_plotdata$VarGroup == 'NoVar'])


mycomparisons <- list(c('DownVar', 'UpVar'), c('DownVar', 'NoVar'), c('DownVar', 'LowVar'), c('DownVar', 'HighVar'), c('UpVar', 'NoVar'), c('UpVar', 'LowVar'), c('UpVar', 'HighVar'), c('LowVar', 'HighVar'), c('LowVar', 'NoVar'), c('HighVar', 'NoVar'))
#Make plot
ggplot(na.omit(vargroup_plotdata), aes(x = VarGroup, y = log2fold.change, color = VarGroup)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom='point', shape = 0, size = 1,  color = 'black') +
  geom_jitter(shape = 16, size = 1.2, position=position_jitter(0.2)) +
  geom_signif(comparisons = mycomparisons, map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05), step_increase = 0.1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size = 15),
        legend.text = element_text(size=10), legend.title = element_text(size = 15)) +
  labs(y = 'log2-fold change')


#plot by var
mycomparisons <- list(c('NoVar', 'Var'))

ggplot(var_plotdata[!is.na(var_plotdata$Var),], aes(x = Var, y = log2fold.change, color = Var, na.rm = TRUE)) +
  geom_boxplot(na.rm = T) +
  stat_summary(fun.y=mean, geom='point', shape = 5, size = 2,  color = 'black') +
  geom_jitter(shape = 16, position=position_jitter(0.2)) +
  geom_signif(comparisons = list(c("Var", "NoVar")), map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size = 20),
        legend.text = element_text(size=15), legend.title = element_text(size = 15))


#find sig loci
sigloci <- c()
for(i in 1:nrow(rna)){
  
  if(is.na(rna$`log2(EXP/Min)`[i]) | is.na(rna$X13[i])){
    next
  }
  
  if(as.numeric(rna$`log2(EXP/Min)`[i]) >= log2(1.5/1) && as.numeric(rna$X13[i]) < 0.05){
    sigloci <- c(sigloci, rna$`1717.genes`[i])
  }
}


sigvar_plotdata <- var_plotdata[which(var_plotdata$Locus %in% sigloci), ]
sigvargroup_plotdata <- vargroup_plotdata[which(vargroup_plotdata$Locus %in% sigloci), ]
sigvar_plotdata <- sigvar_plotdata[which(sigvar_plotdata$log2fold.change >= log2(1.5/1)), ]
sigvargroup_plotdata <- sigvargroup_plotdata[which(sigvargroup_plotdata$log2fold.change >= log2(1.5/1)), ]

removable <- c()
for(i in 1:nrow(sigvargroup_plotdata)){
  if(is.na(sigvargroup_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}

sigvargroup_plotdata <- sigvargroup_plotdata[-removable,]


removable <- c()
for(i in 1:nrow(sigvar_plotdata)){
  if(is.na(sigvar_plotdata$log2fold.change[i])){
    removable <- c(removable, i)
  }
}
sigvar_plotdata <- sigvar_plotdata[-removable,]





#sigplot data
wilcox.test((sigvar_plotdata$log2fold.change[sigvar_plotdata$Var == 'Var']), sigvar_plotdata$log2fold.change[sigvar_plotdata$Var == 'NoVar'])
wilcox.test((sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$Var == 'DownVar']), sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$VarGroup == 'UpVar'])
wilcox.test((sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$Var == 'DownVar']), sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$VarGroup == 'NoVar'])
wilcox.test((sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$Var == 'UpVar']), sigvargroup_plotdata$log2fold.change[sigvargroup_plotdata$VarGroup == 'NoVar'])




mycomparisons <- list(c('DownVar', 'UpVar'), c('DownVar', 'NoVar'), c('DownVar', 'LowVar'), c('DownVar', 'HighVar'), c('UpVar', 'NoVar'), c('UpVar', 'LowVar'), c('UpVar', 'HighVar'), c('LowVar', 'HighVar'), c('LowVar', 'NoVar'), c('HighVar', 'NoVar'))
#Make plot
ggplot(na.omit(sigvargroup_plotdata), aes(x = VarGroup, y = log2fold.change, color = VarGroup)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom='point', shape = 0, size = 1,  color = 'black') +
  geom_jitter(shape = 16, size = 1.2, position=position_jitter(0.2)) +
  geom_signif(comparisons = mycomparisons, map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05), step_increase = 0.1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size = 15),
        legend.text = element_text(size=10), legend.title = element_text(size = 15)) +
  labs(y = 'log2-fold change')


#plot by var
mycomparisons <- list(c('NoVar', 'Var'))

ggplot(sigvar_plotdata[!is.na(sigvar_plotdata$Var),], aes(x = Var, y = log2fold.change, color = Var, na.rm = TRUE)) +
  geom_boxplot(na.rm = T) +
  stat_summary(fun.y=mean, geom='point', shape = 5, size = 2,  color = 'black') +
  geom_jitter(shape = 16, position=position_jitter(0.2)) +
  geom_signif(comparisons = list(c("Var", "NoVar")), map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size = 20),
        legend.text = element_text(size=15), legend.title = element_text(size = 15))


#plot sig loci only
mean(sigvar_plotdata$log2fold.change[which(sigvar_plotdata$Var == 'Var')])
mean(sigvar_plotdata$log2fold.change[which(sigvar_plotdata$Var == 'NoVar')])
mean(vargroup_plotdata$log2fold.change[which(sigvar_plotdata$Var == 'NoVar')])







