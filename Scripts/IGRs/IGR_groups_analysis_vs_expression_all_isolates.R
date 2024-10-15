#Script to investigate how LFC is effected by igr variation.
###Error - need to account for min q value potentially being different!
upvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/upvar_igrs.xlsx')  #334
downvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/downvar_igrs.xlsx') #373
var <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/var_igrs.xlsx') #1067
novar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/novar_igrs.xlsx') #717
lowvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/lowvar_igrs.xlsx') #333
highvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/highvar_igrs.xlsx') #27 loci



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

#check which entries in rna have a qval of 1.
rna1 <- rna[which(as.numeric(rna$X13) == 1),]
rna1 <- rna1[which(as.numeric(rna1$`log2(EXP/Min)`) >= 1),]
