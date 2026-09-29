#run chisquared on top IGR hits

library(openxlsx)
tophits <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/IGR_modelling/unordered_fix/important_feats_annotated.xlsx')
dat_init <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/IGR_modelling/igr_annos_synthres-0.4_flank-1000_id-90_cov-70.xlsx')

x <- 2
locus <- paste('IGR_', substr(tophits$Locus[x], 4, nchar(tophits$Locus[x])), sep ='')
dat <- dat_init[,which(colnames(dat_init) %in% c('disease', locus))]

removables <- c()
for(i in 1:nrow(dat)){
  if(is.na(dat$disease[i])){
    removables <- c(removables, i)
  } else if(dat$disease[i] == 'carrier'){
    dat$disease[i] <- 'carriage'
  } else if(dat$disease[i] == 'other'){
    removables <- c(removables, i)
  } else{
    dat$disease[i] <- 'disease'
  }
}
dat <- dat[-removables,]
dat <- dat[which(!is.na(dat[[locus]])),]

dat$disease <- as.factor(dat$disease)
dat[[locus]] <- as.factor(dat[[locus]])

taballeles <- table(dat[[locus]])
filters <- names(taballeles[which(taballeles > 100)])
#Only most common alleles
dat_common <- dat[which(dat[[locus]] %in% filters),]


dattab <- table(dat_common$disease, dat_common[[locus]])
dattab

library(ggplot2)
library(ggokabeito)
ggplot(data = dat_common, aes(x=.data[[locus]], fill=disease)) +
  geom_bar() +
  theme_classic(base_size = 16) +
  ggokabeito::scale_fill_okabe_ito() +
  xlab('aniA IGR allele') +
  ylab('N isolates')

