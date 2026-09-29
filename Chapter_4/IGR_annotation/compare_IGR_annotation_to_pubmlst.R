#script to evaluate the outputs of the IGR annotation pipeline when run on CDS regions
#compares those annotations to ones output by PubMLST 
#Either genome comparator output or just exported from the existing data on PubMLST

library(openxlsx)
library(ggplot2)
library(ggokabeito)
library(magrittr)
library(stringr)
library(ggpubr)
library(dplyr)

#load datasets of 'ground truths'
pubmlstdat <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/datasetexport_synthres_compare.xlsx')
genomecomp <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/genomecomparator_synthres_compare.xlsx')


#how many isolates correctly have presence or absence of a gene matching between datasets?
find_match <- function(dat, pubmlstdat){
  counter <- 0
  matches <- 0
  missing <- 0
  present <- 0
  
  for(i in 4:ncol(dat)){
    columnname <- colnames(dat)[i]
    if(columnname == 'capsule') next
    column <- dat[,columnname]
    if(!(columnname %in% colnames(pubmlstdat))) next
    refcolumn <- pubmlstdat[,columnname]
    
    counter <- counter + 1
    if(length(column) != length(refcolumn)){
      print(i)
      print(columnname)
      break
    }
    
    for(j in 1:length(column)){
      if(is.na(column[j]) & is.na(refcolumn[j])){
        matches <- matches + 1
      } else if(!is.na(column[j]) & !is.na(refcolumn[j])){
        matches <- matches + 1
      } else if(is.na(column[j]) & !is.na(refcolumn[j])){
        missing <- missing + 1
      } else if(!is.na(column[j]) & is.na(refcolumn[j])){
        present <- present + 1
      }
    }
    
  }
  
  comparisons <- counter*nrow(dat)
  match <- round(100*matches/comparisons, 1)
  
  print(c(match, comparisons-matches, missing, present))
  
  return(c(match, comparisons-matches ,missing, present))
}


#How many isolates have the same alleles as each other between datasets?
find_match_alleles <- function(dat, pubmlstdat){
  counter <- 0
  correlations <- c()
  
  for(i in 3:ncol(dat)){
    columnname <- colnames(dat)[i]
    if(columnname == 'capsule') next
    column <- dat[,columnname]
    refcolumn <- pubmlstdat[,columnname]
    
    counter <- counter + 1
    if(length(column) != length(refcolumn)){
      print(i)
      print(columnname)
      break
    }
    
    #isolates with semicolons in their allele ids should be made NA in both column and refcolumn
    for(j in 1:length(column)){
      if(is.na(refcolumn[j])) next
      if(str_detect(refcolumn[j], ';')){
        refcolumn[j] <- NA
        column[j] <- NA
      }
    }
    
    colvals <- which(!is.na(column))
    refcolvals  <- which(!is.na(refcolumn))
    values <- which(colvals %in% refcolvals)
    
    if(length(values) < 3) next
    
    #correlate the two columns, excluding NA values
    correlation <- cor.test(as.numeric(column), as.numeric(refcolumn))$estimate
    
    correlations <- c(correlations, correlation)
    
  }
  return(mean(na.omit(correlations)))
}




#check CDS annos against pubmlst and GC data 
for(i in 1:nrow(genomecomp)){
  for(j in 1:ncol(genomecomp)){
    if(genomecomp[i,j] == 'X' | genomecomp[i,j] == 'I'){
      genomecomp[i,j] <- NA
    }
  }
}

synteny_thresholds <- seq(0.1, 0.9, 0.1)
out <- data.frame('synteny_threshold' = rep(synteny_thresholds, 15),
                  'flanklength' = NA,
                  'id_cov' = NA
)
flanklength <- c(rep(500,9), rep(1000,9), rep(2000,9))
out$flanklength <- as.factor(flanklength)
id_cov <- c(rep('id = 90%, cov = 90%',27),
            rep('id = 90%, cov = 70%',27),
            rep('id = 70%, cov = 90%',27),
            rep('id = 70%, cov = 70%',27),
            rep('id = 97%, cov = 70%',27)
)
out$id_cov <- id_cov


base <- '~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/'
s1 <- 'skipIGR/id'

for(i in 1:nrow(out)){
  
  size <- out$flanklength[i]
  id <- out$id_cov[i]
  cov <- as.numeric(substr(id, 17,18))
  id <- as.numeric(substr(id, 6,7))
  syn <- out$synteny_threshold[i]
  
  filename <- paste('igr_annos_synthres-', syn, '_flank-', size, '_id-', id, '_cov-', cov, '.xlsx', sep='')
  
  print(filename)
  root <- paste(base, s1, id, 'cov', cov, '/', sep='')
  dat <- read.xlsx(paste(root, filename, sep=''))
  dat %<>% arrange(id, pubmlstdat$id)
  out$pubmlst_export[i] <- find_match(dat, pubmlstdat)[1]
  out$genome_comparator[i] <- find_match(dat, genomecomp)[1]
}


p1 <- ggplot(data=out, aes(x=as.factor(synteny_threshold), y=pubmlst_export, color = id_cov, shape = flanklength)) +
  geom_point(stat = 'identity') +
  theme_classic(base_size=15) +
  xlab('Synteny threshold') +
  ylab('Similarity (%)') +
  ylim(c(50,100)) +
  ggokabeito::scale_color_okabe_ito() +
  guides(shape=guide_legend(title="Flank sequence length (nt)", title.position='top')) +
  guides(color=guide_legend(title="BLAST settings", title.position='top', nrow = 2)) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.box = 'horizontal',
        legend.text = element_text(size=12), legend.title=element_text(size=14),
        legend.spacing.x = unit(20, 'pt'), legend.margin = margin(5,0,2,0))

p2 <- ggplot(data=out, aes(x=as.factor(synteny_threshold), y=genome_comparator, color = id_cov, shape = flanklength)) +
  geom_point(stat = 'identity') +
  theme_classic(base_size=15) +
  ggokabeito::scale_color_okabe_ito() +
  ylim(c(50,100)) +
  xlab('Synteny threshold') +
  ylab('Similarity (%)') +
  theme(legend.position = 'none')

pdf('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Figures/syn_thres_optimisation_CDS_final.pdf',
  width = 10, height = 7)
p <- ggarrange(p1, p2, labels = c('A)', 'B)'), ncol=2, legend.grob=get_legend(p1), common.legend = T, legend = 'bottom')
print(p)
dev.off()


# #check out the correlations among isolates belonging to specific serogroups
# dat <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/igr_annos_synthres_10.xlsx')
# dat <- dat[which(dat$capsule_group =='Y'),]
# pubmlstdat_sub <- pubmlstdat[which(pubmlstdat$capsule_group =='B'),]
# genomecomp_sub <- genomecomp[which(genomecomp$id %in% dat$id),]
# 
# find_match(dat, pubmlstdat_sub)
# find_match(dat, genomecomp_sub)



# 
# #Repeat for IGR data and just see completeness of the excel file
# #initialise reporting object
# synteny_thresholds <- seq(0.1, 0.9, 0.1)
# out <- data.frame('synteny_threshold' = rep(synteny_thresholds, 15),
#                   'flanklength' = NA,
#                   'id_cov' = NA
#                   )
# flanklength <- c(rep(500,9), rep(1000,9), rep(2000,9))
# out$flanklength <- as.factor(flanklength)
# id_cov <- c(rep('id = 90%, cov = 90%',27),
#             rep('id = 90%, cov = 70%',27),
#             rep('id = 70%, cov = 90%',27),
#             rep('id = 70%, cov = 70%',27),
#             rep('id = 97%, cov = 70%',27)
# )
# out$id_cov <- id_cov
# 
# 
# base <- '~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/'
# s1 <- 'noskipIGR/id'
# 
# for(i in 1:nrow(out)){
#   
#   size <- out$flanklength[i]
#   id <- out$id_cov[i]
#   cov <- as.numeric(substr(id, 17,18))
#   id <- as.numeric(substr(id, 6,7))
#   syn <- out$synteny_threshold[i]
#   
#   filename <- paste('igr_annos_synthres-', syn, '_flank-', size, '_id-', id, '_cov-', cov, '.xlsx', sep='')
#   
#   print(filename)
#   root <- paste(base, s1, id, 'cov', cov, '/', sep='')
#   dat <- read.xlsx(paste(root, filename, sep=''))
#   out$completeness[i] <- round(100*sum(!is.na(dat))/(nrow(dat)*ncol(dat)),1)
# }
# 
# ggplot(data=out, aes(x=as.factor(synteny_threshold), y=completeness, color = id_cov, shape = flanklength)) +
#   geom_point(stat = 'identity') +
#   theme_classic(base_size=15) +
#   ggokabeito::scale_color_okabe_ito() +
#   ylim(c(0,100)) +
#   xlab('Synteny threshold') +
#   ylab('Completeness (%)') +
#   guides(color=guide_legend(title="BLAST settings")) +
#   guides(shape=guide_legend(title="Flank sequence length (bp)"))
# 
# 

#compare Genome comparator output that used refernce IGR seqs to annotate each isolate
gc_igr <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/GC_IGRflank_data.xlsx')
for(i in 1:nrow(gc_igr)){
  for(j in 1:ncol(gc_igr)){
    if(gc_igr[i,j] == 'X' | gc_igr[i,j] == 'I'){
      gc_igr[i,j] <- NA
    }
  }
}


synteny_thresholds <- seq(0.1, 0.9, 0.1)
out <- data.frame('synteny_threshold' = rep(synteny_thresholds, 15),
                  'flanklength' = NA,
                  'id_cov' = NA
)
flanklength <- c(rep(500,9), rep(1000,9), rep(2000,9))
out$flanklength <- as.factor(flanklength)
id_cov <- c(rep('id = 90%, cov = 90%',27),
            rep('id = 90%, cov = 70%',27),
            rep('id = 70%, cov = 90%',27),
            rep('id = 70%, cov = 70%',27),
            rep('id = 97%, cov = 70%',27)
)
out$id_cov <- id_cov


base <- '~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/'
s1 <- 'noskipIGR/id'

for(i in 1:nrow(out)){
  
  size <- out$flanklength[i]
  id <- out$id_cov[i]
  cov <- as.numeric(substr(id, 17,18))
  id <- as.numeric(substr(id, 6,7))
  syn <- out$synteny_threshold[i]
  
  filename <- paste('igr_annos_synthres-', syn, '_flank-', size, '_id-', id, '_cov-', cov, '.xlsx', sep='')
  
  print(filename)
  root <- paste(base, s1, id, 'cov', cov, '/', sep='')
  dat <- read.xlsx(paste(root, filename, sep=''))
  dat %<>% arrange(id, gc_igr$id)
  out$gc_igrs[i] <- find_match(dat, gc_igr)[1]
}

p3 <- ggplot(data=out, aes(x=as.factor(synteny_threshold), y=gc_igrs, color = id_cov, shape = flanklength)) +
  geom_point(stat = 'identity') +
  theme_classic(base_size=15) +
  ggokabeito::scale_color_okabe_ito() +
  ylim(c(50,100)) +
  xlab('Synteny threshold') +
  ylab('Similarity (%)') +
  theme(legend.position = 'none')
  # guides(color=guide_legend(title="BLAST settings")) +
  # guides(shape=guide_legend(title="Flank sequence length (bp)"))
  

pdf('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Figures/syn_thres_optimisation_A-CDS_B-IGR_GCdata.pdf',
      width = 10, height = 7)
p <- ggarrange(p2, p3, labels = c('A)', 'B)'), ncol=2, legend.grob=get_legend(p1), common.legend = T, legend = 'bottom')
print(p)
dev.off()



#################################################################
#Repeat everything with correlations
#check CDS annos against pubmlst and GC data 
for(i in 1:nrow(genomecomp)){
  for(j in 1:ncol(genomecomp)){
    if(is.na(genomecomp[i,j])) next
    if(genomecomp[i,j] == 'X' | genomecomp[i,j] == 'I'){
      genomecomp[i,j] <- NA
    }
  }
}
synteny_thresholds <- seq(0.1, 0.9, 0.1)
out <- data.frame('synteny_threshold' = rep(synteny_thresholds, 15),
                  'flanklength' = NA,
                  'id_cov' = NA
)
flanklength <- c(rep(500,9), rep(1000,9), rep(2000,9))
out$flanklength <- as.factor(flanklength)
id_cov <- c(rep('id = 90%, cov = 90%',27),
            rep('id = 90%, cov = 70%',27),
            rep('id = 70%, cov = 90%',27),
            rep('id = 70%, cov = 70%',27),
            rep('id = 97%, cov = 70%',27)
)
out$id_cov <- id_cov


base <- '~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/'
s1 <- 'skipIGR/id'

for(i in 1:nrow(out)){
  
  size <- out$flanklength[i]
  id <- out$id_cov[i]
  cov <- as.numeric(substr(id, 17,18))
  id <- as.numeric(substr(id, 6,7))
  syn <- out$synteny_threshold[i]
  
  filename <- paste('igr_annos_synthres-', syn, '_flank-', size, '_id-', id, '_cov-', cov, '.xlsx', sep='')
  
  print(filename)
  root <- paste(base, s1, id, 'cov', cov, '/', sep='')
  dat <- read.xlsx(paste(root, filename, sep=''))
  dat %<>% arrange(id, pubmlstdat$id)
  out$pubmlst_export[i] <- find_match_alleles(dat, pubmlstdat)
  out$genome_comparator[i] <- find_match_alleles(dat, genomecomp)
}

p1 <- ggplot(data=out, aes(x=as.factor(synteny_threshold), y=pubmlst_export, color = id_cov, shape = flanklength)) +
  geom_point(stat = 'identity') +
  theme_classic(base_size=15) +
  xlab('Synteny threshold') +
  ylab('Correlation') +
  ylim(c(0,1)) +
  ggokabeito::scale_color_okabe_ito() +
  guides(shape=guide_legend(title="Flank sequence length (nt)", title.position='top')) +
  guides(color=guide_legend(title="BLAST settings", title.position='top', nrow = 2)) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.box = 'horizontal',
        legend.text = element_text(size=12), legend.title=element_text(size=14),
        legend.spacing.x = unit(20, 'pt'), legend.margin = margin(5,0,2,0))

p2 <- ggplot(data=out, aes(x=as.factor(synteny_threshold), y=genome_comparator, color = id_cov, shape = flanklength)) +
  geom_point(stat = 'identity') +
  theme_classic(base_size=15) +
  ggokabeito::scale_color_okabe_ito() +
  ylim(c(0,1)) +
  xlab('Synteny threshold') +
  ylab('Correlation') +
  theme(legend.position = 'none')

pdf('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Figures/syn_thres_optimisation_CDS_correlations.pdf',
    width = 10, height = 7)
p <- ggarrange(p1, p2, labels = c('A)', 'B)'), ncol=2, legend.grob=get_legend(p1), common.legend = T, legend = 'bottom')
print(p)
dev.off()



#compare Genome comparator output that used refernce IGR seqs to annotate each isolate
gc_igr <- read.xlsx('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/GC_IGR_data.xlsx')
for(i in 1:nrow(gc_igr)){
  for(j in 1:ncol(gc_igr)){
    if(gc_igr[i,j] == 'X' | gc_igr[i,j] == 'I'){
      gc_igr[i,j] <- NA
    }
  }
}


synteny_thresholds <- seq(0.1, 0.9, 0.1)
out <- data.frame('synteny_threshold' = rep(synteny_thresholds, 15),
                  'flanklength' = NA,
                  'id_cov' = NA
)
flanklength <- c(rep(500,9), rep(1000,9), rep(2000,9))
out$flanklength <- as.factor(flanklength)
id_cov <- c(rep('id = 90%, cov = 90%',27),
            rep('id = 90%, cov = 70%',27),
            rep('id = 70%, cov = 90%',27),
            rep('id = 70%, cov = 70%',27),
            rep('id = 97%, cov = 70%',27)
)
out$id_cov <- id_cov


base <- '~/Documents/PhD/Asides/Disease_carriage_classifier/syn_test_output/'
s1 <- 'noskipIGR/id'

for(i in 1:nrow(out)){
  
  size <- out$flanklength[i]
  id <- out$id_cov[i]
  cov <- as.numeric(substr(id, 17,18))
  id <- as.numeric(substr(id, 6,7))
  syn <- out$synteny_threshold[i]
  
  filename <- paste('igr_annos_synthres-', syn, '_flank-', size, '_id-', id, '_cov-', cov, '.xlsx', sep='')
  
  print(filename)
  root <- paste(base, s1, id, 'cov', cov, '/', sep='')
  dat <- read.xlsx(paste(root, filename, sep=''))
  dat %<>% arrange(id, gc_igr$id)
  out$gc_igrs[i] <- find_match_alleles(dat, gc_igr)
}

p3 <- ggplot(data=out, aes(x=as.factor(synteny_threshold), y=gc_igrs, color = id_cov, shape = flanklength)) +
  geom_point(stat = 'identity') +
  theme_classic(base_size=15) +
  ggokabeito::scale_color_okabe_ito() +
  ylim(c(0,1)) +
  xlab('Synteny threshold') +
  ylab('Correlation') +
  theme(legend.position = 'none')
# guides(color=guide_legend(title="BLAST settings")) +
# guides(shape=guide_legend(title="Flank sequence length (bp)"))


pdf('~/Library/CloudStorage/OneDrive-UniversityofWarwick/Documents/PhD/Asides/Disease_carriage_classifier/Figures/syn_thres_optimisation_A-CDS_B-IGR_GCdata_cors.pdf',
    width = 10, height = 7)
p <- ggarrange(p2, p3, labels = c('A)', 'B)'), ncol=2, legend.grob=get_legend(p1), common.legend = T, legend = 'bottom')
print(p)
dev.off()





#Make plots of CDS similarity and correlation, a
