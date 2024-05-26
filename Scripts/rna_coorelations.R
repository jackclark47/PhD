#Script to correlate rna expression data between pairs of isolates.

#Load libraries
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(stringr)
#read in rna data - expression
rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 6)[, c(8,9,15:22)]

#rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 6)[, c(8,9,25:34)]
#read in rna data - log2(expression/MIN)

###NOTE
#if you use the log2 exp/min, set minimum scale to be slightly above 0 (like 0.000001) to remove the datapoints where the primary isolate is the minimum
#Make sure to check this every time before running the script
#Also, when doing the raw expression, you can adjust the max to be the full range of hits, or set the max to around 5000, as this removes a few ultra highly expressed loci which change the correlation quite a bit




#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

rna <- rna[,c(1:4, 6,5, 8,7, 9,10)]
colnames(rna) <- c('Synonym', 'PubMLST_ID', '27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
rna <- rna[-c(1),-c(1)]

#Plan here is to get correlations for each pairwise isolate-isolate comparison. 
#This can be stored as a simple matrix
#Look at the expression at each locus and see if the second isolate in the comparison has similar expression. 
#i.e. when expression is high at one locus in isolate 1, is it also high at that locus in the second isolate?

#This is best done by making a panel of graphs
#The panel consists of 8 graphs - one for each isolate
#Each graph is a scatter plot where expression is on the x and the y axis
#Each point on the scatter plot is a locus, coloured based on the isolate its from.
#Y axis expression is the unchanging isolate's expression values, x axis is expression for the changing isolates.
#So points for a locus will have the same y axis location, but change on the x axis depending on the exp of the seconday isolate at that locus
#SO there'll be a panel of 8 graphs, each consisting of 7 sets of coloured points. 


#Make a function that performs a pairwise comparison, taking two isolates as input
cor_exp <- function(isolate_1, isolate_2){

}

cor_exp()

#Test case
#Primary isolate will be the first, 
# ggplot(rna, aes(x=`20578_1`, y=`22689_1`)) +
#   geom_point() +
#   geom_smooth(method='lm', se=FALSE)

rna[,2] <- as.numeric(rna[,2])
rna[,3] <- as.numeric(rna[,3])
rna[,4] <- as.numeric(rna[,4])
rna[,5] <- as.numeric(rna[,5])
rna[,6] <- as.numeric(rna[,6])
rna[,7] <- as.numeric(rna[,7])
rna[,8] <- as.numeric(rna[,8])
rna[,9] <- as.numeric(rna[,9])



#Make table just have four columns: Locus, Primary Isolate exp, Secondary Isolate_exp, Secondary Isolate Name
table27509 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "27509"), variable.name="")
colnames(table27509) <- c('Locus', '27509_exp', 'Isolate', 'Secondary_exp')

table27553 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "27553"), variable.name="")
colnames(table27553) <- c('Locus', '27553_exp', 'Isolate', 'Secondary_exp')

table28262 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "28262"), variable.name="")
colnames(table28262) <- c('Locus', '28262_exp', 'Isolate', 'Secondary_exp')

table28269 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "28269"), variable.name="")
colnames(table28269) <- c('Locus', '28269_exp', 'Isolate', 'Secondary_exp')

table28287 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "28287"), variable.name="")
colnames(table28287) <- c('Locus', '28287_exp', 'Isolate', 'Secondary_exp')

table53930 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "53930"), variable.name="")
colnames(table53930) <- c('Locus', '53930_exp', 'Isolate', 'Secondary_exp')

table53948 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "53948"), variable.name="")
colnames(table53948) <- c('Locus', '53948_exp', 'Isolate', 'Secondary_exp')

table53951 <- data.table::melt(rna, id.vars = c("PubMLST_ID", "53951"), variable.name="")
colnames(table53951) <- c('Locus', '53951_exp', 'Isolate', 'Secondary_exp')

#plotting the expression data. For log2/exp comment this segment and uncomment the next one underneath
p1 <- ggplot(table27509, aes(x=`27509_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.5) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table27509$`27509_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))

p2 <- ggplot(table27553, aes(x=`27553_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table27553$`27553_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27509" = "#F8766D", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))

p3 <- ggplot(table28262, aes(x=`28262_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table28262$`28262_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))

p4 <- ggplot(table28269, aes(x=`28269_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table28269$`28269_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))

p5 <- ggplot(table28287, aes(x=`28287_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table28287$`28287_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))

p6 <- ggplot(table53930, aes(x=`53930_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table53930$`53930_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53948" = "#C77CFF", "53951" = "#FF61CC"))

p7 <- ggplot(table53948, aes(x=`53948_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table53948$`53948_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53951" = "#FF61CC"))

p8 <- ggplot(table53951, aes(x=`53951_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
  geom_point(size=0.1) +
  geom_smooth(method='lm', se=FALSE, fullrange=TRUE) +
  labs(y = 'expression') +
  scale_x_continuous(limits=c(0, max(table53951$`53951_exp`))) +
  stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
  scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF"))

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,
                         nrow = 2)

#The plotting for the log2/exp data
# 
# p1 <- ggplot(table27509, aes(x=`27509_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table27509$`27509_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))
# 
# p2 <- ggplot(table27553, aes(x=`27553_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table27553$`27553_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27509" = "#F8766D", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))
# 
# p3 <- ggplot(table28262, aes(x=`28262_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table28262$`28262_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))
# 
# p4 <- ggplot(table28269, aes(x=`28269_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table28269$`28269_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))
# 
# p5 <- ggplot(table28287, aes(x=`28287_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table28287$`28287_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "53930" = "#00A9FF", "53948" = "#C77CFF", "53951" = "#FF61CC"))
# 
# p6 <- ggplot(table53930, aes(x=`53930_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table53930$`53930_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53948" = "#C77CFF", "53951" = "#FF61CC"))
# 
# p7 <- ggplot(table53948, aes(x=`53948_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table53948$`53948_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53951" = "#FF61CC"))
# 
# p8 <- ggplot(table53951, aes(x=`53951_exp`, y=`Secondary_exp`, group=Isolate, color=Isolate)) +
#   geom_point(size=0.1) +
#   geom_smooth(method='lm', se=FALSE, fullrange=TRUE, linewidth=0.6) +
#   labs(y = 'expression') +
#   scale_x_continuous(limits=c(0.000001, max(table53951$`53951_exp`))) +
#   stat_cor(aes(label = after_stat(rr.label)), color = "black", geom="label") +
#   scale_color_manual(values = c("27509" = "#F8766D", "27553" = "#CD9600", "28262" = "#7CAE00", "28269" = "#00BE67", "28287" = "#00BFC4", "53930" = "#00A9FF", "53948" = "#C77CFF"))
# 
# grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,
#              nrow = 2)

