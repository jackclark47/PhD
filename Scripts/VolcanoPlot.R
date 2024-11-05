#Volcano plot making

library(openxlsx)
library(ggplot2)
library(EnhancedVolcano)

#Read in rna dataset
rna <- read.xlsx("~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx", sheet = 6, rowNames = T)[2:2272,c(8, 12:154)]

res = rna[,1:3]
colnames(res) <- c("Locus", "qval", "log2 fold change")

res$qval <- as.numeric(res$qval)

for(i in 1:length(res$qval)){
  num = res$qval[i]
  if(num == 0){
    num = 1*10^-300
    print(num)
    res$qval[i] = num
  }
}

res$Colour <- 1
for(i in 1:length(res$qval)){
  print(i)
  if(is.na(res$qval[i]) | is.na(res$`log2 fold change`[i])){
    next
  }
  if(res$qval[i] < 0.05 && res$`log2 fold change`[i] >= 1){
    res$Colour[i] <- 2
  }
}


res$qval = -log(res$qval)




#volcano plots plot log2fold change of each locus against its p value
alpha = 0.05

plot(res$`log2 fold change`, y=res$qval, main = 'Volcano plot', xlab = "log2-fold change",
     ylab = "-log10(qvalue)", pch=20, cex=0.6,
     col = res$Colour)
abline(v=1, col='brown')
abline(h=-log(alpha), col='brown')
table(res$Colour)




gn.selected <- res$`log2 fold change` >= 2 & res$qval > -log(alpha)
text(res$`log2 fold change`[gn.selected],
     res$qval[gn.selected],
     lab=res$Locus[gn.selected], cex = 0.4)
