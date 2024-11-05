igrs <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs.xlsx')

for(i in 1:nrow(igrs)){
  igrs$`Alleles-Down`[i] <- max(as.numeric(na.omit(as.numeric(igrs[i, 4:11]))))
  igrs$`Alleles-Up`[i] <- max(as.numeric(na.omit(as.numeric(igrs[i,12:19]))))
  igrs$Combined[i] <- igrs$`Alleles-Down`[i] + igrs$`Alleles-Up`[i]
}
write.xlsx(igrs, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs.xlsx')


lowvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/lowvar_igrs.xlsx')
novar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/novar_igrs.xlsx')

removables <- c()
for(i in 1:nrow(lowvar)){
  if(lowvar$Locus[i] %in% novar$Locus){
    removables <- c(removables, i)
  }
}

length(removables)
nrow(lowvar)
1050 - 717
lowvar <- lowvar[-removables,]
write.xlsx(lowvar, '~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/lowvar_igrs.xlsx')
