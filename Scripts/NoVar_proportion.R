#How many loci with log(EXP/MIN) >1 are NoVar?



#Load RNA
rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 9)[, c(8,9,15:22, 25:32)]
#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}
#rearrange columns in the rna dataset so theyre the same order
rna <- rna[,c(1:4, 6,5, 8,7, 9:12, 14,13,16,15, 17, 18)]


#Remove entries which don't have an NEIS code
rna <- rna[!is.na(rna$`1717.genes`),]
rna <- rna[!is.na(rna$`log2(EXP/MIN)`),]
# rna <- rna[rna$`log2(EXP/MIN)` >= 2,]


#Load NoVar list
NoVar <- read.xlsx('igr_new_29112023.xlsx', sheet = 2)

count = 0
for(i in 1:length(NoVar$Locus)){
  if(NoVar$Locus[i] %in% rna$`1717.genes`){
    count = count + 1
  }
}
count

length(NoVar$Locus)
length(unique(rna$`1717.genes`))
table(rna$`1717.genes`)
