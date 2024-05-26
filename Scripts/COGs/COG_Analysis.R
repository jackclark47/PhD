#Script to get a table of COG classes for each locus
#Structure is:
#Locus    #COG class
#NEIS001  #COG54
#NEIS002  #COG19

library(BiocManager)
library(gggenomes)
library(Biostrings)
library(rBLAST)
library(seqinr)
#Read in the annotated .gbk. read_gbk converts the file to a .gff3 and reads it in that format
isolate <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222.2.gbk")
isolateseq <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222.2.fasta", as.string = TRUE)

table(isolate$type)
#First issue is the annotated file does not have NEIS codes, just a numeric order and the protein sequence.


#Load blast database
blast_db <- blast(db="~/Documents/PhD/PhD/N_meningitidis_loci/translations/BLAST_db/coding.fsa", type = 'blastp')

#make empty dataframe
out <- as.data.frame(matrix(data = NA, nrow=nrow(isolate), ncol = 2))
colnames(out) <- c("id", "COG")
#Need to extract that protein sequence and blast it against every locus to find out which locus it is and replace the NEIS code
for(i in 2:nrow(isolate)){
  sequence <- substr(isolateseq, start = isolate$start[i], stop = isolate$end[i])
  #Get the COG class for that sequence, skip any seq that has no COG class as thats all im interested in
  if(is.na(isolate$dbxref[i])){
    next
  }
  
  COG_class <- isolate$dbxref[i]
  
  #translate the dna sequence and convert it to AAStringSet
  sequence <- DNAStringSet(sequence)
  sequence <- Biostrings::translate(sequence)
  
  #Run BLAST
  res <- predict(blast_db, sequence, type = 'blastp')[1,]
  
  #If the percent identity and evalue are low in the top result, discard
  if(is.na(res$pident) | is.na(res$evalue)){
    next
  }
  #Otherwise, assign the sequence the same name as the top blast hit
  if(res$pident > 80 && res$evalue < 0.005){
    pub_id <- res$sseqid
  }
  out$id[i] <- pub_id
  out$COG[i] <- COG_class
  
  
}
out
out2 <- out
out <- out2
length(out$id)
length(unique(out$id))
length(out$COG)
length(unique(out$COG))

#Now process the COG tables downloaded from ncbi to make one table of all cogs
A <- read.csv("~/Downloads/COG/cog_A.tsv" , sep = '\t')
B <- read.csv("~/Downloads/COG/cog_B.tsv" , sep = '\t')
C <- read.csv("~/Downloads/COG/cog_C.tsv" , sep = '\t')
D <- read.csv("~/Downloads/COG/cog_D.tsv" , sep = '\t')
E <- read.csv("~/Downloads/COG/cog_E.tsv" , sep = '\t')
f <- read.csv("~/Downloads/COG/cog_F.tsv" , sep = '\t')
G <- read.csv("~/Downloads/COG/cog_G.tsv" , sep = '\t')
H <- read.csv("~/Downloads/COG/cog_H.tsv" , sep = '\t')
I <- read.csv("~/Downloads/COG/cog_I.tsv" , sep = '\t')
J <- read.csv("~/Downloads/COG/cog_J.tsv" , sep = '\t')
K <- read.csv("~/Downloads/COG/cog_K.tsv" , sep = '\t')
L <- read.csv("~/Downloads/COG/cog_L.tsv" , sep = '\t')
M <- read.csv("~/Downloads/COG/cog_M.tsv" , sep = '\t')
N <- read.csv("~/Downloads/COG/cog_N.tsv" , sep = '\t')
O <- read.csv("~/Downloads/COG/cog_O.tsv" , sep = '\t')
P <- read.csv("~/Downloads/COG/cog_P.tsv" , sep = '\t')
Q <- read.csv("~/Downloads/COG/cog_Q.tsv" , sep = '\t')
R <- read.csv("~/Downloads/COG/cog_R.tsv" , sep = '\t')
S <- read.csv("~/Downloads/COG/cog_S.tsv" , sep = '\t')
t <- read.csv("~/Downloads/COG/cog_T.tsv" , sep = '\t')
U <- read.csv("~/Downloads/COG/cog_U.tsv" , sep = '\t')
V <- read.csv("~/Downloads/COG/cog_V.tsv" , sep = '\t')
W <- read.csv("~/Downloads/COG/cog_W.tsv" , sep = '\t')
X <- read.csv("~/Downloads/COG/cog_X.tsv" , sep = '\t')
Z <- read.csv("~/Downloads/COG/cog_Z.tsv" , sep = '\t')

COGs <- rbind(A, B, C, D, E, f, G, H, I, J, K, L, M, N, O, P, Q, R, S, t, U, V, W, X, Z)

out <- na.omit(out)
out$category <- NA
for(i in 1:nrow(out)){
  out$category[i] <- COGs$Cat[which(COGs$COG == substr(out$COG[i], 5, nchar(out$COG[i]) ) ) ]
}

out$COG[2] %in% COGs$COG
substr(out$COG[i], 5, nchar(out$COG[i])) %in% COGs$COG
COGs$Cat[which(COGs$COG == substr(out$COG[i], 5, nchar(out$COG[i]) ) ) ]

out$COG[1]

#The .gbk contains COG classes for a number of the annotated loci

#How many loci are there?
length(unique(out$id))
#405

#How many COG classes are there?
length(unique(out$COG)) 
#807
length(out$COG)


#How many loci have at least one COG class?


#How often does each COG class appear? Get a table of Class and number of appearances
write.xlsx(table(out$category), file = '~/Downloads/COG/COGs_whole_genome.xlsx')
#Now look at how the COGs are distributed throughout the rna dataset log2() > 1, >2 and the whole set
rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 9)[, c(8,9,14)]
#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

#Remove entries which don't have an NEIS code
rna <- rna[!is.na(rna$`1717.genes`),]
rna <- rna[!is.na(rna$`log2(EXP/Min)`),]
rna <- rna[-1, ]
# rna <- rna[rna$`log2(EXP/MIN)` >= 2,]

for(i in 2:nrow(rna)){
  out_rna$category[i] <- COGs$Cat[which(COGs$COG == substr(out$COG[i], 5, nchar(out_rna$COG[i]) ) ) ]
}

out$log2_exp_change <- NA
for(i in 1:nrow(out)){
  print(i)
  if(length(which(rna$`1717.genes` == out$id[i])) > 0){
    print('o')
    out$log2_exp_change[i] <- rna$`log2(EXP/Min)`[which(rna$`1717.genes` == out$id[i])]
  }
}

#get a table of those with log(EXP/MIN >= 1)
out_1 <- out[which(out$log2_exp_change >=1) ,]
write.xlsx(table(out_1$category), file = '~/Downloads/COG/COGs_above_1.xlsx')


#and a table of those with >= 2
out_2 <- out[which(out$log2_exp_change >=2), ]
write.xlsx(table(out_2$category), file = '~/Downloads/COG/COGs_above_2.xlsx')

#stats tests for testing significantly distributions of COG categories





#Extract a table of loci by COG class


  








  