#This script investigates the operon_mapper output.
#It removes operons containing only one gene and then annotates the remaining hits with PubMLST IDs

library(openxlsx)
#Load in the operon_mapper output from the run on N222.2.fasta
mapperpred <- read.xlsx('~/Downloads/N222.2_operons/list_of_operons_3477439.xlsx')
#Load in my own predicted operons from Operons.R using N222.2.gbk
Allenpred <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/test.xlsx', sheet = 3)

#Fix formatting of the operon-mapper output
for(i in 1:nrow(mapperpred)){
  if(!(is.na(mapperpred$Operon[i]))){
    current_operon <- mapperpred$Operon[i]
  }
  if(is.na(mapperpred$Operon[i])){
    mapperpred$Operon[i] <- current_operon
  }
}

#Remove rows containing only an operon number
removables <- c()
for(i in 1:nrow(mapperpred)){
  if(length(na.omit(unlist(mapperpred[i,]))) == 1){
    removables <- c(removables, i)
  }
}

mapperpred <- mapperpred[-removables,]

#Remove operons containing only a single gene
operon_counts <- count(mapperpred, Operon)
removable_ops <- c()
for(i in 1:nrow(operon_counts)){
  if(operon_counts$n[i] == 1){
    removable_ops <- c(removable_ops, operon_counts$Operon[i])
  }
}

mapperpred <- mapperpred[which(!(mapperpred$Operon %in% removable_ops)),]

length(unique(mapperpred$Operon)) #434 operons identified covering 1254 genes. 2.89 genes per operon
length(unique(Allenpred$operon)) #330 operons identigied covering 817 genes. 2.48 genes per operon

#Check how many Allenpred operons and genes are contained within the mapperpred results
match = 0
for(i in 1:nrow(Allenpred)){
  if(Allenpred$start[i] %in% mapperpred$PosLeft | Allenpred$stop[i] %in% mapperpred$postRight){
    match = match + 1
  }
}
match #813 out of 817 of the genes predicted by Allen's code are contained in the operon-mapper results.
#So this suggests we should use operon-mapper instead as it may contain more genes.
