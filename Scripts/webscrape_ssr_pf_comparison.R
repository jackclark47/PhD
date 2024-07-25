library(openxlsx)
library(rjson)
library(stringr)

##Processing ssr data
long_ssr <- read.xlsx('D:/PhD/pv_webscrape/02412table5.xlsx', sheet = 1)[-379,]
igr_ssr <- read.xlsx('D:/PhD/pv_webscrape/02412table5.xlsx', sheet = 3)[-379,]

#Combine the two datasets
ssrs <- cbind(long_ssr, igr_ssr)
ssrs <- ssrs[, c(1, 13, 26)]
ssrs$total <- ssrs$Sum + ssrs$Sum.1
ssrs <- ssrs[, -c(2,3)]

#Remove rows in the dataset which contain no identified ssr
removables <- c()
for(i in 1:nrow(ssrs)){
  if(ssrs$total[i] == 0){
    removables <- c(removables, i)
  }
}
ssrsfilt <- ssrs[-removables, ]
ssr_species_freq <- as.data.frame(table(ssrsfilt$Genome))


#Processing phasefinder data
invertons <- read.xlsx('D:/PhD/pv_webscrape/invertons.xlsx', sheet = 3)
inverton_genome_ids <- unique(invertons$Genome)
write.csv(inverton_genome_ids, file = 'D:/PhD/pv_webscrape/inverton_genome_idlist.csv', row.names = F, quote = FALSE)

#Get species names for the GCF_XXXXXXX.X codes
species <- c()
for(file in list.files('D:/PhD/pv_webscrape/inverton_genome_metadata', full.names = T)){
  print(file)
  summary <- fromJSON(file=file)
  specie <- summary$reports[[1]]$average_nucleotide_identity$submitted_species
  species <- c(species, specie)
}
nrow(table(species))
#Write initial list to a csv and edit it in notepad to remove strain information
write.csv(species, 'D:/PhD/pv_webscrape/phasefinderspecies.csv')
species <- read.csv('D:/PhD/pv_webscrape/phasefinderspecies.csv', header = F)
pf_species_freq <- as.data.frame(table(species))
#Loading webscraper data for abstracts
webscrape_species <- read.csv('D:/PhD/pv_webscrape/abstract_species.csv', header = F)
webscrape_species_freq <- as.data.frame(table(webscrape_species))

##Compare
#How many unique species in each set?
nrow(ssr_species_freq) #246
nrow(pf_species_freq) #231
nrow(webscrape_species_freq) #76

#How many unique genera in each set?
extract_genera <- function(x, unique = TRUE){
  genera <- c()
  for(specie in x[,1]){
    genus <- str_extract(specie, '([:alpha:]+)(?=[:space:])')
    genera <- c(genera, genus)
  }
  if(unique == TRUE){
    genera <- unique(genera)
  }
  return(genera)
}

ssr_genera <- extract_genera(ssr_species_freq, unique = TRUE)
phasefinder_genera <- extract_genera(pf_species_freq, unique = TRUE)
webscrape_genera <- extract_genera(webscrape_species_freq, unique = TRUE)

length(ssr_genera) #153 genera
length(phasefinder_genera) #133 genera
length(webscrape_genera) #47 genera

#Now look at overlaps on a species level
sum(pf_species_freq$V1 %in% webscrape_species_freq$V1) #21 in common
sum(ssr_species_freq$Var1 %in% webscrape_species_freq$V1) #34 in common
sum(pf_species_freq$V1 %in% ssr_species_freq$Var1) #33

#and on a genus level
sum(phasefinder_genera %in% webscrape_genera) #25 in common
sum(ssr_genera %in% webscrape_genera) #32 in common
sum(phasefinder_genera %in% ssr_genera) #35
