#This script uses an R library of samtools to subset sorted and indexed bam.bai files based on 
#the coordinates of PV genes within the reference used during bam file generation with bwa mem

library(openxlsx)
library(Rsamtools)

#Load in phasomeit.out
PV_data <- read.xlsx('~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out.xlsx', sheet =4)[c(1:79),c(1,40,41)]

#Remove empty rows 
PV_data <- na.omit(PV_data)



#Function to extract reads from a bam file and save them as a new, smaller bam locally
reads <- function(start, stop, isolate, PV_locus, run){
  #Get file path of the bamfile
  bamfile <- paste('~/Documents/PhD/PhD/RNA_IGR/N222.1.2_bams/3_bam/', isolate, '-', run, '.bam', sep ='')
  indexfile <- paste('~/Documents/PhD/PhD/RNA_IGR/N222.1.2_bams/4_index/', isolate, '-', run, sep ='')
  #Create coordinates object for the PV locus
  if(start > stop){
    temp <- stop
    stop <- start
    start <- stop
  }
  which <- GRanges(paste('Neisseria:', start, '-', stop, sep = ''))
  #Define fields to recover from reads in those coordinates
  what <- c('seq')
  #Set which and what as a parameter object
  param <- ScanBamParam(which=which, what = what)
  
  #Subset the bam according to the params and write to a new file
  destination <- paste('~/Documents/PhD/PhD/PhasomeIt_data/PV_bams/', PV_locus, '_', isolate, '.bam', sep ='')
  filterBam(bamfile, param=param, destination = destination, index = indexfile)

}

#Function to loop through all the PC coords for an isolate
main <- function(PV_data, isolates){

  for(isolate in isolates){
    print(isolate)
    #Loop through all the PV loci for this isolate
    isolate_bams(isolate, PV_data)
  }
}

isolate_bams <- function(isolate, PV_data){
  for(locus in 1:length(PV_data$PubMLST_id)){
    reads(start = PV_data$N222.1.2_start[locus], stop = PV_data$N222.1.2_stop[locus], isolate=isolate, PV_locus=PV_data$PubMLST_id[locus], run='1')
    reads(start = PV_data$N222.1.2_start[locus], stop = PV_data$N222.1.2_stop[locus], isolate=isolate, PV_locus=PV_data$PubMLST_id[locus], run='2')
    reads(start = PV_data$N222.1.2_start[locus], stop = PV_data$N222.1.2_stop[locus], isolate=isolate, PV_locus=PV_data$PubMLST_id[locus], run='3')
  }
}

isolates <- c('20578_1', '22689_1', 'N188_1', 'N222_1', 'N222_2', 'N445_1', 'N459_3', 'N459_6')
main(PV_data, isolates)


