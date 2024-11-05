library(cowplot)
library(ggplot2)
library(stringr)
library(rBLAST)
library(ggsignif)
library(openxlsx)

sheetnums <- 2:9
main(sheetnums)

#Operon Numbers
preprocessing <- function(sheetnum){
  
  opdata<-read.xlsx("~/Documents/PhD/PhD/operon_mapper_res/operon_mapper_annos_contigs.xlsx", sheet = sheetnum)
  
  
  new_operons <- c(1) #Stores the new operon numbers for each row
  previous_operon <- opdata$Operon[1] #Initialise previous operon as the value of the first entry in opdata$Operon
  
  #Construct the loop here
  for(i in 2:nrow(opdata)){
    current_operon <- opdata$Operon[i]
    #Check if the current operon number is the same as thopdaae previous. If so, add the last value of new_operons in to new_operons again
    if( current_operon==previous_operon){
      new_operons <-c(new_operons, new_operons[i-1])
      
    }
    
    #if they're different, take the last value of new_operons, increase its value by one, and add it to the end of new_operons
    if(current_operon!=previous_operon ){
      new_operons <- c(new_operons,(new_operons[i-1]+1))
      
    }
    
    previous_operon=current_operon
    
  }
  
  opdata$Operon <- new_operons
  return(opdata)
}

#significant operons
find_sig_operons<-function(opdata){
  
  rna_seq<-read.xlsx("~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx", sheet = 8)[-1,]
  rna_seq$`1717.genes`<-str_trim(rna_seq$`1717.genes`)
  colnames(rna_seq)[13]<-"q_val"
  rna_seq$q_val<-as.numeric(rna_seq$q_val)
  print('file loaded')
  
  sig_loci <- c()
  for (i in 1:nrow(rna_seq)) {
    if (is.na(rna_seq$`log2(EXP/Min)`[i])|is.na(rna_seq$q_val[i])){
      next
    }
    if (rna_seq$`log2(EXP/Min)`[i]>=1 && rna_seq$q_val[i]<0.05){
      sig_loci<-c(sig_loci,rna_seq$`1717.genes`[i])
    }
  }
  print('sig_loci filled')
  
  #Add an empty column to opdata to store whether a given *gene* is significant, as this will allow us to plot later
  opdata$sig_gene <- 0
  #Repeat but to store whether an *operon* contains at least one significant gene
  opdata$sig_operon <- 0
  
  #Construct the for loop
  for (i in 1:nrow(opdata)){
    #Use %in%
    if (opdata$pubmlst_id [i] %in% sig_loci) {
      #Change the value of opdata$sig_operon to 1 for all rows in opdata with the same operon number
      opdata$sig_operon[which(opdata$Operon == opdata$Operon[i])] <- 1
      opdata$sig_gene[i]<-1
    }
  }
  print('sig_operons column filled')
  opdata$log_exp <- NA
  
  for(i in 1:nrow(opdata)){
    
    locus <- opdata$pubmlst_id[i]
    if (length(rna_seq$`log2(EXP/Min)`[which(rna_seq$`1717.genes` == locus)])==0) {
      next
    }
    if (is.na(opdata$pubmlst_id[i])){
      next
    }
    if(locus %in% rna_seq$`1717.genes`){
      #Extract the expression data for the current locus, by checking which row contains the locus name, retrieving the index of that row,
      #and then indexing the expression column with the same index
      opdata$log_exp[i] <- rna_seq$`log2(EXP/Min)`[which(rna_seq$`1717.genes` == locus)]
    }
  }
  return(opdata)
}

sig_opdata <- op_blast(sig_opdata = , out= out, out_row = 2, blastdb = , opseq = , i=1)

sig_plots<-function(opdata,sheetnums){
  O<-ggplot(opdata, aes(x=as.factor(sig_operon), y=as.numeric(log_exp), color = sig_operon)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width=0.1) +
    geom_signif(comparisons = list(c('0', '1')), test = 't.test',y_position=6.5)+
    xlab("sig_operon")+
    ylab("log2 expression change")
  
  G<-ggplot(opdata, aes(x=as.factor(sig_gene), y=as.numeric(log_exp), color = sig_gene)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width=0.1) +
    geom_signif(comparisons = list(c('0', '1')), test = 't.test',y_position=6.5)+
    xlab("sig_gene")+
    ylab("log2 expression change")
  
  plot_grid(O,G)
  
  graph<-plot_grid(O,G)
  
  title<-ggdraw() +
    draw_label(
      "Graph of operon and gene against expression",
      fontface= "bold",
      x=0,
      hjust=0
    )+
    theme(
      plot.margin = margin(0,0,0,7)
    )
  
  pdf(paste("~/Documents/PhD/PhD/operon_mapper_res/Figures/exp_sigop_nsop_siggene_nsgene",as.character(sheetnums),".pdf",sep = ""))
  print(plot_grid(
    title, graph,
    ncol = 1,
    rel_heights = c(0.1,1)
  ))
  dev.off()
}

get_coordinates <- function(operon_number, sig_opdata){
  
  #subset opdata to get all the rows that contain the operon number input to this function. You need to fill in the which and the columns you want to extract
  operon <- sig_opdata[which(sig_opdata$Operon==operon_number),]
 
  #Find the coordinates here. Think about how you might do this - a for loop isnt necessary here
  start<- min(operon$PosLeft)
  stop<- max(operon$postRight)
  contig<-operon$contig[1]
  strand<-operon$Strand[1]

  #Return a vector of the coordinates and contig number
  return(c(start, stop, contig,strand))
}

#Load the seqinr package to read and write fasta files
library(seqinr)
#Load Biostrings to let us reverse complement sequences
library(Biostrings)

extract_seqs <- function(coordinates, opnum, isolate){
  start<-as.numeric(coordinates[1])
  stop<-as.numeric(coordinates[2])
  contig<-coordinates[3]
  strand<-coordinates[4]
  print(coordinates[3])
  filepath <- paste('~/Documents/PhD/RNA_IGR/Isolate_fastas/', isolate, '.fasta', sep='')
  print('opening genome')
  print(filepath)
  genome = read.fasta(file = filepath, as.string = TRUE) #Load the full sequence of a specific isolate here
  print('genome loaded')
  #Write some code to extract the sequence. read.fasta just loads a sequence in as a string so you can manipulate it with the usual functions
  operonseq1 <- substr(genome[contig],start,stop)
  operonseq <- substr(genome[as.numeric(coordinates[3])], start, stop)
  if(nchar(operonseq1) > nchar(operonseq)){
    operonseq <- operonseq1
  }
  print(contig)

  #If our operon is on the reverse strand, we should reverse complement the sequence before saving it. See if you can find the function that does this (they often have intuitive names!)
  if(strand=="-"){
    operonseq<-reverseComplement(DNAString(operonseq)) #already defined operonseq above
  }

  names=paste(isolate,"_",opnum,sep="") 
  write.fasta(file.out=paste("~/Documents/PhD/PhD/operon_mapper_res/operon_fastas/", isolate,"/",isolate, '_operon_', opnum, '.fasta', sep = ''),
              sequences = operonseq, 
              names = names)
}

find_operons<-function(isolate, sig_opdata){
  
  opdata_unique<-unique(sig_opdata$Operon)
  for (i in opdata_unique) {
   coordinates<-get_coordinates(i,sig_opdata)
   print(coordinates)
   extract_seqs(coordinates,i,isolate) #no need to assign since writing file
  }
}

mean_opdata<-function(opdata,isolate){
  for (i in unique(opdata$Operon) ) {
    unique_opdata<-opdata[which(opdata$Operon==i),]
    print(as.numeric(unique_opdata$log_exp))
    avg <-mean(as.numeric(na.omit(unique_opdata$log_exp)))
    opdata$mean_exp_change[which(opdata$Operon==i)] <- avg
  }
  
  ###Need to remove duplicate rows of an operon.
  
  pdf(paste("~/Documents/PhD/PhD/operon_mapper_res/Figures/", isolate,".pdf",sep = ""))
  mean_graph<-ggplot(opdata, aes(x=as.factor(sig_operon), y=(as.numeric(mean_exp_change)), color = sig_operon)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width=0.1) +
    geom_signif(comparisons = list(c('0', '1')), test = 'wilcox.test')+
    xlab("mean_change")+
    ylab("sig_operon")
  
  print(mean_graph)
  dev.off()
}

main <- function(sheetnums){
  isolates<-c("27509","27553","28262","28269","28287","53930","53948","53951")
  for(sheetnum in sheetnums){
    isolate <-isolates[sheetnum-1]
    opdata <- preprocessing(sheetnum)
    print('perprocess done')
    opdata <- find_sig_operons(opdata)
    print('sig operons and genesidentified')
    sig_plots(opdata, sheetnum)
    print('sig plots made')
    mean_opdata(opdata[!duplicated(opdata$Operon),], isolate)
    print('mean plot made')
    sig_opdata<-opdata[which(opdata$sig_operon==1),]
    
    find_operons(isolate,sig_opdata)
    print('operons sequences extracted')

    write.xlsx(opdata, paste("~/Documents/PhD/PhD/operon_mapper_res/", sheetnum, isolate, '.xlsx', sep=''))
  }
}

###Here we will convert your individual fasta files to multifastas assuming that directory structure is present on your computer

#Function to take all the files in a given directory and combine them
multifasta <- function(filepath, isolate){
  #Get a list of the full paths to every file in the current folder
  files <- list.files(filepath, full.names = TRUE)
  #Read every file and assign to a list
  seqs <- lapply(files, readDNAStringSet)
  #run c() on each element of the list, combining them to a single DNAStringSet object
  seqs <- do.call(c, seqs)
  #Write the sequence - feel free to make the filename iterable (maybe include a vector of isolate names in the loop and add them as an argument to this function?)
  writeXStringSet(seqs, paste('~/Documents/PhD/PhD/operon_mapper_res/operon_fastas/multifastas/', isolate, 'multi.fasta', sep =''))
}

#Enter the path to your equivalent of fasta_files/ 
dirpath = '~/Documents/PhD/PhD/operon_mapper_res/operon_fastas'
isolates<-c("27509","27553","28262","28269","28287","53930","53948","53951")

#Loop through every directory and combine the files
i = 0
for(dir in list.dirs(dirpath, recursive=FALSE)[1:8]){
  i = i +1
  isolate <- isolates[i]
  multifasta(dir,isolate)
}


for(isolate in isolates){
  makeblastdb(paste("~/Documents/PhD/PhD/operon_mapper_res/operon_fastas/multifastas/", isolate, "multi.fasta", sep = ''),dbtype = "nucl")
}


#Load one of the databases using the path to it
db <- blast(db = "~/Documents/PhD/PhD/operon_mapper_res/operon_fastas/multifastas/27509multi.fasta", type = 'blastn')
#Obtain an operon sequence to read in 
operons <- readDNAStringSet("~/Documents/PhD/PhD/operon_mapper_res/operon_fastas/multifastas/27509multi.fasta")
opseq <- operons[1]
predict(db, opseq, type = 'blastn') 

#This sigopdata should be the sheet for 53930. We want every row but not all the columns
#Just the opnum, pubmlst id, sig_gene, exp and coords 
for(i in 1:length(isolates)){
  isolate <- isolates[i]
  opdata <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/opdata.xlsx', sheet = i)
  sigopdata <- opdata[which(opdata$sig_operon == 1),]
  write.xlsx(sigopdata, paste('~/Documents/PhD/PhD/operon_mapper_res/only_sigop_data_', isolate,'.xlsx', sep = ''))
}

#Read in sigopdata for 53930
sig_opdata<- read.xlsx("~/Documents/PhD/PhD/operon_mapper_res/only_sigop_data.xlsx", sheet = 6)

sig_opdata<-sig_opdata[which(sig_opdata$sig_operon==1),]

#Initialise out and set column names

out <- as.data.frame(matrix(data = NA, nrow = nrow(sig_opdata), ncol=22))

colnames(out) <- c('Operon', 'pubmlst_id', 'ngenes_27509', 'ngenes_27553', 'ngenes_28262', 'ngenes_28269', 'ngenes_28287', 'ngenes_53930', 'ngenes_53948', 'ngenes_53951', 'coords_27509', 'coords_27553', 'coords_28262', 'coords_28269', 'coords_28287', 'coords_53930', 'coords_53948', 'coords_53951', 'log2_exp(MAX/MIN)', 'q_val', 'up_igr', 'down_igr')



###Now lets add the 53930 data to out. We want one row in out for each operon so need to collapse the opdata rows

#length(unique(sig_opdata$Operon))

for(i in 1:nrow(out)){
  
  #Get all rows with the current operon number, i.
  operon <- sig_opdata[which(sig_opdata$Operon== i),]
  
  #out$opnum is just i thanks to how we define it in the loop
  out$Operon[i] <- i
  
  #ngenes_53930 is column 8 and will just be the number of rows in operon
  out$ngenes_53930[i] <- nrow(operon)

  #pubmlst ids is trickier - we need to combine the ids on multiple rows
  ids <- operon$pubmlst_id
  ids <- paste(ids, collapse=';')
  out$pubmlst_id[i] <- ids
  
  #coords is also less simple, we need to combine the start and stop coords as well as the contig
  #Lets format it like: contig:start-stop
  
  contig <- operon$contig[1]
  start <- operon$PosLeft[1] #Take first start coord
  stop <- operon$postRight[nrow(operon)] #take last stop coord
  
  #join them all up in our format
  coords <- paste(contig, ':', start, '-', stop, sep='')
  out$coords_53930[i] <- coords
  
  #log2_exp(MAX/MIN) is similar to pubmlst ids
  exp <- operon$log_exp
  print(exp)
  exp <- paste(exp, collapse=';')
  print(exp)
  out$`log2_exp(MAX/MIN)`[i] <- exp
  
  #qvals are the same
  #qvals <- operon$q_values
  #qvals <- paste(qvals, collapse=';')
  #out$q_val[i] <- qvals
}

#We also need to remove the empty rows in out
out <- out[which(out$ngenes_53930 != 0),]

#function1
#Define our function
op_format <- function(i,operon, out, out_row){
  
  #ngenes_53930 is column 8 and will just be the number of rows in operon
  out[out_row, (i + 2)] <- nrow(operon)
  
  print('operon rows are:')
  print(nrow(operon))
  print('out row is')
  print(out_row)
  print('new value is:')
  print(out[out_row, (i + 2)] )
  
  #coords is also less simple, we need to combine the start and stop coords as well as the contig
  #Lets format it like: contig:start-stop
  
  contig <- operon$contig[1]
  start <- operon$PosLeft[1] #Take first start coord
  stop <- operon$postRight[nrow(operon)] #take last stop coord
  
  #join them all up in our format
  coords <- paste(contig, ':', start, '-', stop, sep = '')
  out[out_row, (i + 10)] <- coords
  
  #Return out so that it stays updata and can be iterated over with each loop in the higher functions we will make soon
  return(out)
}


#function 2
op_blast <- function(sig_opdata, out, out_row, blastdb, opseq,i){
  
  blastresults<-predict(blastdb, opseq, type = 'blastn')[1,]
  print(blastresults)
  
  if(is.na(blastresults$pident)){
    return(out)
  }
  
  id<-blastresults$pident
  e_score<- as.numeric(blastresults$evalue)
  coverage<-(blastresults$length/nchar(opseq)) * 100
  
  if(id<75| e_score>=0.00001| coverage<85){
    return(out)
  }
  
  opname <- blastresults$sseqid
  opname <- str_extract(string=opname, pattern = '(?<=[:punct:])[:digit:]+')
  print('opname is:')
  print(opname)
  
  operon <- sig_opdata[which(sig_opdata$Operon == as.numeric(opname)),]
  print(operon)
  
  out <- op_format(i, operon, out, out_row)
  print(out[out_row, (i + 2)] )
  return(out)
}


#Function 3
all_op_blast <- function(blastdb, sig_opdata, out, i, isolate){
  
  for(out_row in 1:length(out$Operon)){
    
    ###This assumes list.files() lists operon files in the same order theyre listed in out. Not true!
    
    dirpath <- paste('~/Documents/PhD/PhD/operon_mapper_res/operon_fastas/53930',sep = '')
    files <- list.files(dirpath, full.names = T)
    #May need to check the file match line here to check it is correctly identifying the right file for the row in out
    file <- files[which(str_extract(files, '(?<=_)[:digit:]+') == as.character(out$Operon[out_row]))]
    
    opseq <- readDNAStringSet(file)

    out <- op_blast(sig_opdata, out, out_row, blastdb, opseq, i)
    print(out[out_row, (i + 2)] )
  }

  return(out)
}

#Function 4
main_blast <- function(filepath, out){
  isolates <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
  
  for(i in 1:length(isolates)){
    
    print('loop starting')
    isolate <- isolates[i]
    
    # if(isolate == '53930'){
    #   next
    # }
    
    sig_opdata <- read.xlsx(filepath, sheet = i)
    print('loaded sig_opdata')
    blastdb <- blast(db = paste('~/Documents/PhD/PhD/operon_mapper_res/operon_fastas/multifastas/', isolate, 'multi.fasta', sep = ''), type = 'blastn')
    print('loaded db')
    out <- all_op_blast(blastdb, sig_opdata, out, i, isolate)
    
  }
  return(out)
}

out <- main_blast('~/Documents/PhD/PhD/operon_mapper_res/only_sigop_data.xlsx', out)
igr<-read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx", sheet = 1)

for (i in 1:nrow(out)) {
  outpub<-out$pubmlst_id[i]
  print(outpub)
  library(stringr)
  outpub<-unlist(str_split(outpub,";"))
  first_outpub<-outpub[1]
  last_outpub<-outpub[length(outpub)]
  
  if(first_outpub!="NA"){
    first_igr<-igr[which(igr$Locus==first_outpub),12:19]
    first_igr<-paste(first_igr,collapse = ";")
    out$up_igr[i]<-first_igr
  }
  if(last_outpub!="NA"){
    last_igr<-igr[which(igr$Locus==last_outpub),12:19]
    last_igr<-paste(last_igr,collapse = ";")
    out$down_igr[i]<-last_igr 
  }
  
}
write.xlsx(out,"~/Documents/PhD/PhD/operon_mapper_res/operon_summary.xlsx")


#Load libraries for the script
library(openxlsx)
library(ggplot2)
library(ggsignif)
library(openxlsx)
library(stringr)
#Load data
opdata_53930<-read.xlsx("~/Documents/PhD/PhD/operon_mapper_res/only_sigop_data.xlsx", sheet = 6)
op_summary<-read.xlsx("~/Documents/PhD/PhD/operon_mapper_res/operon_summary.xlsx")


#Initialise an object to store the operon numbers of operons missing one or more igrs 
missing_igr<-c()
up_igr<-c()
#Loop through each row in the data and check to see if at least one of the igr columns is NA for that row. add the operon number to the igr missing object if so

for (i in 1:nrow(op_summary)) {
  
  if(is.na(op_summary$up_igr[i])){ 
    missing_igr<-c(missing_igr,op_summary$up_igr[i])
    up_igr<-c(up_igr,op_summary$up_igr[i])
  }
  if(is.na(op_summary$down_igr[i])){
    missing_igr<-c(missing_igr,op_summary$down_igr[i])
  }
}
#Initialise a new object to store the names of operons where all ngene columns are the same
#Loop through each row again, but this time check for operons where every n_genes column is the same
#I would get the ngenes number for 53930 and check all the others equal that

ngene_values<-c()
for (i in 1:nrow(op_summary)) {
  ngenes<-op_summary[i,c(3:10)]
  target<-op_summary[i,8]
  ngenes<-unlist(as.vector(ngenes))
  ngenes<-as.vector(na.omit(ngenes))
  if(!all(ngenes==target,na.rm=TRUE)){
    ngene_values<-c(ngene_values,op_summary$Operon[i])
    
  }
}



#Initialise three variables that store the operon numbers of operons that contain up, down or no igr variation
up_var<-c()
down_var<-c()
non_var<-c()

#Loop through each row of the table and identify operons containing igr variation based on the numbers in their igr column being different
for (i in 1:nrow(op_summary)) {
  up_variable<- op_summary$up_igr[i]
  up_variable<-str_split(up_variable,";")
  up_variable<-unlist(up_variable)
  down_variable<-op_summary$down_igr[i]
  down_variable<-str_split(down_variable,";" )
  down_variable<-unlist(down_variable)
  
  print(up_variable)
  print(down_variable)
  
  
  if(is.na(op_summary$up_igr[i])&& is.na(op_summary$down_igr[i])){
    next
  }
  if(is.na(op_summary$up_igr[i])){
    if(all(down_variable=='1')){
      non_var<-c(non_var,op_summary$Operon[i])
      next
      
    } else{
      down_var<-c(down_var,op_summary$Operon[i])
      next
    }
    
  }
  if(is.na(op_summary$down_igr[i])){
    if(all(up_variable=='1')){
      non_var<-c(non_var,op_summary$Operon[i])
      next
    } else{
      up_var<-c(up_var,op_summary$Operon[i])
      next
    }
    
  }
  if (all(up_variable=='1' ) && all(down_variable=='1')){
    non_var<-c(non_var,op_summary$Operon[i])
    print("yes")
  } else if(!all(up_variable=='1' ) && !all(down_variable=='1')) {
    down_var<-c(down_var,op_summary$Operon[i])
    up_var<-c(up_var,op_summary$Operon[i])
  }
  else if(all(up_variable=="1")){
    down_var<-c(down_var,op_summary$Operon[i])
    
  }else{
    up_var<-c(up_var,op_summary$Operon[i])
    
  }
  
}
#up and down variable should not be mutually exclusive groups

#Initialise 4 variables to store the categories - these will take the operon numbers
amb_gene<-c()
two_gene<-c()
three_gene<-c()
six_gene<-c()

#Loop through each row and obtain a vector of n_genes for the current row
for (i in 1:nrow(op_summary)) {
  ngenes<-op_summary[i,(3:10)]
  ngenes<-unlist(as.vector(ngenes))
  if(all(ngenes==2|is.na(ngenes))){
    two_gene<-c(two_gene,op_summary$Operon[i])
  } else if(all(ngenes %in% c(3,4,5)|is.na(ngenes))){
    three_gene<-c(three_gene,op_summary$Operon[i])
  }else if(all(ngenes>=6|is.na(ngenes))){
    six_gene<-c(six_gene,op_summary$Operon[i])  
  } else{
    amb_gene<-c(amb_gene,op_summary$Operon[i])
  }
    
}
#Multiple if statements, to check for each category. If an operon belongs to none of the 3 categories, assign it to the ambiguous group

#Loop through each row once more
one_sig<-c()
two_sig<-c()
three_sig<-c()
all_sig<-c()

for (i in 1:nrow(op_summary)) {
  opnum<-op_summary$Operon[i]
  comp_sig<-opdata_53930$sig_gene[which(opdata_53930$Operon==opnum)]
  sum_gene<-sum(comp_sig)
  
  if(sum_gene==length(comp_sig)){
    all_sig<-c(all_sig,op_summary$Operon[i])
  }
  
  if(sum_gene==1){
    one_sig<-c(one_sig,op_summary$Operon[i])
    
  }else if(sum_gene==2){
    two_sig<-c(two_sig,op_summary$Operon[i])
  }else if(sum_gene>=3){
    three_sig<-c(three_sig,op_summary$Operon[i])
  }
  
  
}

one_sig_two_gene<-one_sig[which(one_sig %in% two_gene)]
one_sig_three_gene<-one_sig[which(one_sig %in% three_gene)]
one_sig_six_gene<-one_sig[which(one_sig %in% six_gene)]

two_sig_two_gene<-two_sig[which(two_sig %in% two_gene)]
two_sig_three_gene<-two_sig[which(two_sig %in% three_gene)]
two_sig_six_gene<-two_sig[which(two_sig%in% six_gene)]

three_sig_three_gene<-three_sig[which(three_sig%in%three_gene)]
three_sig_six_gene<-three_sig[which(three_sig%in%six_gene)]


#Get the operon number
#Find rows in opdata that contain this operon number. Extract the sig_gene column only
#Count how many 1s are in the sig_gene column
#Check which group this operon number belongs to, 2 gene, 6+ etc. and assign to a combined group based on that and its sig gene count

op_summary$mean_exp_change<-NA
op_summary$igr<-NA

for (i in 1:nrow(op_summary)) {
  if (op_summary$Operon[i]%in% non_var){
      op_summary$igr[i]<-0
  }
  if(op_summary$Operon[i]%in% up_var| op_summary$Operon[i] %in% down_var){
    op_summary$igr[i]<-1
  }
}


mean_opdata<-function(op_summary,isolate,opdata_53930){
  for (i in unique(opdata_53930$Operon) ) {
    unique_opdata53930<-opdata_53930[which(opdata_53930$Operon==i),]
    print(as.numeric(unique_opdata53930$log_exp))
    avg <-mean(as.numeric(na.omit(unique_opdata53930$log_exp)))
    op_summary$mean_exp_change[which(op_summary$Operon==i)] <- avg
  }
  op_summary<-op_summary[which(!is.na(op_summary$igr)),]
  return(op_summary)
  ###Need to remove duplicate rows of an operon.
  
  pdf(paste("~/Documents/PhD/PhD/operon_mapper_res/Figures/graph", isolate,".pdf",sep = ""))
  mean_graph<-ggplot(op_summary, aes(x=as.factor(igr), y=(as.numeric(mean_exp_change)), color = igr)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width=0.1) +
    geom_signif(comparisons = list(c('0', '1')), test = 'wilcox.test',y_position = 5.5)+
    xlab("igr")+
    ylab("change in expression")
  
  print(mean_graph)
  dev.off()
}

op_summary_graphing <- mean_opdata(op_summary, 'igr_graph', opdata_53930)
temp <- op_summary_graphing[which(op_summary_graphing$igr==0),]
shapiro.test(as.numeric(temp$mean_exp_change))


length(up_var) #70
length(down_var) #70
length(non_var) #77
length(two_gene) #92
length(three_gene) #70
length(six_gene) #20
length(one_sig) #119
length(two_sig) #47
length(three_sig) #22
length(all_sig) #36
length(one_sig_two_gene) #67
length(one_sig_three_gene) #38
length(one_sig_six_gene) #9
length(two_sig_two_gene) #25
length(two_sig_three_gene) #17
length(two_sig_six_gene) #4
length(three_sig_three_gene) #15
length(three_sig_six_gene) #7
length(amb_gene) #6

length(op_summary_graphing$igr[which(op_summary_graphing$igr==1)]) #107

