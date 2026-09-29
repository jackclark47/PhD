#Script to automate downloading isolates from pubmlst
library(httr)
library(seqinr)

GET("https://rest.pubmlst.org/db/pubmlst_neisseria_isolates/isolates/1/contigs_fasta?header=original_designation")

get_isolates <- function(idlist, out_dir){
  url0 <- "https://rest.pubmlst.org/db/pubmlst_neisseria_isolates/isolates/"
  url1 <- "/contigs_fasta?header=original_designation"
  for(id in idlist){
    
    url_final <- paste(url0, id, url1, sep ='')
    outfile <- paste(out_dir, '/', id, '.fasta', sep = '')
    
    file <- GET(url_final)
    cat(content(file, as = 'text'), "\n", file = outfile)
    

    
  }
  return(file)
}

idlist <- read.xlsx('~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/B/cc41_44/MenB_cc41_44_dataset.xlsx')[,1]
#idlist <- (1)
outdir <- '~/Documents/PhD/Asides/Disease_carriage_classifier/Datasets/serogroups/B/cc41_44/IGRs/Annotation/test_fastas'
get_isolates(idlist, outdir)
