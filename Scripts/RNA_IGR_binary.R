#Checking the accuracy of IGR and RNA data
#P VALUES ARENT BEING READ CORRECTLY AS THEY ARE E VALUES. SO P < 0.05 ISNT WORKING. 
#Example data:
#IGR        1     1     0     1     0     0
#RNA        1     0     0     1     1     0
#MATCH?     1     0     1     1     0     1   
#SUM/COMPARISONS = 4/6
#1 in igr means the variation is different between the two isolates
#1 in rna means there is a significant logfold change between the two isolates
#there will be 28 columns, 1 for every pairwise comparison

library(ggplot2)
library(openxlsx)
library(stringr)
library(ggpubr)
library(pheatmap)
library(grid)
isolates <- c("20578", "22689", "N222.1", "N188.1", "N445.1", "N222.2", "N459.3", "N459.6")
#Approach
rna <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx", sheet = 8)
igr_down <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 1)[,c(1,4:11)]
igr_up <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 1)[,c(1,12:19)]

colnames(igr_down) <- c('Locus', '20578', '22689', 'N222.1', 'N188.1', 'N445.1', 'N222.2', 'N459.3', 'N459.6')
colnames(igr_up) <- c('Locus', '20578', '22689', 'N222.1', 'N188.1', 'N445.1', 'N222.2', 'N459.3', 'N459.6')

#rna only needs locus and pairwise expression changes and q values
rna <- rna[-c(1),c(9,13:22, 101:156)]
colnames(rna) <- c("Locus", "Min_q", "Max_change",
                   "20578", "22689", "N188.1", "N222.1", "N222.2", "N445.1", "N459.3", "N459.6",
                   "q20578_22689", "q20578_N188.1", "q20578_N222.1", "q20578_N222.2", "q20578_N445.1", "q20578_N459.3", "q20578_N459.6",
                   "q22689_20578", "q22689_N188.1", "q22689_N222.1", "q22689_N222.2", "q22689_N445.1", "q22689_N459.3", "q22689_N459.6",
                   "qN188.1_20578", "qN188.1_22689", "qN188.1_N222.1", "qN188.1_N222.2", "qN188.1_N445.1", "qN188.1_N459.3", "qN188.1_N459.6",
                   "qN222.1_20578", "qN222.1_22689", "qN222.1_N188.1", "qN222.1_N222.2", "qN222.1_N445.1", "qN222.1_N459.3", "qN222.1_N459.6",
                   "qN222.2_20578", "qN222.2_22689", "qN222.2_N188.1", "qN222.2_N222.1", "qN222.2_N445.1", "qN222.2_N459.3", "qN222.2_N459.6",
                   "qN445.1_20578", "qN445.1_22689", "qN445.1_N188.1", "qN445.1_N222.1", "qN445.1_N222.2", "qN445.1_N459.3", "qN445.1_N459.6",
                   "qN459.3_20578", "qN459.3_22689", "qN459.3_N188.1", "qN459.3_N222.1", "qN459.3_N222.2", "qN459.3_N445.1", "qN459.3_N459.6",
                   "qN459.6_20578", "qN459.6_22689", "qN459.6_N188.1", "qN459.6_N222.1", "qN459.6_N222.2", "qN459.6_N445.1", "qN459.6_N459.3"
)

count = 0
for(locus in na.omit(rna$Locus)){
  if(is.na(rna[which(rna$Locus == locus), "Min_q"][1]) | is.na(rna[which(rna$Locus == locus), "Max_change"][1])){
    next
  }
  
  if(rna[which(rna$Locus == locus), "Min_q"][1] < 0.05 && rna[which(rna$Locus == locus), "Max_change"][1] >=1){
    print(locus)
    count = count + 1
  }
}
count
test <- na.omit(rna)
#remove whitespaces at the end of rna genes
for(i in 1:length(rna$Locus)){
  rna$Locus[i] <- str_trim(rna$Locus[i])
}
locus = 'NEIS0958'
rna[which(rna$Locus == locus), "Min_q"][1] < 0.05
rna[which(rna$Locus == locus), "Max_change"][1] >=1


#Optional - subset rna to remove loci with both no variation and no significant expression change
# removables <- c()
# for(locus in na.omit(rna$Locus)){
#   print(locus)
#   
#   if(is.na(rna[which(rna$Locus == locus), 'Max_change'][1])){
#     next
#   }
#   
#   max_exp <- rna[which(rna$Locus == locus), "Max_change"][1]
#   var_up <- as.vector(igr_up[which(igr_up$Locus == locus), 2:9])
#   var_down <- as.vector(igr_down[which(igr_up$Locus == locus), 2:9])
#   
#   if(length(unique(var_up)) == 1 && length(unique(var_down)) == 1 && max_exp < 1){
#     removables <- c(removables, locus)
#   }
# }
# 
# rna <- rna[which(!(rna$Locus %in% removables)),]



#Function that calculates the log2-fold change in expression between two isolates at a given locus
logchange <- function(isolate1, isolate2, locus){
  exp1 <- rna[which(rna$Locus == locus),][[isolate1]][1]
  exp2 <- rna[which(rna$Locus == locus),][[isolate2]][1]

  exp1 <- as.numeric(exp1)
  exp2 <- as.numeric(exp2)
  
  change <- exp1/exp2
  log_change <- log2(change)
  
  return(abs(log_change))
}

logchange("N188.1", "22689", "NEIS0320")

#Function to check if a pairwise comparison has variation and a significant change in expression
varcheck <- function(isolate1, isolate2, locus){
  
  expression = 0
  var_up = 0
  var_down = 0
  var = 0
  
  check_up = TRUE
  check_down = TRUE
  
  qcol = paste('q', isolate1, '_', isolate2, sep = "")
  log_change <- logchange(isolate1, isolate2, locus)
  exp <- as.numeric(rna[which(rna$Locus == locus),][[qcol]][1])

  vcol1 <- igr_up[which(igr_up$Locus == locus),][[isolate1]]
  vcol2 <- igr_up[which(igr_up$Locus == locus),][[isolate2]]
  vcol3 <- igr_down[which(igr_down$Locus == locus),][[isolate1]]
  vcol4 <- igr_down[which(igr_down$Locus == locus),][[isolate2]]
  
  print(exp)
  print(qcol)
  if(is.na(exp)){
    return(NA)
  }
  
  if(is.na(vcol1) | is.na(vcol2)){
    check_up = FALSE
    var_up = NA
  }
  
  if(is.na(vcol3) | is.na(vcol4)){
    var_down = NA
    check_down = FALSE
  }
  
  if(check_up == FALSE && check_down == FALSE){
    return(NA)
  }
  
  if(exp < 0.05 && log_change >= 1){
    expression = 1
  }
  
  if(check_up == TRUE){
    if(vcol1 != vcol2){
      var_up = 1
    }
  }

  if(check_down == TRUE){
    if(vcol3 != vcol4){
      var_down = 1
    }
  }

  
  if(check_down == TRUE && check_up == TRUE){
    if(var_up == 1 | var_down == 1){
      var = 1
    }
  }
  
  if(is.na(var_up)){
    var = var_down
  }
  
  if(is.na(var_down)){
    var=var_up
  }
  
  return(c(expression, var_up, var_down, var))
}

x <- varcheck("N459.6", "N222.1", "NEIS1330")
x <- as.data.frame(x)
x <- t(x)
colnames(x) <- c("expression", "var_up", "var_down", "var")
rownames(x) <- c("NEIS1330")
x
x <- t(x)
isolate1 = "22689"
isolate2 = "N222.1"
locus = "NEIS1055"



#Function that builds a table based on the igr and rna matches and adds a row for whether the two match
match_check <- function(locus, isolates){
  out <- as.data.frame(matrix(data = NA, nrow = 7, ncol = 57))
  colnames(out) <- c("Category",substr(colnames(rna)[12:67], 2, nchar(colnames(rna)[12:67])))
  out$Category <- c("IGR_up", "IGR_down", "IGR_all", "RNA", "Match_up", "Match_down", "Match_total")
  print(locus)
  for(isolate1 in isolates){
    for(isolate2 in isolates){
      if(isolate1 == isolate2){
        next
      }
      check <- varcheck(isolate1, isolate2, locus)
      
      combo <- paste(isolate1, isolate2, sep='_')
      
      out[[combo]][1] <- check[2]
      out[[combo]][2] <- check[3]
      out[[combo]][3] <- check[4]
      out[[combo]][4] <- check[1]
      
      if(!is.na(out[[combo]][1])){
        out[[combo]][5] = 0
        if(out[[combo]][4] == out[[combo]][1]){
          out[[combo]][5] = 1
        }
      }
      !is.na(out[[combo]][2])
      if(!is.na(out[[combo]][2])){
        out[[combo]][6] = 0
        if(out[[combo]][4] == out[[combo]][2]){
          out[[combo]][6] = 1
        }
      }
      
      !is.na(out[[combo]][3])
      if(!is.na(out[[combo]][3])){
        out[[combo]][7] = 0
        if(out[[combo]][4] == out[[combo]][3]){
          out[[combo]][7] = 1
        }
      }
    }
  }
  
  
  return(out)
}

match_check("NEIS1364", isolates)

#Function to add a column to the output that identifies loci with a significant log2-fold expression change >= 1
sigfinder <- function(rna, out_main){
  out_main$sig_change = 'Nonsignificant'
  out_main$exp = NA
  for(locus in out_main$Locus){
    if(is.na(locus)){
      next
    }
    out_main[which(out_main$Locus == locus), "exp"][1]<- rna[which(rna$Locus==locus),"Max_change"][1]
    if(is.na(rna[which(rna$Locus==locus),3][1]) | is.na(rna[which(rna$Locus==locus),2][1])){
      next
    }
    print(locus)
    if(as.numeric(rna[which(rna$Locus==locus),3][1]) >= 1 && as.numeric(rna[which(rna$Locus==locus),2][1]) < 0.05){
      out_main[which(out_main$Locus==locus),"sig_change"] = 'Significant'
    }
  }
  return(out_main)
}


#Function to create and save a number of plots describing the data
graphing <- function(out_main, run){
  path = "~/Documents/PhD/PhD/RNA_IGR/Figures/Correct_Figs/"
  out_main$exp <- as.numeric(out_main$exp)
  #Boxplot of the data
  pdf(file = paste(path, run, "boxplot.pdf", sep =''), width = 6, height = 5)
  p <- ggplot(out_main, aes(x = sig_change, y = percent_both, color = sig_change)) +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom='point', shape = 5, size = 2,  color = 'black') +
    geom_signif(comparisons = list(c("TRUE", "FALSE")), map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) +
    geom_jitter(shape = 16, position=position_jitter(0.2)) +
    theme(axis.text = element_text(size=15), axis.title = element_text(size = 20),
          legend.text = element_text(size=15), legend.title = element_text(size = 15)) +
    theme(legend.title=element_blank())
  print(p)
  dev.off()
  
  #scatterplot of the data
  pdf(file = paste(path, run, "scatter.pdf", sep =''), width = 8, height = 7)
  p <- ggplot(out_main, aes(x = exp, y = percent_both, color = sig_change)) +
    geom_point() +
    xlab("log2 fold change (MAX/MIN)") +
    ylab("Percentage match") +
    theme(legend.title=element_blank(),
          axis.text = element_text(size=19),
          axis.title=element_text(size=24),
          legend.text = element_text(size=16))
  print(p)
  dev.off()
  
  #scatterplot of just those that are true
  pdf(file = paste(path, run, "scatter_only_sig.pdf", sep =''), width = 8, height = 7)
  p <- ggplot(out_main[which(out_main$sig_change == TRUE),], aes(x = exp, y = percent_both, color = sig_change)) +
    geom_point() +
    geom_smooth(method = lm, se=F, color = 'black') +
    xlab("log2 fold change (MAX/MIN)") +
    ylab("Percentage match") +
    theme(legend.title=element_blank(),
          axis.text = element_text(size=19),
          axis.title=element_text(size=24),
          legend.text = element_text(size=16))
  print(p)
  dev.off()
}

main <- function(rna, igr_up, igr_down, isolates, run){
  out_main <- as.data.frame(matrix(data = NA, nrow = nrow(rna), ncol = 4))
  colnames(out_main) <- c("Locus", "percent_up", "percent_down", "percent_both")
  out_main$Locus <- rna$Locus
  
  comparisons <- length(isolates)*(length(isolates)-1)
  
  for(locus in rna$Locus){
    if(is.na(locus)){
      next
    }
    if(locus %in% igr_up$Locus){
      matches <- match_check(locus, isolates)
      n_up <- length(matches[5,][which(matches[5,] == 0 | matches[5,] == 1)])
      n_down <- length(matches[6,][which(matches[6,] == 0 | matches[6,] == 1)])
      n_both <- length(matches[7,][which(matches[7,] == 0 | matches[7,] == 1)])
      
      out_main[which(out_main$Locus == locus),2] <- 100*(sum(matches[5,2:(comparisons+1)], na.rm = T)/n_up)
      out_main[which(out_main$Locus == locus),3] <- 100*(sum(matches[6,2:(comparisons+1)], na.rm = T)/n_down)
      out_main[which(out_main$Locus == locus),4] <- 100*(sum(matches[7,2:(comparisons+1)], na.rm = T)/n_both)
      
    }
  }
  out_main <- na.omit(sigfinder(rna, out_main))
  graphing(out_main, run)
  
  
  return(out_main)
}

isolates <- c("20578", "22689", "N222.1", "N188.1", "N445.1", "N222.2", "N459.3", "N459.6")
#exlcuding N188.1:
#isolates <- c("20578", "22689", "N222.1", "N445.1", "N222.2", "N459.3", "N459.6")
run <- "test"
y <- main(rna, igr_up, igr_down, isolates, run)

pheatmap(na.omit(y)[,2:4], scale = 'none', cluster_rows = T, cluster_cols = F, show_rownames = F, color= hcl.colors(50, "Blues 3"),
         angle_col = 45, fontsize_col = 13, fontsize_row = 11, display_numbers = F, fontsize_number = 9, number_color = 'black')
grid.text("Comparisons", x=-0.07, rot=90, gp=gpar(fontsize=16))



# 
# locus = "NEIS0412"
# as.numeric(rna[which(rna$Locus==locus),3][1]) >= 1
# as.numeric(rna[which(rna$Locus==locus),2][1]) < 0.05
# 
# 
# count = 0
# for(locus in na.omit(rna$Locus)){
#   if(is.na(rna[which(rna$Locus == locus), "Min_q"][1]) | is.na(rna[which(rna$Locus == locus), "Max_change"][1])){
#     next
#   }
#   
#   if(as.numeric(rna[which(rna$Locus == locus), "Min_q"][1]) < 0.05 && as.numeric(rna[which(rna$Locus == locus), "Max_change"][1]) >=1){
#     count = count + 1
#   }
# }
# count
# 
# 
# missing_loci <- c()
# count = 0
# for(locus in test$Locus){
#   if(locus %in% z$Locus){
#     if(z[which(z$Locus == locus), 5][1] == FALSE && (as.numeric(test[which(test$Locus == locus), "Min_q"][1]) < 0.05 && as.numeric(test[which(test$Locus == locus), "Max_change"][1]) >= 1)){
#       count = count + 1
#       missing_loci <- c(missing_loci, locus)
#     }
#     
#   }
#   if(!(locus %in% z$Locus) && test[which(test$Locus == locus), "Min_q"][1] < 0.05 && test[which(test$Locus == locus), "Max_change"][1] >= 1){
#     missing_loci <- c(missing_loci, locus)
#   }
# }
# length(test$Locus)
# count
# missing_loci
out_main = y


#subset rna dataset to only include loci with no variation
novar_subsettables <- c()
var_subsettables <- c()
for(locus in rna$Locus){
  if(locus %in% igr_down$Locus){
    rowdown <- as.vector(unlist(igr_down[which(igr_down$Locus == locus),2:9]))
    rowup <- as.vector(unlist(igr_up[which(igr_up$Locus == locus),2:9]))
    if(length(unique(na.omit(rowup))) == 1 && length(unique(na.omit(rowdown))) == 1){
      novar_subsettables <- c(novar_subsettables, locus)
    }
    else{
      var_subsettables <- c(var_subsettables, locus)
    }
  }
}

locus = "NEIS0371"
rna_novar <- rna[which(rna$Locus %in% novar_subsettables),]
igr_upnovar <- igr_up[which(igr_up$Locus %in% novar_subsettables),]
igr_downnovar <- igr_down[which(igr_down$Locus %in% novar_subsettables),]
main(rna_novar, igr_upnovar, igr_downnovar, isolates, run = "novar")
#subset rna dataset to only include loci with variation. 
rna_onlyvar <- rna[which(rna$Locus %in% var_subsettables),]
igr_uponlyvar <- igr_up[which(igr_up$Locus %in% var_subsettables),]
igr_downonlyvar <- igr_down[which(igr_down$Locus %in% var_subsettables),]
main(rna_onlyvar, igr_uponlyvar, igr_downonlyvar, isolates, run = "onlyvar")
#subset rna dataset to only include loci with variation only in N188.1
subsettable <- c()
for(locus in rna$Locus){
  print(locus)
  igr_down2 <- igr_down[which(igr_down$N188.1 == 2),]
  igr_up2 <- igr_up[which(igr_up$N188.1 == 2),]
  if(locus %in% igr_down2$Locus){
    print("igrdown2")
    row <- as.vector(unlist(igr_down2[which(igr_down2$Locus == locus),]))
    row_up <- as.vector(unlist(igr_up[which(igr_up$Locus == locus),]))
    if(length(unique(na.omit(row[2:9]))) == 2 && length(row[which(row[2:9] == 2)]) == 1){
      print('2')
      if(length(unique(na.omit(row_up[2:9]))) == 1){
        subsettable <- c(subsettable, locus)
      }
      if(length(unique(na.omit(row_up[2:9]))) == 2 && length(row_up[which(row_up[2:9] == 2)]) == 1){
        subsettable <- c(subsettable, locus)
      }
    }
  }
  
  else if(locus %in% igr_up2$Locus){
    print("igrup2")
    row <- as.vector(unlist(igr_up2[which(igr_up2$Locus == locus),]))
    row_down <- as.vector(unlist(igr_down[which(igr_down$Locus == locus),]))
    if(length(unique(na.omit(row[2:9]))) == 2 && length(row[which(row[2:9] == 2)]) == 1){
      print('1')
      if(length(unique(na.omit(row_down[2:9]))) == 1){
        print('2')
        subsettable <- c(subsettable, locus)
      }
      if(length(unique(na.omit(row_down[2:9]))) == 2 && length(row_down[which(row_down[2:9] == 2)]) == 1){
        print('3')
        subsettable <- c(subsettable, locus)
      }
    }
  }

  
}
subsettable
rna_N188 <- rna[which(rna$Locus %in% subsettable),]
igr_downN188 <- igr_down[which(igr_down$Locus %in% subsettable),]
igr_upN188 <- igr_up[which(igr_up$Locus %in% subsettable), ]
main(rna_N188, igr_upN188, igr_downN188, isolates, run = "N188_var_only")


#finally get results after removing the novar, and the loci that are vairable only in N188
rna_final <- rna_onlyvar[which(!(rna_onlyvar$Locus %in% subsettable)),]
igr_downfinal <- igr_downonlyvar[which(!(igr_downonlyvar$Locus %in% subsettable)),]
igr_upfinal <- igr_uponlyvar[which(!(igr_uponlyvar$Locus %in% subsettable)),]
main(rna_final, igr_upfinal, igr_downfinal, isolates, run = 'onlyvar_noN188onlyvariableloci')

y$Group = NA
for(locus in y$Locus){
  if(locus %in% novar_subsettables && y[which(y$Locus == locus),"sig_change"] == TRUE){
    y[which(y$Locus == locus), "Group"] <- 1
    next
  }
  if(locus %in% novar_subsettables && y[which(y$Locus == locus),"sig_change"] == FALSE){
    y[which(y$Locus == locus), "Group"] <- 4
    next
  }
  
  
  if(locus %in% subsettable && y[which(y$Locus == locus),"sig_change"] == TRUE){
    y[which(y$Locus == locus), "Group"] <- 2
    next
  }
  if(locus %in% subsettable && y[which(y$Locus == locus),"sig_change"] == FALSE){
    y[which(y$Locus == locus), "Group"] <- 5
    next
  }
  
  
  if(locus %in% var_subsettables && y[which(y$Locus == locus),"sig_change"] == TRUE){
    y[which(y$Locus == locus), "Group"] <- 3
    next
  }
  if(locus %in% var_subsettables && y[which(y$Locus == locus),"sig_change"] == FALSE){
    y[which(y$Locus == locus), "Group"] <- 6
    next
  }
  
  
}

y$exp <- as.numeric(y$exp)
y$Group <- as.factor(y$Group)
my_colours <- c("darkblue", "blue", "lightblue", "darkred", "red", "orange")

path = '~/Downloads/'

pdf(file = paste(path, "Summary_scatter.pdf", sep =''), width = 8, height = 5)
ggplot(y, aes(x = exp, y = percent_both, color = Group)) +
  geom_point() + 
  scale_color_manual(labels = c("Significant, No var", "Significant, N188 var only", "Significant, var", "Nonsignificant, No var", "Nonsignificant, N188 var only", "Nonsignificant, var"), values = my_colours) +
  xlab("log2 fold change (MAX/MIN)") +
  ylab("Percentage match")
dev.off()

ggplot(y, aes(x = exp, y= percent_both)) +
  geom_point()

#Extract a list of loci with 100% match and significant logfold expression change
extracts <- c()
for(i in 1:length(y$Locus)){
  if(y[i,4] == 100 && y[i,5] == TRUE){
    if(y[i,7] %in% c(2, 3)){
      extracts <- c(extracts, y$Locus[i])
    }
  }
}
extracts

z <- y[which(y$Group %in% c(1,2,3)),]
ggplot(z, aes(x = exp, y = percent_both, color = Group)) +
  geom_point() + 
  xlab("log2 fold change (MAX/MIN)") +
  ylab("Percentage match")

#Extract a table of significant loci, with their percentage match and variation recorded in columns
sighits <- y[which(y$sig_change == TRUE),]
sighits$Group <- as.numeric(sighits$Group)
sighits[which(sighits$Group == 1),7] <- "NoVar"
sighits[which(sighits$Group == 2),7] <- "N188var"
sighits[which(sighits$Group == 3),7] <- "Var"

#write.xlsx(sighits, file= '~/Documents/PhD/PhD/RNA_IGR/significant_loci_matches.xlsx')

class(sighits$Group)
for(i in nrow(sighits)){
  if(sighits[i,])
}
