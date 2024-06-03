#Load swirl and type swirl() to begin the interactive swirl lessons
library(swirl)


#Read xlsx file
opdata <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/operon_mapper_annos.xlsx', sheet = 2)
read.csv()
#First need to change operon numbers to increment by 1 instead of the slightly random numbering 
previous_operon <- opdata$Operon[1]
new_operons <- c(1)

#For each row, check if the operon number matches that of the previous row. Set new operon numbers accordingly
for(row in 2:nrow(opdata)){
  print(row)
  current_operon <- opdata$Operon[row]
  if(current_operon == previous_operon){
    new_operons <- c(new_operons, new_operons[length(new_operons)])
  }
  if(current_operon != previous_operon){
    new_operons <- c(new_operons, (new_operons[length(new_operons)]+1))
  }
  previous_operon <- current_operon
}


#Now need to assign this to each row in the sheet
opdata$Operon <- new_operons


#Try adapting the code so that it can run for every sheet in the excel file. Can be in a function or not
opdatafunc <- function(sheetnums){
  for(sheetnum in sheetnums){
    #Read the excel file, set sheet to equal the current sheet number
    opdata <- read.xlsx()
    
    #Code to change the operon numbers here
    
    
    #Try printing the output or writing an excel file for each iteration like so
    filepath <- paste('folder1/folder2/folder3/corrected_operons', sheetnum, '.xlsx')
    write.xlsx(opdata, filepath)
  }
  
  return(opdata)
}

#Input for sheetnums if you choose to make a function to do the above
sheetnums <- 2:9
opdata_func(sheetnums)

#Some functions that can be used to check the length, number of rows, and number of columns of an object
length(new_operons)
length(opdata)
nrow(opdata)
ncol(opdata)


#A vector of the numbers 1 to 15 stored in an object called numbers
numbers <- 1:15

#A loop that prints each number in numbers
for(number in numbers){
  if(number == 1 & TRUE){
    print(number)
  }
}

#An if statement
x <- 3
if(x > 4){
  print('TRUE')
}

#A simple function that takes x and y as input and prints their sum
my_func <- function(x, y){
  out <- x+y
  print(out)
}

#Running the function with x=1 and y=2
my_func(1,2)

#For loop structure
for (variable in vector) {
  
}

#A for loop that prints the numbers of each row and the value in that cell
for(row in 1:length(opdata$IdGene)){
  print(row)
  print(opdata$IdGene[row])
}

#Making a vector of numbers 1 to 10, without saving as an object
1:10

#A string, and some manipulations of it
string <- 'hello'
substr(string, start = 1, stop = 1)
paste(string, string, sep = '')

#A numeric and its class
x <- 3
x
class(x)

#A string and its class
x <- 'string'
x
class(x)

#Converting a numeric to a string/character
x <- 3
class(x)
x <- as.character(x)
class(x)
x

#vectors and coercion
x <- c(3, 2, 1)
class(x)
x <- c('hello', 'world')
class(x)
#Here the numeric 3 is coerced to the string class.
x <- c(3, 'hello')

#A list of length 1 containing a vector of the numbers 1,2,3
x <- list(c(1, 2, 3))

#A list of length 3 containing the numbers 1 and 2, as well as the string 'hello'
x <- list(1,2,'hello')
#Lists can contain multiple different classes of object


#Modulus, finds the remainder when the lefthand side is divided by the right hand side of the operator
5 %% 2

#Exponentiation
2**2 


#Indexing vectors by position
my_vec <- c(10,20,30,40)
my_vec
my_vec[1]

#Finding the index of a vector based on the position of a value
my_vec <- c('hello', 'world')
which(my_vec == 'hello')

#Finding the index based on a value is useful for checking across tables when structured like so:
my_vec[which(my_vec == 'hello')]

#For example, here we get all the pubmlst ids of the genes in the first operon
opdata$pubmlst_id[which(opdata$Operon == 1)]



# '=' is an alternate assignment operator
x = 3
x <- 3


#lists and subsetting them
my_list <- list(10, 20, 30)
my_list[1]
my_vec[1]

list_of_lists <- list(list(10,20,30), list(40,50,60))
list_of_lists[[1]][1]

#dataframes and subsetting them
opdata[1,] #row 1
opdata[,1] #column1
opdata$Operon[1] #Entry 1 of the column titled 'Operon'


#Basic if statement structure
if (condition) {
  
}


