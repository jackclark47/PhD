#R Tutorial

#1. Data types

#1.1 Strings
'this is a string'
'11232' 

paste('one', 'two') #Join strings together
substr('abcdef', start = 1, stop = 2) #extract parts of a string

#1.2 Numerics
1
3.142

1 + 4
4/2
5%%2 #Gets the remainder

class(1) #This function, class(), tells us the class of whatever we put inside the brackets
class(43730.01)

#1.3 Vectors and lists
my_numeric_vector <- c(1,2,3) #A vector of numerics, assigned using '<-' to the variable my_numeric_vector
class(my_vector)

my_string_vector <- c('string1', 'string2') #A vector of strings assigned to a new variable
class(my_new_vector)

my_mixed_vector <- c(4, 'string') #A vector of a numeric and a string
class(c(4, 'string')) #the presence of a string in this mixed vector coerces the numeric to take the form of a string.

my_list <- list(1,2,3,4)
class(my_list)

#1.4 Manipulating vectors and lists
my_list <- list(1,2,3,4)
my_list[3] #Select third element of the list
my_list[[3]] #the same but slightly different

my_vector <- c('stringone', '2', 'stringthree')
my_vector[2] #select 2nd element of the vector
my_vector[3] #select 3rd element

my_vector <- c(my_vector, 'stringfour') #We can add elements to the end of vectors like so
my_vector
my_vector <- my_vector[-4] #And remove them 
my_vector

#1.5 Dataframes
#This is the operon data. When you load a table into R it is usually stored as a data frame
library(openxlsx)
operon_data <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/operon_mapper_annos.xlsx', sheet = 2)
class(operon_data)

colnames(operon_data) #This function shows us the column names and allows us to change them
colnames(operon_data)[9] <- "Pubmlst_id"

operon_data$IdGene #You can select one column in the dataframe using its name
operon_data[,2] #Or by using the same indexing method as before
operon_data[1,3] #But remember dataframes are two dimensional. Here we select the first row, third column

class(operon_data$IdGene) #Each column is stored as a vector of a certain class
class(operon_data$Operon)
test <- operon_data$Operon



#Practice
#E1. Create two strings, S1 and S2, taking the values 'Hello' and 'world', and using the paste() function, combine them to create S3


#E2. Using the substr() function, extract the string 'o w' from S3


#E3. Extract the COGgene column from the operon_data object and save it as an object


#E4. Reduce the COGgene object to the first 100 entries


#E4 Reduce the COGgene object to the last 100 entries


#E5 Create a list, L1, of 3 vectors, vec1, vec2 and vec3


#E6 Remove the second vector from L1



#2 Programming structures
#2.1 If/else statements
x <- 4

if(x < 5 ){
  
  print('x is smaller than 5')
  
}

#if (condition) {

  #Code block   

#}
x <- 5
if(x < 4){
  print('x is less than 4')
} else if(x > 4){
  print('x is greater than 4')
} else{
  print('x is 4')
}

#2.2 For loops

# for (variable in vector) {
#   #Code block
# }

my_vector <- c('a', 'b', 'c', 'd')
for(i in my_vector){
  print(i)
}

for(i in 1:length(my_vector)){
  print(paste('i is equal to', i))
  print(paste('index', i, 'of my_vector is equal to', my_vector[i]))
}

#i simply takes the form of an entry in whatever vector is passed to it
1:length(my_vector) #here the vector is not my_vector, but a new one created between 1 and the length of my_vector


#2.3Flow control
for(i in 1:10){
  if(i == 4){
    next #skip this iteration of the loop
  }
  if(i == 8){
    break #stop the loop from running this iteration or any additional iterations.
  }
  print(i)
}


#2.4 Functions
my_function <- function(argument1){
  print(argument1)
}
my_function('Hi')

my_complex_function <- function(start, stop, divisor){
  numbers <- c() #make an empty vector that can be filled as the loop progresses
  for(i in start:stop){
    if(i %% divisor != 0){
      print(paste(i, '/', divisor, 'has a remainder, skipping...'))
      next
    }
    print(paste(divisor, 'is a factor of', i))
    numbers <- c(numbers, i)
  }
  return(numbers) #Return allows us to refernce an object created within the function after it finishes running
}

output <- my_complex_function(1,10,2) #Find the numbers between 1 and 10 that have 2 as a factor
output

