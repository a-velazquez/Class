### Advanced R Ch. 2
# install.packages("lobstr")
library(lobstr)

x <- c(1,2,3)
# 1. creates the object
# 2. binds it to the name "x"

y <- x
# Doesn't create a copy, just creates a second binding

lobstr::obj_addr(x)
lobstr::obj_addr(y)
# change every time you restart R.

# The ability to use strings on the left hand side of the assignment arrow is an historical artefact, used before R supported backticks.

lobstr::obj_addr(mean)
lobstr::obj_addr(base::mean)
lobstr::obj_addr(get("mean"))
lobstr::obj_addr(evalq(mean))
lobstr::obj_addr(match.fun("mean"))

# all have same object address

?read.csv
?make.names

x <- c(1, 2, 3)
cat(tracemem(x), "\n")

y <- x
y[[3]] <- 4L
y[[3]] <- 5L #this DID get copied again, but the book said it wouldn't?

untracemem(x)


x <- c(1L, 2L, 3L)
tracemem(x)

x[[3]] <- 4



# Classes: logical, character

?class
# Can develop our own
# Packages can include their own class definitions, like tibbles
# as.CLASS
# is.CLASS is a test

as.character(c(1,2,3))
# Vectors can only have 1 type of object

## Special character functions
colors<-c("red","yellow","blue")
more_colors<-c(colors, "magenta", "cyan")
substr(colors, start = 1, 2)

paste0(colors,"flowers")
paste("I like",colors,"flowers")

paste("I like ",colors,"flowers", collapse = '')

# nchar
# strsplit
# unlist
# gsub - find element in vector and replace with somehthing else
extreme_statement<-"coding is my life"
gsub("my life", "the bees knees", extreme_statement)

# stringr
# "\n" return

# Logicals
x<-c(0,7,8)
chooser<-c(T,F,T)

# print only elements where chooser is true
x[chooser]



## What does "||"/"&&" do?


state.size<-(voting.age.population > median(voting.age.population))*1
state.size

state.size[state.size==1] <- rep(3,sum(state.size==1))
# Anywhere there is a value of 1, assign it the value of 3 
# repeated some number of times

# Recodes 1s into 3
# repeating is unneeded 


# which(), as in elements of the vector meet the criteria
# any(), T/F
# all(), does every element?


# Adding more complex data structures
## matrix
## list
## dataframes

vap<-voting.age.population

# What are some uses of arrays?

# matrix: collection of row or column vectors

m1<-cbind(vap, tv)

m2<-rbind(vap, tv)

m2[1,2]
m2[,1] #first column

m2[1,1:10]
m2[1:2,1:10]

# objects are pieces of data that contain attributes
# classes define what attributes an object must and can have

dim(m1)
attributes(m1)
dimnames(m1)
colnames(m1)
rownames(m1)

# When do you use "names()" vs "colnames"

colnames(m1)[2]<-"temp"

dimnames(m1)[[2]][1]

# Many R objects are lists
list.a<-list(m1, vap, 3)

# lists are flexible way to store loosely structured data 
# lets you mix objects of different classes


# names(object)<-c(1,2,3)

# 4 ways of accessing elements of a list
# my_list[[1]]
# my_list$NAME
# my_list[1] will still be a list
# my_list["NAME]


# my_list[[3][1,]] first row of matrix in list

# can't do math on lists, but can on matrices
# matrices expect everything in it to be same data type

# mean(voting.data$vap)
# with(voting.data,  mean(vap))

# %*% matrix mult
# solve(matrix, outcome_vector) gives coefficients; solves a system of linear equations





# regex cheat sheet?
# when to use repeat


x <- 10
g01 <- function() {
  x <- 20
  x
}

g01()


install.packages("sloop")

# S3 Class Systems

# What is a class for 
# Everything is an object, but not everything is object oriented
# Real objects have proper class definitions: ex. lm object 
# a class defines an object and what it should look like
# what is stored in each "slot"
# every other function defined for this class will know how to operate 

# R has many class systems
# S3 classes, really loosely structured, enforce very few rules 
# S4 more proper
# S5 is reference class system
# S6 system actually a package, python-like syntax, used in R Shiny


# S3 classes
# Take any object and call it something by assigning class() to a name 
# S3 doesn't impose any rules on who should be the member of the class 


# Methods and generics
methods(print)
methods(plot)
