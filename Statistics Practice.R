
# This is an R script that follows Chapter 1 and 2 from "Introducing Monte Carlo Methods
# with R" by Christian P. Robert and George Casella.  This is meant for student that would
# like to practice statistical methods in R.

#Clears all objects from workspace (i.e. your global environment)
rm(list = ls())

#
# Need to have > R version 3.5 to use library(rtweet)
#

# Load Required Packages for Text
library(combinat) # combinatorics utilities
library(datasets) # The R Datasets Package

#Load packages
library(tidyverse)
library(stringi)
library(stringr)
library(caret)
library(tm)
library(MASS)
library(dplyr)
library(readxl)
library(e1071)
library(ranger)
library(RcppEigen)
library(ggplot2)
library(devtools)
library(varhandle)
library(rtweet)
library(tidytext)
library(textdata)

### Chapter 1: Basic R Programming

## 1.1 Introduction

# Note: help(), help.search() and help.start() are the basic functions need for general inquirary

## 1.2 Getting Started

# Note: > is called the "prompt"
# Note: Control-C will stop any action done in R

# Exercise 1.1
demo(image())
demo(graphics)

# Note: The basic packages in R that are already installed are base, stats, graphics, nmle, and lattice
# Note: Installing packages, such as the ones from this book, is easy as

# install.packages("mcsm") # already installed
# or
# download.packages("mcsm")

library(coda)
library(mcsm)

## 1.3 R Objects

# As with many advanced programming languages, 
# R distinguishes between several types of objects(scalar, vector, matrix,
# time series, data frames, functions, or graphics).
# An R object is mostly characterized by a "mode" that describes its contents
# and a "class" that describes its structure.
# The R function "str" applied to any R object, including R functions, will show its structure.

# The different modes are:
  # null (an empty object)
  # logical (TRUE or FALSE)
  # numeric (such as 3, 0.1349, 2+sqrt(3))
  # complex (such as 3-2i or complex(1,4,-2))
  # character "Blue", "binomial", "male", or "y=a+bx"
  
# The main classes are:
  # vector, matrix, array, factor, time-series, data.frame, and list

## 1.3.1 The vector class

# As indicated logically by its name, the vector object corresponds to a mathematical 
# vector of elements of the same type, such as (TRUE,TRUE,FALSE,TRUE)
# or (1,2,3,5,7,11).

# Creating a vector can be done using the R command c() (combines or concatenates terms together)
# Note: When creating a vector, by default it is a column vector!
a = c(2,6,-4,9,18)
d = c(a,b)

x <- c(3,6,9)

# Exercise 1.2 

# Instead of if (x[1]<-2), use
# if (x[1<-2])
(x[1]< -2)

# Fig. 1.1 Illustrations of the processing of vectors in R

a = c(5,5.6,1,4,-5) # build the object a containing a numeric vector of dimension 5 with elements 5, 5.6, 1, 4, –5
a[1] # display the first element of a
b = a[2:4] # build the numeric vector b of dimension 3 with elements 5.6, 1, 4
d = a[c(1,3,5)] # build the numeric vector d of dimension 3 with elements 5, 1, –5
2*a # multiply each element of a by 2 and display the result
b%%3 # provides each element of b modulo 3
d%%2.4 # computes the integer division of each element of d by 2.4
e = 3/d # build the numeric vector e of dimension 3 and elements 3/5, 3, –3/5
log(d*e) # multiply the vectors d and e term by term and transform each term into its natural logarithm
sum(d) # calculate the sum of d
length(d) # display the length of d
t(d) # transpose d, the result is a row vector
t(d)%*%e # scalar product between the row vector t(b) and the column vector e with identical length
t(d)*e # elementwise product between two vectors with identical lengths
g = c(sqrt(2),log(10)) #build the numeric vector g of dimension 2 and elements p2, log(10)
e[d==5] # build the subvector of e that contains the components e[i] such that d[i]=5
a[-3] # create the subvector of a that contains all components of a but the third.
is.vector(d) # display the logical expression TRUE if a vector and FALSE else

remove(e)

e=gamma(e^2)
log(e)

# Exercise 1.3

help(help)

# help(topic, package = NULL, lib.loc = NULL,
#      verbose = getOption("verbose"),
#      try.all.packages = getOption("help.try.all.packages"),
#      help_type = getOption("help_type"))

help("seq") # Generates regular sequences
help("sample") # Takes a sample of the specified size from the elements of x using either with or without replacement
help("order") # Returns a permutation which rearranged its first argument into ascenng order, breaking ties by further arguements. sort.list does the same, using only one arguement.

# Exercise 1.4

# Explain the difference between order and rank. See link for help: "https://towardsdatascience.com/r-rank-vs-order-753cc7665951"
help(rank) # Returns the sample ranks of the values in a vector.

# Sort actully sorts the data in increasing order.
# Rank gives the possition of that value when that vector is a sorted vector.  Note that it just lists it, and sort() actully does it.
# Order returns the position of the original value and is in the order of sorted sequence, that is smallest value to largest value.

a = c(5, 5.6, 1, 4, -5)
b = sort(a) # Actully rearenges the values in increasing order from top to bottom
c = order(a) # 
d = rank(a)
e = cbind(a,b,c,d) # Careful with merge() and join()

View(e)

help("rep") 

f = rep(a) # Replicates elements of vectors and lists
g = rep(a, times = 3) # Replicates the elements of the object the number of times specified
h = rep(a, length.out = 3) # Replicates the object to the specified length starting from the firt element
i = rep(a, each = 3) # Replicates each element in the vector how ever many times specified.

j = cbind(a,f,g,h,i) # Note that you need to be careful with the cbind!  
                     # It will just fill in values for the column it is getting comninded with
k = rbind(a,f,g,h,i) # Row bind: combines each row and puts it into a column.  This is a clever way to transpose.
View(j)
View(k)

help("cbind") # Note that cbind() just automatically fills the empty columns with elements rather than applying NAs.

# The names() function

help("names")
l = names(j)
View(l) # Weird...

# Let
n = 100
m = c(1:n-1) # Interpreted as (1:n)-1 where there are 100 elements in the sequence
o = seq(1, n-1, by = 1) # 

View(m)

p = c(3:1) # is equal to c(3,2,1)
View(p)

## 1.3.2 The Matrix, Array, and Factor classes
# 
# The MATRIX class provides the R representation of matrices.
help("matrix")

n = 100
p = 10
vec = c(NA) # Note that we need to define an object first befor we use it to turn it into a matrix.  Notice that the matrix is full of the logical oparator NA though.
View(vec)
vec = c(1:n) # Note we can create a vector of the same length of the matrix, and then make the matrix using this vector.  Notice we now hava a matrix of integers.
View(vec) 
x = matrix(vec, nrow = n, ncol = p)
# x = matrix(c(NA), nrow = n, ncol = p) # Different way of doing it.  Notice that we still need a vector to form the matrix from.
View(x)

matrix(1:4, nrow = 2, ncol = 3) # Note the warning message received even though the matrix was created.

# How to see different elements in the matrix

x
x[82,5]
vec[82] # Note that they are the same

i = 10
j = 5

x[i,j] = x[i+n*(j-1)]
x[i,j]
y = x[x>5] # Basically, notice how we lose 50 elements, and create and entire new vector of interegers.
y

z = x[1,] # Returns a vector of the first row of the matrix.
View(z)

# Matrix Operations (i.e. multipication, crossproduct, etc.)
# The standard matrix product is denoted by %*%, while * represents the term-by-term product.

A = matrix(1:10, nrow = 5, ncol = 4) # Produces a 5 by 4 matrix
B = matrix(2:11, nrow = 5, ncol = 4) # Produces a 5 by 4 matrix

View(A)
View(B)
# Term by term product
C = A*B # Gives a matrix of integers where each element is multiplied by the coresponding element in that matrix.
A
B
C

B = t(B) # Transpose object B
B
C = A%*%B # The Standard Matrix product.  Note that we have to transpose B to make it a conformable operation.

# Notice that we can also use the crossproduct() function
B = t(B)
# Where the first element will be the on transposed
t(A)%*%B
crossprod(A,B)
# AND
t(B)%*%A
crossprod(B,A) # where these two operations are equivilent s.t. t(A)%*%B = t(B)%*%A or
# crossproduct(A,B) = crossproduct(B,A)
D = crossprod(1:10^6,1:10^6)
E = c(1:10^6) # Makes an increasing integer vector
G = c(1:10^6) 
system.time(t(E)%*%G)
system.time(crossprod(1:10^6,1:10^6))
# Note the difference in time elapsed.

# Fig. 1.2 Illustrations of the processing of matracies in R

x1 = matrix(1:20, nrow = 5) # Build the numeric matrix x1 of dimension 5x4 with first row 1, 6, 11, 16
x2 = matrix(1:20, nrow = 5, byrow = TRUE) # Build the numeric matrix x2 of dimension 5 ⇥ 4 with first row 1, 2, 3, 4... says to build it across rows.
x1
x2
x3 = t(x2) # Transpose matrix x2 and define object as x3
a = x3%*%x2 # Matrix summation of x2.  Note this also equals t(x2)%*%x2, which also equals crossproduct(x2,x2)
c = x1*x2 # Term-by-term product between x1 and x2.
c
d = dim(x1)
View(d)
a[,2] # Select the second column of a
a[c(3,4),] # Select the third and forth rows of a
a[-2,] # Delete the second row of a
rbind(x1,x2) # Verticle merging of x1 and x2
cbind(x1,x2) # Horizontal merging of x1 and x2
apply(x1,1,sum) # Calculate the sum of each row of x1.  Essentilly, apply() applies the funtion you want to the specified rows that you want
help(apply)
apply(x1,2, sum) # Calculate the sum of each column of x1.

as.matrix(1:10) # Turn the vector 1:10 into a 10x1 matrix.  This is usuful for formating the data.

# Let,
w1 = c(1:10)
w1 # Notice this is a vector.
w2 = as.matrix(w)
w2 # Notice this becomes a matrix.
View(w1)
View(w2) # Note how the vector names change, and that there is now a V1 intead of the variable name being the vector dimensions.
# Changing is back
w3 = w2[,1]
w3

diag(a) # Returns the diagonal elements of matrix.  VERY USEFUL!
diag(1:10) # Creates a diagonal matrix with given elements
## IMPORTANT! ##
diag(x = 1, nrow = 100, ncol = 100) # Creates a 100x100 identity matrix.

library(base)

m = matrix(1:10, nrow = 10, ncol = 10)
View(m)

chol(m) # Returns the upper triangulat factor of the Choleski decomposition of m
eigen(m) # Returns the eigen values and eigen vectors of a square matrix
solve(m) # Need to understand more...

# The ARRAY class provides the R representation of arrays.

x = array(1:50, dim = c(2,5,5)) # Uses the data  to create 5 matracies, each with the dimension of 2x5
View(x)
x # Notice that the array has three different dimiensions now.  Rows, Columns, and now 

# The apply function used in Figure 1.2 is a very powerful device that operates 
# on arrays and, in particular, matrices. Since it can return arrays, it bypasses 
# calls to multiple loops and makes for (sometimes) quicker and (always) cleaner 
# programs. It should not be considered as a panacea, however, as apply hides 
# calls to loops inside a single command. For instance, a comparison of 
# apply(A, 1, mean) with rowMeans(A) shows the second version is about 200 times 
# faster. Using linear algebra whenever possible is therefore a more effcient solution.

# The FACTOR class provides the R representation of arrays.

# A factor is a vector of characters or integers used to specify a 
# discrete classification of the components of other vectors with the same length. 
# Its main di↵erence from a standard vector is that it comes with a level attribute 
# used to specify the possible values of the factor.

# Good for qualatative variables.

state = c("tas","tas","sa","sa","wa") # Create a vector with five values
state
statef = factor(state) # Distinguish entries by group.
statef
levels(statef) # Give the groups.
incomes=c(60, 59, 40, 42, 23) # Create a vector of incomes.
tapply(incomes, statef, mean) # Average the incomes for each group.
statef = factor(state, levels=c("tas","sa","wa","yo")) # Define a new level with one more group than observed.
table(statef) # Return statistics for all levels.

## 1.3.3 The List and Data.frame classes

# A list in R is a rather loose object made of a collection of other arbitrary objects
# known as "components".  For instance, a list can be derived from n existing onjects using the
# function:

a = list(name_1=object_1, ... , name_n=onject_n)

# Lists are very usedful in preserving information abot the values of variables used within
# R funtions in the sense that all relvent values can be put within a list that is the output
# of the corresponding function.

li = list(num = 1:5, y = "color", a = TRUE) # Create a list with three arguments
li

# Fig. 1.4 Chosen features of the list class

a = matrix(c(6, 2, 0 , 2, 6, 0, 0, 0, 36), nrow = 3) # Create a (3,3) matrix
a
res = eigen(a, symmetric = TRUE) # diagonalize a and produce a list wth two argumnts: vectors and values
res
names(res) # gives the names of the arguemnts in the list

res$vectors # returns the called vector on the list res

diag(res$values) # create a diagonal matrix using the eigenvalues

res$vec%*%diag(res$val)%*%t(res$vec) # Returns our initial matrix a using the eigenvalues and eigenvectors of a

# The local version of apply is lapply, which computes a function for each argument of the list

x = list(a = 1:10, beta = exp(-3:3), logic = c(TRUE, FALSE, FALSE, TRUE))
lapply(x,mean) # compute the empirical means
# A user friendly version of lapply is

sapply(x,mean)

# A data frame is a list whose elements are possibly made of differeing modes and attributes,
# but have the same length. as in the example provided in Figure 1.5.  A data frame can be displayed in matrix form, 
# and its rows and columns can be extracted using matrix indexing conventions. A list whose components 
# satisfy the restrictions imposed on a data frame can be coerced into a data frame using the function as.data.frame. 
# The main purpose of this object is to import data from an external file by using the read.table function.

# Fig. 1.5 Definition of a data frame

v1 = sample(1:12, 30, rep = TRUE) # Simulate 30 independent unifrom random variables on {1, 2, ... , 12}
v1                                # Note that this is a vector!
v2 = sample(LETTERS[1:10], 30, rep = TRUE) # Simulate 30 independent uniform random variables on {a, b, ..., j}
v2
v3 = runif(30) # Simulate 30 indep. uniform RVs on [0,1]
v3
v4 = rnorm(30, mean = 0, sd = 1) # Simulate 30 indep. realizations from a standard normal distribution
v4

xx = data.frame(v1, v2, v3, v4) # CREATE DATA FRAME! :}
xx
View(xx)


## 1.4 Probability Distributions in R

# A “core” name, such as norm, is associated with each distribution, 
# and the four basic associated functions, namely the cdf, the pdf, 
# the quantile function, and the simulation procedure, 
# are defined by appending the prefixes d, p, q, and r to the core name, 
# such as dnorm, pnorm, qnorm, and rnorm.

# Exercise 1.6 - Study the properties of the reslm() function using simulated data

n = 20
x = rnorm(n) # Simulate 30 indep. realizations from a std. norm. distribution
x_bar = sum(x)/n
x_bar
y = 3*x+5+rnorm(n, sd = 0.3)
y_bar = sum(y)/n
y_bar
mean(y)

# Regaression
reslm = lm(y~x)
summary(reslm)

# Basic Descriptive Statistics

mean(x)
var(x) # Note that (n-1) is used in the denominator in order to ensure unbiasedness
sd(x)
median(x)
quantile(x)

help("var")

# When the basic descriptive statistic functions are used on a matrix, the output from var(x)
# differs from the output of sd(x)^2

b = matrix(1:9, ncol= 3)
b
var(b)
sd(b)^2

# Classical hypothesis tests, such as the equality of two means or the equality of two
# variances, can be conducted using standard functions.

help.search("test")

######
# Hypothesis Test Example - Checking that the mean of a normal sample with unknown variance is zero
# can be conducted using a t-test (Casella and Berger, 2001)
######

# H_0 : mu = 0
# H_A : mu != 0

x = rnorm(25) # Produced a N(0,1) sample size of 25
t=t.test(x)
print(t)
mean(x)
t
# Note that since the p-value is not less than .05, that we can reject the null and conclude that
# the mu is significantly different than zero.

help("t.test")

# The presence of corralation between two variables can be tested
help("cor.test")

attach(faithful)
View(faithful)

cor.test(faithful[,1],faithful[,2])

# Non-parametric tests
# 1) One-sample and two-sample Kolmogorov–Smirnov adequation tests - Did the data come from a normal distribution?
help("ks.test")
ks.test(jitter(faithful[,1]),pnorm)

# 2) Shapiro’s normality test
help("shapiro.test")
shapiro.test(faithful[,2])

# 3) Kruskall–Wallis homogeneity test
help("kruskal.test")
kruskal.test(faithful[,1], faithful[,2])

# 4) Wilcoxon rank tests
help("wilcox.test")
wilcox.test(faithful[,1])

# Non-parametric Kernal Density Estimates - Natural Splines and Loess

Nit = c(0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,6,6,6)
AOB =c(4.26,4.15,4.68,6.08,5.87,6.92,6.87,6.25,6.84,6.34,6.56,6.52,7.39,7.38,7.74,7.76,8.14,7.22)
print(AOB)

AOBm=tapply(AOB,Nit,mean) #means of AOB
AOBm
Nitm=tapply(Nit,Nit,mean) #means of Nit
Nitm

library(splines)

plot(Nit,AOB,xlim=c(0,6),ylim=c(min(AOB),max(AOB)),pch=19)
fitAOB=lm(AOBm~ns(Nitm,df=2)) #natural spline
xmin=min(Nit)
xmax=max(Nit)
lines(seq(xmin,xmax,.5), predict(fitAOB,data.frame(Nitm=seq(xmin,xmax,.5)))) #fit to means
fitAOB2=loess(AOBm~Nitm,span = 1.25) #loess
lines(seq(xmin,xmax,.5), predict(fitAOB2,data.frame(Nitm=seq(xmin,xmax,.5)))) #fit to means 
fitAOB2=loess(AOBm~Nitm,span = 1.5) #loess
lines(seq(xmin,xmax,.5), predict(fitAOB2,data.frame(Nitm=seq(xmin,xmax,.5)))) #fit to means 
fitAOB2=loess(AOBm~Nitm,span = 2) #loess
lines(seq(xmin,xmax,.5), predict(fitAOB2,data.frame(Nitm=seq(xmin,xmax,.5)))) #fit to means 


x=seq(-3,3,le=5) # equidispersed regressor
y=2+4*x+rnorm(5) # simulated variable
lm(y~x)
summary = summary(lm(y~x))

plot(x,y)
lines(x,y)


