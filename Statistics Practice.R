
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
t(d)




