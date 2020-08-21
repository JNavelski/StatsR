#Load packages


library(tidyverse)
library(stringr)
library(caret)
library(tm)
library(dplyr)
library(readxl)
library(e1071)
library(ranger)
library(RcppEigen)
library(ggplot2)
library(devtools)
library(varhandle)
library(MASS)

# Example with Gaussian Noise

set.seed(123)
#
N <- 100

#Assign Initial Values
beta_0 = 1
beta_1 = 1
#
X <- rnorm(n = N, mean = 10, sd = 2)
e <- rnorm(n = length(X), mean = 0, sd = 2)
Y <- beta_0 + beta_1 * X + e

log_lik <- function(par_vec, y, x) {
  # If the standard deviation prameter is negative, return a large value:
  if(par_vec[3] < 0) return(1e8)
  
  # The likelihood function values:
  # lik <- dnorm(y, mean = par_vec[1] + par_vec[2] * x, sd = par_vec[3])
  
  #This is similar to calculating the likelihood for Y - XB
  res <- y - par_vec[1] - par_vec[2] * x
  lik <- dnorm(res, mean = 0, sd = par_vec[3])
  
  # If all logarithms are zero, return a large value
  if(all(lik == 0)) return(1e8)
  
  # Logarithm of zero = -Inf
  return(-sum(log(lik[lik != 0])))
}

#
#
coef_est <- optim(par = c(0, 0, 10), fn = log_lik, hessian = T, y = Y, x = X)
print(coef_est)









# My practice of Optim Function


x=rnorm(100, mean = 0, sd = 1)
y=rnorm(100,0,1)
e=rnorm(100, mean = 0, sd = 1)
beta=1 
# beta = matrix(1,100)


SlopeInt<-function(beta){
    y = beta*x+e
    return(G)
  }

est<-optim(beta,SlopeInt,method = "L-BFGS-B")
est














