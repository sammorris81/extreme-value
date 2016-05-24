rm(list = ls())

library(fields)
library(SpatialTools)

# library(compiler)
# enableJIT(3)

#### Load simdata
load(file='./simdata.RData')
source('../../R/mcmc_cont_lambda.R', chdir=T)
source('../../R/auxfunctions.R')
source('./max-stab/Bayes_GEV.R')
source('./max-stab/MCMC4MaxStable.R', chdir=T)

options(warn=2)