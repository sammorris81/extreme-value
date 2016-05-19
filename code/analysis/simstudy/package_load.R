library(fields)
library(SpatialTools)
options(warn=2)

library(compiler)
enableJIT(3)

#### Load simdata
rm(list = ls())
load(file='./simdata.RData')
source('../../R/mcmc_cont_lambda.R', chdir=T)
source('../../R/auxfunctions.R')
source('./max-stab/Bayes_GEV.R')
source('./max-stab/MCMC4MaxStable.R', chdir=T)