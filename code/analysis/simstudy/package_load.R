rm(list = ls())

library(fields)
library(SpatialTools)
library(Rcpp)

library(compiler)
enableJIT(3)

#### Load simdata
load(file='./simdata.RData')
source('../../../../usefulR/usefulfunctions.R', chdir = TRUE)
source('../../R/mcmc_cont_lambda.R', chdir = TRUE)
source('../../R/auxfunctions.R')
source('./max-stab/Bayes_GEV.R')
source('./max-stab/MCMC4MaxStable.R', chdir = TRUE)

if (Sys.info()["nodename"] == "cwl-mth-sam-001") {
  openblas.set.num.threads(1)
}

options(warn=2)