rm(list = ls())

library(fields)
library(SpatialTools)
library(mvtnorm)
library(Rcpp)
library(compiler)
enableJIT(3)

load('us-all-setup.RData')
source('~/repos-git/usefulR/usefulfunctions.R', chdir = T)
source('../../../R/mcmc_cont_lambda.R', chdir=T)
source('../../../R/auxfunctions.R')
source('../max-stab/MCMC4MaxStable.R', chdir=T)

if (Sys.info()["nodename"] == "cwl-mth-sam-001") {
  # setMKLthreads(1)
  openblas.set.num.threads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
  # set number of threads to use
  openblas.set.num.threads(1)
}
