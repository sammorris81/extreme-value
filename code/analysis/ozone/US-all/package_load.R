rm(list = ls())

library(fields)
library(SpatialTools)
library(mvtnorm)
library(compiler)
enableJIT(3)

load('us-all-setup.RData')
source('../../../R/mcmc_cont_lambda.R', chdir=T)
source('../../../R/auxfunctions.R')
source('../max-stab/MCMC4MaxStable.R', chdir=T)