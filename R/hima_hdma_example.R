#HIMA HDMA example
rm(list=ls())
load("data/med_dat.rda")
source("R/hima_hdma_utils.R")
source("R/hima.R")
source("R/hdma.R")
str(med_dat)

A <- with(med_dat, as.numeric(scale(A)))
M <- with(med_dat, as.matrix(scale(M)))
Y <- with(med_dat, as.numeric(scale(Y)))

#HIMA
library(ncvreg)
library(MASS)
library(foreach)
out1 <- mediate_hima(A, M, Y); out1

#HDMA
library(hdi)
library(scalreg)
library(lars)
out2 <- mediate_hdma(A, M, Y); out2
