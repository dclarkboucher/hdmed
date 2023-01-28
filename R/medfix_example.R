# MedFix
rm(list=ls())
load("data/med_dat.rda")
source("R/medfix.R")
source("R/medfix_utils.R")

str(med_dat)

A <- with(med_dat, as.numeric(scale(A)))
M <- with(med_dat, as.matrix(scale(M)))
Y <- with(med_dat, as.numeric(scale(Y)))

out <- mediate_medfix()

