#LVMA example
library(psych)
library(Matrix)
rm(list=ls())
load("data/med_dat.rda")
source("R/lvma.R")
source("R/lvma_utils.R")
str(med_dat)

A <- with(med_dat, as.numeric(scale(A)))
M <- with(med_dat, as.matrix(scale(M)))
Y <- with(med_dat, as.numeric(scale(Y)))

# inputs
rhoLM = 2
rhoEL = 2
rhoLY = 2
q <- 4

# more iterations should be used in practice; see default value
out1 <- lvma(A, M, Y, q = q, rhoL = rhoLM, rhoE = rhoEL, rhoY = rhoLY, 50)
str(out1)

# Process out using EBIC
EBIC <- out1$val_EBIC
which_chosen <- out1$model_selectionEBIC
penalty <- c(rhoLM[which_chosen[1]], rhoEL[which_chosen[2]],
             rhoLY[which_chosen[3]])
names(penalty) <- c("rhoLM", "rhoEL", "rhoLY")
LM_effects <- out1$EBIC_sol$Lambda
colnames(LM_effects) <- paste0("L", 1:q)
rownames(LM_effects) <- colnames(M)

LY_effects <- out1$EBIC_sol$betaFY
names(LY_effects) <- paste0("L", 1:q)

AL_effects <- as.vector(out1$EBIC_sol$betaEF)
names(AL_effects) <- paste0("L", 1:q)

AY_direct_effect <- out1$EBIC_sol$betaEY

mediator_active <-
  suppressWarnings(apply(t(LM_effects) * (AL_effects * LY_effects), 2, any)) |>
  as.numeric()
table(mediator_active)

output <-
  list(
    EBIC = EBIC,
    penalty = penalty,
    AL_effects = AL_effects,
    LM_effects = LM_effects,
    AY_direct_effect = AY_direct_effect,
    LY_effects = LY_effects,
    mediator_active
  )

# Binary example
Covars <- matrix(rnorm(200),ncol = 2)
Y1 <- ifelse(Y > mean(Y), 1, 0)
out2 <- lvma(A, M, Y = Y1, C = Covars, q = q, binary_y = T, rhoL = 4,
             rhoE = 4, rhoY = 4, 5000)






