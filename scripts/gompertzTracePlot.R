library(tidyverse)
library(panelPomp)
library(foreach)
library(doParallel)
library(doRNG)

registerDoParallel(25)
registerDoRNG(123456)

gomp <- panelPomp::panelGompertz(N = 50, U = 50, seed = 257297)

lpp <- c(0.01, 0.01, rep(0.02, 50))
upp <- c(0.5, 0.5, rep(0.4,  50))

names(lpp) <- c('r', 'sigma', paste0('tau[unit', 1:50, ']'))
names(upp) <- c('r', 'sigma', paste0('tau[unit', 1:50, ']'))

guesses <- pomp::runif_design(
  lower = lpp,
  upper = upp,
  nseq = 25
)

all_params <- coef(gomp)
fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]

fixed_mat <- matrix(
  rep(fixed_params, nrow(guesses)),
  byrow = TRUE, nrow = nrow(guesses)
)

colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])

# Combine estimated and fixed parameters, and reorder based on original order.
guesses_all <- cbind(guesses, fixed_mat)[, names(coef(gomp))]
gomp_rw.sd <- rw_sd(r = 0.02, sigma = 0.02, tau = 0.02)

mpif_out <- bake(
  file="../data/traceMPIF.rds", {
    foreach(i=1:25, .combine=c) %dopar% {
      mif2(
        gomp,
        Nmif = 50,
        start = unlist(guesses_all[i, ]),
        Np = 1000,
        rw.sd = gomp_rw.sd,
        cooling.type = 'geometric', 
        cooling.fraction.50 = 0.5,
        block = TRUE
      )
    }
  })