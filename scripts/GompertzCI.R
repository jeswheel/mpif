# This script tests MPIF vs PIF on the Gompertz population model. 
# It also attempts to maximize it using a Kalman Filter. 

library(batchtools)
library(panelPomp)
library(tidyverse)
library(data.table)

# This run level controls how long the computations will take. 
RUN_LEVEL <- 3

nprof <- switch(RUN_LEVEL, 2,  15, 30)
nseq  <- switch(RUN_LEVEL, 10, 20, 30)

reg <- makeExperimentRegistry(
  file.dir = paste0("Gompertz-CI-", Sys.Date(), '-RL', RUN_LEVEL),
  seed = 123456,
  packages = c("panelPomp")
)

# Create Gompertz model using defaults, comparable to Breto20.
# Seed was chosen to match the U = 50, N = 50 model from Figure 3 of the
# manuscript.
gomp <- panelPomp::panelGompertz(N = 50, U = 50, seed = 257297)

final_pars <- coef(gomp)

prof_params <- c("r", "sigma")
prof_vars <- c()
for (pp in prof_params) {
  if (pp == 'r') {
    prof_values <- seq(0.025, 0.18, length.out = nseq)
  } else if (pp == 'sigma') {
    prof_values <- seq(0.065, 0.14, length.out = nseq)
  }
  
  prof_cols <- matrix(rep(prof_values, each = nprof), ncol = 1)
  colnames(prof_cols) <- pp
  
  lpp <- rep(0.02, 50)
  upp <- rep(0.4,  50)
  
  names(lpp) <- paste0('tau[unit', 1:50, ']')
  names(upp) <- paste0('tau[unit', 1:50, ']')
  
  if (pp == 'r') {
    lpp <- c('sigma' = 0.05, lpp)
    upp <- c('sigma' = 0.2, upp)
  } else {
    lpp <- c('r' = 0.01, lpp)
    upp <- c('r' = 0.5, upp)
  }
  
  guesses_tmp <- pomp::runif_design(
    lower = lpp,
    upper = upp,
    nseq = nprof * length(prof_values)
  )
  
  guesses <- dplyr::bind_cols(prof_cols, guesses_tmp)
  all_params <- coef(gomp)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]
  
  fixed_mat <- matrix(
    rep(fixed_params, nprof * length(prof_values)),
    byrow = TRUE, nrow = nprof * length(prof_values)
  )
  
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])
  
  # Combine estimated and fixed parameters, and reorder based on original order.
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(gomp))]
  final_pars <- rbind(final_pars, guesses_all)
  prof_vars <- c(prof_vars, rep(pp, nrow(guesses_all)))
}

# First row is true generating parameter, remove these.
final_pars <- final_pars[-1, names(coef(gomp))]
data_obj <- list(starts = final_pars, prof_vars = prof_vars)


#' Function that creates a model and starting parameters for that model. 
#'
#' @param data NULL
#' @param job Job ID
#' @param U Number of units
#' @param N Number of observations per unit
#' @param data_seed random seed to generate the data 
#'    (needed to allow for fixed model + data but varying starting values. )
#' @param ... Additional parameters required by the batchtools package
#'
#' @return a list containing the model (and data), and starting parameter values. 
create_model <- function(data, job, i, U = 50, N = 50, data_seed = 257297, ...) {
  
  # Create Gompertz model using defaults, comparable to Breto20.
  # Seed was chosen to match the U = 50, N = 50 model from Figure 3 of the
  # manuscript.
  gomp <- panelGompertz(N = N, U = U, seed = data_seed)
  
  

  list(
    model = gomp,
    prof_var = data$prof_vars[i],
    starts = unlist(data$starts[i, ])
  )
}

addProblem(name = 'createGomp', fun = create_model, seed = 555555, data = data_obj)

#' Fit Gompertz model using Kalman Filter
#'
#' @param data NULL
#' @param job Job ID
#' @param instance a list containing model (and data), and starting point.
#' @param ... Additional arguments, required by batchtools package
#'
#' @return a 2-d numeric vector of initial likelihood and final likelihood.
fit_gomp_Gaussion <- function(data, job, instance, ...) {

  # Get the starting parameters from the instance
  start <- instance$starts

  # Get the model from the instance
  gomp <- instance$model

  # Set which parameters to be estimated
  est_pars <- c(
    'r', 'sigma', paste0('tau[unit', 1:length(gomp), ']')
  )

  # Use optim to fit the model to data.
  out <- optim(
    par = start[est_pars],
    fn = panelGompertzLikelihood,
    panelPompObject = gomp,
    params = start,
    hessian = FALSE,
    control = list(trace = 0, fnscale = -1, maxit = 1000),
    method = 'L-BFGS-B',
    lower = c(1e-32, 1e-32, rep(1e-32, length(gomp)))
  )

  # Return the likelihood values from start and finish.
  ll_0 <- panelGompertzLikelihood(start, gomp, start)
  ll_1 <- out$value
  c(ll_0 = ll_0, ll_1 = ll_1, convergence = out$convergence)
}

fit_gomp_IF <- function(data, job, instance, J, M, BLOCK, COOLING, COOLING_TYPE, ...) {
  
  if (instance$prof_var == 'r') {
    gomp_rw.sd <- rw_sd(sigma = 0.015, tau = 0.02)
  } else {
    gomp_rw.sd <- rw_sd(r = 0.02, tau = 0.02)
  }
  
  start <- instance$starts
  
  all_pars <- names(start)
  
  # ll0 <- panelGompertzLikelihood(
  #   start, instance$model, start
  # )
  
  mf_all <- mif2(
    instance$model,
    Nmif = M,
    start = start,
    Np = J,
    rw.sd = gomp_rw.sd,
    cooling.type = COOLING_TYPE, 
    cooling.fraction.50 = COOLING,
    block = TRUE
  )
  
  pf_results <- replicate(
    3,
    unitLogLik(pfilter(mf_all, Np = 5000))
  ) |> 
    panel_logmeanexp(
      MARGIN = 1, se = TRUE
    )
 
  
  # all_traces <- apply(
  #   traces(mf_all), 
  #   1, 
  #   function(x) panelGompertzLikelihood(x[all_pars], instance$model, x[all_pars])
  # )
  
  est_pars <- c(
    'r', 'sigma', paste0('tau[unit', 1:length(instance$model), ']')
  )
  
  c(
    coef(mf_all)[est_pars],
    # r = unname(coef(mf_all)['r']),
    # sigma = unname(coef(mf_all)['sigma']),
    pf_ll = unname(pf_results[1]),
    pf_se = unname(pf_results[2]),
    logLik = panelGompertzLikelihood(coef(mf_all), instance$model, coef(mf_all))
  )
}

addAlgorithm(name = "IF", fun = fit_gomp_IF)

if (RUN_LEVEL == 0) {
  U <- c(5, 10)
  N <- c(100, 500)
  J <- c(100, 500)
  M <- 20
  REPS <- 3
} else {
  U <- c(50)
  N <- c(50)
  J <- c(1000)
  M <- 50
  REPS <- 1
}

pdes <- list(
  'createGomp' = CJ(U = U, N = N, i = 1:nrow(final_pars))
)

ades <- list(
  'IF' = CJ(J = J, M = M, BLOCK = c(TRUE), COOLING = 0.5, COOLING_TYPE = 'geometric')
)

addExperiments(
  prob.designs = pdes, algo.designs = ades, repls = REPS
)

getJobPars() |> unwrap() -> job_pars

submitJobs(
  data.table(job.id = 1:(2*nprof*nseq))
)

waitForJobs()

# Once all of the jobs have been submitted, one simply needs to wait until
# each job completes. Once this is done, the final results can be loaded
# as follows:

results <- unwrap(reduceResultsDataTable())
pars <- unwrap(getJobPars())
tab <- ijoin(pars, results)
# 
saveRDS(tab, paste0("data/Gompertz-CI-", Sys.Date(), '-RL', RUN_LEVEL, ".rds"))

instance <- create_model(data = data_obj, job = NULL, i = 1)

est_pars <- c(
  'r', 'sigma', paste0('tau[unit', 1:length(gomp), ']')
)

gomp <- instance$model

# Use optim to fit the model to data. 
out <- optim(
  par = coef(instance$model)[est_pars],
  fn = panelGompertzLikelihood,
  panelPompObject = gomp,
  params = coef(instance$model),
  hessian = FALSE,
  control = list(trace = 0, fnscale = -1, maxit = 1000),
  method = 'L-BFGS-B', 
  lower = c(1e-32, 1e-32, rep(1e-32, length(gomp)))
)

saveRDS(out, '../data/Gompterz50-mle.rds')

