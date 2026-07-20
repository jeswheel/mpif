# This script tests MPIF vs PIF on the Gompertz population model. 
# It also attempts to maximize it using a Kalman Filter. 

library(batchtools)
library(panelPomp)
library(tidyverse)
library(data.table)

# This run level controls how long the computations will take. 
RUN_LEVEL <- 3

reg <- makeExperimentRegistry(
  file.dir = paste0("Gompertz-", Sys.Date(), '-RL', RUN_LEVEL),
  seed = 123456,
  packages = c("panelPomp")
)

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
create_model <- function(data, job, U, N, data_seed, ...) {
  
  # Create Gompertz model using defaults, comparable to Breto20. 
  gomp <- panelGompertz(N = N, U = U, seed = data_seed)
  
  gomp_start <- coef(gomp)
  
  est_pars <- c(
    'r', 'sigma', paste0('tau[unit', 1:U, ']')
  )
  
  gomp_start[est_pars] <- runif(
    length(est_pars), 
    min = gomp_start[est_pars] / 2,
    max = gomp_start[est_pars] * 2
  )
  
  list(
    model = gomp, 
    starts = gomp_start
  )
}

addProblem(name = 'createGomp', fun = create_model, seed = 555555)

#' Fit Gompertz model using Kalman Filter
#'
#' @param data NULL
#' @param job Job ID
#' @param instance a list containing model (and data), and starting point. 
#' @param ... Additional arguments, required by batchtools package 
#'
#' @return a 2-d numeric vector of initial likelihood and final likelihood. 
fit_gomp_Gaussian <- function(data, job, instance, ...) {
  
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

addAlgorithm(name = 'KF', fun = fit_gomp_Gaussian)

fit_gomp_IF <- function(data, job, instance, J, M, BLOCK, COOLING, COOLING_TYPE, ...) {
  gomp_rw.sd <- rw_sd(r = 0.02, sigma = 0.02, tau = 0.02)
  
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
    block = BLOCK
  )
  
  all_traces <- apply(
    traces(mf_all), 
    1, 
    function(x) panelGompertzLikelihood(x[all_pars], instance$model, x[all_pars])
  )
  
  names(all_traces) <- paste0("ll_", 0:M)
  all_traces
}

addAlgorithm(name = "IF", fun = fit_gomp_IF)

if (RUN_LEVEL == 0) {
  U <- c(5, 10)
  N <- c(100, 500)
  J <- c(100, 500)
  M <- 20
  REPS <- 3
} else {
  U <- c(5, 15, 50, 100)
  N <- c(20, 50, 100, 200)
  J <- c(1000)
  M <- 50
  REPS <- 50
}

seeds <- sample(100000:999999, size = length(U) * length(N))

pdes <- list(
  'createGomp' = CJ(U = U, N = N) %>% mutate(data_seed = seeds)
)

ades <- list(
  'IF' = CJ(J = J, M = M, BLOCK = c(TRUE, FALSE), COOLING = 0.5, COOLING_TYPE = 'geometric'),
  'KF' = CJ()
)

addExperiments(
  prob.designs = pdes, algo.designs = ades, repls = REPS
)

getJobPars() |> unwrap() -> job_pars

tot_jobs <- (pdes$createGomp |> nrow()) * REPS * (ades$IF |> nrow()) + (pdes$createGomp |> nrow()) * REPS

speed1 <- job_pars %>% filter(U <= 15) %>% pull(job.id)
speed2 <- job_pars %>% filter(U == 50) %>% pull(job.id)
speed3 <- job_pars %>% filter(U == 100, algorithm == 'IF') %>% pull(job.id)
speed4 <- job_pars %>% filter(U == 100, algorithm != 'IF') %>% pull(job.id)

if (length(setdiff(c(speed1, speed2, speed3, speed4), 1:tot_jobs)) != 0) stop("Missing Jobs")

resources1 <- list(account = 'ACCOUNT', walltime = '1:00:00', memory = '1000m', ncpus = 1)
resources2 <- list(account = 'ACCOUNT', walltime = '2:00:00', memory = '1000m', ncpus = 1)
resources3 <- list(account = 'ACCOUNT', walltime = '3:00:00', memory = '1000m', ncpus = 1)
resources4 <- list(account = 'ACCOUNT', walltime = '10:00:00', memory = '1000m', ncpus = 1)

submitJobs(
  data.table(job.id = speed1, chunk = 1:200), resources = resources1
)

submitJobs(
  data.table(job.id = speed2, chunk = 1:200), resources = resources2
)

submitJobs(
  data.table(job.id = speed3, chunk = 1:100), resources = resources3
)

submitJobs(
  ids = speed4, resources = resources4
)

# waitForJobs()

# Once all of the jobs have been submitted, one simply needs to wait until
# each job completes. Once this is done, the final results can be loaded
# as follows:

# results <- unwrap(reduceResultsDataTable())
# pars <- unwrap(getJobPars())
# tab <- ijoin(pars, results)
# 
# saveRDS(tab, paste0("../data/Gompertz-", Sys.Date(), '-RL', RUN_LEVEL, ".rds"))
