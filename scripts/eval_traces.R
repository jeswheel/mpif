Sys.time()
devtools::session_info()
######### load-packages ###############################
library(foreach)
######## Source functions ############################
source("scripts/functions.R")

######## Get arguments from command line #############
(out_dir = as.character(Sys.getenv("out_dir", unset = NA)))
(array_job_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)))
(parent_dir = as.character(Sys.getenv("parent_dir", unset = NA)))
## ############# OPTIONS #############################
# Set number of cores
ncores = as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE", unset = NA))
if(is.na(ncores)) ncores = 2
print(ncores)

# Set filter parameters
RUN_LEVEL = 3
NP_EVAL      = switch(RUN_LEVEL, 2, 3000, 10000)
NREPS_EVAL   = switch(RUN_LEVEL, ncores, ncores, ncores)
eval_iteration = array_job_id

MAIN_SEED = 169566665
# Add to MAIN_SEED if running array job
if(!is.na(array_job_id)){
  MAIN_SEED = MAIN_SEED + array_job_id
  print(MAIN_SEED)
}

# Set PREVIOUS_FIT_PATH to NULL to choose starting parameters from a box
# instead of from a previous fit. Setting this equal to a path will nullify
# the portion of code which chooses starting parameters from a box.
(PREVIOUS_FIT_PATH = file.path(parent_dir,"fit_results_out.rds"))

# Specify names of output file
RESULTS_FILE = "fit_results_out.rds"

################## SETUP ###########################################
set.seed(MAIN_SEED)
# Create directory for output if it does not exist
if(RUN_LEVEL == 1 & is.na(out_dir)){
  write_path = "./scripts/DEFAULT_OUT/"
} else {
  write_path = out_dir
}
if(!dir.exists(write_path)) dir.create(write_path)
write_results_to = file.path(write_path, RESULTS_FILE)

fit_results_in = readRDS(PREVIOUS_FIT_PATH)

# selected trace parameters
stp = lapply(fit_results_in$fitr_out, function(x){
  z = panelPomp::traces(x) |>
    as.data.frame() |>
    dplyr::select(-loglik, -dplyr::contains("unitLoglik")) |>
    dplyr::slice(1 + eval_iteration) |>
    coef_to_pparams()
})

# prepared panelPomp list
pppl = lapply(seq_along(fit_results_in$fitr_out), function(x){
  panelPomp::panelPomp(
    object = fit_results_in$fitr_out[[x]],
    shared = stp[[x]]$shared,
    specific = as.matrix(stp[[x]]$specific)
  )
})

EL_final = pomp::bake(file = write_results_to, {
  eval_logLik(
    model_obj_list = pppl,
    ncores = ncores,
    np_pf = NP_EVAL,
    nreps = NREPS_EVAL,
    seed = NULL
  )
})
print(as.data.frame(dplyr::arrange(EL_final$fits[,1:2], dplyr::desc(logLik))))
