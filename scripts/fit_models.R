Sys.time()
devtools::session_info()
######### load-packages ###############################
library(foreach)

######## Source functions ############################
source("scripts/functions.R")

######## Get arguments from command line #############
(out_dir = as.character(Sys.getenv("out_dir", unset = NA)))
(array_job_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)))
(mod = Sys.getenv("mod", unset = "cohort"))
(np = as.numeric(Sys.getenv("np", unset = 4)))
(fitr = Sys.getenv("fitr", unset = "mpif"))

## ############# OPTIONS #############################
# Set number of cores
ncores = as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE", unset = NA))
if(is.na(ncores)) ncores = 2
print(ncores)

# Set fitting and filter parameters
RUN_LEVEL = 3
NP_FITR      = switch(RUN_LEVEL, 2, np/2,  np)
NFITR        = switch(RUN_LEVEL, 200,  200,  200)
NREPS_FITR   = switch(RUN_LEVEL, ncores, ncores, ncores)
NP_EVAL      = switch(RUN_LEVEL, 2, 3000,10000)
NREPS_EVAL   = switch(RUN_LEVEL, ncores, ncores, ncores)

DATA = clean_twentycities()
# Units to select from data
UNITS =  unique(panelPomp::twentycities$measles$unit)
MODEL = switch(mod,
	cohort = model_mechanics_001(shared_params = c("mu", "cohort", "alpha")),
	iota = model_mechanics_011(shared_params = "mu"),
	many = model_mechanics_011(
		shared_params = c(
			"mu", "cohort", "amplitude", "sigmaSE", "R0", "gamma", "sigma"
		)
	)
)
# Time step
DT = 1/365.25
BLOCK_MIF2 = switch(fitr, mpif = TRUE, pif = FALSE)
INTERP_METHOD = "shifted_splines"
# Cooling fraction for rw_sd.
COOLING_FRAC = 0.5

MAIN_SEED = 169566665
SIM_MODEL_SEED = array_job_id
# Add to MAIN_SEED if running array job
if(!is.na(array_job_id)){
  MAIN_SEED = MAIN_SEED + array_job_id
  print(MAIN_SEED)
}

# Use INITIAL_RW_SD to set random walk standard deviations for parameters.
DEFAULT_SD = 0.02
IVP_DEFAULT_SD = DEFAULT_SD*12
INITIAL_RW_SD = switch(mod,
	cohort = c(
  	S_0 = IVP_DEFAULT_SD,
  	E_0 = IVP_DEFAULT_SD,
  	I_0 = IVP_DEFAULT_SD,
  	R_0 = IVP_DEFAULT_SD,
  	R0 = DEFAULT_SD*0.25,
  	sigmaSE = DEFAULT_SD,
  	amplitude = DEFAULT_SD*0.5,
  	gamma = DEFAULT_SD*0.5,
  	rho = DEFAULT_SD*0.5,
  	psi = DEFAULT_SD*0.25,
  	iota = DEFAULT_SD,
  	sigma = DEFAULT_SD*0.25,
  	cohort = DEFAULT_SD*0.5,
  	alpha = DEFAULT_SD*10^(-2)*0,
  	mu = 0
	),
	iota = c(
  	S_0 = IVP_DEFAULT_SD,
  	E_0 = IVP_DEFAULT_SD,
  	I_0 = IVP_DEFAULT_SD,
  	R_0 = IVP_DEFAULT_SD,
  	R0 = DEFAULT_SD*0.25,
 		sigmaSE = DEFAULT_SD,
  	amplitude = DEFAULT_SD*0.5,
  	gamma = DEFAULT_SD*0.5,
  	rho = DEFAULT_SD*0.5,
  	psi = DEFAULT_SD*0.25,
  	iota_1 = DEFAULT_SD*0.5,
  	iota_2 = DEFAULT_SD*0.5,
  	sigma = DEFAULT_SD*0.25,
  	cohort = DEFAULT_SD*0.5,
  	mu = 0
	),
	many = c(                                                                      
    S_0 = IVP_DEFAULT_SD,                                                        
    E_0 = IVP_DEFAULT_SD,                                                        
    I_0 = IVP_DEFAULT_SD,                                                        
    R_0 = IVP_DEFAULT_SD,                                                        
    R0 = DEFAULT_SD*0.25,                                                        
    sigmaSE = DEFAULT_SD,                                                        
    amplitude = DEFAULT_SD*0.5,                                                  
    gamma = DEFAULT_SD*0.5,                                                      
    rho = DEFAULT_SD*0.5,                                                        
    psi = DEFAULT_SD*0.25,                                                       
    iota_1 = DEFAULT_SD*0.5,                                                     
    iota_2 = DEFAULT_SD*0.5,                                                     
    sigma = DEFAULT_SD*0.25,                                                     
    cohort = DEFAULT_SD*0.5,                                                     
    mu = 0                                                                       
  )
)
stopifnot(
  all(names(INITIAL_RW_SD) %in% MODEL$paramnames,
  MODEL$paramnames %in% names(INITIAL_RW_SD))
)

# Specify names of output files
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

###### Starting parameters #############################

bounds_tbl = switch(mod,
	cohort = tibble::tribble(
  	~param,       ~lower,        ~upper,
  	"R0",             10,            60,
  	"rho",           0.1,           0.9,
  	"sigmaSE",      0.04,           0.1,
  	"amplitude",     0.1,           0.6,
  	"S_0",          0.01,          0.07,
  	"E_0",      0.000004,        0.0001,
  	"I_0",      0.000003,         0.001,
  	"R_0",           0.9,          0.99,
  	"sigma",          25,           100,
  	"iota",        0.004,             3,
  	"alpha",           1,             1,
  	"psi",          0.05,             3,
  	"cohort",        0.1,           0.7,
  	"gamma",          25,           320,
  	"mu",           0.02,          0.02
	),
	iota = tibble::tribble(
  	~param,       ~lower,        ~upper,
  	"R0",             10,            60,
  	"rho",           0.1,           0.9,
  	"sigmaSE",      0.04,           0.1,
  	"amplitude",     0.1,           0.6,
  	"S_0",          0.01,          0.07,
  	"E_0",      0.000004,        0.0001,
  	"I_0",      0.000003,         0.001,
  	"R_0",           0.9,          0.99,
  	"sigma",          25,           100,
  	"iota_1",        -10,            -2,
  	"iota_2",        0.4,           1.2,
  	"psi",          0.05,             3,
  	"cohort",        0.1,           0.7,
  	"gamma",          25,           320,
  	"mu",           0.02,          0.02
	),
	many = tibble::tribble(                                                        
    ~param,       ~lower,        ~upper,                                         
    "R0",             10,            60,                                         
    "rho",           0.1,           0.9,                                         
    "sigmaSE",      0.04,           0.1,                                         
    "amplitude",     0.1,           0.6,                                         
    "S_0",          0.01,          0.07,                                         
    "E_0",      0.000004,        0.0001,                                         
    "I_0",      0.000003,         0.001,                                         
    "R_0",           0.9,          0.99,                                         
    "sigma",          25,           100,                                         
    "iota_1",        -10,            -2,                                         
    "iota_2",        0.4,           1.2,                                         
    "psi",          0.05,             3,                                         
    "cohort",        0.1,           0.7,                                         
    "gamma",          25,           320,                                         
    "mu",           0.02,          0.02                                          
  )
) 
stopifnot(setequal(bounds_tbl$param, MODEL$paramnames))
bounds_tbl = bounds_tbl |>
  dplyr::mutate(shared = param %in% MODEL$shared_params)

# Sample initial parameters and place into lists
initial_pparams_list = sample_initial_pparams_ul(
  sh_ul = dplyr::filter(bounds_tbl, shared == TRUE),
  sp_ul = dplyr::filter(bounds_tbl, shared == FALSE),
  units = UNITS,
  n_draws = NREPS_FITR
)

################## Construct panelPomp object ##########################
measlesPomp_mod = make_measlesPomp(
  data = DATA,
  model = MODEL,
  interp_method = INTERP_METHOD,
  dt = DT
)

###### MODEL FITTING #####################################
round_out = run_round(
  measlesPomp_mod,
  initial_pparams_list = initial_pparams_list,
  rw_sd_obj = make_rw_sd(INITIAL_RW_SD),
  cooling_frac = COOLING_FRAC,
  N_fitr = NFITR,
  np_fitr = NP_FITR,
  np_eval = NP_EVAL,
  nreps_eval = NREPS_EVAL,
  panel_block = BLOCK_MIF2,
  ncores = ncores,
  write_results_to = write_results_to,
  print_times = TRUE
)

EL_final = round_out$EL_out
print(as.data.frame(dplyr::arrange(EL_final$fits[,1:2], dplyr::desc(logLik))))

