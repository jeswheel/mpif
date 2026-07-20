# Code included here is copied over from a personal package that isn't ready for
# release.

#' He10 POMP model for UK measles
#'
#' @return List of objects required for `make_measlesPomp`.
#' @param shared_params Character vector of parameters to be treated as shared.
#'
model_mechanics_001 = function(
    shared_params = "mu"
){
  rproc <- pomp::Csnippet("
    double beta, br, seas, foi, dw, births;
    double rate[6], trans[6];

    // cohort effect
    if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
      br = cohort*birthrate/dt + (1-cohort)*birthrate;
    else
      br = (1.0-cohort)*birthrate;

    // term-time seasonality
    t = (t-floor(t))*365.25;
    if ((t>=7 && t<=100) ||
        (t>=115 && t<=199) ||
        (t>=252 && t<=300) ||
        (t>=308 && t<=356))
        seas = 1.0+amplitude*0.2411/0.7589;
    else
        seas = 1.0-amplitude;

    // transmission rate
    beta = R0*seas*(1.0-exp(-(gamma+mu)*dt))/dt;

    // expected force of infection
    foi = beta*pow(I+iota,alpha)/pop;

    // white noise (extrademographic stochasticity)
    dw = rgammawn(sigmaSE,dt);

    rate[0] = foi*dw/dt;  // stochastic force of infection
    rate[1] = mu;         // natural S death
    rate[2] = sigma;      // rate of ending of latent stage
    rate[3] = mu;         // natural E death
    rate[4] = gamma;      // recovery
    rate[5] = mu;         // natural I death

    // Poisson births
    births = rpois(br*dt);

    // transitions between classes
    reulermultinom(2,S,&rate[0],dt,&trans[0]);
    reulermultinom(2,E,&rate[2],dt,&trans[2]);
    reulermultinom(2,I,&rate[4],dt,&trans[4]);

    S += births   - trans[0] - trans[1];
    E += trans[0] - trans[2] - trans[3];
    I += trans[2] - trans[4] - trans[5];
    R = pop - S - E - I;
    W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
    C += trans[4];           // true incidence
  ")
  
  dmeas <- pomp::Csnippet("
    double m = rho*C;
    double v = m*(1.0 - rho + psi*psi*m);
    double tol = 1.0e-18; // 1.0e-18 in He10 model; 0.0 is 'correct'
    if(ISNA(cases)) {lik = 1;} else {
        if (C < 0) {lik = 0;} else {
          if (cases > tol) {
            lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) -
              pnorm(cases - 0.5 , m, sqrt(v) + tol, 1, 0) + tol;
          } else {
            lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) + tol;
          }
        }
      }
    if (give_log) lik = log(lik);
  ")
  
  rmeas <- pomp::Csnippet("
    double m = rho*C;
    double v = m*(1.0-rho+psi*psi*m);
    double tol = 1.0e-18; // 1.0e-18 in He10 model; 0.0 is 'correct'
    cases = rnorm(m,sqrt(v)+tol);
    if (cases > 0.0) {
      cases = nearbyint(cases);
    } else {
      cases = 0.0;
    }
  ")
  
  rinit <- pomp::Csnippet("
    double m = pop/(S_0+E_0+I_0+R_0);
    S = nearbyint(m*S_0);
    E = nearbyint(m*E_0);
    I = nearbyint(m*I_0);
    R = nearbyint(m*R_0);
    W = 0;
    C = 0;
  ")
  
  pt <- pomp::parameter_trans(
    log = c("sigma","gamma","sigmaSE","psi","R0", "mu", "alpha", "iota"),
    logit = c("cohort","amplitude", "rho"),
    barycentric = c("S_0","E_0","I_0","R_0")
  )
  
  paramnames = c("R0","mu","sigma","gamma","alpha","iota", "rho",
                 "sigmaSE","psi","cohort","amplitude",
                 "S_0","E_0","I_0","R_0")
  states = c("S", "E", "I", "R", "W", "C")
  
  if(!all(shared_params %in% paramnames)){
    stop(
      "At least one parameter name given to shared_params is not in the model.",
      call. = FALSE
    )
  }
  panel_mechanics(
    rproc = rproc,
    dmeas = dmeas,
    rmeas = rmeas,
    rinit = rinit,
    pt = pt,
    shared_params = shared_params,
    specific_params = setdiff(paramnames, shared_params),
    states = states
  )
}

#' panelPOMP model with log-log relationship between iota and the
#' standardized 1950 population
#'
#' @name model_mechanics_011
#' @param shared_params Character vector of parameters to be treated as shared.
#'
#' @return List of objects required for `make_measlesPomp()`.
#'
model_mechanics_011 = function(shared_params = "mu"){
  rproc <- pomp::Csnippet("
    double beta, br, seas, foi, dw, births;
    double rate[6], trans[6];

    // Population-varying parameters
    double iota = exp(iota_2*std_log_pop_1950 + iota_1);

    // cohort effect
    if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
      br = cohort*birthrate/dt + (1-cohort)*birthrate;
    else
      br = (1.0-cohort)*birthrate;

    // term-time seasonality
    t = (t-floor(t))*365.25;
    if ((t>=7 && t<=100) ||
        (t>=115 && t<=199) ||
        (t>=252 && t<=300) ||
        (t>=308 && t<=356))
        seas = 1.0+amplitude*0.2411/0.7589;
    else
        seas = 1.0-amplitude;

    // transmission rate
    beta = R0*seas*(1.0-exp(-(gamma+mu)*dt))/dt;

    // expected force of infection
    foi = beta*(I+iota)/pop;

    // white noise (extrademographic stochasticity)
    dw = rgammawn(sigmaSE,dt);

    rate[0] = foi*dw/dt;  // stochastic force of infection
    rate[1] = mu;         // natural S death
    rate[2] = sigma;      // rate of ending of latent stage
    rate[3] = mu;         // natural E death
    rate[4] = gamma;      // recovery
    rate[5] = mu;         // natural I death

    // Poisson births
    births = rpois(br*dt);

    // transitions between classes
    reulermultinom(2,S,&rate[0],dt,&trans[0]);
    reulermultinom(2,E,&rate[2],dt,&trans[2]);
    reulermultinom(2,I,&rate[4],dt,&trans[4]);

    S += births   - trans[0] - trans[1];
    E += trans[0] - trans[2] - trans[3];
    I += trans[2] - trans[4] - trans[5];
    R = pop - S - E - I;
    W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
    C += trans[4];           // true incidence
  ")
  
  dmeas <- pomp::Csnippet("
    double m = rho*C;
    double v = m*(1.0-rho+psi*psi*m);
    double tol = 1.0e-18; // 1.0e-18 in He10 model; 0.0 is 'correct'
    if(ISNA(cases)) {lik=1;} else {
        if (C < 0) {lik = 0;} else {
          if (cases > tol) {
            lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-
              pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
          } else {
            lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
          }
        }
      }
    if (give_log) lik = log(lik);
  ")
  
  rmeas <- pomp::Csnippet("
    double m = rho*C;
    double v = m*(1.0-rho+psi*psi*m);
    double tol = 1.0e-18; // 1.0e-18 in He10 model; 0.0 is 'correct'
    cases = rnorm(m,sqrt(v)+tol);
    if (cases > 0.0) {
      cases = nearbyint(cases);
    } else {
      cases = 0.0;
    }
  ")
  
  rinit <- pomp::Csnippet("
    double m = pop/(S_0+E_0+I_0+R_0);
    S = nearbyint(m*S_0);
    E = nearbyint(m*E_0);
    I = nearbyint(m*I_0);
    R = nearbyint(m*R_0);
    W = 0;
    C = 0;
  ")
  
  pt <- pomp::parameter_trans(
    log = c("sigmaSE", "R0", "mu", "psi", "sigma", "gamma"),
    logit = c("cohort", "amplitude", "rho"),
    barycentric = c("S_0", "E_0", "I_0", "R_0")
  )
  
  paramnames = c("R0","mu","rho","sigmaSE","cohort","amplitude",
                 "S_0","E_0","I_0","R_0", "gamma", "psi", "iota_2", "iota_1",
                 "sigma")
  full_shared_params = union(
    shared_params, c("iota_2", "iota_1")
  )
  states = c("S", "E", "I", "R", "W", "C")
  
  if(!all(shared_params %in% paramnames)){
    stop(
      "At least one parameter name given to shared_params is not in the model.",
      call. = FALSE
    )
  }
  
  panel_mechanics(
    rproc = rproc,
    dmeas = dmeas,
    rmeas = rmeas,
    rinit = rinit,
    pt = pt,
    shared_params = full_shared_params,
    specific_params = setdiff(paramnames, full_shared_params),
    states = states
  )
}

#' Make a list containing necessary specifications for a panelPomp model
#'
#' @name panel_mechanics
#' @param rproc Csnippet that simulates process step.
#' @param dmeas Csnippet that calculates conditional log likelihood.
#' @param rmeas Csnippet that simulates measurement.
#' @param rinit Csnippet that simulates initial states.
#' @param pt Csnippet specifiying how to transform parameters onto the
#'   estimation scale.
#' @param shared_params Character vector of shared parameter names.
#' @param specific_params Character vector of unit-specific parameter names.
#' @param states Character vector of unobserved states.
#'
#' @return The arguments in list form with class `panel_mechanics`.
#'

new_panel_mechanics = function(
    rproc,
    dmeas,
    rmeas,
    rinit,
    pt,
    shared_params,
    specific_params,
    states
){
  out = list(
    rproc = rproc,
    dmeas = dmeas,
    rmeas = rmeas,
    rinit = rinit,
    pt = pt,
    shared_params = shared_params,
    specific_params = specific_params,
    paramnames = c(shared_params, specific_params),
    states = states
  )
  structure(out, class = "panel_mechanics")
}

validate_panel_mechanics = function(x){
  if(
    !all(c("rproc", "dmeas", "rmeas", "rinit", "pt", "paramnames", "states") %in% names(x))
  ){
    stop(
      "`x` must have names 'rproc', 'dmeas', 'rmeas', 'rinit', 'pt',
      'shared_params', 'specific_params', 'paramnames', and 'states'",
      call. = FALSE
    )
  }
  combined_params = c(x$shared_params, x$specific_params)
  if(!setequal(combined_params, x$paramnames)){
    stop(
      "The union of shared_params and specific_params is not paramnames.",
      call. = FALSE
    )
  }
  if(any(x$specific_params %in% x$shared_params)){
    stop(
      "The intersection of specific_params and shared_params is not empty.",
      call. = FALSE
    )
  }
  x
}

#' @rdname panel_mechanics
panel_mechanics = function(
    rproc,
    dmeas,
    rmeas,
    rinit,
    pt,
    shared_params,
    specific_params,
    states
){
  x = new_panel_mechanics(
    rproc,
    dmeas,
    rmeas,
    rinit,
    pt,
    shared_params,
    specific_params,
    states
  )
  validate_panel_mechanics(x)
}


#' Clean the `twentycities` data from panelPomp
#'
#' @return Returns `twentycities` but with tweaked observations.
clean_twentycities = function(){
  ukm = panelPomp::twentycities
  
  #            unit       date cases
  # 13770 Liverpool 1955-11-04    10
  # 13771 Liverpool 1955-11-11    25
  # 13772 Liverpool 1955-11-18   116
  # 13773 Liverpool 1955-11-25    17
  # 13774 Liverpool 1955-12-02    18
  
  ukm$measles[
    ukm$measles$unit == "Liverpool" &
      ukm$measles$date == "1955-11-18",
    "cases"
  ] = NA
  
  # 13950 Liverpool 1959-04-17   143
  # 13951 Liverpool 1959-04-24   115
  # 13952 Liverpool 1959-05-01   450
  # 13953 Liverpool 1959-05-08    96
  # 13954 Liverpool 1959-05-15   157
  
  ukm$measles[
    ukm$measles$unit == "Liverpool" &
      ukm$measles$date == "1959-05-01",
    "cases"
  ] = NA
  
  # 19552 Nottingham 1961-08-18     6
  # 19553 Nottingham 1961-08-25     7
  # 19554 Nottingham 1961-09-01    66
  # 19555 Nottingham 1961-09-08     8
  # 19556 Nottingham 1961-09-15     7
  
  ukm$measles[
    ukm$measles$unit == "Nottingham" &
      ukm$measles$date == "1961-09-01",
    "cases"
  ] = NA
  
  # London 1955-08-12   124
  # London 1955-08-19    82
  # London 1955-08-26     0
  # London 1955-09-02    58
  # London 1955-09-09    38
  
  ukm$measles[
    ukm$measles$unit == "London" &
      ukm$measles$date == "1955-08-26",
    "cases"
  ] = NA
  # The value 76 was used in He10, but it seems safer to use NA.
  
  # Sheffield 1961-05-05   266
  # Sheffield 1961-05-12   346
  # Sheffield 1961-05-19     0
  # Sheffield 1961-05-26   314
  # Sheffield 1961-06-02   297
  
  ukm$measles[
    ukm$measles$unit == "Sheffield" &
      ukm$measles$date == "1961-05-19",
    "cases"
  ] = NA
  
  # Hull 1956-06-22    72
  # Hull 1956-06-29    94
  # Hull 1956-07-06     0
  # Hull 1956-07-13    91
  # Hull 1956-07-20    87
  
  ukm$measles[
    ukm$measles$unit == "Hull" &
      ukm$measles$date == "1956-07-06",
    "cases"
  ] = NA
  
  # 1 Hornsey 1957-01-27    51
  # 2 Hornsey 1957-02-03    70
  # 3 Hornsey 1957-02-10    88
  # 4 Hornsey 1957-02-17     8
  # 5 Hornsey 1957-02-24    99
  # 6 Hornsey 1957-03-03    64
  # 7 Hornsey 1957-03-10    87
  
  ukm
}

#' Use upper and lower bounds to sample initial parameters from box
#'
#' @param sh_ul `tbl` with `param`, `lower`, and `upper` columns.
#' @param sp_ul `tbl` with `param`, `lower`, and `upper` columns.
#' format of pparams.
#' @param units Character vector of unit names.
#' @param n_draws Number of initial parameter sets to draw.
#'
#' @return A list of parameters sets in the `pparams()` format.
#'
sample_initial_pparams_ul = function(
    sh_ul,
    sp_ul,
    units,
    n_draws
){
  helper_df = tidyr::expand_grid(
    x = sp_ul$param,
    y = units
  ) |>
    dplyr::rename(param = "x", unit = "y")
  
  expanded_specific = sp_ul |>
    dplyr::right_join(helper_df, by = "param") |>
    dplyr::mutate(`param[unit]` = paste0(.data$param,"[",.data$unit,"]"))
  
  to_named_vec = function(x, name_col, val_col){
    named_vec = x[[val_col]]
    names(named_vec) = x[[name_col]]
    named_vec
  }
  
  initial_parameters_tbl = dplyr::bind_cols(
    pomp::runif_design(
      lower = to_named_vec(sh_ul, "param", "lower"),
      upper = to_named_vec(sh_ul, "param", "upper"),
      nseq = n_draws
    ),
    pomp::runif_design(
      lower = to_named_vec(expanded_specific, "param[unit]", "lower"),
      upper = to_named_vec(expanded_specific, "param[unit]", "upper"),
      nseq = n_draws
    )
  )
  lapply(1:nrow(initial_parameters_tbl), function(z)
    coef_to_pparams(initial_parameters_tbl[z,])
  )
}

#' Convert coef-style object to pparams-style object.
#'
#' Vector entries with name in the style of `param[unit]` are assumed to be
#' unit-specific whereas those with name in the style of `param` are assumed to
#' be shared.
#'
#' @param coef Vector in the style of `coef(panelPomp_obj)`. That is, a numeric
#'   vector with names styled as "`shared_parameter`" or
#'   "`specific_parameter[unit]`".
#'
#' @return A list of length 2 in the style of `pparams(panelPomp_obj)`. That is,
#'   a numeric vector with shared parameter names, and a matrix with specific
#'   parameters as row names and units as column names.
#'
coef_to_pparams = function(coef){
  coef_tibble = tibble::tibble(
    full_param_name = names(coef),
    value = as.numeric(coef)
  ) |>
    tidyr::separate(
      .data$full_param_name,
      into = c("param", "unit"),
      sep = "\\[",
      fill = "right"
    ) |>
    dplyr::mutate(unit = gsub(pattern = "\\]", "", x = .data$unit))
  shared_tibble = dplyr::filter(coef_tibble, is.na(.data$unit))
  specific_tibble = dplyr::filter(coef_tibble, !is.na(.data$unit))
  shared_params = shared_tibble$value
  names(shared_params) = shared_tibble$param
  
  specific_params = specific_tibble |>
    tidyr::pivot_wider(names_from = "unit", values_from = "value") |>
    tibble::column_to_rownames(var = "param")
  
  list(shared = shared_params, specific = specific_params)
}

#' Turn named numeric vector into `rw_sd`
#'
#' @param rw_sd_vec Named numerical vector. Names are parameters, values are
#'   desired random walk standard deviations. Names of initial value parameters
#'   should end in _0.
#' @param weighted_param Name of the weighted parameter in a weighted model. Set
#'   to NULL (the default) if model does not contain a weighted parameter.
#'
#' @return Returns `rw_sd` object.
#'
make_rw_sd = function(rw_sd_vec, weighted_param = NULL){
  special = !is.null(weighted_param)
  ivp_indices = grep("_0", x = names(rw_sd_vec))
  if(special){
    special_indices = grep(
      paste0("^",weighted_param,"[1-9][0-9]*$"),
      x = names(rw_sd_vec)
    )
    si_ordered = sapply(seq_along(special_indices), function(x){
      grep(paste0("^",weighted_param, x,"$"), names(rw_sd_vec))
    })
  } else {
    special_indices = NULL
  }
  if(length(intersect(special_indices, ivp_indices)) > 0){
    stop(
      "IVPs cannot overlap with weighted parameters.",
      call. = FALSE
    )
  }
  reg_indices = setdiff(seq_along(rw_sd_vec), c(special_indices, ivp_indices))
  if(length(ivp_indices) > 0){
    ivp_rw_sd = lapply(names(rw_sd_vec[ivp_indices]), function(x){
      eval(bquote(expression(ivp(rw_sd_vec[[.(x)]]))))
    })
    names(ivp_rw_sd) = names(rw_sd_vec[ivp_indices])
    reg_rw_sd = as.list(rw_sd_vec[reg_indices])
  } else {
    ivp_rw_sd = NULL
    reg_rw_sd = as.list(rw_sd_vec)
  }
  if(special){
    special_rw_sd = as.list(rw_sd_vec[si_ordered])
    out = lapply(seq_along(special_indices), function(x){
      special_rw_sd[-x] = 0
      do.call(pomp::rw_sd, c(reg_rw_sd, ivp_rw_sd, special_rw_sd))
    })
  } else {
    out = do.call(pomp::rw_sd, c(reg_rw_sd, ivp_rw_sd))
  }
  out
}

#' Make a panelPomp or spatPomp model using measles data
#'
#' @param model `panel_mechanics` object.
#' @param data List in the format of `twentycities`.
#' @param starting_pparams Parameters in the format of `pparams()` output. Set
#'   to NULL to assign NA values. Only for panelPomp models currently.
#' @param interp_method Method used to interpolate population and births.
#'   Possible options are `"shifted_splines"` and `"linear"`.
#' @param first_year Integer for the first full year of data desired.
#' @param last_year Integer for the last full year of data desired.
#' @param custom_obs_list List of observations where each element supplies
#'   observations for a different city. Useful when using simulated
#'   observations. Set to `NULL` to use real observations. Only for panelPomp
#'   models currently.
#' @param dt Size of the time step.
#'
#' @return A panelPomp or spatPomp object using the model and data supplied.
#'
make_measlesPomp = function(
    model,
    data,
    starting_pparams = NULL,
    interp_method = c("shifted_splines", "linear"),
    first_year = 1950,
    last_year = 1963,
    custom_obs_list = NULL,
    dt = 1/365.25
){
  rproc = model$rproc
  dmeas = model$dmeas
  rmeas = model$rmeas
  rinit = model$rinit
  pt = model$pt
  paramnames = model$paramnames
  states = model$states
  measles = data$measles
  demog = data$demog
  
  ## ----prep-data-------------------------------------------------
  units = unique(measles$unit)
  # Obs list
  dat_list = vector("list", length(units))
  # Population list
  demog_list = vector("list", length(units))
  for(i in seq_along(units)){
    dat_list[[i]] = measles |>
      dplyr::mutate(year = as.integer(format(date,"%Y"))) |>
      dplyr::filter(
        .data$unit == units[[i]] & .data$year >= first_year &
          .data$year < (last_year + 1)
      ) |>
      dplyr::mutate(
        time = julian(
          .data$date,
          origin = as.Date(paste0(first_year, "-01-01"))
        )/365.25 + first_year
      ) |>
      dplyr::filter(.data$time > first_year & .data$time < (last_year + 1)) |>
      dplyr::select("time", "cases")
    if(!is.null(custom_obs_list)) dat_list[[i]]$cases = custom_obs_list[[i]]
    
    demog_list[[i]] = demog |>
      dplyr::filter(.data$unit == units[[i]]) |>
      dplyr::select(-"unit")
  }
  ## ----prep-covariates-------------------------------------------------
  covar_list = vector("list", length(units))
  for(i in seq_along(units)){
    dmgi = demog_list[[i]]
    times = seq(from = min(dmgi$year), to = max(dmgi$year), by = 1/12)
    switch(interp_method[[1]],
           shifted_splines = {
             pop_interp = stats::predict(
               stats::smooth.spline(x = dmgi$year, y = dmgi$pop),
               x = times
             )$y
             births_interp = stats::predict(
               stats::smooth.spline(x = dmgi$year + 0.5, y = dmgi$births),
               x = times - 4
             )$y
           },
           linear = {
             pop_interp = stats::approx(
               x = dmgi$year,
               y = dmgi$pop,
               xout = times
             )$y
             births_interp = stats::approx(
               x = dmgi$year,
               y = dmgi$births,
               xout = times - 4
             )$y
           }
    )
    covar_list[[i]] = dmgi |>
      dplyr::reframe(
        time = times,
        pop = pop_interp,
        birthrate = births_interp
      )
    covar_list[[i]] = covar_list[[i]] |>
      dplyr::mutate(
        pop_1950 = dplyr::filter(
          covar_list[[i]], covar_list[[i]]$time == 1950
        )$pop
      )
  }
  for(i in seq_along(units)){
    log_pop_1950 = sapply(seq_along(units), function(x)
      log(covar_list[[x]][["pop_1950"]][[1]])
    )
    covar_list[[i]] = covar_list[[i]] |>
      dplyr::mutate(
        std_log_pop_1950 = (log(.data$pop_1950) - mean(log_pop_1950))/
          stats::sd(log_pop_1950),
        unit_num = i
      )
  }
  
  ## ----pomp-construction-----------------------------------------------
  lapply(seq_along(units), function(i){
    time = covar_list[[i]]$time
    dat_list[[i]] |>
      pomp::pomp(
        t0 = with(dat_list[[i]], 2*time[1] - time[2]),
        times = "time",
        rprocess = pomp::euler(rproc, delta.t = dt),
        rinit = rinit,
        dmeasure = dmeas,
        rmeasure = rmeas,
        covar = pomp::covariate_table(covar_list[[i]], times = "time"),
        accumvars = c("C","W"),
        partrans = pt,
        statenames = states,
        paramnames = paramnames
      )
  }) -> pomp_list
  names(pomp_list) = units
  
  ## ----panelPomp-construction-----------------------------------------------
  if(is.null(starting_pparams)){
    shared = as.numeric(rep(NA, length(model$shared_params)))
    specific = matrix(
      NA,
      nrow = length(model$specific_params),
      ncol = length(units)
    )
    class(specific) = "numeric"
    storage.mode(specific) = "numeric"
    rownames(specific) = model$specific_params
    colnames(specific) = units
    names(shared) = model$shared_params
  } else {
    shared = starting_pparams$shared
    specific = as.matrix(starting_pparams$specific)
    if(!setequal(names(shared), model$shared_params)){
      stop(
        "Starting shared parameters do not match parameters in model mechanics.",
        call. = FALSE
      )
    }
    if(!setequal(rownames(specific), model$specific_params)){
      stop(
        "Starting unit-specific parameters do not match parameters in model mechanics.",
        call. = FALSE
      )
    }
  }
  panelPomp::panelPomp(
    pomp_list,
    shared = shared,
    specific = specific
  )
}

#' Perform one round of model fitting
#'
#' Depending on the model used, `run_round` will use either [panelPomp::mif2] or
#' [spatPomp::ibpf] to fit it, then use `eval_logLik()` to estimate the log
#' likelihood of the fit. The results are saved using [pomp::bake].
#'
#' @param x A `panelPomp` object.
#' @param initial_pparams_list List of initial parameters in the format of
#'   `panelPomp::coef()` with `format = "list`. Each entry in the list specifies
#'   the initial parameters for one replication of the fitting algorithm.
#' @param write_results_to File path to save Rds file containing results to.
#' @param ncores Number of cores to use.
#' @param np_fitr Number of particles to use when running the fitting algorithm.
#' @param cooling_frac Cooling fraction to use when running fitting algorithm.
#' @param rw_sd_obj Object of class `rw_sd` specifying random walk standard
#'   deviations to use when running the fitting algorithm.
#' @param N_fitr Number of iterations to use when running the fitting algorithm.
#' @param np_eval Number of particles to use when running `eval_logLik()`.
#' @param nreps_eval Number of replications to use when running `eval_logLik()`.
#' @param print_times Boolean for whether times to run the fitting algorithm and
#'   `eval_logLik()` should be printed.
#' @param panel_block Boolean specifying whether to perform block resampling of
#'   specific parameters for a `panelPomp` model. Only used when `x` is a
#'   `panelPomp` object.
#'
#' @return Object of class `fit_results` containing a list of `mif2d.ppomp` or
#'   `ibpfd_spatPomp` objects and a list of `EL_list` objects.
#'
run_round = function(
    x,
    initial_pparams_list,
    write_results_to,
    ncores,
    np_fitr,
    cooling_frac,
    rw_sd_obj,
    N_fitr,
    panel_block = FALSE,
    np_eval,
    nreps_eval,
    print_times = FALSE
){
  doParallel::registerDoParallel(cores = ncores)
  RNGkind("L'Ecuyer-CMRG")
  doRNG::registerDoRNG()
  fit_results = pomp::bake(file = write_results_to, {
    if(print_times) start_t = Sys.time()
    mif2_out = foreach::foreach(
      z = initial_pparams_list,
      .packages = "panelPomp"
    ) %dopar% {
      use_shared = !is.null(z$shared)
      use_specific = !is.null(z$specific)
      args = c(
        list(
          x,
          Np = np_fitr,
          cooling.fraction.50 = cooling_frac,
          rw.sd = rw_sd_obj,
          cooling.type = "geometric",
          Nmif = N_fitr,
          block = panel_block
        ),
        list(shared.start = z$shared)[use_shared],
        list(specific.start = z$specific)[use_specific]
      )
      do.call(panelPomp::mif2, args = args)
    }
    fitr_out = mif2_out
    if(print_times) print(Sys.time() - start_t)
    if(print_times) start_t = Sys.time()
    EL_out = eval_logLik(
      model_obj_list = fitr_out,
      ncores = ncores,
      np_pf = np_eval,
      nreps = nreps_eval,
      seed = NULL
    )
    if(print_times) print(Sys.time() - start_t)
    new_fit_results(fitr_out = fitr_out, EL_out = EL_out)
  })
  fit_results
}

#' Evaluate the log likelihood of a model using `pfilter()`
#'
#' @param model_obj_list List of `panelPomp` objects to evaluate
#'   the log likelihood of.
#' @param block_size The number of spatial units per block. Only used when
#'   evaluating are `spatPomp` models. (NOTE: function will break for any block
#'   size other than 1; this will be fixed when choosing other block sizes seems
#'   worthwhile.)
#' @param ncores Number of cores to use for parallel computing.
#' @param np_pf Number of particles to use.
#' @param nreps Number of particle filter repetitions.
#' @param seed Seed for particle filter. If NULL, does not set new seed.
#'
#'
#' @return Object of type `EL_list`, a list of data frames containing log
#'   likelihood and standard error estimates.
#'
#'   For `panelPomp`, log likelihood estimates are obtained by using
#'   [panelPomp::panel_logmeanexp()] over the `pfilter()` log likelihood
#'   replications. For `spatPomp`, log likelihood estimates are obtained by
#'   using [pomp::logmeanexp()] over the `pfilter()` log likelihood
#'   replications. For both model types, unit log likelihood estimates are
#'   obrained by using `logmeanexp()` over the unit log likelihood replications
#'   for each unit, and conditional log likelihood estimates are obtained by
#'   using `logmeanexp()` over the condtional log likelihood replications for
#'   each unit and time point.
#'
eval_logLik = function(
    model_obj_list,
    block_size = 1,
    ncores,
    np_pf,
    nreps,
    seed = NULL
){
  N_models = length(model_obj_list)
  units = names(model_obj_list[[1]])
  
  pf_logLik_frame = data.frame(
    logLik = rep(0, N_models),
    se = rep(0, N_models)
  ) |>
    cbind(
      rbind(t(sapply(model_obj_list, panelPomp::coef)))
    )
  
  pf_unitlogLik_list = vector("list", N_models)
  pf_unitSE_list = vector("list", N_models)
  pf_cll_list = vector("list", N_models)
  
  doParallel::registerDoParallel(cores = ncores)
  RNGkind("L'Ecuyer-CMRG")
  doRNG::registerDoRNG(seed)
  
  foreach_out = foreach::foreach(
    i = 1:nreps,
    .packages = "panelPomp"
  ) %dopar% {
    lapply(model_obj_list, function(x){
      out = panelPomp::pfilter(x, Np = np_pf)
      out = list(
        ull = panelPomp::unitLogLik(out),
        cll = sapply(
          tryCatch(out@unit_objects, error = function(x) out@unit.objects),
          function(u){
            u@cond.logLik
          }) |> t() |> `rownames<-`(units)
      )
      out
    })
  }
  
  ull_matrices = lapply(1:N_models, function(i){
    lapply(1:nreps, function(j){
      foreach_out[[j]][[i]]$ull
    }) |> dplyr::bind_rows() |> as.matrix()
  })
  
  cllse_matrices = lapply(1:N_models, function(i){
    lapply(units, function(u){
      sapply(1:nreps, function(j){
        foreach_out[[j]][[i]]$cll[u,]
      }) |> apply(MARGIN = 1, FUN = pomp::logmeanexp, se = TRUE)
    }) |> `names<-`(units)
  })
  
  llse = sapply(1:N_models, function(i){
    panelPomp::panel_logmeanexp(ull_matrices[[i]], MARGIN = 2, se = TRUE)
  }) |> t()
  pf_logLik_frame[,1:2] = llse
  
  ullse = lapply(1:N_models, function(i){
    apply(ull_matrices[[i]], MARGIN = 2, FUN = pomp::logmeanexp, se = TRUE)
  }) |> t()
  
  ull = lapply(1:N_models, function(i){
    out = as.data.frame(ullse[[i]][1,, drop = FALSE])
    rownames(out) = NULL
    out
  }) |> dplyr::bind_rows()
  
  se = lapply(1:N_models, function(i){
    out = as.data.frame(ullse[[i]][2,, drop = FALSE])
    rownames(out) = NULL
    out
  }) |> dplyr::bind_rows()
  
  cll = lapply(units, function(u){
    sapply(1:N_models, function(i){
      cllse_matrices[[i]][[u]]["est",]
    }) |> t()
  })
  names(cll) = units
  cll_se = lapply(units, function(u){
    sapply(1:N_models, function(i){
      cllse_matrices[[i]][[u]]["se",]
    }) |> t()
  })
  names(cll_se) = units
  
  new_EL_list(
    fits = pf_logLik_frame,
    ull = ull,
    se = se,
    cll = cll,
    cll_se = cll_se,
    np_pf = np_pf,
    nreps = nreps
  )
}

#' Make a list containing the results of `eval_logLik()`
#'
#' @name EL_list
#' @param fits Data frame of fit results. Columns should be named `logLik`,
#'   `se`, followed by the parameter names.
#' @param ull Data frame of unit log likelihoods. Column names are unit names.
#' @param se Data frame of unit standard errors. Column names are unit names.
#' @param cll List of matrices containing estimated conditional log likelihoods
#'   for each unit.
#' @param cll_se List of matrices containing standard errors for estimated
#'   conditional log likelihoods for each unit.
#' @param np_pf Number of particles used by `eval_logLik()`.
#' @param nreps Number of replications used by `eval_logLik()`.
#'
#' @return The arguments in list form with class `EL_list`.
#'

new_EL_list = function(
    fits,
    ull,
    se,
    cll,
    cll_se,
    np_pf,
    nreps
){
  stopifnot(is.data.frame(fits))
  stopifnot(is.data.frame(ull))
  stopifnot(is.data.frame(se))
  stopifnot(is.list(cll))
  stopifnot(is.list(cll_se))
  stopifnot(is.numeric(np_pf))
  stopifnot(is.numeric(nreps))
  out = list(
    fits = dplyr::as_tibble(fits),
    ull = dplyr::as_tibble(ull),
    se = dplyr::as_tibble(se),
    cll = cll,
    cll_se = cll_se,
    np_pf = np_pf,
    nreps = nreps
  )
  structure(out, class = "EL_list")
}

validate_EL_list = function(x){
  if(!all(c("fits", "ull", "se", "np_pf", "nreps") %in% names(x))){
    stop(
      "`x` must have names 'fits', 'ull', 'se', 'np_pf', and 'nreps'",
      call. = FALSE
    )
  }
  if(!all(is.data.frame(x$fits), is.data.frame(x$ull), is.data.frame(x$se))){
    stop(
      "'fits', 'ull', and 'se' must be of class data.frame"
    )
  }
  if(!all(is.numeric(x$np_pf), is.numeric(x$nreps))){
    stop(
      "'np_pf' and 'nreps' must be of class numeric"
    )
  }
  if(!all(c("logLik", "se") %in% colnames(x$fits))){
    stop(
      "'fits' must have 'logLik' and 'se' columns",
      call. = FALSE
    )
  }
  x
}

EL_list = function(
    fits,
    ull,
    se,
    cll,
    cll_se,
    np_pf,
    nreps
){
  x = new_EL_list(fits, ull, se, cll, cll_se, np_pf, nreps)
  validate_EL_list(x)
}

#' Make new fit_results object
#'
#' @param fitr_out List of `mif2d.ppomp` objects.
#' @param EL_out `EL_list` object.
#'
#' @return List of class `fit_results` with first entry `fitr_out` and second
#'  entry `EL_out`. `fitr_out` is a list of `mif2d.ppomp` objects, and `EL_out`
#'  is an `EL_list`.
#'
new_fit_results = function(fitr_out, EL_out){
  stopifnot(is.list(fitr_out))
  stopifnot(class(EL_out) == "EL_list")
  out = list(
    fitr_out = fitr_out,
    EL_out = EL_out
  )
  structure(out, class = "fit_results")
}

#' Gather information on fit results from subdirectories and tidy them
#'
#' Given a parent directory and a vector of rds file names, this function will
#' attempt to load any rds files with those names under the parent directory. If
#' the rds file contains a `fit_results` object or an `EL_list` object, it will
#' call `tidy_results()` on it and row bind the resulting data frames. If an
#' error occurs at any point when trying to load and tidy a given file, the file
#' will be skipped over.
#'
#' @param file_names Character vector of names of rds files to load.
#' @param parent_dir Path to the directory that we want to search the
#'   subdirectories of.
#'
#' @return `data.frame` where each row contains information for a unit from a
#'   model replication. Columns include information such as the path to the Rds
#'   file, the number of particles used for the fitting algorithm and likelihood
#'   evaluation, and the parameters that the likelihood was evaluated at.
#'
gather_results = function(file_names, parent_dir = "."){
  out_paths = c()
  for(fn in file_names){
    out_paths = c(
      out_paths,
      list.files(
        path = parent_dir,
        pattern = paste0("*", fn),
        recursive = TRUE
      )
    )
  }
  out_paths = paste0(parent_dir,"/",out_paths)
  lapply(out_paths, function(path){
    out = try({
      loaded_obj = readRDS(path)
      tidy_results(loaded_obj, path = path)
    })
    if(inherits(out, "try-error")){
      NULL
    } else {
      out
    }
  }) |>
    dplyr::bind_rows()
}

#' Convert `fit_results` or `EL_list` output into tidy data frame
#'
#' @param x Object of type `fit_results` or `EL_list`.
#' @param path Character specifying where `x` was obtained from. Mostly exists
#'   for use by [measlespkg::gather_results()]. Very optional.
#'
#' @return `data.frame` where each row contains information for a unit from a
#'   model replication. Columns include information such as the number of
#'   particles used for the fitting algorithm and likelihood evaluation, and the
#'   parameters that the likelihood was evaluated at.
#'
tidy_results = function(x, path = NA){
  UseMethod("tidy_results")
}

tidy_results_helper = function(
    ELL,
    path,
    N_fitr,
    Np,
    cf,
    np_eval,
    nreps_eval,
    block,
    obs_hash
){
  ELL |>
    tidy_pfilter_dfs() |>
    dplyr::mutate(
      path = path,
      N_fitr = N_fitr,
      np_fitr = Np,
      cooling_frac = cf,
      np_eval = np_eval,
      nreps_eval = nreps_eval,
      block = block,
      obs_hash = obs_hash
    ) |>
    dplyr::select(
      "path", "N_fitr", "np_fitr", "cooling_frac", "block", "np_eval",
      "nreps_eval", dplyr::everything()
    )
}

#' @rdname tidy_results
tidy_results.EL_list = function(x, path = NA){
  ELL = x
  N_fitr = NA
  cf = NA
  block = NA
  Np = NA
  obs_hash = NA
  np_eval = ELL$np_pf
  nreps_eval = ELL$nreps
  tidy_results_helper(
    ELL = ELL,
    path = path,
    N_fitr = N_fitr,
    Np = Np,
    cf = cf,
    np_eval = np_eval,
    nreps_eval = nreps_eval,
    block = block,
    obs_hash = obs_hash
  )
}

#' @rdname tidy_results
tidy_results.fit_results = function(x, path = NA){
  FITR_O = x$fitr_out[[1]]
  ELL = x$EL_out
  ###
  if(inherits(FITR_O, "mif2d.ppomp")){
    N_fitr = FITR_O@Nmif
    block = FITR_O@block
    Np = FITR_O@Np
  }
  obs_hash = rlang::hash(obs2(FITR_O))
  ###############################
  cf = FITR_O@cooling.fraction.50
  np_eval = ELL$np_pf
  nreps_eval = ELL$nreps
  tidy_results_helper(
    ELL = ELL,
    path = path,
    N_fitr = N_fitr,
    Np = Np,
    cf = cf,
    np_eval = np_eval,
    nreps_eval = nreps_eval,
    block = block,
    obs_hash = obs_hash
  )
}

#' Extract the data matrix from a pomp/panelPomp object.
#'
#' @param x An object of class `pomp` or `panelPomp``.
#'
#' @return A matrix of observations with rows corresponding to units and columns
#'   corresponding to time points.
#' @export
#'
#' @examples
#' \dontrun{
#' obs2(AK_model())
#' }
obs2 = function(x){
  if(inherits(x, "panelPomp")){
    out = tryCatch(x@unit_objects, error = function(z) x@unit.objects) |>
      sapply(pomp::obs) |>
      t()
  } else {
    out = pomp::obs(x)
  }
  out
}

#' Create tidy data frame from `eval_logLik()` output.
#'
#' @param x Object of class `EL_list`.
#'
#' @return Data frame composed from information in `x` data in tidy format.
#'
tidy_pfilter_dfs = function(x){
  stopifnot(class(x) == "EL_list")
  lapply(1:nrow(x$fits), function(z){
    z_pparams = coef_to_pparams(x$fits[z,-c(1, 2)])
    tidy_LL_df = z_pparams$specific |>
      t() |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "unit") |>
      dplyr::mutate(rep = z)
    if(length(z_pparams$shared) > 0){
      tidy_LL_df = z_pparams$shared |>
        tibble::enframe() |>
        tidyr::pivot_wider() |>
        dplyr::bind_cols(tidy_LL_df)
    }
    tidy_ull_df = subset(x$ull, subset = rownames(x$ull) == z) |>
      t() |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "unit") |>
      dplyr::rename(ull = 2)
    tidy_se_df = subset(x$se, subset = rownames(x$se) == z) |>
      t() |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "unit") |>
      dplyr::rename(se = 2)
    tidy_df = dplyr::left_join(tidy_ull_df, tidy_LL_df, by = "unit") |>
      dplyr::left_join(tidy_se_df, by = "unit") |>
      dplyr::mutate(total_ll = x$fits$logLik[[z]], total_se = x$fits$se[[z]]) |>
      dplyr::select(
        "rep", "total_ll", "total_se", "unit", "ull", "se", dplyr::everything()
      )
    tidy_df
  }) |>
    dplyr::bind_rows() -> tidy_coef_df
  tidy_coef_df
}
