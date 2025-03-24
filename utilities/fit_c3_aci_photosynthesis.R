library(photosynthesis)

# Here we assume that `defaults.R` has already been run

# This function is based on the example from the "Combining PhotoGEA With Other
# Packages" vignette from the PhotoGEA R package, which is available online at
# https://eloch216.github.io/PhotoGEA/articles/combining_with_other_packages.html

# Make a wrapper for `photosynthesis::fit_aci_response`
fit_c3_aci_photosynthesis <- function(
    replicate_exdf, # an `exdf` object representing a single A-Ci curve
    a_column_name = 'A',                     # micromol / m^2 / s
    ci_column_name = 'Ci',                   # micromol / mol
    oxygen_column_name = 'oxygen',           # percent
    pressure_column_name = 'total_pressure', # bar
    qin_column_name = 'Qin',                 # micromol / m^2 / s
    tleaf_column_name = 'TleafCnd',          # degrees C
    ... # Not used
)
{
    ### Check inputs

    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_aci_photosynthesis requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]        <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]]       <- 'micromol mol^(-1)'
    required_variables[[oxygen_column_name]]   <- 'percent'
    required_variables[[pressure_column_name]] <- 'bar'
    required_variables[[qin_column_name]]      <- 'micromol m^(-2) s^(-1)'
    required_variables[[tleaf_column_name]]    <- 'degrees C'

    check_required_variables(replicate_exdf, required_variables)

    # Convert leaf temperature to Kelvin (??)
    replicate_exdf[, tleaf_column_name] <- replicate_exdf[, tleaf_column_name] + 273.15

    # Just keep the minimum necessary columns
    dataf <- replicate_exdf[, c(a_column_name, tleaf_column_name, ci_column_name, qin_column_name)]

    # Get the pressure in bar and kPa
    pressure_bar <- mean(replicate_exdf[, pressure_column_name])
    pressure_kPa <- pressure_bar * 100

    ### Call function from external packge with appropriate units
    fit_res <- photosynthesis::fit_aci_response(
        dataf,
        varnames = list(
            A_net = a_column_name,
            T_leaf = tleaf_column_name,
            C_i = ci_column_name,
            PPFD = qin_column_name
        ),
        P = pressure_kPa, # kPa
        fitTPU = TRUE,
        alpha_g = 0,
        R_d_meas = NULL,
        useR_d = FALSE,
        useg_mc = FALSE,
        useg_mct = FALSE,
        usegamma_star = FALSE,
        useK_M = FALSE,
        useK_C_K_O = TRUE, # changed from default
        alpha = ALPHA_J,
        theta_J = THETA_J,
        gamma_star25 = with(TEMPERATURE_PARAM$Gamma_star, {exp(c - Ea / PhotoGEA:::f)}), # micromol / mol
        Ea_gamma_star = TEMPERATURE_PARAM$Gamma_star$Ea * 1e3,                           # J / mol
        K_M25 = NULL,   #718.4
        Ea_K_M = NULL,  #65508.28
        g_mc25 = NULL,  #0.08701
        Ea_g_mc = NULL, #0
        K_C25 = with(TEMPERATURE_PARAM$Kc, {exp(c - Ea / PhotoGEA:::f)}) * pressure_bar, # microbar
        Ea_K_C = TEMPERATURE_PARAM$Kc$Ea * 1e3,                                          # J / mol
        K_O25 = with(TEMPERATURE_PARAM$Ko, {exp(c - Ea / PhotoGEA:::f)}) * 1e-3 * pressure_kPa, # kPa
        Ea_K_O = TEMPERATURE_PARAM$Ko$Ea * 1e3,                                                 # J / mol
        Oconc = mean(replicate_exdf[, oxygen_column_name]), # percent
        gamma_star_set = NULL,
        K_M_set = NULL
    )

  ### Collect outputs

  # Get the identifier columns from the original exdf object
  replicate_identifiers <- identifier_columns(replicate_exdf)

  # Get a data frame with the fitted values of assimilation and convert it to an
  # `exdf` object, setting the category to `fit_c3_aci_photosynthesis` and
  # specifying the units for each column
  fits <- exdf(fit_res[[3]])

  # Append the identifier columns to the fits
  fits <- cbind(replicate_identifiers, fits)

  # Add a column for the residuals
  fits <- set_variable(
    fits,
    'A_residuals',
    'micromol m^(-2) s^(-1)',
    'fit_c3_aci_photosynthesis',
    fits[, 'A'] - fits[, 'A_model']
  )

  # Rename a few columns
  colnames(fits)[colnames(fits) == 'A_carbox'] <- 'Ac'
  colnames(fits)[colnames(fits) == 'A_regen']  <- 'Aj'
  colnames(fits)[colnames(fits) == 'A_tpu']    <- 'Ap'

  # Document the fit units
  fits <- document_variables(
    fits,
    c('fit_c3_aci_photosynthesis', 'A_model', 'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'Ac',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'Aj',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'Ap',      'micromol m^(-2) s^(-1)')
  )

  # Create an exdf object with the parameter values that are included in the
  # fitting result.
  parameters <- exdf(data.frame(
    RSS_reported  = as.numeric(fit_res[[1]]$cost) * 2,
    citransition1 = as.numeric(fit_res[[1]]$citransition1),
    citransition2 = as.numeric(fit_res[[1]]$citransition2),
    alpha         = ALPHA_J,
    theta_J       = THETA_J,
    V_cmax        = as.numeric(fit_res[[1]]$V_cmax),
    V_cmax_se     = as.numeric(fit_res[[1]]$V_cmax_se),
    J_max         = as.numeric(fit_res[[1]]$J_max),
    J             = as.numeric(fit_res[[1]]$J),
    J_se          = as.numeric(fit_res[[1]]$J_se),
    V_TPU         = as.numeric(fit_res[[1]]$V_TPU),
    V_TPU_se      = as.numeric(fit_res[[1]]$V_TPU_se),
    R_d           = as.numeric(fit_res[[1]]$R_d),
    R_d_se        = as.numeric(fit_res[[1]]$R_d_se),
    V_cmax_pts    = as.numeric(fit_res[[1]]$V_cmax_pts),
    J_max_pts     = as.numeric(fit_res[[1]]$J_max_pts),
    V_TPU_pts     = as.numeric(fit_res[[1]]$V_TPU_pts)
  ))

  # `photosynthesis::fit_aci_response` returns absurdly high values of `V_TPU`
  # when it cannot get a good estimate for `TPU`. These can cause problems when
  # plotting the results, so we replace any `V_TPU` values above 500 micromol /
  # m^2 / s with NA.
  if (parameters[1, 'V_TPU'] > 500) {
    parameters[1, 'V_TPU'] <- NA
  }

  fits[fits[, 'Ap'] > 500, 'Ap'] <- NA

  # Identify limiting process
  fits <- identify_c3_limiting_processes(
    fits,
    'A_model',
    'Ac',
    'Aj',
    'Ap'
  )

  # The Rd values as returned are negative
  parameters[1, 'R_d'] <- parameters[1, 'R_d'] * -1

  # Calculate upper and lower confidence interval limits for key parameters.
  # photosynthesis reports standard errors (SE). The limits for a 95% confidence
  # interval are mean +/- 1.96 * SE.
  ci_f <- 1.96

  parameters[, 'J_lower']      <- parameters[, 'J']      - parameters[, 'J_se']      * ci_f
  parameters[, 'J_upper']      <- parameters[, 'J']      + parameters[, 'J_se']      * ci_f
  parameters[, 'R_d_lower']    <- parameters[, 'R_d']    - parameters[, 'R_d_se']    * ci_f
  parameters[, 'R_d_upper']    <- parameters[, 'R_d']    + parameters[, 'R_d_se']    * ci_f
  parameters[, 'V_TPU_lower']  <- parameters[, 'V_TPU']  - parameters[, 'V_TPU_se']  * ci_f
  parameters[, 'V_TPU_upper']  <- parameters[, 'V_TPU']  + parameters[, 'V_TPU_se']  * ci_f
  parameters[, 'V_cmax_lower'] <- parameters[, 'V_cmax'] - parameters[, 'V_cmax_se'] * ci_f
  parameters[, 'V_cmax_upper'] <- parameters[, 'V_cmax'] + parameters[, 'V_cmax_se'] * ci_f

  # Determine whether this is an "inadmissible fit".
  #
  # A note about the Ci transitions: In the fit_aci_response code, we have:
  #
  #  datc <- data[data$C_i < citransdf$ci1[i], ]
  #  datj <- data[data$C_i > citransdf$ci1[i] & data$C_i < citransdf$ci2[i], ]
  #  datp <- data[data$C_i > citransdf$ci2[i], ]
  #
  ci1 <- parameters[, 'citransition1']
  ci2 <- parameters[, 'citransition2']

  parameters[, 'admissible'] <-
    all(fits[fits[, 'Ci'] < ci1, 'Ac_limiting']) &&
    all(fits[fits[, 'Ci'] > ci1 & fits[, 'Ci'] < ci2, 'Aj_limiting']) &&
    all(fits[fits[, 'Ci'] > ci2, 'Ap_limiting'])

  # Document the parameter units
  parameters <- document_variables(
    parameters,
    c('fit_c3_aci_photosynthesis', 'RSS_reported',   ''),
    c('fit_c3_aci_photosynthesis', 'citransition1',  'micromol mol^(-1)'),
    c('fit_c3_aci_photosynthesis', 'citransition2',  'micromol mol^(-1)'),
    c('fit_c3_aci_photosynthesis', 'alpha',          'dimensionless'),
    c('fit_c3_aci_photosynthesis', 'theta_J',        'dimensionless'),
    c('fit_c3_aci_photosynthesis', 'V_cmax',         'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_cmax_se',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'J_max',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'J',              'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'J_se',           'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_TPU',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_TPU_se',       'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'R_d',            'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'R_d_se',         'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_cmax_pts',     ''),
    c('fit_c3_aci_photosynthesis', 'J_max_pts',      ''),
    c('fit_c3_aci_photosynthesis', 'V_TPU_pts',      ''),
    c('fit_c3_aci_photosynthesis', 'J_lower',        'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'J_upper',        'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'R_d_lower',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'R_d_upper',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_TPU_lower',    'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_TPU_upper',    'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_cmax_lower',   'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'V_cmax_upper',   'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_photosynthesis', 'admissible',     'micromol m^(-2) s^(-1)')
  )

  # Append the identifier columns to the parameters
  parameters <- cbind(replicate_identifiers, parameters)

  # Attach the residual stats to the parameters
  nparam <- 4 # J, Vcmax, Rd, TPU

  parameters <- cbind(
      parameters,
      residual_stats(
          fits[, 'A_residuals'],
          fits$units[['A_residuals']],
          nparam
      )
  )

  # Return a list of two data frames: `fits` and `parameters`
  return(list(
      fits = fits,
      parameters = parameters
  ))
}
