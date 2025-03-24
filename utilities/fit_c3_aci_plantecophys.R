library(plantecophys)

# Here we assume that `defaults.R` has already been run

# This function is based on the example from the "Combining PhotoGEA With Other
# Packages" vignette from the PhotoGEA R package, which is available online at
# https://eloch216.github.io/PhotoGEA/articles/combining_with_other_packages.html

# Make a wrapper for `plantecophys::fitaci`
fit_c3_aci_plantecophys <- function(
  replicate_exdf, # an `exdf` object representing a single A-Ci curve
  a_column_name = 'A',                     # micromol / m^2 / s
  ci_column_name = 'Ci',                   # micromol / mol
  gamma_star_column_name = 'Gamma_star',   # micromol / mol
  kc_column_name = 'Kc',                   # micromol / mol
  ko_column_name = 'Ko',                   # mmol / mol
  oxygen_column_name = 'oxygen',           # percent
  pressure_column_name = 'total_pressure', # bar
  qin_column_name = 'Qin',                 # micromol / m^2 / s
  tleaf_column_name = 'TleafCnd',          # degrees C
  ... # Not used
)
{
  ### Check inputs

  if (!is.exdf(replicate_exdf)) {
      stop('fit_c3_aci_plantecophys requires an exdf object')
  }

  # Make sure the required variables are defined and have the correct units
  required_variables <- list()
  required_variables[[a_column_name]]          <- 'micromol m^(-2) s^(-1)'
  required_variables[[ci_column_name]]         <- 'micromol mol^(-1)'
  required_variables[[gamma_star_column_name]] <- 'micromol mol^(-1)'
  required_variables[[kc_column_name]]         <- 'micromol mol^(-1)'
  required_variables[[ko_column_name]]         <- 'mmol mol^(-1)'
  required_variables[[oxygen_column_name]]     <- 'percent'
  required_variables[[pressure_column_name]]   <- 'bar'
  required_variables[[qin_column_name]]        <- 'micromol m^(-2) s^(-1)'
  required_variables[[tleaf_column_name]]      <- 'degrees C'

  check_required_variables(replicate_exdf, required_variables)

  # Specify or extract some values
  PPFD <- mean(replicate_exdf[, qin_column_name]) # micromol / m^2 / s

  Kc <- mean(replicate_exdf[, kc_column_name] * replicate_exdf[, pressure_column_name]) # microbar
  Ko <- mean(replicate_exdf[, ko_column_name] * replicate_exdf[, pressure_column_name]) # mbar
  PO <- mean(replicate_exdf[, oxygen_column_name] * 1e-2 * replicate_exdf[, pressure_column_name] * 1e3) # mbar

  Km <- Kc * (1.0 + PO / Ko)

  dataf <- replicate_exdf[, c(a_column_name, tleaf_column_name, qin_column_name, ci_column_name)]

  ### Call function from external packge with appropriate units
  fit_res <- plantecophys::fitaci(
    dataf,
    varnames = list(
        ALEAF = a_column_name,
        Tleaf = tleaf_column_name,
        PPFD = qin_column_name,
        Ci = ci_column_name
    ),
    Tcorrect = FALSE, # changed from default
    Patm = mean(replicate_exdf[, pressure_column_name]) * 100, # kPa
    citransition = NULL,
    quiet = FALSE,
    startValgrid = TRUE,
    fitmethod = 'bilinear', # must be bilinear when fitTPU is TRUE
    algorithm = 'default',
    fitTPU = TRUE, # changed from default
    alphag = 0,
    useRd = FALSE,
    PPFD = NULL,
    Tleaf = NULL,
    alpha = ALPHA_J,
    theta = THETA_J,
    gmeso = NULL,
    EaV = TEMPERATURE_PARAM$Vcmax_norm$Ea * 1e3,
    EdVC = 0,
    delsC = 0,
    EaJ = TEMPERATURE_PARAM$J_norm$Ea * 1e3,
    EdVJ = 0,
    delsJ = 0,
    GammaStar = mean(replicate_exdf[, gamma_star_column_name] * replicate_exdf[, pressure_column_name]), # microbar
    Km = Km, # microbar
    id = NULL
  )

  ### Collect outputs

  # Get the identifier columns from the original exdf object
  replicate_identifiers <- identifier_columns(replicate_exdf)

  # Get a data frame with the fitted values of assimilation and convert it to an
  # `exdf` object, setting the category to `fit_c3_aci_plantecophys` and
  # specifying the units for each column
  fits <- exdf(fit_res$df)

  fits <- document_variables(
    fits,
    c('fit_c3_aci_plantecophys', 'Ci',          'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'Ameas',       'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Amodel',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Ac',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Aj',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Ap',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Rd',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'VPD',         'kPa'),
    c('fit_c3_aci_plantecophys', 'Tleaf',       'degrees C'),
    c('fit_c3_aci_plantecophys', 'Cc',          'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'PPFD',        'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Patm',        'kPa'),
    c('fit_c3_aci_plantecophys', 'Ci_original', 'micromol mol^(-1)')
  )

  # Append the identifier columns to the fits
  fits <- cbind(replicate_identifiers, fits)

  # `plantecophys::fitaci` returns `Ap = 1000` when it cannot get a good
  # estimate for `TPU`. These can cause problems when plotting the results, so
  # we replace any `Ap` values above 500 micromol / m^2 / s with NA.
  fits[fits[, 'Ap'] > 500, 'Ap'] <- NA

  # Convert `plantecophys::fitaci` outputs to net CO2 assimilation rates for
  # consistency with `fit_c3_aci`.
  fits[, 'Ac'] <- fits[, 'Ac'] - fits[, 'Rd']
  fits[, 'Aj'] <- fits[, 'Aj'] - fits[, 'Rd']
  fits[, 'Ap'] <- fits[, 'Ap'] - fits[, 'Rd']

  # Add a column for the residuals
  fits <- set_variable(
    fits,
    'A_residuals',
    'micromol m^(-2) s^(-1)',
    'fit_c3_aci_plantecophys',
    fits[, 'Ameas'] - fits[, 'Amodel']
  )

  # Identify limiting process
  fits <- identify_c3_limiting_processes(
    fits,
    'Amodel',
    'Ac',
    'Aj',
    'Ap'
  )

  # Create an exdf object with the parameter values that are included in the
  # fitting result. Note that we do not include RMSE because it is not
  # calculated correctly; the returned value of "RMSE" is actually the root of
  # the sum of squared errors, not the root of the mean squared error.
  parameters <- exdf(data.frame(
    RSS_reported   = as.numeric(fit_res$RMSE)^2,
    Ci_transition  = as.numeric(fit_res$Ci_transition),
    Ci_transition2 = as.numeric(fit_res$Ci_transition2),
    Tcorrect       = fit_res$Tcorrect,
    Rd_measured    = fit_res$Rd_measured,
    GammaStar      = fit_res$GammaStar,
    Km             = fit_res$Km,
    kminput        = fit_res$kminput,
    gstarinput     = fit_res$gstarinput,
    fitmethod      = fit_res$fitmethod,
    citransition   = fit_res$citransition,
    gmeso          = fit_res$gmeso,
    fitTPU         = fit_res$fitTPU,
    alphag         = fit_res$alphag,
    Vcmax          = fit_res$par[1, 1],
    Vcmax_err      = fit_res$par[1, 2],
    Qin_avg        = mean(replicate_exdf[, qin_column_name]),
    alpha          = ALPHA_J,
    theta          = THETA_J,
    Rd             = fit_res$par[3, 1],
    Rd_err         = fit_res$par[3, 2],
    ci_star        = fit_res$Ci(0),
    A_transition   = fit_res$Photosyn(Ci=fit_res$Ci_transition)$ALEAF
  ))

  # If no points are J-limited, the return value for Jmax at leaf temperature
  # will be 1e6. Temperature corrections may lower this a bit.
  parameters[, 'Jmax'] <- if (fit_res$par[2, 1] > 1e5) {
    NA
  } else {
    fit_res$par[2, 1]
  }

  parameters[, 'Jmax_err'] <- if (fit_res$par[2, 1] > 1e5) {
    NA
  } else {
    fit_res$par[2, 2]
  }

  # Calculate J from J_max using the plantecophys function
  parameters[, 'J'] <-
    plantecophys:::Jfun(PPFD, ALPHA_J, parameters[, 'Jmax'], THETA_J)

  # Estimate standard error in J using the error for J_max
  parameters[, 'J_err'] <- 0.5 *
    (plantecophys:::Jfun(PPFD, ALPHA_J, parameters[, 'Jmax'] + parameters[, 'Jmax_err'], THETA_J) -
    plantecophys:::Jfun(PPFD, ALPHA_J, parameters[, 'Jmax'] - parameters[, 'Jmax_err'], THETA_J))

  # The value of `TPU` and its error will depend on `fitTPU`
  parameters[, 'TPU'] <- if (fit_res$fitTPU) {
    fit_res$par[4, 1]
  } else {
    NA
  }

  parameters[, 'TPU_err'] <- if (fit_res$fitTPU) {
    fit_res$par[4, 2]
  } else {
    NA
  }

  # Calculate upper and lower confidence interval limits for key parameters.
  # plantecophys reports standard errors (SE). The limits for a 95% confidence
  # interval are mean +/- 1.96 * SE.
  ci_f <- 1.96

  parameters[, 'J_lower']     <- parameters[, 'J']     - parameters[, 'J_err']     * ci_f
  parameters[, 'J_upper']     <- parameters[, 'J']     + parameters[, 'J_err']     * ci_f
  parameters[, 'Rd_lower']    <- parameters[, 'Rd']    - parameters[, 'Rd_err']    * ci_f
  parameters[, 'Rd_upper']    <- parameters[, 'Rd']    + parameters[, 'Rd_err']    * ci_f
  parameters[, 'TPU_lower']   <- parameters[, 'TPU']   - parameters[, 'TPU_err']   * ci_f
  parameters[, 'TPU_upper']   <- parameters[, 'TPU']   + parameters[, 'TPU_err']   * ci_f
  parameters[, 'Vcmax_lower'] <- parameters[, 'Vcmax'] - parameters[, 'Vcmax_err'] * ci_f
  parameters[, 'Vcmax_upper'] <- parameters[, 'Vcmax'] + parameters[, 'Vcmax_err'] * ci_f

  # Determine whether this is an "inadmissible fit".
  #
  # A note about the Ci transitions: In the fitaci code, we have:
  #
  #  datv <- data[data$Ci < citransition & data$Ci < citransition2,]
  #  datj <- data[data$Ci >= citransition & data$Ci < citransition2,]
  #  datp <- data[data$Ci >= citransition2,]
  #
  # So if citransition2 < citransition, no points are J-limited, and the switch
  # from Rubisco limitations to TPU limitations occurs at citransition2. If no
  # points are TPU limited, citransition2 will be NA; here we change it to 1e6
  # to avoid calculation errors.
  ci1 <- parameters[, 'Ci_transition']

  ci2 <- if (is.na(parameters[, 'Ci_transition2'])) {
    1e6
  } else {
    parameters[, 'Ci_transition2']
  }

  parameters[, 'admissible'] <-
    all(fits[fits[, 'Ci'] < ci1 & fits[, 'Ci'] < ci2, 'Ac_limiting']) &&
    all(fits[fits[, 'Ci'] >= ci1 & fits[, 'Ci'] < ci2, 'Aj_limiting']) &&
    all(fits[fits[, 'Ci'] >= ci2, 'Ap_limiting'])

  # Document the parameter units
  parameters <- document_variables(
    parameters,
    c('fit_c3_aci_plantecophys', 'Ci_transition',  'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'Ci_transition2', 'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'alpha',          'dimensionless'),
    c('fit_c3_aci_plantecophys', 'theta',          'dimensionless'),
    c('fit_c3_aci_plantecophys', 'Tcorrect',       ''),
    c('fit_c3_aci_plantecophys', 'Rd_measured',    ''),
    c('fit_c3_aci_plantecophys', 'GammaStar',      'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'Km',             'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'kminput',        ''),
    c('fit_c3_aci_plantecophys', 'gstarinput',     ''),
    c('fit_c3_aci_plantecophys', 'fitmethod',      ''),
    c('fit_c3_aci_plantecophys', 'citransition',   'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'gmeso',          'mol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'fitTPU',         ''),
    c('fit_c3_aci_plantecophys', 'alphag',         'dimensionless'),
    c('fit_c3_aci_plantecophys', 'Vcmax',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Vcmax_err',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'J',              'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'J_err',          'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Jmax',           'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Jmax_err',       'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Rd',             'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Rd_err',         'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'TPU',            'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'TPU_err',        'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'ci_star',        'micromol mol^(-1)'),
    c('fit_c3_aci_plantecophys', 'A_transition',   'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'J_lower',        'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'J_upper',        'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Rd_lower',       'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Rd_upper',       'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'TPU_lower',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'TPU_upper',      'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Vcmax_lower',    'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'Vcmax_upper',    'micromol m^(-2) s^(-1)'),
    c('fit_c3_aci_plantecophys', 'admissible',     'micromol m^(-2) s^(-1)')
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
