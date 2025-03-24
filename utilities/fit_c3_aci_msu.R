library(msuRACiFit)

# Here we assume that `defaults.R` has already been run

# This function is based on the example from the "Combining PhotoGEA With Other
# Packages" vignette from the PhotoGEA R package, which is available online at
# https://eloch216.github.io/PhotoGEA/articles/combining_with_other_packages.html

# Make a wrapper for `msuRACiFit::fitComplete`
fit_c3_aci_msu <- function(
    replicate_exdf,
    a_column_name = 'A',                     # micromol / m^2 / s
    gamma_star_column_name = 'Gamma_star',   # micromol / mol
    kc_column_name = 'Kc',                   # micromol / mol
    ko_column_name = 'Ko',                   # mmol / mol
    oxygen_column_name = 'oxygen',           # percent
    pci_column_name = 'Pci',                 # Pa
    pressure_column_name = 'total_pressure', # bar
    tleaf_column_name = 'TleafCnd',          # degrees C
    ... # Not used
)
{
    ### Check inputs
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_aci_msu requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]          <- 'micromol m^(-2) s^(-1)'
    required_variables[[gamma_star_column_name]] <- 'micromol mol^(-1)'
    required_variables[[kc_column_name]]         <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]]         <- 'mmol mol^(-1)'
    required_variables[[oxygen_column_name]]     <- 'percent'
    required_variables[[pci_column_name]]        <- 'Pa'
    required_variables[[pressure_column_name]]   <- 'bar'
    required_variables[[tleaf_column_name]]      <- 'degrees C'

    check_required_variables(replicate_exdf, required_variables)

    # Just keep the minimum necessary columns
    dataf <- replicate_exdf[, c(a_column_name, pci_column_name)]

    fit_res <- msuRACiFit::fitComplete(
        tibble::as_tibble(dataf, .name_repair = 'unique'),
        name_assimilation = a_column_name,
        name_ci = pci_column_name,
        gammastar = mean(replicate_exdf[, gamma_star_column_name] * 1e-6 * replicate_exdf[, pressure_column_name] * 1e5), # Pa
        O2 = mean(replicate_exdf[, oxygen_column_name]),               # percent
        tleaf = mean(replicate_exdf[, tleaf_column_name]),             # degrees C
        pressure = mean(replicate_exdf[, pressure_column_name]) * 100, # kPa
        forceValues = c(NA, NA, NA, Inf, NA, NA, NA),                  # Don't fit gm; just set to Inf
        Kc = mean(replicate_exdf[, kc_column_name] * 1e-6 * replicate_exdf[, pressure_column_name] * 1e5), # Pa
        Ko = mean(replicate_exdf[, ko_column_name] * 1e-3 * replicate_exdf[, pressure_column_name] * 1e2)  # kPa
    )

    ## Try to calculate confidence intervals
    #tryCatch(
    #    {str(confint(fit_res[[2]]))},
    #    condition = function(c) {
    #        # do nothing
    #    }
    #)

    # Get the identifier columns from the original exdf object
    replicate_identifiers <- identifier_columns(replicate_exdf)

    # Get a data frame with the fitted values of assimilation and convert it to an
    # `exdf` object, setting the category to `fit_c3_aci_photosynthesis` and
    # specifying the units for each column
    fits <- exdf(as.data.frame(fit_res[[3]]))

    # Append the identifier columns to the fits
    fits <- cbind(replicate_identifiers, fits)

    # Rename a few columns
    colnames(fits)[colnames(fits) == 'Rubisco Limited'] <- 'Ac'
    colnames(fits)[colnames(fits) == 'ET Limited']      <- 'Aj'
    colnames(fits)[colnames(fits) == 'TPU Limited']     <- 'Ap'

    # We set gmc = Inf, so Cc = Ci. The results are in Pa by default, so we need
    # to convert back to micromol / mol for consistency when plotting.
    fits[, 'Ci'] <- fits[, 'Cc'] / (mean(replicate_exdf[, pressure_column_name]) * 1e5) * 1e6

    # Add an A_fit column
    fits[, 'A_fit'] <- NA

    for (i in seq_len(nrow(fits))) {
        lp <- fits[i, 'Limiting process']
        fits[i, 'A_fit'] <- if (lp == 1) {
            fits[i, 'Ac']
        } else if (lp == 2) {
            fits[i, 'Aj']
        } else {
            fits[i, 'Ap']
        }
    }

    # Document the fit units
    fits <- document_variables(
        fits,
        c('fit_c3_aci_msu', 'A_fit', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'Ac',     'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'Aj',     'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'Ap',     'micromol m^(-2) s^(-1)')
    )

    # Create an exdf object with the parameter values that are included in the
    # fitting result.
    parameters <- exdf(as.data.frame(fit_res[[1]]))

    # Append the identifier columns to the parameters
    parameters <- cbind(replicate_identifiers, parameters)

    # Calculate upper and lower confidence interval limits for key parameters.
    # msuRACiFit does not return any uncertainty estimate, so these will all be
    # NA.
    parameters[, 'J_lower']     <- NA
    parameters[, 'J_upper']     <- NA
    parameters[, 'rL_lower']    <- NA
    parameters[, 'rL_upper']    <- NA
    parameters[, 'TPU_lower']   <- NA
    parameters[, 'TPU_upper']   <- NA
    parameters[, 'VcMax_lower'] <- NA
    parameters[, 'VcMax_upper'] <- NA

    # Document the parameter units
    parameters <- document_variables(
        parameters,
        c('fit_c3_aci_msu', 'VcMax',       'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'J',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'TPU',         'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'gm',          ''),
        c('fit_c3_aci_msu', 'rL',          'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'ag',          'dimensionless'),
        c('fit_c3_aci_msu', 'as',          'dimensionless'),
        c('fit_c3_aci_msu', 'J_lower',     'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'J_upper',     'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'rL_lower',    'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'rL_upper',    'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'TPU_lower',   'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'TPU_upper',   'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'VcMax_lower', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci_msu', 'VcMax_upper', 'micromol m^(-2) s^(-1)')
    )

    # `msuRACiFit::fitComplete` returns absurdly high values of Ap
    # when it cannot be the limiting rate. These can cause problems when
    # plotting the results, so we replace any values above 500 micromol /
    # m^2 / s with NA.
    fits[fits[, 'Ap'] > 500, 'Ap'] <- NA

    # Identify limiting process, reporting it in a standardized way
    fits <- identify_c3_limiting_processes(
        fits,
        'A_fit',
        'Ac',
        'Aj',
        'Ap'
    )

    # Attach the residual stats to the parameters
    nparam <- 6 # VcMax, J, TPU, (not gm), Rd, aG, aS

    parameters <- cbind(
        parameters,
        residual_stats(
            fits[, 'residual'],
            fits$units[['A_fit']],
            nparam
        )
    )

    # Return a list of two data frames: `fits` and `parameters`
    return(list(
      fits = fits,
      parameters = parameters
    ))
}
