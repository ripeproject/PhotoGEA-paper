###               ###
### PRELIMINARIES ###
###               ###

# Load libraries
library(PhotoGEA)
library(lattice)
library(lhs)

# Clear workspace
to_remove <- ls()
to_remove <- to_remove[to_remove != 'ALPHA_TYPE']
rm(list=to_remove)

# Load some helping functions
source(file.path('utilities', 'output_tools.R'))
source(file.path('utilities', 'defaults.R'))

source(file.path('utilities', 'plot_synthetic_curves.R'))
source(file.path('utilities', 'extract_fit_param.R'))
source(file.path('utilities', 'fit_c3_aci_plantecophys.R'))
source(file.path('utilities', 'fit_c3_aci_photosynthesis.R'))
source(file.path('utilities', 'fit_c3_aci_msu.R'))

# Load data
load(file.path(OUTPUT_DIR, RDATA_DIR, 'tobacco_aci_raw_data.RData')) # load licor_data; created by `load_tobacco_data.R`

# Choose some settings
SAVE_TO_PDF <- TRUE
MAKE_NEW_CALCULATIONS <- TRUE

REMOVE_BAD_CURVES <- TRUE

# Related to curve generation
NCURVES <- 200

J_RANGE     <- c(20, 200)
RL_RANGE    <- c(0.5, 2.5)
TP_RANGE    <- c(2, 18)
VCMAX_RANGE <- c(20, 120)

# Check ALPHA_TYPE, which must be defined before calling this script
if (!tolower(ALPHA_TYPE) %in% c('none', 'alpha_old', 'alpha_gs')) {
    stop('The value of ALPHA_TYPE must be none, alpha_old, or alpha_gs')
}

BASE_NAME <- paste0('synthetic_aci_', ALPHA_TYPE)

# Set a seed number because this script uses randomness
set.seed(123)

###
### CHOOSE CI VALUES BASED ON MEASURED TOBACCO A-CI CURVES
###

# Get average Ci sequence from measured tobacco A-Ci curves
Ci_seq <- round(as.numeric(tapply(licor_data[, 'Ci'], licor_data[, 'CO2_r_sp'], mean)))

###
### GENERATE SYNTHETIC A-CI CURVES USING THE FVCB MODEL
###

# Define a helping function for generating a synthetic A-Ci curve
generate_curve <- function(
    Ci_vals,
    alpha_g,
    alpha_old,
    alpha_s,
    J_at_25,
    RL_at_25,
    Tp_at_25,
    Vcmax_at_25,
    cid
)
{
    # Generate an initial exdf object
    curve_exdf <- exdf(
        data.frame(
            Ci = Ci_vals,
            Cc = Ci_vals,
            Ca = 420,
            gmc = Inf,
            TleafCnd = 25,
            total_pressure = 1,
            oxygen = 21,
            Qin = 1800, #plantecophys requires this
            curve_identifier = cid
        ),
        units = data.frame(
            Ci = 'micromol mol^(-1)',
            Cc = 'micromol mol^(-1)',
            Ca = 'micromol mol^(-1)',
            gmc = 'mol m^(-2) s^(-1) bar^(-1)',
            TleafCnd = 'degrees C',
            total_pressure = 'bar',
            oxygen = 'percent',
            Qin = 'micromol m^(-2) s^(-1)',
            curve_identifier = ''
        )
    )

    # Convert Ci to Pa (msuRACiFit requires this). Ci * total_pressure is in
    # microbar. 1 microbar = 0.1 Pa.
    curve_exdf[, 'Pci'] <- curve_exdf[, 'Ci'] * curve_exdf[, 'total_pressure'] * 0.1

    # Calculate temperature-dependent values of key photosynthetic parameters
    curve_exdf <-
        calculate_temperature_response(curve_exdf, TEMPERATURE_PARAM)

    # Calculate assimilation rates
    assim <- calculate_c3_assimilation(
        curve_exdf,
        alpha_g = alpha_g,
        alpha_old = alpha_old,
        alpha_s = alpha_s,
        alpha_t = 0,
        Gamma_star = '', # use value from column
        J_at_25 = J_at_25,
        RL_at_25 = RL_at_25,
        Tp_at_25 = Tp_at_25,
        Vcmax_at_25 = Vcmax_at_25
    )

    # Add new columns to the curve:
    # - The net CO2 assimilation rate
    # - The number of points where each potential process is limiting
    # - The actual parameter value used to generate the curve
    curve_exdf[, 'A'] <- assim[, 'An']
    curve_exdf[, 'true_n_Ac_limiting'] <- PhotoGEA:::n_C3_A_limiting(assim, 'An', 'Ac')
    curve_exdf[, 'true_n_Aj_limiting'] <- PhotoGEA:::n_C3_A_limiting(assim, 'An', 'Aj')
    curve_exdf[, 'true_n_Ap_limiting'] <- PhotoGEA:::n_C3_A_limiting(assim, 'An', 'Ap')
    curve_exdf[, 'true_alpha_g']       <- alpha_g
    curve_exdf[, 'true_alpha_old']     <- alpha_old
    curve_exdf[, 'true_alpha_s']       <- alpha_s
    curve_exdf[, 'true_J_at_25']       <- J_at_25
    curve_exdf[, 'true_RL_at_25']      <- RL_at_25
    curve_exdf[, 'true_Tp_at_25']      <- Tp_at_25
    curve_exdf[, 'true_Vcmax_at_25']   <- Vcmax_at_25

    # Return the exdf
    document_variables(
        curve_exdf,
        c('', 'A',                assim$units[['An']]),
        c('', 'true_alpha_g',     'dimensionless'),
        c('', 'true_alpha_old',   'dimensionless'),
        c('', 'true_alpha_s',     'dimensionless'),
        c('', 'true_J_at_25',     'micromol m^(-2) s^(-1)'),
        c('', 'true_RL_at_25',    'micromol m^(-2) s^(-1)'),
        c('', 'true_Tp_at_25',    'micromol m^(-2) s^(-1)'),
        c('', 'true_Vcmax_at_25', 'micromol m^(-2) s^(-1)'),
        c('', 'Pci',              'Pa')
    )
}

# Use Latin hypercube sampling to choose parameter sets
hypercube <- randomLHS(NCURVES, 7)

j_delta     <- max(J_RANGE) - min(J_RANGE)
j_min       <- min(J_RANGE)
rl_delta    <- max(RL_RANGE) - min(RL_RANGE)
rl_min      <- min(RL_RANGE)
tp_delta    <- max(TP_RANGE) - min(TP_RANGE)
tp_min      <- min(TP_RANGE)
vcmax_delta <- max(VCMAX_RANGE) - min(VCMAX_RANGE)
vcmax_min   <- min(VCMAX_RANGE)

curve_data <- do.call(rbind, lapply(seq_len(NCURVES), function(i) {
    # Get alpha values as necessary
    alpha_old <- switch(
        ALPHA_TYPE,
        none = 0,
        alpha_old = hypercube[i, 1],
        alpha_gs = 0
    )

    alpha_g <- switch(
        ALPHA_TYPE,
        none = 0,
        alpha_old = 0,
        alpha_gs = hypercube[i, 2]
    )

    alpha_s <- switch(
        ALPHA_TYPE,
        none = 0,
        alpha_old = 0,
        alpha_gs = hypercube[i, 3] * alpha_g * 3 / 4
    )

    # Generate a curve
    generate_curve(
        Ci_seq,
        alpha_g     = alpha_s,
        alpha_old   = alpha_old,
        alpha_s     = alpha_g,
        J_at_25     = j_min     + j_delta     * hypercube[i, 4],
        RL_at_25    = rl_min    + rl_delta    * hypercube[i, 5],
        Tp_at_25    = tp_min    + tp_delta    * hypercube[i, 6],
        Vcmax_at_25 = vcmax_min + vcmax_delta * hypercube[i, 7],
        paste('Curve', i)
    )
}))

curve_data[, 'curve_identifier'] <- factor(
    curve_data[, 'curve_identifier'],
    levels = paste('Curve', seq_len(NCURVES))
)

# If ALPHA_TYPE is alpha_old, make an image showing selected curves with
# various characteristics
if (ALPHA_TYPE == 'alpha_old') {
    plot_synthetic_curves(curve_data, BASE_NAME)
}

###
### FIT THE CURVES USING EACH R PACKAGE
###
if (MAKE_NEW_CALCULATIONS) {
    # Fit the A-Ci curves using PhotoGEA
    FIT_OPTIONS <- if (ALPHA_TYPE == 'alpha_gs') {
        list(alpha_g = 'fit', alpha_s = 'fit', alpha_old = 0)
    } else {
        list()
    }

    OPTIMIZER <- if (ALPHA_TYPE == 'alpha_gs') {
        optimizer_deoptim(C3_ACI_ITERMAX + 300)
    } else {
        optimizer_deoptim(C3_ACI_ITERMAX)
    }

    c3_aci_results <- consolidate(by(
        curve_data,                       # The `exdf` object containing the curves
        curve_data[, 'curve_identifier'], # A factor used to split `curve_data` into chunks
        fit_c3_aci,                       # The function to apply to each chunk of `curve_data`
        optim_fun = OPTIMIZER,
        fit_options = FIT_OPTIONS,
        relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
        hard_constraints = HARD_CONSTRAINTS,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    ))

    # Fit the A-Ci curves using plantecophys
    plantecophys_results <- consolidate(by(
        curve_data,                       # The `exdf` object containing the curves
        curve_data[, 'curve_identifier'], # A factor used to split `curve_data` into chunks
        fit_c3_aci_plantecophys           # The function to apply to each chunk of `curve_data`
    ))

    # Fit the A-Ci curves using photosynthesis
    photosynthesis_results <- consolidate(by(
        curve_data,                       # The `exdf` object containing the curves
        curve_data[, 'curve_identifier'], # A factor used to split `curve_data` into chunks
        fit_c3_aci_photosynthesis         # The function to apply to each chunk of `curve_data`
    ))

    # Fit the A-Ci curves using msuRACiFit
    msu_results <- consolidate(by(
        curve_data,                       # The `exdf` object containing the curves
        curve_data[, 'curve_identifier'], # A factor used to split `curve_data` into chunks
        fit_c3_aci_msu                    # The function to apply to each chunk of `curve_data`
    ))

    # Save the results
    write.csv.exdf(
        c3_aci_results$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_photogea_parameters.csv'))
    )
    write.csv.exdf(
        c3_aci_results$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_photogea_fits.csv'))
    )

    write.csv.exdf(
        plantecophys_results$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_plantecophys_parameters.csv'))
    )

    write.csv.exdf(
        plantecophys_results$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_plantecophys_fits.csv'))
    )

    write.csv.exdf(
        photosynthesis_results$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_photosynthesis_parameters.csv'))
    )

    write.csv.exdf(
        photosynthesis_results$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_photosynthesis_fits.csv'))
    )

    write.csv.exdf(
        msu_results$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_msu_parameters.csv'))
    )

    write.csv.exdf(
        msu_results$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_msu_fits.csv'))
    )

    save(
        file_paths, c3_aci_results, plantecophys_results,
        photosynthesis_results, msu_results,
        file = file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData'))
    )
} else {
    load(file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData')))
}

###                                 ###
### COMPARE RESULTS ACROSS PACKAGES ###
###                                 ###

# Extract some of the results from the different fitting methods and save them
fit_param <- extract_fit_param(
    c3_aci_results$parameters,
    NULL,
    plantecophys_results$parameters,
    photosynthesis_results$parameters,
    msu_results$parameters,
    NULL,
    include_true_values = TRUE,
    include_ppfd = FALSE
)

# Identify false positives (a parameter is estimated when it shouldn't be) and
# false negatives (a parameter is not estimated when it should be)
fit_param <- within(fit_param, {
    false_positive_J     = true_n_Aj_limiting == 0 & !is.na(J)
    false_positive_Tp    = true_n_Ap_limiting == 0 & !is.na(Tp)
    false_positive_Vcmax = true_n_Ac_limiting == 0 & !is.na(Vcmax)

    false_negative_J     = true_n_Aj_limiting > 0 & is.na(J)
    false_negative_Tp    = true_n_Ap_limiting > 0 & is.na(Tp)
    false_negative_Vcmax = true_n_Ac_limiting > 0 & is.na(Vcmax)
})

# Calculate absolute and relative errors
fit_param <- within(fit_param, {
    J_diff_abs     = J - true_J
    RL_diff_abs    = RL - true_RL
    Tp_diff_abs    = Tp - true_Tp
    Vcmax_diff_abs = Vcmax - true_Vcmax

    J_diff_abs_lower     = J_lower - true_J
    RL_diff_abs_lower    = RL_lower - true_RL
    Tp_diff_abs_lower    = Tp_lower - true_Tp
    Vcmax_diff_abs_lower = Vcmax_lower - true_Vcmax

    J_diff_abs_upper     = J_upper - true_J
    RL_diff_abs_upper    = RL_upper - true_RL
    Tp_diff_abs_upper    = Tp_upper - true_Tp
    Vcmax_diff_abs_upper = Vcmax_upper - true_Vcmax

    J_diff_rel     = 100 * J_diff_abs / true_J
    RL_diff_rel    = 100 * RL_diff_abs / true_RL
    Tp_diff_rel    = 100 * Tp_diff_abs / true_Tp
    Vcmax_diff_rel = 100 * Vcmax_diff_abs / true_Vcmax
})

# Indicate whether the confidence interval includes the true value (to within
# a small tolerance, allowing for roundoff errors)
tol <- 1e-2
fit_param <- within(fit_param, {
    true_J_in_confidence_interval     = !is.na(J_lower)     & !is.na(J_upper)     & J_lower     <= true_J     * (1.0 + tol) & J_upper     >= true_J     * (1.0 - tol)
    true_RL_in_confidence_interval    = !is.na(RL_lower)    & !is.na(RL_upper)    & RL_lower    <= true_RL    * (1.0 + tol) & RL_upper    >= true_RL    * (1.0 - tol)
    true_Tp_in_confidence_interval    = !is.na(Tp_lower)    & !is.na(Tp_upper)    & Tp_lower    <= true_Tp    * (1.0 + tol) & Tp_upper    >= true_Tp    * (1.0 - tol)
    true_Vcmax_in_confidence_interval = !is.na(Vcmax_lower) & !is.na(Vcmax_upper) & Vcmax_lower <= true_Vcmax * (1.0 + tol) & Vcmax_upper >= true_Vcmax * (1.0 - tol)
})

write.csv(
    fit_param,
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_fit_param.csv')),
    row.names = FALSE
)

# Create and save summary table
packages <- c(
    'PhotoGEA',
    'plantecophys',
    'photosynthesis',
    'msuRACiFit'
)

pinfo <- list(
    c('J',     'true_n_Aj_limiting'),
    c('Tp',    'true_n_Ap_limiting'),
    c('Vcmax', 'true_n_Ac_limiting'),
    c('RL',    '')
)

get_param_info <- function(
    pname,  # parameter name
    ct,     # curve type
    gtz,    # greater than zero?
    pkg,    # package
    en,     # error name
    total_n # total number of curves of this type
)
{
    conf_i_name <- paste0('true_', pname, '_in_confidence_interval')

    if (ct == '') {
        dataf <- fit_param[fit_param$package == pkg & !is.na(fit_param[, en]), ]

        n_est = nrow(dataf[!is.na(dataf[, en]), ])

        n_gci <- sum(dataf[, conf_i_name])

        data.frame(
            parameter = pname,
            curve_type = 'all',
            n_curves = total_n,
            package = pkg,
            n_estimates = n_est,
            error_percent = (total_n - n_est) / total_n * 100,
            n_good_CI = n_gci,
            good_CI_percent = n_gci / total_n * 100,
            RMSE = sqrt(mean(dataf[, en]^2, na.rm = TRUE))
        )
    } else {
        if (gtz) {
            dataf <- fit_param[fit_param$package == pkg & fit_param[, ct] > 0, ]

            n_est = nrow(dataf[!is.na(dataf[, en]), ])

            n_gci <- sum(dataf[, conf_i_name])

            data.frame(
                parameter = pname,
                curve_type = paste(ct, '> 0'),
                n_curves = total_n,
                package = pkg,
                n_estimates = n_est,
                error_percent = (total_n - n_est) / total_n * 100,
                n_good_CI = n_gci,
                good_CI_percent = n_gci / total_n * 100,
                RMSE = sqrt(mean(dataf[, en]^2, na.rm = TRUE))
            )
        } else {
            dataf <- fit_param[fit_param$package == pkg & fit_param[, ct] == 0, ]

            n_est = nrow(dataf[!is.na(dataf[, en]), ])

            n_gci <- sum(dataf[, conf_i_name])

            data.frame(
                parameter = pname,
                curve_type = paste(ct, '== 0'),
                n_curves = total_n,
                package = pkg,
                n_estimates = n_est,
                error_percent = n_est / total_n * 100,
                n_good_CI = n_gci,
                good_CI_percent = n_gci / total_n * 100,
                RMSE = NA
            )
        }
    }
}

fit_summary <- do.call(rbind, lapply(pinfo, function(x) {
    pname <- x[1]
    ctype <- x[2]
    err_name <- paste0(pname, '_diff_abs')

    if (ctype == '') {

        do.call(rbind, lapply(packages, function(pkg) {
            get_param_info(
                pname,                                           # parameter name
                ctype,                                           # curve type
                TRUE,                                            # greater than zero?
                pkg,                                             # package
                err_name,                                        # error name
                length(unique(curve_data[, 'curve_identifier'])) # total number of curves of this type
            )
        }))

    } else {

        rbind(
            do.call(rbind, lapply(packages, function(pkg) {
                get_param_info(
                    pname,                                                                  # parameter name
                    ctype,                                                                  # curve type
                    TRUE,                                                                   # greater than zero?
                    pkg,                                                                    # package
                    err_name,                                                               # error name
                    length(unique(curve_data[curve_data[, ctype] > 0, 'curve_identifier'])) # total number of curves of this type
                )
            })),
            do.call(rbind, lapply(packages, function(pkg) {
                get_param_info(
                    pname,                                                                   # parameter name
                    ctype,                                                                   # curve type
                    FALSE,                                                                   # greater than zero?
                    pkg,                                                                     # package
                    err_name,                                                                # error name
                    length(unique(curve_data[curve_data[, ctype] == 0, 'curve_identifier'])) # total number of curves of this type
                )
            }))
        )
    }
}))

write.csv(
    fit_summary,
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_fit_summary.csv')),
    row.names = FALSE
)

# Ensure packages are plotted in the order we want
fit_param$package <- factor(fit_param$package, levels = rev(packages))

# Replace Infinite upper limits by finite values so they can be plotted
fit_param[is.infinite(fit_param$J_diff_abs_upper), 'J_diff_abs_upper'] <- 1e6
fit_param[is.infinite(fit_param$RL_diff_abs_upper), 'RL_diff_abs_upper'] <- 1e6
fit_param[is.infinite(fit_param$Tp_diff_abs_upper), 'Tp_diff_abs_upper'] <- 1e6
fit_param[is.infinite(fit_param$Vcmax_diff_abs_upper), 'Vcmax_diff_abs_upper'] <- 1e6

# Plot counts of false positives
pdf_print(
    barchart(
        with(fit_param, {tapply(false_positive_J, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of false positive J estimates'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_false_positive_J.pdf'))
)

pdf_print(
    barchart(
        with(fit_param, {tapply(false_positive_Tp, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of false positive Tp estimates'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_false_positive_Tp.pdf'))
)

pdf_print(
    barchart(
        with(fit_param, {tapply(false_positive_Vcmax, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of false positive Vcmax estimates'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_false_positive_Vcmax.pdf'))
)

# Plot counts of false negatives
pdf_print(
    barchart(
        with(fit_param, {tapply(false_negative_J, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of false negative J estimates'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_false_negative_J.pdf'))
)

pdf_print(
    barchart(
        with(fit_param, {tapply(false_negative_Tp, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of false negative Tp estimates'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_false_negative_Tp.pdf'))
)

pdf_print(
    barchart(
        with(fit_param, {tapply(false_negative_Vcmax, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of false negative Vcmax estimates'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_false_negative_Vcmax.pdf'))
)

# Plot counts of confidence interval that include the true value
pdf_print(
    barchart(
        with(fit_param, {tapply(true_J_in_confidence_interval, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of J confidence intervals that include the true value'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_good_confidence_J.pdf'))
)

pdf_print(
    barchart(
        with(fit_param, {tapply(true_RL_in_confidence_interval, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of RL confidence intervals that include the true value'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_good_confidence_RL.pdf'))
)

pdf_print(
    barchart(
        with(fit_param, {tapply(true_Tp_in_confidence_interval, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of Tp confidence intervals that include the true value'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_good_confidence_Tp.pdf'))
)

pdf_print(
    barchart(
        with(fit_param, {tapply(true_Vcmax_in_confidence_interval, package, sum)}),
        horizontal = FALSE,
        ylab = 'Number of Vcmax confidence intervals that include the true value'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_good_confidence_Vcmax.pdf'))
)

# Plot relative error distributions
pdf_print(
    bwplot(
        J_diff_rel ~ factor(true_n_Aj_limiting > 0) | package,
        data = fit_param,
        ylab = 'Relative error of J estimate (percent)',
        xlab = 'true_n_Aj_limiting > 0'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_distribution_J.pdf'))
)

pdf_print(
    bwplot(
        RL_diff_rel ~ package,
        data = fit_param,
        ylab = 'Relative error of RL estimate (percent)'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_distribution_RL.pdf'))
)

pdf_print(
    bwplot(
        Tp_diff_rel ~ factor(true_n_Ap_limiting > 0) | package,
        data = fit_param,
        ylab = 'Relative error of Tp estimate (percent)',
        xlab = 'true_n_Ap_limiting > 0'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_distribution_Tp.pdf'))
)

pdf_print(
    bwplot(
        Vcmax_diff_rel ~ factor(true_n_Ac_limiting > 0) | package,
        data = fit_param,
        ylab = 'Relative error of Vcmax estimate (percent)',
        xlab = 'true_n_Ac_limiting > 0'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_distribution_Vcmax.pdf'))
)

# Plot absolute parameter estimate errors
pdf_print(
    xyplot(
        J_diff_abs ~ true_J | package,
        group = true_n_Aj_limiting > 0,
        data = fit_param,
        type = 'p',
        pch = 16,
        auto = TRUE,
        grid = TRUE,
        xlab = 'True J [ micromol m^(-2) s^(-1) ]',
        ylab = 'Estimated J - True J [ micromol m^(-2) s^(-1) ]'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_J.pdf'))
)

pdf_print(
    xyplot(
        RL_diff_abs ~ true_RL | package,
        data = fit_param,
        type = 'p',
        pch = 16,
        auto = TRUE,
        grid = TRUE,
        xlab = 'True RL [ micromol m^(-2) s^(-1) ]',
        ylab = 'Estimated RL - True RL [ micromol m^(-2) s^(-1) ]'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_RL.pdf'))
)

pdf_print(
    xyplot(
        Tp_diff_abs ~ true_Tp | package,
        group = true_n_Ap_limiting > 0,
        data = fit_param,
        type = 'p',
        pch = 16,
        auto = TRUE,
        grid = TRUE,
        xlab = 'True Tp [ micromol m^(-2) s^(-1) ]',
        ylab = 'Estimated Tp - True Tp [ micromol m^(-2) s^(-1) ]'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_Tp.pdf'))
)

pdf_print(
    xyplot(
        Vcmax_diff_abs ~ true_Vcmax | package,
        group = true_n_Ac_limiting > 0,
        data = fit_param,
        type = 'p',
        pch = 16,
        auto = TRUE,
        grid = TRUE,
        xlab = 'True Vcmax [ micromol m^(-2) s^(-1) ]',
        ylab = 'Estimated Vcmax - True Vcmax [ micromol m^(-2) s^(-1) ]'
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_Vcmax.pdf'))
)

# Plot parameter estimates (detailed)
pdf_print(
    xyplot(
        true_J + J ~ true_J | factor(true_n_Aj_limiting) * package,
        data = fit_param,
        type = 'p',
        par.settings = list(
            superpose.symbol = list(pch = c(16, 1), col = c('gray', 'black'))
        ),
        scales = list(alternating = TRUE),
        auto.key = list(space = 'right', text = c('True', 'Estimated')),
        xlab = 'True J [ micromol m^(-2) s^(-1) ]',
        ylab = 'J [ micromol m^(-2) s^(-1) ]',
        main = 'Grouped by number of Aj-limited points'
    ),
    width = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_J_detailed.pdf'))
)

pdf_print(
    xyplot(
        true_RL + RL ~ true_RL | factor(true_n_Ac_limiting) * package,
        data = fit_param,
        type = 'p',
        par.settings = list(
            superpose.symbol = list(pch = c(16, 1), col = c('gray', 'black'))
        ),
        scales = list(alternating = TRUE),
        auto.key = list(space = 'right', text = c('True', 'Estimated')),
        xlab = 'True RL [ micromol m^(-2) s^(-1) ]',
        ylab = 'RL [ micromol m^(-2) s^(-1) ]',
        main = 'Grouped by number of Ac-limited points'
    ),
    width = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_RL_detailed.pdf'))
)

pdf_print(
    xyplot(
        true_Tp + Tp ~ true_Tp | factor(true_n_Ap_limiting) * package,
        data = fit_param,
        type = 'p',
        par.settings = list(
            superpose.symbol = list(pch = c(16, 1), col = c('gray', 'black'))
        ),
        scales = list(alternating = TRUE),
        auto.key = list(space = 'right', text = c('True', 'Estimated')),
        xlab = 'True Tp [ micromol m^(-2) s^(-1) ]',
        ylab = 'Tp [ micromol m^(-2) s^(-1) ]',
        main = 'Grouped by number of Ap-limited points'
    ),
    width = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_Tp_detailed.pdf'))
)

pdf_print(
    xyplot(
        true_Vcmax + Vcmax ~ true_Vcmax | factor(true_n_Ac_limiting) * package,
        data = fit_param,
        type = 'p',
        par.settings = list(
            superpose.symbol = list(pch = c(16, 1), col = c('gray', 'black'))
        ),
        scales = list(alternating = TRUE),
        auto.key = list(space = 'right', text = c('True', 'Estimated')),
        xlab = 'True Vcmax [ micromol m^(-2) s^(-1) ]',
        ylab = 'Vcmax [ micromol m^(-2) s^(-1) ]',
        main = 'Grouped by number of Ac-limited points'
    ),
    width = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_error_Vcmax_detailed.pdf'))
)

# Plot fit quality indicators
pdf_print(
    bwplot(
        RMSE ~ package,
        data = fit_param
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_RMSE.pdf'))
)

pdf_print(
    bwplot(
        AIC ~ package,
        data = fit_param
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_AIC.pdf'))
)
