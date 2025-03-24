###               ###
### PRELIMINARIES ###
###               ###

# Load libraries
library(PhotoGEA)
library(lattice)
library(latticeExtra)
library(RColorBrewer)

# Clear workspace
rm(list=ls())

# Load some helping functions
source(file.path('utilities', 'output_tools.R'))
source(file.path('utilities', 'defaults.R'))

source(file.path('utilities', 'extract_fit_param.R'))
source(file.path('utilities', 'fit_c3_aci_plantecophys.R'))
source(file.path('utilities', 'fit_c3_aci_photosynthesis.R'))
source(file.path('utilities', 'fit_c3_aci_msu.R'))
source(file.path('utilities', 'plot_aci_errors.R'))
source(file.path('utilities', 'plot_c3_aci_fits.R'))
source(file.path('utilities', 'plot_c3_aci_fit_comparison.R'))
source(file.path('utilities', 'param_lim_mean.R'))
source(file.path('utilities', 'plot_parameters.R'))
source(file.path('utilities', 'plot_parameters_avg.R'))
source(file.path('utilities', 'plot_parameter_comparisons.R'))
source(file.path('utilities', 'plot_c3_aci_param_comparison.R'))

# Load data
load(file.path(OUTPUT_DIR, RDATA_DIR, 'tobacco_aci_raw_data.RData')) # load licor_data; created by `load_tobacco_data.R`

# Choose some settings
SAVE_TO_PDF <- TRUE
MAKE_NEW_CALCULATIONS <- TRUE

MAKE_VALIDATION_PLOTS <- TRUE
MAKE_ERROR_PLOTS <- TRUE
MAKE_PARAMETER_PLOTS <- TRUE

RESTRICT_PPFD_AVG <- TRUE

INDIVIDUAL_CURVE_ID <- '800 - wt-5 - mcgrath1'

BASE_NAME <- 'tobacco_aci'

A_LIM     <- c(-5, 50)
AIC_LIM   <- c(-50, 100)
ALPHA_LIM <- c(0, 2)
CI_LIM    <- c(-100, 1500)
GSW_LIM   <- c(0, 1)
J_LIM     <- c(0, 250)
QIN_LIM   <- c(-50, 1550)
RL_LIM    <- c(0, 5)
RMSE_LIM  <- 10^c(-1.5, 2)
TP_LIM    <- c(0, 15)
VCMAX_LIM <- c(0, 130)

if (MAKE_NEW_CALCULATIONS) {
    ###                                ###
    ### FIT A-CI CURVES USING PHOTOGEA ###
    ###                                ###

    # Calculate total pressure in the Licor chamber
    licor_data <- calculate_total_pressure(licor_data)

    # Calculate temperature-dependent values of C3 photosynthetic parameters
    licor_data <- calculate_temperature_response(licor_data, TEMPERATURE_PARAM)

    # Set a seed number because the deoptiom algorithm uses randomness
    set.seed(123)

    # Fit the A-Ci curves. Normally we would set `remove_unreliable_param` to
    # 2 (the default), but we need to have all of the best-fit values
    # (including unreliable ones) to generate the likelihood plots below.
    c3_aci_results <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci,                       # The function to apply to each chunk of `licor_data`
        optim_fun = optimizer_deoptim(C3_ACI_ITERMAX),
        fit_options = FIT_OPTIONS,
        Ca_atmospheric = CA_ATMOSPHERIC,
        atp_use = ATP_USE,
        nadph_use = NADPH_USE,
        relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
        hard_constraints = HARD_CONSTRAINTS,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 0
    ))

    # Fit the A-Ci curves using the new TPU model
    c3_aci_results_new_alpha <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci,                       # The function to apply to each chunk of `licor_data`
        optim_fun = optimizer_deoptim(C3_ACI_ITERMAX + 300),
        fit_options = list(alpha_old = 0, alpha_g = 'fit', alpha_s = 'fit'),
        Ca_atmospheric = CA_ATMOSPHERIC,
        atp_use = 4,
        nadph_use = 8,
        relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
        hard_constraints = HARD_CONSTRAINTS,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    ))

    ###                                    ###
    ### FIT A-CI CURVES USING PLANTECOPHYS ###
    ###                                    ###

    # Fit the A-Ci curves using plantecophys::fitaci
    plantecophys_results <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci_plantecophys           # The function to apply to each chunk of `licor_data`
    ))

    ###                                      ###
    ### FIT A-CI CURVES USING PHOTOSYNTHESIS ###
    ###                                      ###

    # Fit the A-Ci curves using plantecophys::fitaci
    photosynthesis_results <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci_photosynthesis         # The function to apply to each chunk of `licor_data`
    ))

    ###                                  ###
    ### FIT A-CI CURVES USING MSURACIFIT ###
    ###                                  ###

    # Fit the A-Ci curves using plantecophys::fitaci
    msu_results <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci_msu                    # The function to apply to each chunk of `licor_data`
    ))

    ###                                            ###
    ### FIT A-CI CURVES WITH A DIFFERENT OPTIMIZER ###
    ###                                            ###

    # Fit the A-Ci curves using stats::nlminb
    c3_aci_results_deriv <- consolidate(by(
        licor_data,                       # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
        fit_c3_aci,                       # The function to apply to each chunk of `licor_data`
        optim_fun = optimizer_nlminb(1e-10),
        fit_options = FIT_OPTIONS,
        Ca_atmospheric = CA_ATMOSPHERIC,
        atp_use = ATP_USE,
        nadph_use = NADPH_USE,
        relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
        hard_constraints = HARD_CONSTRAINTS,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    ))

    # Save results
    by(licor_data, licor_data[, 'curve_identifier'], function(x) {
        cid <- x[1, 'curve_identifier']

        col_to_keep <- c(
            'curve_identifier', 'replicate', 'instrument', 'ppfd',
            'A', 'Ci', 'CO2_r_sp', 'TleafCnd', 'total_pressure'
        )

        x_sub <- x[, col_to_keep, TRUE]

        x_sub[, 'total_pressure'] <- x_sub[, 'total_pressure'] * 100
        x_sub$units$total_pressure <- 'kPa'

        write.csv(
            x_sub,
            file = file.path(OUTPUT_DIR, SPREADSHEET_DIR, paste0(cid, '.csv')),
            row.names = FALSE
        )
    })

    write.csv.exdf(
        c3_aci_results$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters.csv'))
    )

    write.csv.exdf(
        c3_aci_results$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits.csv'))
    )

    write.csv.exdf(
        c3_aci_results_new_alpha$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_new_alpha_aci_parameters.csv'))
    )

    write.csv.exdf(
        c3_aci_results_new_alpha$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_new_alpha_aci_fits.csv'))
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

    write.csv.exdf(
        c3_aci_results_deriv$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_deriv_aci_parameters.csv'))
    )

    write.csv.exdf(
        c3_aci_results_deriv$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_deriv_aci_fits.csv'))
    )

    save(
        c3_aci_results, c3_aci_results_new_alpha, plantecophys_results,
        photosynthesis_results, msu_results, c3_aci_results_deriv,
        file = file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData'))
    )
} else {
    load(file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData')))
}

# Get a subset of light levels to use for average plots
rows_to_keep_for_avg <- if (RESTRICT_PPFD_AVG) {
  licor_data[, 'ppfd'] %in% NICE_PPFD
} else {
  rep_len(TRUE, nrow(licor_data))
}

licor_data_subset <- licor_data[rows_to_keep_for_avg, , TRUE]

# Make validation plots if necessary
if (MAKE_VALIDATION_PLOTS) {
    # Plot average curves at each ppfd
    pdf_print(
        xyplot_avg_rc(
            licor_data_subset[, 'A'],
            licor_data_subset[, 'Ci'],
            licor_data_subset[, 'CO2_r_sp'],
            licor_data_subset[, 'ppfd'],
            type = 'b',
            pch = 16,
            auto.key = list(space = 'right'),
            grid = TRUE,
            xlim = CI_LIM,
            ylim = c(min(A_LIM), 40), # Use a different upper limit for A in this plot
            xlab = paste('Ci [', licor_data_subset$units$Ci, ']'),
            ylab = paste('A [', licor_data_subset$units$A, ']'),
            cols = if (RESTRICT_PPFD_AVG) {NICE_PPFD_COLORS} else {multi_curve_colors()}
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_avg_acis.pdf'))
    )

    pdf_print(
        xyplot_avg_rc(
            licor_data_subset[, 'gsw'],
            licor_data_subset[, 'Ci'],
            licor_data_subset[, 'CO2_r_sp'],
            licor_data_subset[, 'ppfd'],
            type = 'b',
            pch = 16,
            auto.key = list(space = 'right'),
            grid = TRUE,
            xlim = CI_LIM,
            ylim = GSW_LIM,
            xlab = paste('Ci [', licor_data_subset$units$Ci, ']'),
            ylab = paste('gsw [', licor_data_subset$units$gsw, ']'),
            cols = if (RESTRICT_PPFD_AVG) {NICE_PPFD_COLORS} else {multi_curve_colors()}
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_avg_gswcis.pdf'))
    )
}

# Make error plots for selected curves
if (MAKE_ERROR_PLOTS) {
    curves_for_error_plots <- list(
        list(id = INDIVIDUAL_CURVE_ID,      j_zoom_lim = c(128, 155), ej_thresh = 0.995),
        list(id = '1000 - wt-1 - mcgrath1', j_zoom_lim = c(110, 140), ej_thresh = 1.0)
    )

    lapply(curves_for_error_plots, function(x) {
        plot_aci_errors(
            c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == x$id, , TRUE],
            c3_aci_results$fits_interpolated[c3_aci_results$fits_interpolated[, 'curve_identifier'] == x$id, , TRUE],
            c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == x$id, , TRUE],
            CI_LIM,
            x$j_zoom_lim,
            x$ej_thresh,
            paste(BASE_NAME, x$id, sep = '_'),
            SAVE_TO_PDF
        )
    })
}

# Remove unreliable parameter estimates. Typically this is done automatically
# when calling `fit_c3_aci`, but we did disabled this feature (by setting
# `remove_unreliable_param = 0`) so we could make the error plots above.
fit_parameters_clean <- c3_aci_results$parameters

for (cid in unique(fit_parameters_clean[, 'curve_identifier'])) {
    if (fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Vcmax_trust'] < 2) {
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Vcmax_at_25']  <- NA
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Vcmax_tl_avg'] <- NA

        if (fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Vcmax_trust'] < 1) {
            c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == cid, 'Wc'] <- NA
            c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == cid, 'Ac'] <- NA
        }
    }

    if (fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'J_trust'] < 2) {
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'J_at_25']  <- NA
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'J_tl_avg'] <- NA

        if (fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'J_trust'] < 1) {
            c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == cid, 'Wj'] <- NA
            c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == cid, 'Aj'] <- NA
        }
    }

    if (fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Tp_trust'] < 2) {
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Tp_at_25']  <- NA
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Tp_tl_avg'] <- NA
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'alpha_g']   <- NA
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'alpha_old'] <- NA
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'alpha_s']   <- NA
        fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'alpha_t']   <- NA

        if (fit_parameters_clean[fit_parameters_clean[, 'curve_identifier'] == cid, 'Tp_trust'] < 1) {
            c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == cid, 'Wp'] <- NA
            c3_aci_results$fits[c3_aci_results$fits[, 'curve_identifier'] == cid, 'Ap'] <- NA
        }
    }
}

# Plot the C3 A-Ci fits (including limiting rates and operating point)
plot_c3_aci_fits(
    c3_aci_results$fits,
    'A',
    'A_fit',
    BASE_NAME,
    CI_LIM,
    A_LIM,
    plot_operating_point = TRUE,
    parameter_exdf = c3_aci_results$parameters
)

plot_c3_aci_fits(
    c3_aci_results$fits,
    'A',
    'A_fit',
    paste0(BASE_NAME, '_zoom'),
    c(0, 100),
    c(-5, 10),
    plot_operating_point = TRUE,
    parameter_exdf = c3_aci_results$parameters
)

# Load parameter values fit using PCE calculator
spreadsheet_parameters <-
    read.csv(file = file.path('data', 'pce_spreadsheet_fits', 'spreadsheet_fit_parameters.csv'))

spreadsheet_parameters$RMSE <- NA
spreadsheet_parameters$AIC  <- NA

spreadsheet_parameters$alpha_old_lower <- NA
spreadsheet_parameters$alpha_old_upper <- NA

spreadsheet_parameters$J_at_25_lower <- NA
spreadsheet_parameters$J_at_25_upper <- NA
spreadsheet_parameters$RL_at_25_lower <- NA
spreadsheet_parameters$RL_at_25_upper <- NA
spreadsheet_parameters$Tp_at_25_lower <- NA
spreadsheet_parameters$Tp_at_25_upper <- NA
spreadsheet_parameters$Vcmax_at_25_lower <- NA
spreadsheet_parameters$Vcmax_at_25_upper <- NA

spreadsheet_parameters$J_tl_lower <- NA
spreadsheet_parameters$J_tl_upper <- NA
spreadsheet_parameters$RL_tl_lower <- NA
spreadsheet_parameters$RL_tl_upper <- NA
spreadsheet_parameters$Tp_tl_lower <- NA
spreadsheet_parameters$Tp_tl_upper <- NA
spreadsheet_parameters$Vcmax_tl_lower <- NA
spreadsheet_parameters$Vcmax_tl_upper <- NA

# Load a single curve that was fit using PCE calculator
spreadsheet_fit_800 <-
    read.csv.exdf(file = file.path('data', 'pce_spreadsheet_fits', '800 - wt-5 - mcgrath1.csv'))

spreadsheet_fit_800 <- set_variable(
    spreadsheet_fit_800,
    'A_fit',
    spreadsheet_fit_800$units$A,
    'pce_calculator'
)

for (i in seq_len(nrow(spreadsheet_fit_800))) {
    if (spreadsheet_fit_800[i, 'limiting_rate'] == 1) {
        spreadsheet_fit_800[i, 'A_fit'] <- spreadsheet_fit_800[i, 'Ac']
    } else if (spreadsheet_fit_800[i, 'limiting_rate'] == 2) {
        spreadsheet_fit_800[i, 'A_fit'] <- spreadsheet_fit_800[i, 'Aj']
    } else if (spreadsheet_fit_800[i, 'limiting_rate'] == 3) {
        spreadsheet_fit_800[i, 'A_fit'] <- spreadsheet_fit_800[i, 'Ap']
    }
}

# Save spreadsheet results since we have added some new columns
write.csv(
    spreadsheet_parameters,
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_spreadsheet_parameters.csv')),
    row.names = FALSE
)

write.csv.exdf(
    spreadsheet_fit_800,
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_spreadsheet_fit_800.csv'))
)

# Compare results from each package for an individual curve
plot_c3_aci_fit_comparison(
    INDIVIDUAL_CURVE_ID,
    c3_aci_results$fits,
    plantecophys_results$fits,
    photosynthesis_results$fits,
    msu_results$fits,
    CI_LIM,
    c(min(A_LIM), 40),
    BASE_NAME,
    spreadsheet_fit = spreadsheet_fit_800$main_data
)

# Get average parameter values at each PPFD
avg_param_by_ppfd <- do.call(rbind, by(
    fit_parameters_clean,
    fit_parameters_clean[, 'ppfd'],
    function(x) {
        res <- data.frame(
            ppfd = x[1, 'ppfd']
        )

        for (pn in c('alpha_old', 'J_at_25', 'RL_at_25', 'Tp_at_25', 'Vcmax_at_25')) {
            res <- cbind(res, param_lim_mean(x, pn))
        }

        res
    }
))

# Save clean and average parameters
write.csv.exdf(
    fit_parameters_clean,
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_clean.csv'))
)

write.csv(
    avg_param_by_ppfd,
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_avg.csv')),
    row.names = FALSE
)

avg_param_by_ppfd_for_fitting <- avg_param_by_ppfd[avg_param_by_ppfd[, 'ppfd'] <= 500, ]

# Reorder parameters by incident PPFD
fit_parameters_clean <- fit_parameters_clean[order(fit_parameters_clean[, 'ppfd']), , TRUE]

# Make plots of each parameter vs. Qin
if (MAKE_PARAMETER_PLOTS) {
    plot_parameters(
        fit_parameters_clean$main_data,
        QIN_LIM,
        J_LIM,
        'p',
        paste0(BASE_NAME, '_param_raw'),
        SAVE_TO_PDF,
        RL_LIM,
        TP_LIM,
        VCMAX_LIM
    )

    plot_parameters(
        avg_param_by_ppfd,
        QIN_LIM,
        J_LIM,
        'b',
        paste0(BASE_NAME, '_param_avg'),
        SAVE_TO_PDF,
        RL_LIM,
        TP_LIM,
        VCMAX_LIM
    )

    plot_parameters_avg(
        fit_parameters_clean$main_data,
        QIN_LIM,
        J_LIM,
        'b',
        paste0(BASE_NAME, '_param_avg_eb'),
        SAVE_TO_PDF,
        RL_LIM,
        TP_LIM,
        VCMAX_LIM
    )
}

# Plot the C3 A-Ci fits from the new TPU model (including limiting rates and operating point)
plot_c3_aci_fits(
    c3_aci_results_new_alpha$fits,
    'A',
    'A_fit',
    paste0(BASE_NAME, '_new_alpha'),
    CI_LIM,
    A_LIM,
    plot_operating_point = TRUE,
    parameter_exdf = c3_aci_results_new_alpha$parameters
)


# Plot the plantecophys C3 A-Ci fits (including limiting rates)
plot_c3_aci_fits(
    plantecophys_results$fits,
    'Ameas',
    'Amodel',
    paste0(BASE_NAME, '_plantecophys'),
    CI_LIM,
    A_LIM,
    plot_operating_point = FALSE
)

# Plot the photosynthesis C3 A-Ci fits (including limiting rates)
plot_c3_aci_fits(
    photosynthesis_results$fits,
    'A_net',
    'A_model',
    paste0(BASE_NAME, '_photosynthesis'),
    CI_LIM,
    c(-30, max(A_LIM)),
    plot_operating_point = FALSE
)

# Plot the msuRACiFit C3 A-Ci fits (including limiting rates)
plot_c3_aci_fits(
    msu_results$fits,
    'A',
    'A_model',
    paste0(BASE_NAME, '_msu'),
    CI_LIM,
    A_LIM,
    plot_operating_point = FALSE
)

###                                 ###
### COMPARE RESULTS ACROSS PACKAGES ###
###                                 ###

# Extract some of the results from the different fitting methods
fit_param <- extract_fit_param(
    fit_parameters_clean,
    c3_aci_results_new_alpha$parameters,
    plantecophys_results$parameters,
    photosynthesis_results$parameters,
    msu_results$parameters,
    spreadsheet_parameters,
    include_true_values = FALSE,
    include_ppfd = TRUE
)

# Make sure curve_identifier is treated as a factor and ordered by PPFD
tmp <- unique(fit_param[, c('ppfd', 'curve_identifier')])

tmp <- tmp[order(tmp$ppfd), ]

fit_param$curve_identifier <-
    factor(fit_param$curve_identifier, levels = tmp$curve_identifier)

# Compare results from each package
plot_parameter_comparisons(
    fit_param,
    QIN_LIM,
    RL_LIM,
    J_LIM,
    VCMAX_LIM,
    TP_LIM,
    RMSE_LIM,
    AIC_LIM,
    ALPHA_LIM,
    seq(-1.625, 2.625, by = 0.25), # RMSE_breaks
    BASE_NAME,
    SAVE_TO_PDF
)

# Compare parameters from each package for an individual curve
plot_c3_aci_param_comparison(
    INDIVIDUAL_CURVE_ID,
    fit_parameters_clean,
    c3_aci_results_new_alpha$parameters,
    plantecophys_results$parameters,
    photosynthesis_results$parameters,
    msu_results$parameters,
    spreadsheet_parameters,
    BASE_NAME
)

# Print the number of fit failures when using nlminb
print('Number of fit failures when using nlminb:')
print(sum(c3_aci_results_deriv$parameters[, 'convergence']))

# Compare RMSE when using nlminb
print('Mean RMSE of PhotoGEA (alpha_old, DEoptim) fits:')
print(mean(c3_aci_results$parameters[, 'RMSE']))

print('Mean RMSE of PhotoGEA (alpha_old, nlminb) fits:')
print(mean(c3_aci_results_deriv$parameters[c3_aci_results_deriv$parameters[, 'convergence'] == 0, 'RMSE']))

# Print info about inadmissible fits
print('Number of inadmissible fits returned by photosynthesis')
print(sum(!photosynthesis_results$parameters[, 'admissible']))

print('Number of inadmissible fits returned by plantecophys')
print(sum(!plantecophys_results$parameters[, 'admissible']))
