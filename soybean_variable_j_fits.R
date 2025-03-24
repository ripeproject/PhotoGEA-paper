###               ###
### PRELIMINARIES ###
###               ###

# Load libraries
library(PhotoGEA)
library(lattice)
library(RColorBrewer)
library(latticeExtra)

# Clear workspace
rm(list=ls())

# Load some helping functions and data
source(file.path('utilities', 'output_tools.R'))
source(file.path('utilities', 'defaults.R'))

source(file.path('utilities', 'validation_plots.R'))
source(file.path('utilities', 'plot_aci_errors.R'))
source(file.path('utilities', 'plot_variable_j_fits_one_curve.R'))

load(file.path(OUTPUT_DIR, RDATA_DIR, 'rubisco_avg_arrhenius_parameters.RData')) # load avg_arrhenius_parameters; created by `soybean_rubisco.R`

# Choose some settings
SAVE_TO_PDF <- TRUE
MAKE_NEW_CALCULATIONS <- TRUE

USE_SOYBEAN_ARRHENIUS <- TRUE

MAKE_VALIDATION_PLOTS  <- TRUE
MAKE_ERROR_PLOTS       <- TRUE

REMOVE_SPECIFIC_POINTS <- TRUE

BASE_NAME <- 'soybean_variable_j'

INDIVIDUAL_CURVE <- '2022 - ripe2 - 4'

A_LIM       <- c(-10, 55)
CI_LIM      <- c(0, 1500)
CC_LIM      <- c(0, 600)
GMC_LIM     <- c(-0.2, 0.8)
GMC_AVG_LIM <- c(-0.4, 0.3)
J_LIM       <- c(0, 450)

# Specialized for `2022 - ripe2 - 4`
#GMC_LIM <- GMC_AVG_LIM
#A_LIM   <- c(0, 45)
#CC_LIM  <- c(0, 350)
#J_LIM   <- c(0, 350)

# Specialized for `2022 − ripe2 − 1`
#GMC_LIM <- GMC_AVG_LIM
#CC_LIM  <- c(0, 400)

if (MAKE_NEW_CALCULATIONS) {
    ###                                  ###
    ### LOAD, VALIDATE, AND PROCESS DATA ###
    ###                                  ###

    # Load the Licor data, which has already been consolidated into a single
    # CSV file
    licor_data <- read.csv.exdf(file.path('data', 'soybean', 'ld11_aci.csv'))

    # Get the time as a POSIX vector
    timestamps <- as.POSIXlt(licor_data[, 'time'])

    # Add columns for year and day of year
    licor_data[, 'year'] <- format(timestamps, '%Y')
    licor_data[, 'doy'] <- format(timestamps, '%j')

    # Create a new identifier column
    licor_data[ , 'curve_identifier'] <- paste(
        licor_data[, 'year'],
        licor_data[, 'instrument'],
        licor_data[, 'plot'],
        sep = ' - '
    )

    # Make sure the data meets basic requirements
    check_response_curve_data(licor_data, 'curve_identifier', 16, 'CO2_r_sp')

    # Remove points with duplicate `CO2_r_sp` and order by Ci
    licor_data <- organize_response_curve_data(
        licor_data,
        'curve_identifier',
        c(1, 9),
        'Ci'
    )

    # Make sure curve_identifier is treated as a factor and ordered by year
    tmp <- unique(licor_data[, c('year', 'curve_identifier')])

    tmp <- tmp[order(tmp$year), ]

    licor_data[, 'curve_identifier'] <-
        factor(licor_data[, 'curve_identifier'], levels = tmp$curve_identifier)

    # Make some plots before removing any points
    pdf_print(
        xyplot(
            A ~ seq_num | curve_identifier,
            group = Stable,
            data = licor_data$main_data,
            type = 'p',
            pch = 16,
            auto = TRUE,
            grid = TRUE
        ),
        width = 20,
        height = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_seqnum_a.pdf'))
    )

    pdf_print(
        xyplot(
            Ci ~ seq_num | curve_identifier,
            group = Stable,
            data = licor_data$main_data,
            type = 'p',
            pch = 16,
            auto = TRUE,
            grid = TRUE
        ),
        width = 20,
        height = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_seqnum_ci.pdf'))
    )

    pdf_print(
        xyplot(
            gsw ~ seq_num | curve_identifier,
            group = Stable,
            data = licor_data$main_data,
            type = 'p',
            pch = 16,
            auto = TRUE,
            grid = TRUE
        ),
        width = 20,
        height = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_seqnum_gsw.pdf'))
    )

    pdf_print(
        xyplot(
            PhiPS2 ~ seq_num | curve_identifier,
            group = Stable,
            data = licor_data$main_data,
            type = 'p',
            pch = 16,
            auto = TRUE,
            grid = TRUE
        ),
        width = 20,
        height = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_seqnum_phips2.pdf'))
    )

    # Remove points if desired
    if (REMOVE_SPECIFIC_POINTS) {
        licor_data <- remove_points(
            licor_data,
            list(curve_identifier = '2021 - ripe1 - 5', seq_num = c(11, 12)),
            list(curve_identifier = '2021 - ripe2 - 1', seq_num = c(11, 12))
        )
    }

    ###                                                 ###
    ### FIT CURVES USING PHOTOGEA's VARIABLE J FUNCTION ###
    ###                                                 ###

    # Calculate total pressure in the Licor chamber
    licor_data <- calculate_total_pressure(licor_data)

    # Calculate temperature-dependent values of C3 photosynthetic parameters
    arrhenius_param <- TEMPERATURE_PARAM

    if (USE_SOYBEAN_ARRHENIUS) {
        # Override the default tobacco Gamma_star parameters with soybean values
        arrhenius_param$Gamma_star$c  <- avg_arrhenius_parameters[avg_arrhenius_parameters$parameter == 'Gamma_star', 'c']
        arrhenius_param$Gamma_star$Ea <- avg_arrhenius_parameters[avg_arrhenius_parameters$parameter == 'Gamma_star', 'Ea']
        arrhenius_param$Kc$c  <- avg_arrhenius_parameters[avg_arrhenius_parameters$parameter == 'Kc', 'c']
        arrhenius_param$Kc$Ea <- avg_arrhenius_parameters[avg_arrhenius_parameters$parameter == 'Kc', 'Ea']
        arrhenius_param$Ko$c  <- avg_arrhenius_parameters[avg_arrhenius_parameters$parameter == 'Ko', 'c']
        arrhenius_param$Ko$Ea <- avg_arrhenius_parameters[avg_arrhenius_parameters$parameter == 'Ko', 'Ea']
    }

    licor_data <- calculate_temperature_response(licor_data, arrhenius_param)

    # Set a seed number because the deoptiom algorithm uses randomness
    set.seed(123)

    # Fit the curves
    c3_variable_j_results <- consolidate(by(
        licor_data,                                           # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'],                     # A factor used to split `licor_data` into chunks
        fit_c3_variable_j,                                    # The function to apply to each chunk of `licor_data`
        optim_fun = optimizer_deoptim(itermax = VJ_ITERMAX),  # Use the slow but reliable optimizer
        fit_options = FIT_OPTIONS,
        lower = list(RL_at_25 = -1),
        Ca_atmospheric = CA_ATMOSPHERIC,
        atp_use = ATP_USE,
        nadph_use = NADPH_USE,
        relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2,
        hard_constraints = HARD_CONSTRAINTS,
        require_positive_gmc = REQUIRE_POSITIVE_GMC,
        gmc_max = GMC_MAX,
        check_j = TRUE
    ))

    ###                                       ###
    ### FIT CURVES WITH A DIFFERENT OPTIMIZER ###
    ###                                       ###

    # Fit the A-Ci curves using stats::nlminb
    c3_variable_j_results_deriv <- consolidate(by(
        licor_data,                          # The `exdf` object containing the curves
        licor_data[, 'curve_identifier'],    # A factor used to split `licor_data` into chunks
        fit_c3_variable_j,                   # The function to apply to each chunk of `licor_data`
        optim_fun = optimizer_nlminb(1e-10), # Use the slow but reliable optimizer
        fit_options = FIT_OPTIONS,
        lower = list(RL_at_25 = -1),
        Ca_atmospheric = CA_ATMOSPHERIC,
        atp_use = ATP_USE,
        nadph_use = NADPH_USE,
        relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2,
        hard_constraints = HARD_CONSTRAINTS,
        require_positive_gmc = REQUIRE_POSITIVE_GMC,
        gmc_max = GMC_MAX,
        check_j = TRUE
    ))

    # Save fit results
    write.csv(
        c3_variable_j_results$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_variable_j_parameters_clean.csv')),
        row.names = FALSE
    )

    write.csv(
        c3_variable_j_results$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_variable_j_fits_clean.csv')),
        row.names = FALSE
    )

    write.csv(
        c3_variable_j_results_deriv$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_variable_j_deriv_parameters_clean.csv')),
        row.names = FALSE
    )

    write.csv(
        c3_variable_j_results_deriv$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_variable_j_deriv_fits_clean.csv')),
        row.names = FALSE
    )

    save(
        licor_data, c3_variable_j_results, c3_variable_j_results_deriv,
        file = file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData'))
    )
} else {
    load(file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData')))
}

# Make validation plots if necessary
if (MAKE_VALIDATION_PLOTS) {
    validation_plots(licor_data, SAVE_TO_PDF, 'Ci', CI_LIM, BASE_NAME)
}

# Plot the fits
pdf_print(
    plot_c3_aci_fit(
        c3_variable_j_results,
        'curve_identifier',
        'Ci',
        xlim = CI_LIM,
        ylim = A_LIM
    ),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_fits_ci.pdf'))
)

pdf_print(
    plot_c3_aci_fit(
        c3_variable_j_results,
        'curve_identifier',
        'Ci',
        xlim = c(-10, 150),
        ylim = c(-10, 15)
    ),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_fits_ci_zoom.pdf'))
)

pdf_print(
    plot_c3_aci_fit(
        c3_variable_j_results,
        'curve_identifier',
        'Cc',
        xlim = c(-10, 150),
        ylim = c(-10, 15)
    ),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_fits_cc_zoom.pdf'))
)

pdf_print(
    plot_c3_aci_fit(
        c3_variable_j_results,
        'curve_identifier',
        'Cc',
        xlim = CC_LIM,
        ylim = A_LIM
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_fits_cc.pdf'))
)


# Plot each Cc-Ci curve
pdf_print(
    xyplot(
        Cc ~ Ci | curve_identifier,
        data = c3_variable_j_results$fits$main_data,
        type = 'b',
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CI_LIM,
        #ylim = c(0, 1100),
        xlab = paste('Ci [', c3_variable_j_results$fits$units$Ci, ']'),
        ylab = paste('Cc [', c3_variable_j_results$fits$units$Cc, ']')
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_separate_cccis.pdf'))
)

# Plot each J-Ci curve
pdf_print(
    xyplot(
        J_F + J_tl ~ Ci | curve_identifier,
        data = c3_variable_j_results$fits$main_data,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CI_LIM,
        ylim = J_LIM,
        xlab = paste('Ci [', c3_variable_j_results$fits$units$Ci, ']'),
        ylab = paste('J [', c3_variable_j_results$fits$units$J_F, ']'),
    ),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_separate_js.pdf'))
)

# Plot each gmc-Ci curve
pdf_print(
    xyplot(
        gmc ~ Ci | curve_identifier,
        data = c3_variable_j_results$fits$main_data,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CI_LIM,
        ylim = GMC_LIM,
        xlab = paste('Ci [', c3_variable_j_results$fits$units$Ci, ']'),
        ylab = paste('gmc [', c3_variable_j_results$fits$units$gmc, ']')
    ),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_separate_gmcs.pdf'))
)

# Plot each gmc-Ci curve on one set of axes
pdf_print(
    xyplot(
        gmc ~ Ci,
        group = curve_identifier,
        data = c3_variable_j_results$fits$main_data,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CI_LIM,
        ylim = GMC_LIM,
        xlab = paste('Ci [', c3_variable_j_results$fits$units$Ci, ']'),
        ylab = paste('gmc [', c3_variable_j_results$fits$units$gmc, ']')
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_combined_gmcs.pdf'))
)

# Plot each gmc-Cc curve on one set of axes
pdf_print(
    xyplot(
        gmc ~ Cc,
        group = curve_identifier,
        data = c3_variable_j_results$fits$main_data,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CC_LIM,
        ylim = GMC_LIM,
        xlab = paste('Cc [', c3_variable_j_results$fits$units$Cc, ']'),
        ylab = paste('gmc [', c3_variable_j_results$fits$units$gmc, ']')
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_combined_gmcccs.pdf'))
)

# Plot each A-Cc curve on one set of axes
pdf_print(
    xyplot(
        A ~ Cc,
        group = curve_identifier,
        data = c3_variable_j_results$fits$main_data,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CC_LIM,
        ylim = A_LIM,
        xlab = paste('Cc [', c3_variable_j_results$fits$units$Cc, ']'),
        ylab = paste('A [', c3_variable_j_results$fits$units$A, ']')
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_combined_accs.pdf'))
)

# Plot some results for just one curve
plot_variable_j_fits_one_curve(
    c3_variable_j_results$fits,
    INDIVIDUAL_CURVE,
    CI_LIM,
    c(0, 0.3),
    c(0, 300),
    BASE_NAME
)

###
### PLOTTING AVERAGE CURVES
###

# Here we want to exlcude some curves
exclude_for_average <- c(
    '2021 - ripe3 - 4', # humidity was not properly controlled so the Ci values for this curve are much different than the others
    '2022 - ripe2 - 5'  # gmc values for this curve are clearly outliers
)

fits_for_average <-
    c3_variable_j_results$fits[!c3_variable_j_results$fits[, 'curve_identifier'] %in% exclude_for_average, , TRUE]

# Plot average A-Ci curves
pdf_print(
    xyplot_avg_rc(
        fits_for_average[, 'A'],
        fits_for_average[, 'Ci'],
        fits_for_average[, 'CO2_r_sp'],
        fits_for_average[, 'species'],
        x_error_bars = TRUE,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CI_LIM,
        ylim = A_LIM,
        xlab = paste('Ci [', fits_for_average$units$Ci, ']'),
        ylab = paste('A [', fits_for_average$units$A, ']'),
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_avg_acis.pdf'))
)

# Plot average gmc curves for each year
pdf_print(
    xyplot_avg_rc(
        fits_for_average[, 'gmc'],
        fits_for_average[, 'Ci'],
        fits_for_average[, 'CO2_r_sp'],
        fits_for_average[, 'species'],
        x_error_bars = TRUE,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        grid = TRUE,
        xlim = CI_LIM,
        ylim = GMC_AVG_LIM,
        xlab = paste('Ci [', fits_for_average$units$Ci, ']'),
        ylab = paste('gmc [', fits_for_average$units$gmc, ']'),
    ),
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_avg_gmcs.pdf'))
)

###
### PLOTTING PARAMETER VALUES
###

# Plot values for each curve
for (pn in c('RL_at_25', 'Vcmax_at_25', 'J_at_25', 'Tp_at_25', 'alpha_old', 'tau')) {
    pn_lower <- paste0(pn, '_lower')
    pn_upper <- paste0(pn, '_upper')

    p_lower <- c3_variable_j_results$parameters[, pn_lower]
    p       <- c3_variable_j_results$parameters[, pn]
    p_upper <- c3_variable_j_results$parameters[, pn_upper]
    cid     <- c3_variable_j_results$parameters[, 'curve_identifier']

    pdf_print(
        xyplot(
            p_lower + p + p_upper ~ cid,
            type = 'p',
            pch = 16,
            scales = list(x = list(rot = 90)),
            ylab = pn
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_compare_', pn, '.pdf'))
    )
}

###
### PLOTTING ERROR VALUES
###

if (MAKE_ERROR_PLOTS) {
    curves_for_error_plots <- list(
        list(id = INDIVIDUAL_CURVE, j_zoom_lim = c(157, 162), ej_thresh = 0.7),
        list(id = '2022 - ripe2 - 5', j_zoom_lim = c(157, 162), ej_thresh = 0.7)
    )

    lapply(curves_for_error_plots, function(x) {
        plot_aci_errors(
            c3_variable_j_results$fits[c3_variable_j_results$fits[, 'curve_identifier'] == x$id, , TRUE],
            c3_variable_j_results$fits_interpolated[c3_variable_j_results$fits_interpolated[, 'curve_identifier'] == x$id, , TRUE],
            c3_variable_j_results$parameters[c3_variable_j_results$parameters[, 'curve_identifier'] == x$id, , TRUE],
            CI_LIM,
            x$j_zoom_lim,
            x$ej_thresh,
            paste(BASE_NAME, x$id, sep = '_'),
            SAVE_TO_PDF
        )
    })
}

# Print the number of fit failures when using nlminb
print('Number of fit failures when using nlminb:')
print(sum(c3_variable_j_results_deriv$parameters[, 'convergence']))

# Compare RMSE when using nlminb
print('Mean RMSE of PhotoGEA (alpha_old, DEoptim) fits:')
print(mean(c3_variable_j_results$parameters[, 'RMSE']))

print('Mean RMSE of PhotoGEA (alpha_old, nlminb) fits:')
print(mean(c3_variable_j_results_deriv$parameters[c3_variable_j_results_deriv$parameters[, 'convergence'] == 0, 'RMSE']))
