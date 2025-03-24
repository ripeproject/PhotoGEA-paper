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

# Load some helping functions
source(file.path('utilities', 'output_tools.R'))
source(file.path('utilities', 'defaults.R'))

source(file.path('utilities', 'validation_plots.R'))

# Choose some settings
SAVE_TO_PDF <- TRUE
MAKE_NEW_CALCULATIONS <- TRUE

REMOVE_BAD_POINTS <- TRUE

MAKE_VALIDATION_PLOTS  <- TRUE

EXCLUDE_CURVES_WITH_OUTLIER_VPMAX <- FALSE

LOW_CI_THRESHOLD <- 60

FIT_OPTIONS_1 <- list(Vcmax_at_25 = 'fit', J_at_25 = 1000,  Vpr = 1000) # Fit Vcmax
FIT_OPTIONS_2 <- list(Vcmax_at_25 = 1000,  J_at_25 = 'fit', Vpr = 1000) # Fit J
FIT_OPTIONS_3 <- list(Vcmax_at_25 = 'fit', J_at_25 = 'fit', Vpr = 1000) # Fit Vcmax and J
FIT_OPTIONS_4 <- list(Vcmax_at_25 = 1000,  J_at_25 = 1000,  Vpr = 1000, gbs = 0, gmc_at_25 = Inf) # Only fit Vpmax

BASE_NAME <- 'c4_aci'

A_LIM  <- c(-5, 75)
CI_LIM <- c(-100, 1500)

CI_LIM_LOW <- c(-10, LOW_CI_THRESHOLD + 5)
A_LIM_LOW  <- c(-5, 50)

if (MAKE_NEW_CALCULATIONS) {
    ###                                  ###
    ### LOAD, VALIDATE, AND PROCESS DATA ###
    ###                                  ###

    # Load the Licor data, which has already been consolidated into a single
    # CSV file
    licor_data <- read.csv.exdf(file.path('data', 'c4', 'c4_aci.csv'))

    # Get the time as a POSIX vector
    timestamps <- as.POSIXlt(licor_data[, 'time'])

    # Add columns for year and day of year
    licor_data[, 'year'] <- format(timestamps, '%Y')
    licor_data[, 'doy'] <- format(timestamps, '%j')

    # Create a new identifier column
    licor_data[ , 'curve_identifier'] <- paste(
        licor_data[, 'species'],
        licor_data[, 'instrument'],
        licor_data[, 'plot'],
        licor_data[, 'year'],
        sep = ' - '
    )

    # Make sure the data meets basic requirements, remove points with duplicated
    # `CO2_r_sp`, and order by `Ci`
    licor_data <- do.call(rbind, lapply(
        list(
            list(year = '2021', npts = 16, to_remove = c(9, 10)),
            list(year = '2022', npts = 14, to_remove = c(8, 9))
        ),
        function(x) {
            check_response_curve_data(
                licor_data[licor_data[, 'year'] == x$year, , TRUE],
                'curve_identifier',
                x$npts,
                'CO2_r_sp'
            )

            organize_response_curve_data(
                licor_data[licor_data[, 'year'] == x$year, , TRUE],
                'curve_identifier',
                x$to_remove,
                'Ci'
            )
        }
    ))

    # Make some plots before removing any points
    pdf_print(
        xyplot(
            A ~ seq_num | curve_identifier,
            group = Stable,
            data = licor_data$main_data,
            type = 'p',
            pch = 16
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
            pch = 16
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
            pch = 16
        ),
        width = 20,
        height = 10,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_seqnum_gsw.pdf'))
    )

    # Some points are not reasonable
    if (REMOVE_BAD_POINTS) {
        licor_data <- remove_points(
            licor_data,
            list(curve_identifier = 'maize - ripe3 - 2 - 2021',    seq_num = c(14, 15)), # A low here
            list(curve_identifier = 'maize - ripe4 - 5 - 2021',    seq_num = 16),        # Ci sharply decreases here
            list(curve_identifier = 'sorghum - ripe1 - 3 - 2021',  seq_num = c(15, 16)), # A decreases here
            list(curve_identifier = 'sorghum - ripe14 - 4 - 2022', seq_num = c(3, 4)),   # Ci sharply decreases here
            list(curve_identifier = 'sorghum - ripe15 - 5 - 2022', seq_num = 16),        # A decreases here
            list(curve_identifier = 'sorghum - ripe15 - 5 - 2022', seq_num = 3)          # Ci sharply decreases here
        )
    }

    # Make sure curve_identifier is treated as a factor and ordered by species / year
    tmp <- unique(licor_data[, c('species', 'year', 'curve_identifier')])

    tmp <- tmp[order(tmp$species, tmp$year), ]

    licor_data[, 'curve_identifier'] <-
        factor(licor_data[, 'curve_identifier'], levels = tmp$curve_identifier)

    ###                                ###
    ### FIT A-CI CURVES USING PHOTOGEA ###
    ###                                ###

    # Calculate the total pressure in the Licor chamber
    licor_data <- calculate_total_pressure(licor_data)

    # Calculate temperature-dependent values of C4 photosynthetic parameters
    licor_data <- calculate_temperature_response(licor_data, c4_temperature_param_vc)

    # Get a truncated version only including points where Ci <= 50 ppm
    licor_data_low_ci <- licor_data[licor_data[, 'Ci'] <= 50, , TRUE]

    # Set a seed number because the deoptiom algorithm uses randomness
    set.seed(123)

    # Fit the C4 A-Ci curves
    c4_aci_results_v1 <- consolidate(by(
      licor_data,
      licor_data[, 'curve_identifier'],
      fit_c4_aci,
      optim_fun = optimizer_deoptim(C4_ACI_ITERMAX),
      fit_options = FIT_OPTIONS_1,
      Ca_atmospheric = CA_ATMOSPHERIC,
      relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
      hard_constraints = HARD_CONSTRAINTS,
      calculate_confidence_intervals = TRUE,
      remove_unreliable_param = 2
    ))

    # Fit the C4 A-Ci curves
    c4_aci_results_v2 <- consolidate(by(
      licor_data,
      licor_data[, 'curve_identifier'],
      fit_c4_aci,
      optim_fun = optimizer_deoptim(C4_ACI_ITERMAX),
      fit_options = FIT_OPTIONS_2,
      Ca_atmospheric = CA_ATMOSPHERIC,
      relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
      hard_constraints = HARD_CONSTRAINTS,
      calculate_confidence_intervals = TRUE,
      remove_unreliable_param = 2
    ))

    # Fit the C4 A-Ci curves
    c4_aci_results_v3 <- consolidate(by(
      licor_data,
      licor_data[, 'curve_identifier'],
      fit_c4_aci,
      optim_fun = optimizer_deoptim(1),
      fit_options = FIT_OPTIONS_3,
      Ca_atmospheric = CA_ATMOSPHERIC,
      relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
      hard_constraints = HARD_CONSTRAINTS,
      calculate_confidence_intervals = TRUE,
      remove_unreliable_param = 2
    ))

    # Fit the C4 A-Ci curves
    c4_aci_results_v4 <- consolidate(by(
      licor_data_low_ci,
      licor_data_low_ci[, 'curve_identifier'],
      fit_c4_aci,
      optim_fun = optimizer_deoptim(C4_ACI_ITERMAX),
      fit_options = FIT_OPTIONS_4,
      relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
      hard_constraints = HARD_CONSTRAINTS,
      calculate_confidence_intervals = TRUE,
      remove_unreliable_param = 2
    ))

    # Fit the C4 A-Ci curves
    c4_aci_results_v5 <- consolidate(by(
      licor_data,
      licor_data[, 'curve_identifier'],
      fit_c4_aci_hyperbola,
      relative_likelihood_threshold = ERROR_THRESHOLD_FACTOR,
      hard_constraints = HARD_CONSTRAINTS,
      calculate_confidence_intervals = TRUE
    ))

    # Save results
    write.csv.exdf(
        c4_aci_results_v1$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_v1.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v1$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_v1.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v2$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_v2.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v2$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_v2.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v3$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_v3.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v3$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_v3.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v4$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_v4.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v4$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_v4.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v5$parameters,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_parameters_v5.csv'))
    )

    write.csv.exdf(
        c4_aci_results_v5$fits,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_v5.csv'))
    )

    save(
        licor_data, licor_data_low_ci,
        c4_aci_results_v1, c4_aci_results_v2, c4_aci_results_v3,
        c4_aci_results_v4, c4_aci_results_v5,
        file = file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData'))
    )
} else {
    load(file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '.RData')))
}

# Make validation plots if necessary
if (MAKE_VALIDATION_PLOTS) {
    validation_plots(licor_data, SAVE_TO_PDF, 'Ci', CI_LIM, BASE_NAME)

    # Plot average curves for each species
    pdf_print(
        xyplot_avg_rc(
            licor_data[, 'A'],
            licor_data[, 'Ci'],
            licor_data[, 'CO2_r_sp'],
            licor_data[, 'species'],
            type = 'b',
            pch = 16,
            auto.key = list(space = 'right'),
            grid = TRUE,
            xlim = CI_LIM,
            ylim = A_LIM,
            xlab = paste('Ci [', licor_data$units$Ci, ']'),
            ylab = paste('A [', licor_data$units$A, ']')
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_avg_acis.pdf'))
    )
}

fit_parameters_clean <- c4_aci_results_v1$parameters

# Plot the C4 A-Ci fits (including limiting rates)
pdf_print(
    plot_c4_aci_fit(c4_aci_results_v1, 'curve_identifier', 'Ci', xlim = CI_LIM, ylim = A_LIM),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fits_v1.pdf'))
)

pdf_print(
    plot_c4_aci_fit(c4_aci_results_v2, 'curve_identifier', 'Ci', xlim = CI_LIM, ylim = A_LIM),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fits_v2.pdf'))
)

pdf_print(
    plot_c4_aci_fit(c4_aci_results_v3, 'curve_identifier', 'Ci', xlim = CI_LIM, ylim = A_LIM),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fits_v3.pdf'))
)

pdf_print(
    plot_c4_aci_fit(c4_aci_results_v4, 'curve_identifier', 'Ci', xlim = CI_LIM_LOW, ylim = A_LIM_LOW),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fits_v4.pdf'))
)

pdf_print(
    plot_c4_aci_hyperbola_fit(c4_aci_results_v5, 'curve_identifier', xlim = CI_LIM, ylim = A_LIM),
    width = 15,
    height = 15,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fits_v5.pdf'))
)

# Plot the C4 A-Ci fit residuals
pdf_print(
    xyplot(
      A_residuals ~ Ci | curve_identifier,
      data = c4_aci_results_v1$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right', lines = TRUE, points = TRUE),
      grid = TRUE,
      xlim = CI_LIM,
      #ylim = A_LIM,
      xlab = paste('Intercellular CO2 concentration [', c4_aci_results_v1$fits$units$Ci, ']'),
      ylab = paste('Net CO2 assimilation rate residual [', c4_aci_results_v1$fits$units$A, ']')
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fit_residuals_v1.pdf'))
)

pdf_print(
    xyplot(
      A_residuals ~ Ci | curve_identifier,
      data = c4_aci_results_v2$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right', lines = TRUE, points = TRUE),
      grid = TRUE,
      xlim = CI_LIM,
      #ylim = A_LIM,
      xlab = paste('Intercellular CO2 concentration [', c4_aci_results_v2$fits$units$Ci, ']'),
      ylab = paste('Net CO2 assimilation rate residual [', c4_aci_results_v2$fits$units$A, ']')
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fit_residuals_v2.pdf'))
)

pdf_print(
    xyplot(
      A_residuals ~ Ci | curve_identifier,
      data = c4_aci_results_v4$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right', lines = TRUE, points = TRUE),
      grid = TRUE,
      xlim = CI_LIM,
      #ylim = A_LIM,
      xlab = paste('Intercellular CO2 concentration [', c4_aci_results_v4$fits$units$Ci, ']'),
      ylab = paste('Net CO2 assimilation rate residual [', c4_aci_results_v4$fits$units$A, ']')
    ),
    width = 20,
    height = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_aci_fit_residuals_v4.pdf'))
)


# Make parameter plots
for (pn in c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'alpha_psii', 'gbs', 'Rm_frac', 'Vpr')) {
    pn_lower <- paste0(pn, '_lower')
    pn_upper <- paste0(pn, '_upper')

    p_lower <- fit_parameters_clean[, pn_lower]
    p       <- fit_parameters_clean[, pn]
    p_upper <- fit_parameters_clean[, pn_upper]
    cid     <- fit_parameters_clean[, 'curve_identifier']

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

# Compare key parameters across methods (excluding Vpmax_at_25 outliers)
param_colnames <- c(
    'curve_identifier', 'species', 'plot', 'instrument', 'year', 'RMSE',
    'RL_tl_avg', 'RL_at_25_lower', 'RL_at_25', 'RL_at_25_upper',
    'Vcmax_tl_avg', 'Vcmax_at_25_lower', 'Vcmax_at_25', 'Vcmax_at_25_upper',
    'Vpmax_tl_avg', 'Vpmax_at_25_lower', 'Vpmax_at_25', 'Vpmax_at_25_upper',
    'J_tl_avg', 'J_at_25_lower', 'J_at_25', 'J_at_25_upper'
)

all_param_df <- rbind(
    within(c4_aci_results_v1$parameters[, param_colnames], {method = 'Vcmax'}),
    within(c4_aci_results_v2$parameters[, param_colnames], {method = 'J'}),
    within(c4_aci_results_v3$parameters[, param_colnames], {method = 'Vcmax & J'}),
    within(c4_aci_results_v4$parameters[, param_colnames], {method = 'Vpmax'})
)

all_param_df_clean <- if (EXCLUDE_CURVES_WITH_OUTLIER_VPMAX) {
    exclude_outliers(
        all_param_df,
        'Vpmax_at_25',
        list(all_param_df$species, all_param_df$method)
    )
} else {
    all_param_df
}

all_param_df_clean$species_method <-
    paste(all_param_df_clean$species, all_param_df_clean$method)

all_param_df_clean$species_method <- factor(
    all_param_df_clean$species_method,
    levels = c(
        'maize Vcmax',
        'maize J',
        'maize Vcmax & J',
        'maize Vpmax',
        'sorghum Vcmax',
        'sorghum J',
        'sorghum Vcmax & J',
        'sorghum Vpmax'
    )
)

# Get average parameter values and save them to a CSV file
avg_param <- basic_stats(exdf(all_param_df_clean), 'species_method')

write.csv(
    avg_param$main_data,
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_aci_fits_avg_param.csv'))
)

# Plot parameters as box-whisker plots
pdf_print(
    bwplot_wrapper(
        all_param_df_clean$J_tl_avg,
        all_param_df_clean$species_method,
        ylim = c(0, 500),
        ylab = 'J at leaf temperature [ micromol m^(-2) s^(-1) ]'
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_boxplot_J.pdf'))
)

pdf_print(
    bwplot_wrapper(
        all_param_df_clean$Vcmax_tl_avg,
        all_param_df_clean$species_method,
        ylim = c(0, 100),
        ylab = 'Vcmax at leaf temperature [ micromol m^(-2) s^(-1) ]'
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_boxplot_Vcmax.pdf'))
)

pdf_print(
    bwplot_wrapper(
        all_param_df_clean$Vpmax_tl_avg,
        all_param_df_clean$species_method,
        ylim = c(0, 450),
        ylab = 'Vpmax at leaf temperature [ micromol m^(-2) s^(-1) ]'
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_boxplot_Vpmax.pdf'))
)

pdf_print(
    bwplot_wrapper(
        all_param_df_clean$RMSE,
        all_param_df_clean$species_method,
        ylim = c(0, 12),
        ylab = 'RMSE [ micromol m^(-2) s^(-1) ]'
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_boxplot_RMSE.pdf'))
)


# Get some correlation coefficients
correlation_info <- function(y, x, species, yname, xname) {
    maize_fit   <- lm(y ~ x, subset = species == 'maize' & x > 0)
    sorghum_fit <- lm(y ~ x, subset = species == 'sorghum' & x > 0)

    maize_f   <- summary(maize_fit)$fstatistic
    sorghum_f <- summary(sorghum_fit)$fstatistic

    data.frame(
        species = c('maize', 'sorghum'),
        yname = yname,
        xname = xname,
        intercept     = c(summary(maize_fit)$coefficients[1, 1], summary(sorghum_fit)$coefficients[1, 1]),
        intercept_err = c(summary(maize_fit)$coefficients[1, 2], summary(sorghum_fit)$coefficients[1, 2]),
        slope         = c(summary(maize_fit)$coefficients[2, 1], summary(sorghum_fit)$coefficients[2, 1]),
        slope_err     = c(summary(maize_fit)$coefficients[2, 2], summary(sorghum_fit)$coefficients[2, 2]),
        r_squared     = c(summary(maize_fit)[['r.squared']],  summary(sorghum_fit)[['r.squared']]),
        p_value = c(
            pf(maize_f[1], maize_f[2], maize_f[3], lower.tail = FALSE),
            pf(sorghum_f[1], sorghum_f[2], sorghum_f[3], lower.tail = FALSE)
        )
    )
}

c4_aci_results_v4$parameters[c4_aci_results_v4$parameters[, 'curve_identifier'] == 'sorghum - ripe14 - 4 - 2022', 'Vpmax_at_25'] <- NA # remove an outlier

correlations <- rbind(
    correlation_info(
        c4_aci_results_v1$parameters[, 'Vcmax_at_25'],
        c4_aci_results_v5$parameters[, 'Vmax'],
        c4_aci_results_v1$parameters[, 'species'],
        'Vcmax_at_25',
        'Vmax'
    ),
    correlation_info(
        c4_aci_results_v2$parameters[, 'J_at_25'],
        c4_aci_results_v5$parameters[, 'Vmax'],
        c4_aci_results_v2$parameters[, 'species'],
        'J_at_25',
        'Vmax'
    ),
    correlation_info(
        c4_aci_results_v1$parameters[, 'Vpmax_at_25'],
        c4_aci_results_v4$parameters[, 'Vpmax_at_25'],
        c4_aci_results_v1$parameters[, 'species'],
        'Vpmax (assuming Rubisco)',
        'Vpmax (low Ci)'
    ),
    correlation_info(
        c4_aci_results_v2$parameters[, 'Vpmax_at_25'],
        c4_aci_results_v4$parameters[, 'Vpmax_at_25'],
        c4_aci_results_v2$parameters[, 'species'],
        'Vpmax (assuming light)',
        'Vpmax (low Ci)'
    )
)

write.csv(
    correlations[correlations$species == 'maize', ],
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_correlations_maize.csv')),
    row.names = FALSE
)

write.csv(
    correlations[correlations$species == 'sorghum', ],
    file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASE_NAME, '_correlations_sorghum.csv')),
    row.names = FALSE
)

# Compare hyperbola to mechanistic fits

pdf_print(
    xyplot(
        c4_aci_results_v1$parameters[, 'Vcmax_at_25'] ~ c4_aci_results_v5$parameters[, 'Vmax'] | c4_aci_results_v1$parameters[, 'species'],
        type = 'p',
        pch = 16,
        auto = TRUE,
        xlim = c(40, 65),
        ylim = c(20, 45),
        xlab = 'Vmax (hyperbola) [ micromol m^(-2) s^(-1) ]',
        ylab = 'Vcmax at 25 C (mechanistic) [ micromol m^(-2) s^(-1) ]',
        panel = function(...) {
            panel.xyplot(...)

            corr_x <- c(40, 65)
            corr_y1 <- correlations$intercept[1] + correlations$slope[1] * corr_x # maize
            corr_y2 <- correlations$intercept[2] + correlations$slope[2] * corr_x # sorghum

            panel.lines(corr_y1 ~ corr_x, col = 'black')
            panel.lines(corr_y2 ~ corr_x, col = 'red')
        }
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_compare_fits_vcmax_25.pdf'))
)

pdf_print(
    xyplot(
        c4_aci_results_v2$parameters[, 'J_at_25'] ~ c4_aci_results_v5$parameters[, 'Vmax'] | c4_aci_results_v1$parameters[, 'species'],
        type = 'p',
        pch = 16,
        auto = TRUE,
        xlim = c(40, 65),
        ylim = c(160, 280),
        xlab = 'Vmax (hyperbola) [ micromol m^(-2) s^(-1) ]',
        ylab = 'J at 25 C (mechanistic) [ micromol m^(-2) s^(-1) ]',
        panel = function(...) {
            panel.xyplot(...)

            corr_x <- c(40, 65)
            corr_y1 <- correlations$intercept[3] + correlations$slope[3] * corr_x # maize
            corr_y2 <- correlations$intercept[4] + correlations$slope[4] * corr_x # sorghum

            panel.lines(corr_y1 ~ corr_x, col = 'black')
            panel.lines(corr_y2 ~ corr_x, col = 'red')
        }
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_compare_fits_j_25.pdf'))
)

pdf_print(
    xyplot(
        c4_aci_results_v1$parameters[, 'Vpmax_at_25'] + c4_aci_results_v2$parameters[, 'Vpmax_at_25'] ~ c4_aci_results_v4$parameters[, 'Vpmax_at_25'] | c4_aci_results_v1$parameters[, 'species'],
        type = 'p',
        pch = 16,
        auto.key = list(space = 'top'),
        xlim = c(25, 125),
        ylim = c(50, 300),
        xlab = 'Vpmax at 25 (low Ci) [ micromol m^(-2) s^(-1) ]',
        ylab = 'Vpmax at 25 (whole curve) [ micromol m^(-2) s^(-1) ]',
        panel = function(...) {
            panel.xyplot(...)

            corr_x <- c(25, 125)
            corr_y1 <- correlations$intercept[5] + correlations$slope[5] * corr_x # maize
            corr_y2 <- correlations$intercept[6] + correlations$slope[6] * corr_x # sorghum
            corr_y3 <- correlations$intercept[7] + correlations$slope[7] * corr_x # maize
            corr_y4 <- correlations$intercept[8] + correlations$slope[8] * corr_x # sorghum

            panel.lines(corr_y1 ~ corr_x, col = 'black')
            panel.lines(corr_y2 ~ corr_x, col = 'red')
            panel.lines(corr_y3 ~ corr_x, col = 'black', lty = 2)
            panel.lines(corr_y4 ~ corr_x, col = 'red', lty = 2)
        }
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_compare_fits_vpmax_25.pdf'))
)

pdf_print(
    xyplot(
        c4_aci_results_v1$parameters[, 'Vcmax_tl_avg'] ~ c4_aci_results_v5$parameters[, 'Vmax'] | c4_aci_results_v1$parameters[, 'species'],
        type = 'p',
        pch = 16,
        auto = TRUE,
        xlim = c(40, 65),
        ylim = c(45, 70),
        xlab = 'Vmax (hyperbola) [ micromol m^(-2) s^(-1) ]',
        ylab = 'Vcmax at Tleaf (mechanistic) [ micromol m^(-2) s^(-1) ]'
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_compare_fits_vcmax_tleaf.pdf'))
)

pdf_print(
    xyplot(
        c4_aci_results_v2$parameters[, 'J_tl_avg'] ~ c4_aci_results_v5$parameters[, 'Vmax'] | c4_aci_results_v1$parameters[, 'species'],
        type = 'p',
        pch = 16,
        auto = TRUE,
        xlim = c(40, 65),
        ylim = c(250, 450),
        xlab = 'Vmax (hyperbola) [ micromol m^(-2) s^(-1) ]',
        ylab = 'J at Tleaf (mechanistic) [ micromol m^(-2) s^(-1) ]'
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_compare_fits_j_tleaf.pdf'))
)

pdf_print(
    xyplot(
        c4_aci_results_v1$parameters[, 'Vpmax_tl_avg'] + c4_aci_results_v2$parameters[, 'Vpmax_tl_avg'] ~ c4_aci_results_v4$parameters[, 'Vpmax_tl_avg'] | c4_aci_results_v1$parameters[, 'species'],
        type = 'p',
        pch = 16,
        auto.key = list(space = 'top'),
        xlim = c(50, 350),
        ylim = c(50, 350),
        xlab = 'Vpmax at Tleaf (low Ci) [ micromol m^(-2) s^(-1) ]',
        ylab = 'Vpmax at Tleaf (whole curve) [ micromol m^(-2) s^(-1) ]'
    ),
    width = 10,
    save_to_pdf = SAVE_TO_PDF,
    file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASE_NAME, '_compare_fits_vpmax_tleaf.pdf'))
)
