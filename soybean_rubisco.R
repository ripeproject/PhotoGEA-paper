# Here we use Glycine max values of Rubisco Kc, Ko, and specificity, reported in
# Orr et al. (2016) to estimate Arrhenius temperature response parameters

# Soybean Rubisco catalytic values at three temperatures are available in the
# supplemental information from Orr et al. 2016. “Surveying Rubisco Diversity
# and Temperature Response to Improve Crop Photosynthetic Efficiency.” Plant
# Physiology 172 (2): 707–17. https://doi.org/10.1104/pp.16.00750.

###
### PRELIMINARIES
###

# Load required packages
library(PhotoGEA)
library(lattice)

# Clear workspace
rm(list=ls())

# Load some helping functions
source(file.path('utilities', 'output_tools.R'))

# Specify some settings
BASENAME <- 'rubisco'

DEFAULT_WIDTH <- 10
SAVE_TO_PDF <- TRUE

# Specify Rubisco catalytic properties
rubisco_catalytic_values <- exdf(
    data.frame(
        temperature    = c(20,    25,   30),
        Kc_aq          = c(7.5,   13.2, 19.6),
        Ko_aq          = c(518,   660,  549),
        specificity_aq = c(116.6, 99.6, 90.6),
        oxygen         = 21,
        species        = 'Glycine max'
    ),
    units = data.frame(
        temperature = 'degrees C',
        Kc_aq = 'microM',
        Ko_aq = 'microM',
        specificity_aq = 'M / M',
        oxygen = 'percent',
        species = ''
    )
)

# Calculate Gamma_star and Henry's constants
rubisco_catalytic_values <- calculate_gamma_star(
    rubisco_catalytic_values,
    specificity_at_tleaf_column_name = 'specificity_aq',
    tleaf_column_name = 'temperature'
)

# Convert Kc and Ko to a gas basis:
#
# * We have Kc_aq and Ko_aq in microM = micromol / L, and Henry's constants in
#   micromol / kg / Pa
#
# * The density of water across 20 - 30 degrees C is 0.997 kg / L (to around 1%)
#
# So the partial pressure of CO2 (in Pa) corresponding to Kc_aq can be found
# using Kc_pp = Kc_aq * density / Henry_constant_CO2. Likewise for O2 and Ko.
#
# Then we can convert to a dimensionless mole fraction assuming a standard
# atmosphere (101325 Pa): Kc = Kc_pp / 101325. Likewise for O2 and Ko.
#
# Finally we can convert to micromol / mol or mmol / mol using a standard SI
# prefix conversion.
density_of_water     <- 0.997  # kg / L
atmospheric_pressure <- 101325 # Pa

rubisco_catalytic_values <- set_variable(
    rubisco_catalytic_values,
    'Kc',
    'micromol mol^(-1)',
    value = rubisco_catalytic_values[, 'Kc_aq'] * density_of_water / rubisco_catalytic_values[, 'H_CO2'] / atmospheric_pressure * 1e6
)

rubisco_catalytic_values <- set_variable(
    rubisco_catalytic_values,
    'Ko',
    'mmol mol^(-1)',
    value = rubisco_catalytic_values[, 'Ko_aq'] * density_of_water / rubisco_catalytic_values[, 'H_O2'] / atmospheric_pressure * 1e3
)

# Add a new column for the inverse thermal energy
rubisco_catalytic_values <- set_variable(
    rubisco_catalytic_values,
    name = 'arrhenius_index',
    units = 'mol kJ^(-1)',
    value = 1 / (PhotoGEA:::ideal_gas_constant * (rubisco_catalytic_values[, 'temperature'] - PhotoGEA:::absolute_zero))
)

# Save the catalytic information
write.csv.exdf(rubisco_catalytic_values, file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASENAME, '_catalytic_values.csv')))

###
### VALIDATION
###

plot_rubisco_t_response <- function(yname, xname, groupname, exdf_obj) {
    pdf_print(
        xyplot(
            exdf_obj[, yname] ~ exdf_obj[, xname],
            group = exdf_obj[, groupname],
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            xlab = paste(xname, '[', exdf_obj$units[[xname]], ']'),
            ylab = paste(yname, '[', exdf_obj$units[[yname]], ']')
        ),
        width = DEFAULT_WIDTH,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(BASENAME, '_validation_', yname, '_', xname, '.pdf'))
    )
}

to_plot <- list(
    list(yname = 'Kc_aq',          xname = 'temperature'),
    list(yname = 'Kc',             xname = 'temperature'),
    list(yname = 'Kc',             xname = 'arrhenius_index'),
    list(yname = 'Ko_aq',          xname = 'temperature'),
    list(yname = 'Ko',             xname = 'temperature'),
    list(yname = 'Ko',             xname = 'arrhenius_index'),
    list(yname = 'specificity_aq', xname = 'temperature'),
    list(yname = 'Gamma_star',     xname = 'temperature'),
    list(yname = 'Gamma_star',     xname = 'arrhenius_index')
)

lapply(to_plot, function(x) {
    plot_rubisco_t_response(x$yname, x$xname, 'species', rubisco_catalytic_values)
})

###
### PROCESSING
###

# Define helping functions for estimating Arrhenius parameters
fit_arrhenius <- function(yvals, xvals) {
    lfit <- lm(log(yvals) ~ xvals)

    list(
        c = as.numeric(lfit$coefficients[1]),
        Ea = -as.numeric(lfit$coefficients[2])
    )
}

fit_arrhenius_exdf <- function(exdf_obj, yname) {
    arrhenius_param <-
        fit_arrhenius(exdf_obj[, yname], exdf_obj[, 'arrhenius_index'])

    data.frame(
        species = exdf_obj[1, 'species'],
        parameter = yname,
        units = exdf_obj$units[[yname]],
        c = arrhenius_param$c,
        Ea = arrhenius_param$Ea
    )
}

fit_results_list <- by(rubisco_catalytic_values, rubisco_catalytic_values[, 'species'], function(x) {
    rbind(
        fit_arrhenius_exdf(x, 'Kc'),
        fit_arrhenius_exdf(x, 'Ko'),
        fit_arrhenius_exdf(x, 'Gamma_star')
    )
})

fit_results <- do.call(rbind, fit_results_list)

fit_results <- fit_results[order(fit_results$parameter, fit_results$species), ]

# Save the fit information
write.csv(fit_results, file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASENAME, '_arrhenius_parameters.csv')), row.names = FALSE)

# Get average Arrhenius parameters and save them
avg_fit_results_list <- by(fit_results, fit_results$parameter, function(x) {
    data.frame(
        parameter = x$parameter[1],
        units = x$units[1],
        c = mean(x$c),
        Ea = mean(x$Ea)
    )
})

avg_arrhenius_parameters <- do.call(rbind, avg_fit_results_list)

write.csv(avg_arrhenius_parameters, file.path(OUTPUT_DIR, TABLE_DIR, paste0(BASENAME, '_avg_arrhenius_parameters.csv')), row.names = FALSE)
save(avg_arrhenius_parameters, file = file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASENAME, '_avg_arrhenius_parameters.RData')))
