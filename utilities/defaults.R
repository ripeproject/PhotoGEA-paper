# Set some default values to use in multiple places
CA_ATMOSPHERIC <- 420

ATP_USE <- 4
NADPH_USE <- 8

TEMPERATURE_PARAM <- c3_temperature_param_sharkey

HARD_CONSTRAINTS <- 0

REQUIRE_POSITIVE_GMC <- 'positive_a'

GMC_MAX <- Inf

FIT_OPTIONS <- list(alpha_old = 'fit', alpha_g = 0, alpha_s = 0, alpha_t = 0)

C3_ACI_ITERMAX <- 200

VJ_ITERMAX <- 800

C4_ACI_ITERMAX <- 200

ERROR_THRESHOLD_FACTOR <- 0.147

NICE_PPFD <- c(100, 150, 200, 250, 300, 400, 500, 1000, 1500)

library(RColorBrewer)

NICE_PPFD_COLORS <- colorRampPalette(rev(brewer.pal(11, "Spectral")[8:11]))(length(NICE_PPFD))

# Specify parameters related to the response of J to Qin. alpha is the initial
# slope and theta is the curvature. These parameters are only used by the
# plantecophys and photosynthesis tools. In principle they should only influence
# the values of Jmax, as determined from the "raw" estimates of J. Yet, in
# practice, they can drastically alter the fit quality, and the estimated values
# of all parameters. Both packages use alpha = 0.24 and theta = 0.85 as the
# default values, which were found to produce some very bad fits. Here we use
# values found in von Caemmerer (2000).
ALPHA_J <- 0.85 * (1 - 0.15) / 2 # von Caemmerer (2000) (Equation 2.14)
THETA_J <- 0.7                   # von Caemmerer (2000)
