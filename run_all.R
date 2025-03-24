# Clear the workspace
rm(list=ls())

# Check the packages
print('running: check_packages.R')
source('check_packages.R')

# Load the tobacco data
print('running: load_tobacco_data.R')
source('load_tobacco_data.R')

# Analyze synthetic C3 A-Ci curves
print('running: synthetic_c3_aci_fits.R with alpha = 0')
ALPHA_TYPE <- 'none'; source('synthetic_c3_aci_fits.R')

print('running: synthetic_c3_aci_fits.R with variable alpha_old')
ALPHA_TYPE <- 'alpha_old'; source('synthetic_c3_aci_fits.R')

print('running: synthetic_c3_aci_fits.R with variable alpha_g and alpha_s')
ALPHA_TYPE <- 'alpha_gs'; source('synthetic_c3_aci_fits.R')

# Analyze tobacco A-Ci curves
print('running: tobacco_aci_fits.R')
source('tobacco_aci_fits.R')

# Get soybean Rubisco Arrhenius parameters
print('running: soybean_rubisco.R')
source('soybean_rubisco.R')

# Analyze soybean Variable J curves
print('running: soybean_variable_j_fits.R')
source('soybean_variable_j_fits.R')

# Analyze C4 A-Ci curves
print('running: c4_aci_fits.R')
source('c4_aci_fits.R')
