###               ###
### PRELIMINARIES ###
###               ###

# Load libraries
library(PhotoGEA)
library(lattice)

# Clear workspace
rm(list=ls())

# Load some helping functions
source(file.path('utilities', 'output_tools.R'))
source(file.path('utilities', 'defaults.R'))

source(file.path('utilities', 'validation_plots.R'))

# Choose some settings
SAVE_TO_PDF <- TRUE

REMOVE_BAD_CURVES <- TRUE

MAKE_VALIDATION_PLOTS <- TRUE

BASE_NAME <- 'tobacco_aci'

CI_LIM <- c(-100, 1500)

###                                  ###
### LOAD, VALIDATE, AND PROCESS DATA ###
###                                  ###

# Define a vector of paths to the files we wish to load
file_paths <- c(
    file.path('data', 'tobacco', '2023-03-21 ed aci mcgrath1'),
    file.path('data', 'tobacco', '2023-03-21 ed aci mcgrath2'),
    file.path('data', 'tobacco', '2023-03-23 ed aci mcgrath1'),
    file.path('data', 'tobacco', '2023-03-23 ed aci ripe5'),
    file.path('data', 'tobacco', '2023-03-24 ed aci mcgrath1'),
    file.path('data', 'tobacco', '2023-03-24 ed aci ripe5')
)

# Load each file, storing the result in a list. Add an empty PhiPS2 column if
# the file does not have this column.
licor_exdf_list <- lapply(file_paths, function(fpath) {
    lf <- read_gasex_file(fpath)
    if (!'PhiPS2' %in% colnames(lf)) {
        print(fpath)
        lf <- set_variable(lf, 'PhiPS2', category = 'FLR')
    }
    lf
})

# Get the names of all columns that are present in all of the Licor files
columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)

# Extract just these columns
licor_exdf_list <- lapply(licor_exdf_list, function(x) {
    x[ , columns_to_keep, TRUE]
})

# Use `rbind` to combine all the data
licor_data <- do.call(rbind, licor_exdf_list)

# Create a new column for the date as YYYY-MM-DD
licor_data[, 'date_ymd'] <- paste(
    substring(licor_data[, 'date'], 1, 4),
    substring(licor_data[, 'date'], 5, 6),
    substring(licor_data[, 'date'], 7, 8),
    sep = '-'
)

# Create a new identifier column formatted like `ppfd - replicate - instrument`
licor_data[ , 'curve_identifier'] <- paste(
    licor_data[ , 'ppfd'],
    licor_data[ , 'replicate'],
    licor_data[, 'instrument'],
    sep = ' - '
)

# Create a new identifier column formatted like `ppfd - replicate`
licor_data[ , 'ppfd_replicate'] <- paste(
    licor_data[ , 'ppfd'],
    licor_data[ , 'replicate'],
    sep = ' - '
)

# Set the units for the ppfd column
licor_data <- document_variables(
    licor_data,
    c('UserDefCon', 'ppfd', 'micromol m^(-2) s^(-1)')
)

# Make sure the data meets basic requirements, remove points with duplicated
# `CO2_r_sp`, and order by `Ci`
licor_data <- do.call(rbind, lapply(
    list(
        list(date_ymd = '2023-03-21', npts = 18, to_remove = c(12, 13)),
        list(date_ymd = '2023-03-23', npts = 19, to_remove = c(12, 13)),
        list(date_ymd = '2023-03-24', npts = 19, to_remove = c(12, 13))
    ),
    function(x) {
        check_response_curve_data(
            licor_data[licor_data[, 'date_ymd'] == x$date_ymd, , TRUE],
            'curve_identifier',
            x$npts,
            'CO2_r_sp'
        )

        organize_response_curve_data(
            licor_data[licor_data[, 'date_ymd'] == x$date_ymd, , TRUE],
            'curve_identifier',
            x$to_remove,
            'Ci'
        )
    }
))

# Order by ppfd and make sure PPFD will always be plotted in the same order
licor_data <- licor_data[order(as.numeric(licor_data[, 'ppfd'])), , TRUE]

# Remove bad curves, if necessary
if (REMOVE_BAD_CURVES) {
    licor_data <- remove_points(
        licor_data,
        list(curve_identifier = '100 - wt-3 - mcgrath1') # CO2 cylinder ran out during this curve
    )
}

# Save results for spreadsheet fits, including the total chamber pressure
licor_data <- calculate_total_pressure(licor_data)

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

save(
    file_paths, licor_exdf_list, licor_data,
    file = file.path(OUTPUT_DIR, RDATA_DIR, paste0(BASE_NAME, '_raw_data.RData'))
)

# Make validation plots if necessary
if (MAKE_VALIDATION_PLOTS) {
    validation_plots(licor_data, SAVE_TO_PDF, 'Ci', CI_LIM, BASE_NAME)
}
