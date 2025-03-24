# Specify the packages used by the scripts
required_packages <- c(
    'lattice',
    'lhs',
    'msuRACiFit',
    'PhotoGEA',
    'photosynthesis',
    'plantecophys',
    'RColorBrewer',
    'latticeExtra'
)

# Check to make sure these packages can be loaded
package_loaded <- sapply(required_packages, function(pkgname) {
    require(pkgname, character.only = TRUE)
})

# Send an error message if any were not found
if (any(!package_loaded)) {
    msg <- paste0(
        'The following required R packages could not be loaded: ',
        paste(required_packages[!package_loaded], collapse = ', '),
        '. See README.md for installation instructions.'
    )
    stop(msg)
}

# Specify specific versions for some of the packages
required_versions <- list(
    msuRACiFit = '1.2.0',
    PhotoGEA = '1.1.0',
    photosynthesis = '2.1.4',
    plantecophys = '1.4-6'
)

correct_version <- sapply(seq_along(required_versions), function(i) {
    required_versions[[i]] == packageVersion(names(required_versions)[i])
})

# Send an error message if any versions are incorrect
if (any(!correct_version)) {
    msg <- paste0(
        'Some installed packages do not have the required version: ',
        paste(names(required_versions)[!correct_version], collapse = ', '),
        '. See README.md for installation instructions.'
    )
    stop(msg)
}
