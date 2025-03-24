# Define a helping function for getting average parameter values at each PPFD.
# Exclude any infinite or NA values from the averages.
param_lim_mean <- function(x, pn) {
    pn_lower <- paste0(pn, '_lower')
    pn_upper <- paste0(pn, '_upper')

    lower <- x[, pn_lower]
    best  <- x[, pn]
    upper <- x[, pn_upper]

    lower[is.infinite(lower)] <- NA
    best[is.infinite(best)]   <- NA
    upper[is.infinite(upper)] <- NA

    tmp <- data.frame(
        lower_avg = mean(lower,   na.rm = TRUE),
        best_avg =  mean(x[, pn], na.rm = TRUE),
        upper_avg = mean(upper,   na.rm = TRUE)
    )

    colnames(tmp) <- c(pn_lower, pn, pn_upper)
    tmp
}
