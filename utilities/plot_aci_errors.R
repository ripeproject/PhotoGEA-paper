# A helping function for plotting lkelihood vs. parameter values. Used for
# Figure 3.
plot_aci_errors <- function(
    curve_for_errors, # an exdf object
    fits_interpolated,
    best_fit_param,   # data frame of parameters found for the curve by fit_c3_aci
    ci_lim,
    j_zoom_lim,
    ej_thresh,
    base_name,
    save_to_pdf,
    ej_thresh_upper = Inf
)
{
    use_vj <- 'tau' %in% colnames(best_fit_param)
    use_c4 <- 'Vpmax_at_25' %in% colnames(best_fit_param)

    fit_options <- if (use_vj) {
        c(FIT_OPTIONS, list(J_at_25 = 'fit', RL_at_25 = 'fit', tau = 'fit', Tp_at_25 = 'fit', Vcmax_at_25 = 'fit'))
    } else if (use_c4) {
        list(alpha_psii = 'fit', gbs = 'fit', RL_at_25 = 'fit', Rm_frac = 'fit', Vcmax_at_25 = 'fit', Vpmax_at_25 = 'fit', Vpr = 'fit')
    }else {
        c(FIT_OPTIONS, list(J_at_25 = 'fit', RL_at_25 = 'fit', Tp_at_25 = 'fit', Vcmax_at_25 = 'fit'))
    }

    fit_options <- fit_options[fit_options == 'fit']

    best_fit_param_df <- best_fit_param$main_data

    sd_A <- if (is.na(best_fit_param_df$RMSE)) {
        0.0
    } else {
        best_fit_param_df$RMSE
    }

    # Get the associated error function
    erf <- if (use_vj) {
        # Here we need the "v3" settings, which are being used as the defaults
        error_function_c3_variable_j(
            curve_for_errors,
            fit_options = fit_options,
            atp_use = ATP_USE,
            nadph_use = NADPH_USE,
            sd_A = sd_A,
            hard_constraints = HARD_CONSTRAINTS,
            require_positive_gmc = REQUIRE_POSITIVE_GMC,
            gmc_max = GMC_MAX
        )
    } else if (use_c4) {
        error_function_c4_aci(
            curve_for_errors,
            fit_options = fit_options,
            sd_A = sd_A,
            hard_constraints = HARD_CONSTRAINTS
        )
    } else {
        error_function_c3_aci(
            curve_for_errors,
            fit_options = fit_options,
            atp_use = ATP_USE,
            nadph_use = NADPH_USE,
            sd_A = sd_A,
            hard_constraints = HARD_CONSTRAINTS
        )
    }

    # Specify parameters to vary and the number of values to use for each
    param <- names(fit_options)

    npts <- 5001

    # Get the best-fit values for each of these parameters.
    base_inputs <- sapply(param, function(pn) {best_fit_param_df[[pn]]})

    # Get the smallest value of the error function
    best_val <- best_fit_param_df$optimum_val

    # Make a vector of the error threshold to use for plotting
    ethresh <- rep_len(ERROR_THRESHOLD_FACTOR, npts)

    # Set error limits to use for plotting
    elim <- c(0, 1.1)

    # Set padding to use for plots
    pad <- 0.05

    # Make a line plot for each parameter
    for (i in seq_along(param)) {
        # Get the current parameter name
        pn <- param[i]

        # Get its best-fit info
        lower_val <- best_fit_param_df[[paste0(pn, '_lower')]]
        base_val  <- best_fit_param_df[[pn]]
        upper_val <- best_fit_param_df[[paste0(pn, '_upper')]]

        # Skip this parameter if it was not actually fit
        if (is.na(lower_val) || is.na(upper_val)) {
            next
        }

        # Create a sequence of parameter values to use
        pseq <- if (pn %in% c('Rm_frac')) {
            # We can always test the full range of some parameters
            seq(-0.05, 1.05, length.out = npts)
        } else {
            width <- (upper_val - lower_val) * (1 + pad)
            if (abs(base_val) > 1e-10) {
                seq(
                    max(base_val - width, 0),
                    min(base_val + width, base_val * 2.0),
                    length.out = npts
                )
            } else {
                seq(
                    base_val - width,
                    base_val + width,
                    length.out = npts
                )
            }
        }

        # Calculate the error along this sequence
        eseq <- sapply(pseq, function(x) {
            inputs <- base_inputs
            inputs[i] <- x
            exp(-erf(inputs) + best_val)
        })

        # Plot the results
        pdf_print(
            xyplot(
                eseq + ethresh ~ pseq,
                type = 'l',
                xlab = pn,
                ylab = 'Likelihood',
                main = base_name,
                xlim = c(min(pseq), max(pseq)),
                ylim = elim,
                panel = function(...) {
                    panel.xyplot(...)
                    panel.points(
                        c(ERROR_THRESHOLD_FACTOR, 1, ERROR_THRESHOLD_FACTOR) ~ c(lower_val, base_val, upper_val),
                        type = 'p',
                        col = 'red',
                        pch = 16
                    )
                }
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', pn, '_error.pdf'))
        )

        # Make a zoomed-in version for Tp_at_25
        if (pn == 'Tp_at_25') {
            pseq <- seq(2.4, 3.0, length.out = npts)

            eseq <- sapply(pseq, function(x) {
                inputs <- base_inputs
                inputs[i] <- x
                exp(-erf(inputs) + best_val)
            })

            pdf_print(
                xyplot(
                    eseq + ethresh ~ pseq,
                    type = 'l',
                    xlab = pn,
                    ylab = 'Likelihood',
                    main = base_name,
                    xlim = c(min(pseq), 3),
                    ylim = elim,
                    panel = function(...) {
                        panel.xyplot(...)
                        panel.points(
                            c(ERROR_THRESHOLD_FACTOR, 1, ERROR_THRESHOLD_FACTOR) ~ c(lower_val, base_val, upper_val),
                            type = 'p',
                            col = 'red',
                            pch = 16
                        )
                    }
                ),
                save_to_pdf = save_to_pdf,
                file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', pn, '_zoom_error.pdf'))
            )
        }


        # Make a zoomed-in version for J
        if (pn == 'J_at_25') {
            pseq <- seq(min(j_zoom_lim), max(j_zoom_lim), length.out = npts)

            eseq <- sapply(pseq, function(x) {
                inputs <- base_inputs
                inputs[i] <- x
                exp(-erf(inputs) + best_val)
            })

            # Do a really quick gaussian fit over a range of the J values
            pseq_sub <- pseq[eseq >= ej_thresh & eseq <= ej_thresh_upper]
            eseq_sub <- eseq[eseq >= ej_thresh & eseq <= ej_thresh_upper]

            f <- function(p) {
                p[3] * dnorm(pseq_sub, mean = p[1], sd = p[2])
            }

            pbest <- optim(
                c(min(base_val, 190), 4, 10),
                function(p) {sum((eseq_sub - f(p))^2)},
                method = 'L-BFGS-B',
                lower = c(0, 0, 0),
                upper = c(1e6, 1e6, 1e6)
            )

            eseq_gauss <- sapply(pseq, function(x) {
                pbest$par[3] * dnorm(x, mean = pbest$par[1], sd = pbest$par[2])
            })

            pdf_print(
                xyplot(
                    eseq + ethresh + eseq_gauss ~ pseq,
                    type = 'l',
                    xlab = pn,
                    ylab = 'Likelihood',
                    main = base_name,
                    xlim = c(min(pseq), max(pseq)),
                    ylim = elim,
                    panel = function(...) {
                        panel.xyplot(...)
                        panel.points(
                            c(ERROR_THRESHOLD_FACTOR, 1, ERROR_THRESHOLD_FACTOR) ~ c(lower_val, base_val, upper_val),
                            type = 'p',
                            col = 'red',
                            pch = 16
                        )
                    }
                ),
                save_to_pdf = save_to_pdf,
                file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', pn, '_zoom_error.pdf'))
            )
        }
    }

    # Also plot A-Ci and PhiPS2-Ci curves
    if (!use_c4) {
        pdf_print(
            PhotoGEA::plot_c3_aci_fit(
                list(fits = curve_for_errors, fits_interpolated = fits_interpolated, parameters = best_fit_param),
                'curve_identifier',
                'Ci',
                xlim = ci_lim,
                ylim = c(
                    min(curve_for_errors[, 'A']) * (1 + 9 * pad),
                    max(curve_for_errors[, 'A']) * (1 + 3 * pad)
                )
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_aci.pdf'))
        )
    }

    pdf_print(
        xyplot(
            PhiPS2 ~ Ci,
            data = curve_for_errors$main_data,
            type = 'b',
            pch = 16,
            auto.key = list(space = 'right', lines = TRUE, points = TRUE),
            grid = TRUE,
            xlab = paste('Ci [', curve_for_errors$units$Ci, ']'),
            ylab = 'PhiPS2 [ dimensionless ]',
            xlim = ci_lim,
            ylim = c(0, 0.8)
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_phips2ci.pdf'))
    )

    return(invisible())
}
