# Helping function for plotting each individual C3 A-Ci parameter estimate vs.
# Qin, grouped by fitting tool. Used for Figure 2c-f.
plot_parameter_comparisons <- function(
    dataf,
    qin_lim,
    rl_lim,
    j_lim,
    vcmax_lim,
    tp_lim,
    rmse_lim,
    aic_lim,
    alpha_lim,
    RMSE_breaks,
    base_name,
    save_to_pdf
)
{
    # Convert `package` to a factor so we can control the order of traces when
    # plotting
    packages <- c(
        'msuRACiFit',
        'photosynthesis',
        'plantecophys',
        'PCE calculator',
        'PhotoGEA (new alpha)',
        'PhotoGEA'
    )

    dataf$package <- factor(
        dataf$package,
        levels = packages
    )

    # Make a set of colors
    param_col <- c(
        '#00798C', # msuRACiFit
        '#D1495B', # photosynthesis
        '#EDAE49', # plantecophys
        '#BD93DA', # PCE calculator
        '#999999', # PhotoGEA (new alpha)
        '#000000'  # PhotoGEA
    )

    # Make a set of symbols
    param_pch <- c(
        1,  # msuRACiFit
        1,  # photosynthesis
        1,  # plantecophys
        1,  # PCE calculator
        16, # PhotoGEA (new alpha)
        16  # PhotoGEA
    )

    # Plots of each parameter value vs. Q
    pdf_print(
        histogram(
            ~ RMSE | package,
            data = dataf,
            scales = list(x = list(log = TRUE)),
            ylim = c(0, 50),
            layout = c(1, length(levels(dataf$package))),
            breaks = RMSE_breaks
        ),
        width = 5,
        height = 12,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_hist_RMSE.pdf'))
    )

    pdf_print(
        xyplot(
            RMSE ~ ppfd,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            scales = list(y = list(log = TRUE)),
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'RMSE',
            xlim = qin_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_RMSE.pdf'))
    )

    pdf_print(
        xyplot(
            AIC ~ ppfd,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            scales = list(y = list(log = FALSE)),
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'AIC',
            xlim = qin_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_AIC.pdf'))
    )

    pdf_print(
        xyplot(
            J ~ ppfd,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'J [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_j.pdf'))
    )

    pdf_print(
        xyplot(
            RL ~ ppfd,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'RL [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_rl.pdf'))
    )

    pdf_print(
        xyplot(
            Vcmax ~ ppfd,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'Vcmax [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_vcmax.pdf'))
    )

    # Plots of each parameter value from each curve
    pdf_print(
        xyplot(
            J ~ curve_identifier,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            scales = list(x = list(rot = 90)),
            ylab = 'J [ micromol m^(-2) s^(-1) ]',
            ylim = j_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_curve_j.pdf'))
    )

    pdf_print(
        xyplot(
            RL ~ curve_identifier,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            scales = list(x = list(rot = 90)),
            ylab = 'RL [ micromol m^(-2) s^(-1) ]',
            ylim = rl_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_curve_rl.pdf'))
    )

    pdf_print(
        xyplot(
            Vcmax ~ curve_identifier,
            group = package,
            data = dataf,
            type = 'p',
            auto = list(space = 'top'),
            scales = list(x = list(rot = 90)),
            ylab = 'Vcmax [ micromol m^(-2) s^(-1) ]',
            ylim = vcmax_lim,
            par.settings = list(
                superpose.line = list(col = param_col),
                superpose.symbol = list(col = param_col, pch = param_pch)
            )
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_curve_vcmax.pdf'))
    )

    # Plots of average parameter values vs Q
    error_bars <- FALSE

    rmse_order <- c(1, 2, 3, 5, 6)

    packages_rmse  <- packages[rmse_order]
    param_pch_rmse <- param_pch[rmse_order]
    param_col_rmse <- param_col[rmse_order]

    dataf_rmse <- dataf[dataf$package %in% packages_rmse, ]

    dataf_rmse$package <- factor(
        as.character(dataf_rmse$package),
        levels = packages_rmse
    )

    pdf_print(
        xyplot_avg_rc(
            dataf_rmse$RMSE,
            dataf_rmse$ppfd,
            dataf_rmse$ppfd,
            dataf_rmse$package,
            type = 'b',
            pch = param_pch_rmse,
            auto = list(space = 'top'),
            scales = list(y = list(log = TRUE)),
            y_error_bars = FALSE,
            cols = param_col_rmse,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'RMSE [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = rmse_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_rmse_avg.pdf'))
    )

    pdf_print(
        xyplot_avg_rc(
            dataf_rmse$AIC,
            dataf_rmse$ppfd,
            dataf_rmse$ppfd,
            dataf_rmse$package,
            type = 'b',
            pch = param_pch_rmse,
            auto = list(space = 'top'),
            scales = list(y = list(log = FALSE)),
            y_error_bars = FALSE,
            cols = param_col_rmse,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'AIC',
            xlim = qin_lim,
            ylim = aic_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_aic_avg.pdf'))
    )

    pdf_print(
        xyplot_avg_rc(
            dataf$J,
            dataf$ppfd,
            dataf$ppfd,
            dataf$package,
            type = 'b',
            pch = param_pch,
            auto = list(space = 'top'),
            y_error_bars = error_bars,
            cols = param_col,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'J [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = j_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_j_avg.pdf'))
    )

    pdf_print(
        xyplot_avg_rc(
            dataf$RL,
            dataf$ppfd,
            dataf$ppfd,
            dataf$package,
            type = 'b',
            pch = param_pch,
            auto = list(space = 'top'),
            y_error_bars = error_bars,
            cols = param_col,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'RL [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = rl_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_rl_avg.pdf'))
    )

    pdf_print(
        xyplot_avg_rc(
            dataf$Tp,
            dataf$ppfd,
            dataf$ppfd,
            dataf$package,
            type = 'b',
            pch = param_pch,
            auto = list(space = 'top'),
            y_error_bars = error_bars,
            cols = param_col,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'Tp [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = tp_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_tp_avg.pdf'))
    )

    pdf_print(
        xyplot_avg_rc(
            dataf$Vcmax,
            dataf$ppfd,
            dataf$ppfd,
            dataf$package,
            type = 'b',
            pch = param_pch,
            auto = list(space = 'top'),
            y_error_bars = error_bars,
            cols = param_col,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'Vcmax [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = vcmax_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_vcmax_avg.pdf'))
    )

    alpha_old_indx <- c(2, 6)

    alpha_old_pkg        <- packages[alpha_old_indx]
    alpha_old_df         <- dataf[dataf$package %in% alpha_old_pkg, ]
    alpha_old_df$package <- factor(alpha_old_df$package, levels = alpha_old_pkg)
    alpha_old_col        <- param_col[alpha_old_indx]
    alpha_old_pch        <- param_pch[alpha_old_indx]

    pdf_print(
        xyplot_avg_rc(
            alpha_old_df$alpha_old,
            alpha_old_df$ppfd,
            alpha_old_df$ppfd,
            alpha_old_df$package,
            type = 'b',
            pch = alpha_old_pch,
            auto = list(space = 'top'),
            y_error_bars = error_bars,
            cols = alpha_old_col,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'alpha_old [ dimensionless ]',
            xlim = qin_lim,
            ylim = alpha_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_alpha_old_avg.pdf'))
    )

    alpha_new_indx <- c(1, 3)

    alpha_new_pkg        <- packages[alpha_new_indx]
    alpha_new_df         <- dataf[dataf$package %in% alpha_new_pkg, ]
    alpha_new_df$package <- factor(alpha_new_df$package, levels = alpha_new_pkg)
    alpha_new_col        <- param_col[alpha_new_indx]
    alpha_new_pch        <- param_pch[alpha_new_indx]

    pdf_print(
        xyplot_avg_rc(
            alpha_new_df$alpha_g,
            alpha_new_df$ppfd,
            alpha_new_df$ppfd,
            alpha_new_df$package,
            type = 'b',
            pch = alpha_new_pch,
            auto = list(space = 'top'),
            y_error_bars = error_bars,
            cols = alpha_new_col,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'alpha_g [ dimensionless ]',
            xlim = qin_lim,
            ylim = alpha_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_alpha_g_avg.pdf'))
    )

    pdf_print(
        xyplot_avg_rc(
            alpha_new_df$alpha_s,
            alpha_new_df$ppfd,
            alpha_new_df$ppfd,
            alpha_new_df$package,
            type = 'b',
            pch = alpha_new_pch,
            auto = list(space = 'top'),
            y_error_bars = error_bars,
            cols = alpha_new_col,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'alpha_s [ dimensionless ]',
            xlim = qin_lim,
            ylim = alpha_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_alpha_s_avg.pdf'))
    )

    # Plot individual estimates, separated by packages
    pdf_print(
        xyplot(
            J ~ ppfd | package,
            data = fit_param,
            type = 'p',
            pch = 16,
            auto = TRUE,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'J [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = j_lim,
            panel = function(...) {
                for (q in unique(fit_param$ppfd)) {
                    panel.lines(
                        j_lim ~ c(q, q),
                        col = 'gray',
                        lty = 2
                    )
                }
                panel.xyplot(...)
            }
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_j_pkg.pdf'))
    )

    pdf_print(
        xyplot(
            Tp ~ ppfd | package,
            data = fit_param,
            type = 'p',
            pch = 16,
            auto = TRUE,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'Tp [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = tp_lim,
            panel = function(...) {
                for (q in unique(fit_param$ppfd)) {
                    panel.lines(
                        tp_lim ~ c(q, q),
                        col = 'gray',
                        lty = 2
                    )
                }
                panel.xyplot(...)
            }
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_tp_pkg.pdf'))
    )

    pdf_print(
        xyplot(
            Vcmax ~ ppfd | package,
            data = fit_param,
            type = 'p',
            pch = 16,
            auto = TRUE,
            xlab = 'Qin [ micromol m^(-2) s^(-1) ]',
            ylab = 'Vcmax [ micromol m^(-2) s^(-1) ]',
            xlim = qin_lim,
            ylim = vcmax_lim,
            panel = function(...) {
                for (q in unique(fit_param$ppfd)) {
                    panel.lines(
                        vcmax_lim ~ c(q, q),
                        col = 'gray',
                        lty = 2
                    )
                }
                panel.xyplot(...)
            }
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_compare_vcmax_pkg.pdf'))
    )
}
