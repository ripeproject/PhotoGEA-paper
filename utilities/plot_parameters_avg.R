# Helping function used internally by plot_parameters_avg
get_param_df <- function(dataf, pname) {
    lname <- paste0(pname, '_lower')
    uname <- paste0(pname, '_upper')

    bdf <- dataf[, c('ppfd', pname)]
    bdf$ptype <- 'base'

    ldf <- dataf[, c('ppfd', lname)]
    ldf$ptype <- 'lower'
    colnames(ldf) <- c('ppfd', pname, 'ptype')

    udf <- dataf[, c('ppfd', uname)]
    udf$ptype <- 'upper'
    colnames(udf) <- c('ppfd', pname, 'ptype')

    rbind(bdf, ldf, udf)
}

# Helping function for plotting average best-fit values and upper/lower
# confindence limits vs. Qin. Used for Figure 2c-f.
plot_parameters_avg <- function(
    dataf,
    qin_lim,
    j_lim,
    line_type,
    base_name,
    save_to_pdf,
    rl_lim,
    tp_lim,
    vcmax_lim,
    gstar_lim = c(0, 100)
)
{
    # Replace any infinities with NA before calculating averages. Then the
    # upper limits will only be Inf if all of them at that Qin are Inf.
    dataf <- data.frame(lapply(dataf, function(x) {
        x[!is.finite(x)] <- NA
        x
    }))

    # Plot alpha_old vs. Q
    if ('alpha_old_lower' %in% colnames(dataf)) {
        param_df <- get_param_df(dataf, 'alpha_old')

        pdf_print(
            with(param_df, {xyplot_avg_rc(
                alpha_old,
                ppfd,
                ppfd,
                ptype,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.1, 2.1),
                xlab = 'Qin',
                ylab = 'alpha_old'
            )}),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_alpha_old_q.pdf'))
        )
    }

    # Plot alpha_g vs. Q
    if ('alpha_g_lower' %in% colnames(dataf)) {
        param_df <- get_param_df(dataf, 'alpha_g')

        pdf_print(
            with(param_df, {xyplot_avg_rc(
                alpha_g,
                ppfd,
                ppfd,
                ptype,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.1, 1.1),
                xlab = 'Qin',
                ylab = 'alpha_g'
            )}),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_alpha_g_q.pdf'))
        )
    }

    # Plot alpha_s vs. Q
    if ('alpha_s_lower' %in% colnames(dataf)) {
        param_df <- get_param_df(dataf, 'alpha_s')

        pdf_print(
            with(param_df, {xyplot_avg_rc(
                alpha_s,
                ppfd,
                ppfd,
                ptype,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.1, 1.1),
                xlab = 'Qin',
                ylab = 'alpha_s'
            )}),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_alpha_s_q.pdf'))
        )
    }

    # Plot Gamma_star vs. Q
    if ('Gamma_star_lower' %in% colnames(dataf) && !is.na(dataf[1, 'Gamma_star'])) {
        param_df <- get_param_df(dataf, 'Gamma_star')

        pdf_print(
            with(param_df, {xyplot_avg_rc(
                Gamma_star,
                ppfd,
                ppfd,
                ptype,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = gstar_lim,
                xlab = 'Qin',
                ylab = 'Gamma_star'
            )}),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_gstar_q.pdf'))
        )
    }

    # Plot RL vs. Q
    param_df <- get_param_df(dataf, 'RL_tl_avg')

    pdf_print(
        with(param_df, {xyplot_avg_rc(
            RL_tl_avg,
            ppfd,
            ppfd,
            ptype,
            type = line_type,
            pch = 16,
            grid = TRUE,
            auto = list(space = 'top'),
            xlim = qin_lim,
            ylim = rl_lim,
            xlab = 'Qin',
            ylab = 'RL_tl_avg'
        )}),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_rl_q.pdf'))
    )

    # Plot just the confidence intervals
    tmp <- as.data.frame(tapply(
        param_df$RL_tl_avg,
        list(param_df$ppfd, param_df$ptype),
        function(x) {mean(x, na.rm = TRUE)}
    ))

    tmp$ppfd <- as.numeric(rownames(tmp))
    tmp$upper[is.infinite(tmp$upper)] <- 1e6
    tmp$base[is.na(tmp$base)] <- 1e4

    pdf_print(
        xyplot(
            base ~ ppfd,
            data = tmp,
            type = 'l',
            xlim = qin_lim,
            ylim = rl_lim,
            xlab = 'Qin',
            ylab = 'RL_tl_avg',
            panel = function(...) {
                panel.xyplot(...)

                for (i in seq_len(nrow(tmp))) {
                    panel.arrows(
                        tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'upper'],
                        angle = 90,
                        col = 'black',
                        length = 0.05,
                        lwd = 1
                    )

                    panel.arrows(
                        tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'lower'],
                        angle = 90,
                        col = 'black',
                        length = 0.05,
                        lwd = 1
                    )
                }
            }
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_rl_q_ci.pdf'))
    )

    # Plot J vs. Q
    param_df <- get_param_df(dataf, 'J_tl_avg')

    pdf_print(
        with(param_df, {xyplot_avg_rc(
            J_tl_avg,
            ppfd,
            ppfd,
            ptype,
            type = line_type,
            pch = 16,
            grid = TRUE,
            auto = list(space = 'top'),
            xlim = qin_lim,
            ylim = j_lim,
            xlab = 'Qin',
            ylab = 'J_tl_avg'
        )}),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_j_q.pdf'))
    )

    # Plot just the confidence intervals
    tmp <- as.data.frame(tapply(
        param_df$J_tl_avg,
        list(param_df$ppfd, param_df$ptype),
        function(x) {mean(x, na.rm = TRUE)}
    ))

    tmp$ppfd <- as.numeric(rownames(tmp))
    tmp$upper[is.infinite(tmp$upper)] <- 1e6
    tmp$base[is.na(tmp$base)] <- 1e4

    pdf_print(
        xyplot(
            base ~ ppfd,
            data = tmp,
            type = 'l',
            xlim = qin_lim,
            ylim = j_lim,
            xlab = 'Qin',
            ylab = 'J_tl_avg',
            panel = function(...) {
                panel.xyplot(...)

                for (i in seq_len(nrow(tmp))) {
                    panel.arrows(
                        tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'upper'],
                        angle = 90,
                        col = 'black',
                        length = 0.05,
                        lwd = 1
                    )

                    panel.arrows(
                        tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'lower'],
                        angle = 90,
                        col = 'black',
                        length = 0.05,
                        lwd = 1
                    )
                }
            }
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_j_q_ci.pdf'))
    )

    # Plot tau vs. Q
    if ('tau_lower' %in% colnames(dataf)) {
        param_df <- get_param_df(dataf, 'tau')

        pdf_print(
            with(param_df, {xyplot_avg_rc(
                tau,
                ppfd,
                ppfd,
                ptype,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.05, 1.05),
                xlab = 'Qin',
                ylab = 'tau'
            )}),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_tau_q.pdf'))
        )
    }

    # Plot Tp vs. Q
    if ('Tp_tl_avg_lower' %in% colnames(dataf)) {
        param_df <- get_param_df(dataf, 'Tp_tl_avg')

        pdf_print(
            with(param_df, {xyplot_avg_rc(
                Tp_tl_avg,
                ppfd,
                ppfd,
                ptype,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = tp_lim,
                xlab = 'Qin',
                ylab = 'Tp_tl_avg'
            )}),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_tp_q.pdf'))
        )

        # Plot just the confidence intervals
        tmp <- as.data.frame(tapply(
            param_df$Tp_tl_avg,
            list(param_df$ppfd, param_df$ptype),
            function(x) {mean(x, na.rm = TRUE)}
        ))

        tmp$ppfd <- as.numeric(rownames(tmp))
        tmp$upper[is.infinite(tmp$upper)] <- 1e6
        tmp$base[is.na(tmp$base)] <- 1e4

        pdf_print(
            xyplot(
                base ~ ppfd,
                data = tmp,
                type = 'l',
                xlim = qin_lim,
                ylim = tp_lim,
                xlab = 'Qin',
                ylab = 'Tp_tl_avg',
                panel = function(...) {
                    panel.xyplot(...)

                    for (i in seq_len(nrow(tmp))) {
                        panel.arrows(
                            tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'upper'],
                            angle = 90,
                            col = 'black',
                            length = 0.05,
                            lwd = 1
                        )

                        panel.arrows(
                            tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'lower'],
                            angle = 90,
                            col = 'black',
                            length = 0.05,
                            lwd = 1
                        )
                    }
                }
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_tp_q_ci.pdf'))
        )
    }


    # Plot Vcmax vs. Q
    param_df <- get_param_df(dataf, 'Vcmax_tl_avg')

    pdf_print(
        with(param_df, {xyplot_avg_rc(
            Vcmax_tl_avg,
            ppfd,
            ppfd,
            ptype,
            type = line_type,
            pch = 16,
            grid = TRUE,
            auto = list(space = 'top'),
            xlim = qin_lim,
            ylim = vcmax_lim,
            xlab = 'Qin',
            ylab = 'Vcmax_tl_avg'
        )}),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_vcmax_q.pdf'))
    )

    # Plot just the confidence intervals
    tmp <- as.data.frame(tapply(
        param_df$Vcmax_tl_avg,
        list(param_df$ppfd, param_df$ptype),
        function(x) {mean(x, na.rm = TRUE)}
    ))

    tmp$ppfd <- as.numeric(rownames(tmp))
    tmp$upper[is.infinite(tmp$upper)] <- 1e6
    tmp$base[is.na(tmp$base)] <- 1e4

    pdf_print(
        xyplot(
            base ~ ppfd,
            data = tmp,
            type = 'l',
            xlim = qin_lim,
            ylim = vcmax_lim,
            xlab = 'Qin',
            ylab = 'Vcmax_tl_avg',
            panel = function(...) {
                panel.xyplot(...)

                for (i in seq_len(nrow(tmp))) {
                    panel.arrows(
                        tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'upper'],
                        angle = 90,
                        col = 'black',
                        length = 0.05,
                        lwd = 1
                    )

                    panel.arrows(
                        tmp[i, 'ppfd'], tmp[i, 'base'], tmp[i, 'ppfd'], tmp[i, 'lower'],
                        angle = 90,
                        col = 'black',
                        length = 0.05,
                        lwd = 1
                    )
                }
            }
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_vcmax_q_ci.pdf'))
    )
}
