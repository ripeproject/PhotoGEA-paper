# Helping function for plotting each individual best-fit value and upper/lower
# confindence limits vs. Qin. Not used for any figures in the main text or
# supplemental information.
plot_parameters <- function(
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
    # Plot n limits vs. Q
    if ('n_Wc_smallest' %in% colnames(dataf)) {
        pdf_print(
            xyplot(
                n_Wc_smallest + n_Wj_smallest + n_Wp_smallest ~ ppfd,
                data = dataf,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-1, 20),
                xlab = 'Qin',
                ylab = 'n lowest'
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_nlowest_q.pdf'))
        )
    }

    # Plot alpha_old vs. Q
    if ('alpha_old_lower' %in% colnames(dataf)) {
        pdf_print(
            xyplot(
                alpha_old_lower + alpha_old + alpha_old_upper ~ ppfd,
                data = dataf,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.1, 1.1),
                xlab = 'Qin',
                ylab = 'alpha_old'
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_alpha_old_q.pdf'))
        )
    }

    # Plot alpha_g vs. Q
    if ('alpha_g_lower' %in% colnames(dataf)) {
        pdf_print(
            xyplot(
                alpha_g_lower + alpha_g + alpha_g_upper ~ ppfd,
                data = dataf,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.1, 1.1),
                xlab = 'Qin',
                ylab = 'alpha_g'
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_alpha_g_q.pdf'))
        )
    }

    # Plot alpha_s vs. Q
    if ('alpha_s_lower' %in% colnames(dataf)) {
        pdf_print(
            xyplot(
                alpha_s_lower + alpha_s + alpha_s_upper ~ ppfd,
                data = dataf,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.1, 1.1),
                xlab = 'Qin',
                ylab = 'alpha_s'
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_alpha_s_q.pdf'))
        )
    }

    # Plot Gamma_star vs. Q
    if ('Gamma_star_lower' %in% colnames(dataf)) {
        pdf_print(
            xyplot(
                Gamma_star_lower + Gamma_star + Gamma_star_upper ~ ppfd,
                data = dataf,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = gstar_lim,
                xlab = 'Qin',
                ylab = 'Gamma_star'
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_gstar_q.pdf'))
        )
    }

    # Plot RL vs. Q
    pdf_print(
        xyplot(
            RL_at_25_lower + RL_at_25 + RL_at_25_upper ~ ppfd,
            data = dataf,
            type = line_type,
            pch = 16,
            grid = TRUE,
            auto = list(space = 'top'),
            xlim = qin_lim,
            ylim = rl_lim,
            xlab = 'Qin',
            ylab = 'RL_at_25'
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_rl_q.pdf'))
    )

    # Plot J vs. Q
    pdf_print(
        xyplot(
            J_at_25_lower + J_at_25 + J_at_25_upper ~ ppfd,
            data = dataf,
            type = line_type,
            pch = 16,
            grid = TRUE,
            auto = list(space = 'top'),
            xlim = qin_lim,
            ylim = j_lim,
            xlab = 'Qin',
            ylab = 'J_at_25'
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_j_q.pdf'))
    )

    # Plot tau vs. Q
    if ('tau_lower' %in% colnames(dataf)) {
        pdf_print(
            xyplot(
                tau_lower + tau + tau_upper ~ ppfd,
                data = dataf,
                type = line_type,
                pch = 16,
                grid = TRUE,
                auto = list(space = 'top'),
                xlim = qin_lim,
                ylim = c(-0.05, 1.05),
                xlab = 'Qin',
                ylab = 'tau'
            ),
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_tau_q.pdf'))
        )
    }

    # Plot Tp vs. Q
    pdf_print(
        xyplot(
            Tp_at_25_lower + Tp_at_25 + Tp_at_25_upper ~ ppfd,
            data = dataf,
            type = line_type,
            pch = 16,
            grid = TRUE,
            auto = list(space = 'top'),
            xlim = qin_lim,
            ylim = tp_lim,
            xlab = 'Qin',
            ylab = 'Tp_at_25'
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_tp_q.pdf'))
    )

    # Plot Vcmax vs. Q
    pdf_print(
        xyplot(
            Vcmax_at_25_lower + Vcmax_at_25 + Vcmax_at_25_upper ~ ppfd,
            data = dataf,
            type = line_type,
            pch = 16,
            grid = TRUE,
            auto = list(space = 'top'),
            xlim = qin_lim,
            ylim = vcmax_lim,
            xlab = 'Qin',
            ylab = 'Vcmax_at_25'
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_vcmax_q.pdf'))
    )
}
