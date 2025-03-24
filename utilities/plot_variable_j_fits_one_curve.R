# Helping function that plots Variable J fitting results for a single curve.
# Used for Figure 4a-d.
plot_variable_j_fits_one_curve <- function(
    all_fits,
    curve_to_plot,
    ci_lim,
    gmc_lim,
    j_lim,
    base_name
)
{
    curve_df <- all_fits[all_fits[, 'curve_identifier'] == curve_to_plot, ]

    pdf_print(
        xyplot(gmc ~ Ci,
            data = curve_df,
            type = 'b',
            pch = 16,
            grid = TRUE,
            xlab = paste('Ci [', all_fits$units$Ci, ']'),
            ylab = paste('gmc [', all_fits$units$gmc, ']'),
            xlim = ci_lim,
            ylim = gmc_lim
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_to_plot, '_gmc_fits.pdf'))
    )

    pdf_print(
        xyplot(gmc ~ Ci,
            group = A > 0,
            data = curve_df,
            type = 'b',
            pch = 16,
            grid = TRUE,
            xlab = paste('Ci [', all_fits$units$Ci, ']'),
            ylab = paste('gmc [', all_fits$units$gmc, ']'),
            xlim = ci_lim,
            ylim = gmc_lim
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_to_plot, '_gmc_fits_a_group.pdf'))
    )

    pdf_print(
        xyplot(J_F + J_tl + ETR ~ Ci,
            data = curve_df,
            type = 'b',
            pch = 16,
            auto.key = list(space = 'right', lines = TRUE, points = TRUE),
            grid = TRUE,
            xlab = paste('Ci [', all_fits$units$Ci, ']'),
            ylab = paste('J [', all_fits$units$J_F, ']'),
            xlim = ci_lim,
            ylim = j_lim
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_to_plot, '_jci_fits.pdf'))
    )
}
