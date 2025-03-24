# Plot the C3 A-Ci fits for all curves, as fit by a single fitting tool
# (including limiting rates and, potentially, the operating point). Used for
# Supplemental Figures S1-4.
plot_c3_aci_fits <- function(
    exdf_obj,
    a_column_name,
    a_fit_column_name,
    base_name,
    ci_lim,
    a_lim,
    plot_operating_point,
    parameter_exdf = NULL
)
{
    tmp_df <- exdf_obj$main_data
    colnames(tmp_df)[colnames(tmp_df) == a_column_name]     <- 'A'
    colnames(tmp_df)[colnames(tmp_df) == a_fit_column_name] <- 'A_fit'

    pdf_print(
        xyplot(
            A + Ac + Aj + Ap + A_fit ~ Ci | curve_identifier,
            data = tmp_df,
            type = 'b',
            auto.key = list(space = 'right', lines = TRUE, points = TRUE),
            grid = TRUE,
            xlab = 'Ci [ micromol mol^(-1) ]',
            ylab = 'An [ micromol m^(-2) s^(-1) ]',
            xlim = ci_lim,
            ylim = a_lim,
            par.settings = list(
                superpose.line = list(col = multi_curve_line_colors(), lty = c(1, 2, 5, 4, 1), lwd = 2),
                superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
            ),
            curve_ids = exdf_obj[, 'curve_identifier'],
            panel = function(...) {
                panel.xyplot(...)
                if (plot_operating_point) {
                    args <- list(...)
                    curve_id <- args$curve_ids[args$subscripts][1]
                    fit_param <-
                        parameter_exdf[parameter_exdf[, 'curve_identifier'] == curve_id, ]
                    panel.points(
                        fit_param$operating_An_model ~ fit_param$operating_Ci,
                        type = 'p',
                        col = 'black',
                        pch = 1
                    )
                }
            }
        ),
        width = 15,
        height = 15,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_aci_fits.pdf'))
    )
}
