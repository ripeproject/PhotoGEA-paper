# A helping function for plotting multiple fits of a single C3 curve (one fit
# from each of the tested fitting tools, where the PCE calculator fit is
# optional). Used for Figure 1b-f.
plot_c3_aci_fit_comparison <- function(
    curve_for_package_comparison,
    all_photogea_fits,
    all_plantecophys_fits,
    all_photosynthesis_fits,
    all_msu_fits,
    ci_lim,
    a_lim,
    base_name,
    spreadsheet_fit = NULL
)
{
    include_spreadsheet_fit <- !is.null(spreadsheet_fit)

    col_to_keep <- c('curve_identifier', 'Ci', 'A', 'A_fit', 'Ac', 'Aj', 'Ap')

    photogea_fit <- all_photogea_fits[all_photogea_fits[, 'curve_identifier'] == curve_for_package_comparison, col_to_keep]
    photogea_fit$package <- 'PhotoGEA'
    colnames(photogea_fit) <- c(col_to_keep, 'package')

    plantecophys_fit <- all_plantecophys_fits[all_plantecophys_fits[, 'curve_identifier'] == curve_for_package_comparison, c('curve_identifier', 'Ci', 'Ameas', 'Amodel', 'Ac', 'Aj', 'Ap')]
    plantecophys_fit$package <- 'plantecophys'
    colnames(plantecophys_fit) <- c(col_to_keep, 'package')

    photosynthesis_fit <- all_photosynthesis_fits[all_photosynthesis_fits[, 'curve_identifier'] == curve_for_package_comparison, c('curve_identifier', 'Ci', 'A_net', 'A_model', 'Ac', 'Aj', 'Ap')]
    photosynthesis_fit$package <- 'photosynthesis'
    colnames(photosynthesis_fit) <- c(col_to_keep, 'package')

    msu_fit <- all_msu_fits[all_msu_fits[, 'curve_identifier'] == curve_for_package_comparison, col_to_keep]
    msu_fit$package <- 'msuRACiFit'
    colnames(msu_fit) <- c(col_to_keep, 'package')

    fit_comparison <- rbind(
        photogea_fit,
        plantecophys_fit,
        photosynthesis_fit,
        msu_fit
    )

    if (include_spreadsheet_fit) {
        spreadsheet_fit$package <- 'PCE calculator'

        fit_comparison <- rbind(
            fit_comparison,
            spreadsheet_fit[, c(col_to_keep, 'package')]
        )
    }

    packages <- if (include_spreadsheet_fit) {
        c('PhotoGEA', 'plantecophys', 'photosynthesis', 'msuRACiFit', 'PCE calculator')
    } else {
        c('PhotoGEA', 'plantecophys', 'photosynthesis', 'msuRACiFit')
    }

    fit_comparison$package <- factor(
        fit_comparison$package,
        levels = packages
    )

    pdf_print(
        xyplot(
            A + Ac + Aj + Ap + A_fit ~ Ci | package,
            data = fit_comparison,
            type = 'b',
            auto.key = list(space = 'right', lines = TRUE, points = TRUE),
            grid = TRUE,
            xlab = paste('Ci [', all_photogea_fits$units$Ci, ']'),
            ylab = paste('A [', all_photogea_fits$units$A, ']'),
            xlim = ci_lim,
            ylim = a_lim,
            layout = c(length(packages), 1),
            par.settings = list(
                superpose.line = list(col = multi_curve_line_colors(), lty = c(1, 2, 5, 4, 1), lwd = 2),
                superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
            )
        ),
        height = 3.5,
        width = 12,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_for_package_comparison, '_aci_fit_comparison.pdf'))
    )

    pdf_print(
        xyplot(
            A_fit ~ Ci,
            group = package,
            data = fit_comparison,
            type = 'l',
            auto.key = list(space = 'right', lines = TRUE, points = FALSE),
            grid = TRUE,
            xlab = paste('Ci [', all_photogea_fits$units$Ci, ']'),
            ylab = paste('A [', all_photogea_fits$units$A, ']'),
            xlim = ci_lim,
            ylim = a_lim,
            par.settings = list(
                superpose.line = list(col = multi_curve_colors())
            )
        ),
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_for_package_comparison, '_aci_fit_comparison2.pdf'))
    )
}
