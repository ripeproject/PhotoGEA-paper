library(latticeExtra)

# Define a helping function that makes a bunch of basic plots showing
# information about all the curves in a data set. Used for Supplemental Figures
# S5, S6, and S14.
validation_plots <- function(exdf_obj, save_to_pdf, x_name, x_lim, base_name) {
    # Make a plot to check humidity control
    pdf_print(
        xyplot(
            RHcham + `Humidifier_%` + `Desiccant_%` ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            ylim = c(0, 100),
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            xlim = x_lim
        ),
        width = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_humidity_control.pdf'))
    )

    # Make a plot to check temperature control
    pdf_print(
        xyplot(
            TleafCnd + Txchg ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('Temperature [', exdf_obj$units$TleafCnd, ']'),
            xlim = x_lim
        ),
        width = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_temperature_control.pdf'))
    )

    # Make a plot to check CO2 control
    pdf_print(
        xyplot(
            CO2_s + CO2_r + CO2_r_sp ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('CO2 concentration [', exdf_obj$units$CO2_r, ']'),
            xlim = x_lim
        ),
        width = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_co2_control.pdf'))
    )

    # Make a plot to check Q control
    pdf_print(
        xyplot(
            Qin ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('Incident PPFD [', exdf_obj$units$Qin, ']'),
            xlim = x_lim
        ),
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_Q_control.pdf'))
    )

    # Make a plot to check stability criteria
    pdf_print(
        xyplot(
            `A:OK` + `gsw:OK` + Stable ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            xlim = x_lim
        ),
        width = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_stability_criteria.pdf'))
    )

    # Plot each A curve in a separate panel
    pdf_print(
        xyplot(
            A ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('Net CO2 assimilation rate [', exdf_obj$units$A, ']'),
            xlim = x_lim
        ),
        width = 20,
        height = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_A_separate.pdf'))
    )

    # Plot all A curves on one set of axes
    pdf_print(
        xyplot(
            A ~ exdf_obj[, x_name],
            group = curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto.key = list(space = 'right'),
            grid = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('Net CO2 assimilation rate [', exdf_obj$units$A, ']'),
            xlim = x_lim
        ),
        width = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_A_together.pdf'))
    )

    # Plot each gsw curve in a separate panel
    pdf_print(
        xyplot(
            gsw ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            grid = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('Stomatal conductance to H2O [', exdf_obj$units$gsw, ']'),
            xlim = x_lim
        ),
        width = 20,
        height = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_gsw_separate.pdf'))
    )

    # Overlay A-Ci and gsw-Ci curves
    plot_obj <- xyplot(
        A ~ exdf_obj[, x_name] | curve_identifier,
        data = exdf_obj$main_data,
        type = 'b',
        pch = 16,
        auto.key = list(space = 'right'),
        xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
        ylab = paste('Net CO2 assimilation rate [', exdf_obj$units$A, ']'),
        xlim = x_lim,
        ylim = c(-5, 45),
        scales = list(alternating = 1, tck = c(1,0)),
        col = c('black')
    ) + as.layer(
        xyplot(
            gsw ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto = TRUE,
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('Stomatal conductance to H2O [', exdf_obj$units$gsw, ']'),
            xlim = x_lim,
            ylim = c(0, 0.7),
            scales = list(alternating = 2, tck = c(0, 1)),
            col = c('deeppink')
        ),
        x.same = TRUE,
        y.same = FALSE,
        outside = TRUE
    )

    pdf_print(
        plot_obj,
        width = 20,
        height = 10,
        save_to_pdf = save_to_pdf,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_a_gsw_overlay.pdf'))
    )

    if ('PhiPS2' %in% colnames(exdf_obj)) {
        # Plot each PhiPS2 curve in a separate panel
        pdf_print(
            xyplot(
                PhiPS2 ~ exdf_obj[, x_name] | curve_identifier,
                data = exdf_obj$main_data,
                type = 'b',
                pch = 16,
                auto = TRUE,
                grid = TRUE,
                xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
                ylab = 'PhiPS2 (dimensionless)',
                xlim = x_lim,
                ylim = c(0, 0.8)
            ),
            width = 20,
            height = 10,
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_PhiPSII_separate.pdf'))
        )

        # Plot all PhiPS2 curves on one set of axes
        pdf_print(
            xyplot(
                PhiPS2 ~ exdf_obj[, x_name],
                group = curve_identifier,
                data = exdf_obj$main_data,
                type = 'b',
                pch = 16,
                auto.key = list(space = 'right'),
                grid = TRUE,
                xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
                ylab = 'PhiPS2 (dimensionless)',
                xlim = x_lim
            ),
            width = 10,
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_PhiPSII_together.pdf'))
        )

        # Overlay A-Ci and PhiPS2-Ci curves
        plot_obj <- xyplot(
            A ~ exdf_obj[, x_name] | curve_identifier,
            data = exdf_obj$main_data,
            type = 'b',
            pch = 16,
            auto.key = list(space = 'right'),
            xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
            ylab = paste('Net CO2 assimilation rate [', exdf_obj$units$A, ']'),
            xlim = x_lim,
            ylim = c(-5, 75),
            scales = list(alternating = 1, tck = c(1,0)),
            col = c('black')
        ) + as.layer(
            xyplot(
                PhiPS2 ~ exdf_obj[, x_name] | curve_identifier,
                data = exdf_obj$main_data,
                type = 'b',
                pch = 16,
                auto = TRUE,
                xlab = paste(x_name, '[', exdf_obj$units[[x_name]], ']'),
                ylab = paste('PhiPS2 [', exdf_obj$units$PhiPS2, ']'),
                xlim = x_lim,
                ylim = c(0, 0.8),
                scales = list(alternating = 2, tck = c(0, 1)),
                col = c('cornflowerblue')
            ),
            x.same = TRUE,
            y.same = FALSE,
            outside = TRUE
        )

        pdf_print(
            plot_obj,
            width = 15,
            height = 15,
            save_to_pdf = save_to_pdf,
            file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_a_phips2_overlay.pdf'))
        )
    }
}
