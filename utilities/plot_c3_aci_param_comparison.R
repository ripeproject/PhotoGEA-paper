# Plot individual parameter estimates (and confidence intervals) for a single C3
# curve, as reported by each fitting tool. Used for Figure 1g-j.
plot_c3_aci_param_comparison <- function(
    curve_for_package_comparison,
    all_photogea_param,
    all_photogea_new_alpha_param,
    all_plantecophys_param,
    all_photosynthesis_param,
    all_msu_param,
    all_spreadsheet_param,
    base_name
)
{
    # Get fit parameters for just this curve
    photogea_param           <- all_photogea_param[all_photogea_param[, 'curve_identifier'] == curve_for_package_comparison, ]
    photogea_new_alpha_param <- all_photogea_new_alpha_param[all_photogea_new_alpha_param[, 'curve_identifier'] == curve_for_package_comparison, ]
    plantecophys_param       <- all_plantecophys_param[all_plantecophys_param[, 'curve_identifier'] == curve_for_package_comparison, ]
    photosynthesis_param     <- all_photosynthesis_param[all_photosynthesis_param[, 'curve_identifier'] == curve_for_package_comparison, ]
    msu_param                <- all_msu_param[all_msu_param[, 'curve_identifier'] == curve_for_package_comparison, ]
    spreadsheet_param        <- all_spreadsheet_param[all_spreadsheet_param[, 'curve_identifier'] == curve_for_package_comparison, ]

    # Extract confidence intervals for each package
    all_param <- rbind(
        data.frame(
            J       = with(photogea_param, {c(J_tl_avg_lower,     J_tl_avg,     J_tl_avg_upper)}),
            RL      = with(photogea_param, {c(RL_tl_avg_lower,    RL_tl_avg,    RL_tl_avg_upper)}),
            Tp      = with(photogea_param, {c(Tp_tl_avg_lower,    Tp_tl_avg,    Tp_tl_avg_upper)}),
            Vcmax   = with(photogea_param, {c(Vcmax_tl_avg_lower, Vcmax_tl_avg, Vcmax_tl_avg_upper)}),
            package = 'PhotoGEA (old alpha)'
        ),
        data.frame(
            J       = with(photogea_new_alpha_param, {c(J_tl_avg_lower,     J_tl_avg,     J_tl_avg_upper)}),
            RL      = with(photogea_new_alpha_param, {c(RL_tl_avg_lower,    RL_tl_avg,    RL_tl_avg_upper)}),
            Tp      = with(photogea_new_alpha_param, {c(Tp_tl_avg_lower,    Tp_tl_avg,    Tp_tl_avg_upper)}),
            Vcmax   = with(photogea_new_alpha_param, {c(Vcmax_tl_avg_lower, Vcmax_tl_avg, Vcmax_tl_avg_upper)}),
            package = 'PhotoGEA (new alpha)'
        ),
        data.frame(
            J       = with(plantecophys_param, {c(J_lower,     J,     J_upper)}),
            RL      = with(plantecophys_param, {c(Rd_lower,    Rd,    Rd_upper)}),
            Tp      = with(plantecophys_param, {c(TPU_lower,   TPU,   TPU_upper)}),
            Vcmax   = with(plantecophys_param, {c(Vcmax_lower, Vcmax, Vcmax_upper)}),
            package = 'plantecophys'
        ),
        data.frame(
            J       = with(photosynthesis_param, {c(J_lower,      J,      J_upper)}),
            RL      = with(photosynthesis_param, {c(R_d_lower,    R_d,    R_d_upper)}),
            Tp      = with(photosynthesis_param, {c(V_TPU_lower,  V_TPU,  V_TPU_upper)}),
            Vcmax   = with(photosynthesis_param, {c(V_cmax_lower, V_cmax, V_cmax_upper)}),
            package = 'photosynthesis'
        ),
        data.frame(
            J       = with(msu_param, {c(J_lower,     J,     J_upper)}),
            RL      = with(msu_param, {c(rL_lower,    rL,    rL_upper)}),
            Tp      = with(msu_param, {c(TPU_lower,   TPU,   TPU_upper)}),
            Vcmax   = with(msu_param, {c(VcMax_lower, VcMax, VcMax_upper)}),
            package = 'msuRACiFit'
        ),
        data.frame(
            J       = with(spreadsheet_param, {c(J_tl_lower,     J_tl,     J_tl_upper)}),
            RL      = with(spreadsheet_param, {c(RL_tl_lower,    RL_tl,    RL_tl_upper)}),
            Tp      = with(spreadsheet_param, {c(Tp_tl_lower,    Tp_tl,    Tp_tl_upper)}),
            Vcmax   = with(spreadsheet_param, {c(Vcmax_tl_lower, Vcmax_tl, Vcmax_tl_upper)}),
            package = 'PCE calculator'
        )
    )

    all_param$package <- factor(
        all_param$package,
        levels = c('PhotoGEA (old alpha)', 'PCE calculator', 'PhotoGEA (new alpha)', 'msuRACiFit', 'plantecophys', 'photosynthesis')
    )

    write.csv(
        all_param,
        file = file.path(OUTPUT_DIR, TABLE_DIR, paste0(base_name, '_', curve_for_package_comparison, '_aci_parameter_comparison.csv')),
        row.names = FALSE
    )

    # In PhotoGEA, upper limits may be Inf, which will interfere with plotting.
    # Here we temporarily set them to NA when determining axis limits.
    all_param_tmp <- all_param

    for (j in seq_len(ncol(all_param_tmp))) {
        if (all_param_tmp[3, j] == Inf) {
            all_param_tmp[3, j] <- NA # PhotoGEA (old alpha)
        }
        if (all_param_tmp[6, j] == Inf) {
            all_param_tmp[6, j] <- NA # PhotoGEA (new alpha)
        }
    }

    # Determine limits for plots
    j_lim <- c(
        min(all_param_tmp$J, na.rm = TRUE),
        max(all_param_tmp$J, na.rm = TRUE)
    )

    j_lim <- j_lim + (j_lim - mean(j_lim))

    rl_lim <- c(
        min(all_param_tmp$RL, na.rm = TRUE),
        max(all_param_tmp$RL, na.rm = TRUE)
    )

    rl_lim <- rl_lim + (rl_lim - mean(rl_lim))

    tp_lim <- c(
        min(all_param_tmp$Tp, na.rm = TRUE),
        max(all_param_tmp$Tp, na.rm = TRUE)
    )

    tp_lim <- tp_lim + (tp_lim - mean(tp_lim))

    vcmax_lim <- c(
        min(all_param_tmp$Vcmax, na.rm = TRUE),
        max(all_param_tmp$Vcmax, na.rm = TRUE)
    )

    vcmax_lim <- vcmax_lim + (vcmax_lim - mean(vcmax_lim))

    # In PhotoGEA, upper limits may be Inf, which will interfere with plotting.
    # Here we set them to a large number to make sure they show up in plots.
    for (j in seq_len(ncol(all_param))) {
        if (all_param[3, j] == Inf) {
            all_param[3, j] <- 1e6 # PhotoGEA (old alpha)
        }
        if (all_param[6, j] == Inf) {
            all_param[6, j] <- 1e6 # PhotoGEA (new alpha)
        }
    }

    # In PhotoGEA, best-fit values may be NA when there is still a lower limit.
    # Here we set them to the lower limit to make sure they show up in plots.
    for (j in seq_len(ncol(all_param))) {
        if (is.na(all_param[2, j])) {
            all_param[2, j] <- all_param[1, j] # PhotoGEA (old alpha)
        }
        if (is.na(all_param[5, j])) {
            all_param[5, j] <- all_param[4, j] # PhotoGEA (new alpha)
        }
    }

    # Make plots
    height = 4

    pdf_print(
        xyplot(
            J ~ package,
            group = package,
            data = all_param,
            type = 'b',
            pch = 16,
            ylim = j_lim,
            ylab = paste('J [', all_photogea_param$units$J_tl_avg, ']'),
            main = curve_for_package_comparison
        ),
        height = height,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_for_package_comparison, '_aci_parameter_comparison_J.pdf'))
    )

    pdf_print(
        xyplot(
            RL ~ package,
            group = package,
            data = all_param,
            type = 'b',
            pch = 16,
            ylim = rl_lim,
            ylab = paste('RL [', all_photogea_param$units$RL_tl_avg, ']'),
            main = curve_for_package_comparison
        ),
        height = height,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_for_package_comparison, '_aci_parameter_comparison_rl.pdf'))
    )

    pdf_print(
        xyplot(
            Tp ~ package,
            group = package,
            data = all_param,
            type = 'b',
            pch = 16,
            ylim = tp_lim,
            ylab = paste('Tp [', all_photogea_param$units$Tp_tl_avg, ']'),
            main = curve_for_package_comparison
        ),
        height = height,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_for_package_comparison, '_aci_parameter_comparison_Tp.pdf'))
    )

    pdf_print(
        xyplot(
            Vcmax ~ package,
            group = package,
            data = all_param,
            type = 'b',
            pch = 16,
            ylim = vcmax_lim,
            ylab = paste('Vcmax [', all_photogea_param$units$Vcmax_tl_avg, ']'),
            main = curve_for_package_comparison
        ),
        height = height,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_', curve_for_package_comparison, '_aci_parameter_comparison_Vcmax.pdf'))
    )
}
