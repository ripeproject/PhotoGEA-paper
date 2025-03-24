plot_synthetic_curves <- function(all_curves, base_name) {
    # Get a table of the main characteristics of each curve
    id_col <- c(
        'curve_identifier',
        'true_alpha_old',
        'true_n_Ac_limiting',
        'true_n_Aj_limiting',
        'true_n_Ap_limiting'
    )

    curve_info <- unique(all_curves[, id_col])

    # Find the curve with the smallest alpha_old among those that exhibit TPU
    # limitations
    tpu_curve_info <- curve_info[curve_info$true_n_Ap_limiting > 0, ]
    min_alpha_old <- min(tpu_curve_info$true_alpha_old)

    small_alpha_old_curve <-
        tpu_curve_info[tpu_curve_info$true_alpha_old == min_alpha_old, 'curve_identifier']

    selected_curves <- within(
        all_curves[all_curves[, 'curve_identifier'] == small_alpha_old_curve, ],
        {curve_type = 'Weak reverse sensitivity to CO2'}
    )

    # Find the largest alpha_old
    max_alpha_old <- max(tpu_curve_info$true_alpha_old)

    large_alpha_old_curve <-
        tpu_curve_info[tpu_curve_info$true_alpha_old == max_alpha_old, 'curve_identifier']

    selected_curves <- rbind(
        selected_curves,
        within(
            all_curves[all_curves[, 'curve_identifier'] == large_alpha_old_curve, ],
            {curve_type = 'Strong reverse sensitivity to CO2'}
        )
    )

    # Find curves with only one limiting process
    only_Ac_curves <- curve_info[curve_info$true_n_Ac_limiting > 0 & curve_info$true_n_Aj_limiting == 0 & curve_info$true_n_Ap_limiting == 0, 'curve_identifier']

    if (length(only_Ac_curves) > 0) {
        only_Ac_curve <- only_Ac_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == only_Ac_curve, ],
                {curve_type = 'Only Rubisco limitations'}
            )
        )
    }

    only_Aj_curves <- curve_info[curve_info$true_n_Ac_limiting == 0 & curve_info$true_n_Aj_limiting > 0 & curve_info$true_n_Ap_limiting == 0, 'curve_identifier']

    if (length(only_Aj_curves) > 0) {
        only_Aj_curve <- only_Aj_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == only_Aj_curve, ],
                {curve_type = 'Only RuBP regeneration limitations'}
            )
        )
    }

    only_Ap_curves <- curve_info[curve_info$true_n_Ac_limiting == 0 & curve_info$true_n_Aj_limiting == 0 & curve_info$true_n_Ap_limiting > 0, 'curve_identifier']

    if (length(only_Ap_curves) > 0) {
        only_Ap_curve <- only_Ap_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == only_Ap_curve, ],
                {curve_type = 'Only TPU limitations'}
            )
        )
    }

    # Find curves missing one limiting process
    missing_Ac_curves <- curve_info[curve_info$true_n_Ac_limiting == 0 & curve_info$true_n_Aj_limiting > 0 & curve_info$true_n_Ap_limiting > 0, 'curve_identifier']

    if (length(missing_Ac_curves) > 0) {
        missing_Ac_curve <- missing_Ac_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == missing_Ac_curve, ],
                {curve_type = 'No Rubisco limitations'}
            )
        )
    }

    missing_Aj_curves <- curve_info[curve_info$true_n_Ac_limiting > 0 & curve_info$true_n_Aj_limiting == 0 & curve_info$true_n_Ap_limiting > 0, 'curve_identifier']

    if (length(missing_Aj_curves) > 0) {
        missing_Aj_curve <- missing_Aj_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == missing_Aj_curve, ],
                {curve_type = 'No RuBP regeneration limitations'}
            )
        )
    }

    missing_Ap_curves <- curve_info[curve_info$true_n_Ac_limiting > 0 & curve_info$true_n_Aj_limiting > 0 & curve_info$true_n_Ap_limiting == 0, 'curve_identifier']

    if (length(missing_Ap_curves) > 0) {
        missing_Ap_curve <- missing_Ap_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == missing_Ap_curve, ],
                {curve_type = 'No TPU limitations'}
            )
        )
    }

    # A curve with all three processes
    full_curves <- curve_info[curve_info$true_n_Ac_limiting > 0 & curve_info$true_n_Aj_limiting > 0 & curve_info$true_n_Ap_limiting > 0, 'curve_identifier']

    if (length(full_curves) > 0) {
        full_curve <- full_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == full_curve, ],
                {curve_type = 'All three processes'}
            )
        )
    }

    # A curve with just one TPU-limited point
    one_tpu_curves <- curve_info[curve_info$true_n_Ap_limiting == 1, 'curve_identifier']

    if (length(one_tpu_curves) > 0) {
        one_tpu_curve <- one_tpu_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == one_tpu_curve, ],
                {curve_type = 'Just one TPU-limited point'}
            )
        )
    }

    # A curve with just one RuBP-regeneration-limited point
    one_rubp_curves <- curve_info[curve_info$true_n_Aj_limiting == 1, 'curve_identifier']

    if (length(one_rubp_curves) > 0) {
        one_rubp_curve <- one_rubp_curves[1]

        selected_curves <- rbind(
            selected_curves,
            within(
                all_curves[all_curves[, 'curve_identifier'] == one_rubp_curve, ],
                {curve_type = 'Just one RuBP-regeneration-limited point'}
            )
        )
    }

    # Plot the representative curves
    pdf_print(
        xyplot(
            A ~ Ci | curve_type,
            data = selected_curves,
            type = 'b',
            pch = 16,
            xlab = 'Ci [ micromol mol^(-1) ]',
            ylab = 'A [ micromol mol^(−2) s^(−1) ]'
        ),
        width = 12,
        height = 9,
        save_to_pdf = SAVE_TO_PDF,
        file = file.path(OUTPUT_DIR, FIGURE_DIR, paste0(base_name, '_selected_curves.pdf'))
    )

}
