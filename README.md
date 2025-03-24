## PhotoGEA-paper

This repository includes R scripts and input data that were used to produce
the figures for the manuscript titled "PhotoGEA: An R Package for Closer Fitting
of Photosynthetic Gas Exchange Data With Non-Gaussian Confidence Interval
Estimation"

These scripts have been tested using the following installations:
- Windows:
  - R version 4.3.2 (2023-10-31 ucrt)
  - Windows 10 Enterprise version 10.0.19045 Build 19045
- Windows:
  - R version 4.4.2 (2024-10-31 ucrt)
  - Microsoft Windows 11 Enterprise version 10.0.26100 Build 26100

All outputs were generated using these scripts and stored in the
`output_archive` directory. For instructions for running these scripts, see the
[Reproducing the Outputs](#reproducing-the-outputs) section below. For more
information about the purpose of each script, see the [Scripts](#scripts)
section.

Figures and tables in the manuscript were produced from the files in
`output_archive/figures` and `output_archive/tables`. For figures, Adobe
Illustrator was used to perform the final combinations, recoloring, annotations,
etc. For more details, see the
[Main Figures and Tables](#main-figures-and-tables) section below.

## Reproducing the Outputs

### Requirements
- The [R environment](https://cran.r-project.org/)
- The `lattice` R package
  - To install this package, type: `install.packages('lattice')`
- The `RColorBrewer` R package
  - To install this package, type: `install.packages('RColorBrewer')`
- The `latticeExtra` R package
  - To install this package, type: `install.packages('latticeExtra')`
- The `lhs` R package
  - To install this package, type: `install.packages('lhs')`
- The `remotes` R package
  - This is used to install some of the other packages
  - To install this package, type: `install.packages('remotes')`
- The `PhotoGEA` R package
  - Version 1.1.0
  - To install this version, type: `remotes::install_github('eloch216/PhotoGEA', ref = 'v1.1.0')`
- The `plantecophys` R package
  - Version 1.4-6
  - To install this version, type: `remotes::install_github('RemkoDuursma/plantecophys', ref = 'c974982')`
- The `photosynthesis` R package
  - Version 2.1.4
  - To install this version, type: `remotes::install_github('cdmuir/photosynthesis', ref = 'be43024')`
- The `msuRACiFit` R package:
  - Version 1.2.0
  - To install this version, type: `remotes::install_github('poales/msuRACiFit', ref = 'a8f2da7')`

These scripts may work with other versions of `PhotoGEA`, `plantecophys`,
`photosynthesis`, and `msuRACiFit`, but it is not guaranteed that the same
results will be obtained when using versions other than the ones specified
above.

### Steps
- Start a fresh R session and set the working directory to this one
- Run the main script: `source('run_all.R')`
- All outputs (will be generated in a new `output` directory, whose content
  should be identical to the `output_archive` directory provided with the
  repository.
- In total, the scripts may take a few hours to run.

### Additional Comments
- The main script (`run_all.R`) simply calls several other scripts. As an
  alternative to sourcing that script, individual analysis scripts can be
  sourced instead to produce a subset of all outputs.
- Many settings can be changed in these scripts, as well as in the settings and
  helping functions defined in the `utilities` directory.

## Scripts
- `run_all.R`
  - Calls the other scripts in an appropriate order.
- `check_packages.R`
  - This script checks to make sure the required packages can be installed, and
    makes sure the correct version of certain packages is installed.
- `load_tobacco_data.R`
  - This script reads the Licor files that contain the tobacco A-Ci curve
    measurements and stores them in a convenient R data file for other scripts
    to use.
  - It also exports each curve to an individual CSV file, to be used when
    fitting the curves using the `PCE Calculator` spreadsheet tool.
  - The names of output files produced by this script are prefixed with
    `tobacco_aci`.
- `synthetic_c3_aci_fits.R`
  - This script generates synthetic C3 A-Ci curves and fits them using several
    tools.
  - It is necessary to run `load_tobacco_data.R` before running this script.
  - It is also necessary to define an R variable called `ALPHA_TYPE` before
    running this script. Its value must be `none`, `alpha_old`, or `alpha_gs`.
  - With default settings, this script generates and fits 200 curves, so it may
    take more than 30 minutes to run, especially when `ALPHA_TYPE` is
    `alpha_gs`.
  - The names of output files produced by this script are prefixed with
    `synthetic_aci_` and `ALPHA_TYPE`. For example, `synthetic_aci_alpha_old`.
- `tobacco_aci_fits.R`
  - This script fits each tobacco A-Ci curve using several fitting tools.
  - It is necessary to run `load_tobacco_data.R` before running this script.
  - The names of output files produced by this script are prefixed with
    `tobacco_aci`.
- `soybean_rubisco.R`
  - This script fits an Arrhenius equation to temperature-dependendent values
    of soybean `Kc`, `Ko`, and `Gamma_star` values as reported in Orr et al.
    (2016).
  - The names of output files produced by this script are prefixed with
    `rubisco`.
- `soybean_variable_j_fits.R`
  - This script fits each soybean A-Ci + CF curve using the Variable J method.
  - It is necessary to run `soybean_rubisco.R` before running this script.
  - The names of output files produced by this script are prefixed with
    `soybean_variable_j`.
- `c4_aci_fits.R`
  - This script fits each maize and sorghum A-Ci curve using several
    approaches.
  - The names of output files produced by this script are prefixed with
    `c4_aci`.

## Main Figures and Tables
- Figure 1 was created from:
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_phips2ci.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_aci_fit_comparison.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_aci_parameter_comparison_J.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_aci_parameter_comparison_rl.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_aci_parameter_comparison_Tp.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_aci_parameter_comparison_Vcmax.pdf`
- Figure 2 was created from:
  - `output_archive/figures/tobacco_aci_avg_acis.pdf`
  - `output_archive/figures/tobacco_aci_compare_rmse_avg.pdf`
  - `output_archive/figures/tobacco_aci_compare_j_avg.pdf`
  - `output_archive/figures/tobacco_aci_param_avg_eb_j_q_ci.pdf`
  - `output_archive/figures/tobacco_aci_compare_vcmax_avg.pdf`
  - `output_archive/figures/tobacco_aci_param_avg_eb_vcmax_q_ci.pdf`
  - `output_archive/figures/tobacco_aci_compare_tp_avg.pdf`
  - `output_archive/figures/tobacco_aci_param_avg_eb_tp_q_ci.pdf`
  - `output_archive/figures/tobacco_aci_compare_rl_avg.pdf`
  - `output_archive/figures/tobacco_aci_compare_rl_laisk.pdf`
  - `output_archive/figures/tobacco_aci_param_avg_eb_rl_q_ci.pdf`
- Figure 3 was created from:
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_Vcmax_at_25_error.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_Tp_at_25_error.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_RL_at_25_error.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_alpha_old_error.pdf`
  - `output_archive/figures/tobacco_aci_800 - wt-5 - mcgrath1_J_at_25_zoom_error.pdf`
  - `output_archive/figures/tobacco_aci_1000 - wt-1 - mcgrath1_J_at_25_zoom_error.pdf`
- Figure 4 was created from:
  - `output_archive/figures/soybean_variable_j_fits_ci.pdf`
  - `output_archive/figures/soybean_variable_j_fits_cc.pdf`
  - `output_archive/figures/soybean_variable_j_2022 - ripe2 - 4_gmc_fits.pdf`
  - `output_archive/figures/soybean_variable_j_2022 - ripe2 - 4_jci_fits.pdf`
  - `output_archive/figures/soybean_variable_j_avg_gmcs.pdf`
- Figure 5 was created from:
  - `output_archive/c4_aci_aci_fits_v4.pdf`
  - `output_archive/c4_aci_aci_fits_v5.pdf`
  - `output_archive/c4_aci_aci_fits_v1.pdf`
  - `output_archive/c4_aci_aci_fits_v2.pdf`
  - `output_archive/c4_aci_compare_fits_vcmax_25.pdf`
  - `output_archive/c4_aci_compare_fits_j_25.pdf`
  - `output_archive/c4_aci_compare_fits_vpmax_25.pdf`

- Table 2 was created from `soybean_variable_j_variable_j_parameters_clean.csv`

## License
This repository is licensed under the CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/)
