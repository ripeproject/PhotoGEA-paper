# Utilities

This directory mostly contains scripts that define helping functions used in the
main analysis scripts. Most are related to plotting. It is helpful to separate
these plotting commands from the main scripts because they are often long and
complicated. In more detail:

- `defaults.R` specifies some global default settings used throughout the
  analysis scripts.

- `output_tools.R` ensures that the output directories exist, and specifies
  their names so they can be more easily referenced elsewhere in a consistent
  way.

- `extract_fit_param.R` is used to combine the results of parameter estimates
  from multiple C3 A-Ci curve fitting tools.

The remaining scripts are all related to creating various types of plots.
