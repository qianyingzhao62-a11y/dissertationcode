# dissertation_figures.R
# Reproduces all 7 figures in "Model and Calibration-Window Risk in UK
# Longevity" from the raw HMD GBR_NP files.
#
#   Fig 3.1  mortality surface ln m(x,t) + period e65          (data only)
#   Fig 3.2  rates by age, selected years                      (data only)
#   Fig 3.3  rates over time at ages 65/75/85                  (data only)
#   Fig 3.4  annual improvement, ages 60-89, with means        (data only)
#   Fig 5.1  out-of-sample validation at age 75 (LC & CBD)     (StMoMo)
#   Fig 5.2  best-estimate annuity by model x window           (StMoMo, ~min)
#   Fig 5.3  longevity capital SCR% heatmap                    (StMoMo, ~min)
#
# Usage: place this file in the same directory as GBR_data/ and run it
# section by section (RStudio: Ctrl+Enter). Part 0 must be run first;
# Part 1 needs only the data; Parts 2-3 need StMoMo (Part 3 takes ~2-5 min).
# Output: 7 PNG files in ./figures/ plus capital_grid_R.csv (use it to
# update Table 5.2).
# Note: the numbers produced by Part 3 are the results of MY run; the
# figures in Tables 5.1/5.2/5.3 and in the text of the dissertation should
# be updated to match them. Small differences from any earlier prototype
# results are normal and genuine.
# ============================================================================

## ---- PART 0 . packages, palette, data -------------------------------------
# install.packages(c("ggplot2", "patchwork", "StMoMo"))   
# install once
