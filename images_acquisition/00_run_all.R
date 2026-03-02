# ============================================================
# run_all.R — Full pipeline runner
# ============================================================

# Packages
pkgs <- c(
  "dplyr","tidyr","stringr","lubridate","forcats","purrr","tibble",
  "glmmTMB","DHARMa","emmeans","ggplot2","performance","readr",
  "car","MuMIn"
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# Run scripts (order matters)
source(file.path("images_acquisition/01_functions.R"))
source(file.path("images_acquisition/02_prepare_fit_models.R"))
source(file.path("images_acquisition/03_tables_figures.R"))
