# DownscalingFABLE/R/run_country.R

# ==================================================================
# USER SETTINGS
# ==================================================================

# Select the configuration file to use for this run.
# The file must be located in the config/ folder.
config_file <- "UZB.yml"


# ------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------

# Restore the project environment once after cloning the repository:
#install.packages("renv")
#renv::restore()

# In normal use, load the installed package.
# In developer mode, devtools::load_all() will already have loaded it.
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all("../FABLEDownscalR", quiet = TRUE)
}

# remotes::install_github("davidecozza99/FABLEDownscalR", dependencies = TRUE, force = TRUE) #Run if FABLEDownscalR has been modified
remotes::install_github("FABLE-consortium/FABLEDownscalR")

if (!"package:FABLEDownscalR" %in% search()) {
  library(FABLEDownscalR)
}

library(dplyr)
library(here)
library(ggnewscale)
library(ggpattern)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)
library(tidyr)
library(stringr)


# ------------------------------------------------------------------
# 2. Read configuration
# ------------------------------------------------------------------

cfg <- fdr_read_config(
  here::here("config", config_file)
)

# Automatically create the date stamp used in output filenames
cfg$stamp <- format(Sys.Date(), "%y%m%d")

# Ensure reproducible MNL estimation
set.seed(cfg$seed)


# ------------------------------------------------------------------
# 3. Prepare country information
# ------------------------------------------------------------------

country <- countrycode::countrycode(
  cfg$country,
  origin = "iso3c",
  destination = "country.name"
)

border_sf <- rnaturalearth::ne_countries(
  country = country,
  scale = "medium",
  returnclass = "sf"
)

# (3) Load raw inputs (geojson + mapping + grid + FABLE)
inputs <- fdr_load_inputs(
  data_root        = cfg$data_root,
  country          = cfg$country,
  start_map_source = cfg$start_map_source,
  pathway          = cfg$pathway
)

# (4) Build land-cover change calibration table (DownscalR format)
luc <- lc_build_country_luc(
  LandCoverChange_df = inputs$spatial$landcoverchange,
  map_LUC            = inputs$mapping$map_LUC,
  Ts                 = 2015
)

# (5) a  Harmonise start map to match FABLE baseline totals
harm_starting <- fdr_harmonize_start_map(
  LandCoverStarting_df  = inputs$spatial$landcoverstarting,
  LandCoverInitial_df   = inputs$spatial$landcoverinitial,
  type          = "starting",
  map_LC        = inputs$mapping$map_LC,
  LC_targets    = inputs$LC_targets
)

# (5) b Harmonise initial map to match FABLE baseline totals
harm_initial <- fdr_harmonize_start_map(
  LandCoverStarting_df  = inputs$spatial$landcoverstarting,
  LandCoverInitial_df   = inputs$spatial$landcoverinitial,
  type = "initial",
  map_LC        = inputs$mapping$map_LC,
  LC_targets    = inputs$LC_targets
)

# (6) Build ns_map + rasterized ID layer (resolution controlled by YAML)
id <- fdr_build_id_maps(
  grid_sp         = inputs$grid_sp,      # e.g. Travel as sp/sf geometry for cells
  ns_map          = harm_starting$ns_map,
  pixel_res_m     = cfg$pixel_res_m
)

# (7) Build priors (X matrix etc.) + drop cells with NA covariates
priors <- fdr_build_priors(
  inputs         = inputs,
  start_map      = harm_initial$start_map_reproj,
  good_ns_only   = TRUE
)

# (8) Fit MNL + downscale
results <- fdr_run_downscaling(
  targets      = fdr_wrangle_fable_targets(inputs$FABLE_targets, min_year = 2020),
  country_luc  = luc$country_luc,
  priors       = priors,
  mnl_niter    = cfg$mnl_niter,
  mnl_nburn    = cfg$mnl_nburn,
  EF_LUC       = inputs$EF_LUC
)

# (9) Save outputs (consistent naming)
fdr_save_outputs(
  country     = cfg$country,
  tag         = fdr_make_tag(cfg),
  output_root = cfg$output_root,
  outputs     = list(
    start_map_reproj   = harm_starting$start_map_reproj,
    ns_map             = harm_starting$ns_map,
    rasterized_layer   = id$rasterized_layer,
    grid_sf            = id$grid_sf,
    X_long             = results$X_long,
    betas              = results$pred_coeff_long,
    country_start_areas= results$country_start_areas,
    downscaled_LUC     = results$downscaled_LUC,
    luc_hist           = luc$country_luc
  )
)

message("✅ Done: ", cfg$country, " (", cfg$pathway, ")")


