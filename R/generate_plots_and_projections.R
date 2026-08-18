## Generate figures and computations on GHG emissions from Downscaled Results

# ============================================================
# 1. Set the envrionment to load outputs from run_country.R
# ============================================================

library(here)
library(terra)
library(dplyr)
library(writexl)

devtools::load_all("../FABLEDownscalR")

cfg <- fdr_read_config(here::here("config", "BRA.yml"))
cfg$stamp <- format(Sys.Date(), "%y%m%d")
set.seed(cfg$seed)

date <- "260810"
path <- paste(date, cfg$country, cfg$pathway, cfg$start_map_source, sep = "_")

# Setting directory for storing figues maps and excel
figure_directory <- here("output", cfg$country)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
print(figure_directory)

country <- countrycode::countrycode(cfg$country, origin = 'iso3c', destination = 'country.name')
border_sf <- rnaturalearth::ne_countries(country = country, scale = "medium", returnclass = "sf")

# Loading the already generated downscaled LUC and necessry data to creat maps
out.res          <- readRDS(here(cfg$country, paste(path, "downscaled_LUC.rds", sep = "_")))
rasterized_layer <- rast(here(cfg$country, paste(path, "rasterized_layer.tif", sep = "_")))
ns_map           <- readRDS(here(cfg$country, paste(path, "ns_map.rds", sep = "_")))
luc              <- readRDS(here(cfg$country, paste(path, "luc_hist.rds", sep = "_")))

# ============================================================
# 2. Generate LU and LUC maps
# ============================================================

# ------------------------------------------------------------
# 2020 Dominant land cover per cell
# /!\ Takes time to run, only generate if needed /!\
# ------------------------------------------------------------

# p_main_LU <- fdr_plot_downscaled_LU_one(
#   out_res          = out.res, 
#   rasterized_layer = rasterized_layer,
#   ns_map           = ns_map,
#   LU               = c("cropland", "forest", "pasture", "otherland"),
#   year             = "2020"
# )
# 
# print(p_main_LU)
# 
# ggplot2::ggsave(
#   filename = here("Output", cfg$country, paste0("Dominant_Land_Cover_", cfg$pathway, ".tiff")),
#   plot = p_main_LU,
#   units = "in",
#   height = 7.45, width = 7.45, dpi = 300)

# ------------------------------------------------------------
#  Land Use Map for each land cover type over 2020-2050
# (1000 ha)
# ------------------------------------------------------------

p_LU <- fdr_plot_downscaled_LU(
  out_res          = out.res, 
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf
)

print(p_LU)

ggplot2::ggsave(
  filename = here("Output", cfg$country, paste0("LU_", cfg$pathway, ".tiff")),
  plot = p_LU,
  units = "in",
  height = 6, width = 9, dpi = 600)

# ------------------------------------------------------------
#  Land Use Change Map for each land cover type over 2020-2050
# (1000 ha)
# ------------------------------------------------------------

p_LUC <- fdr_plot_downscaled_LUC(
  out_res          = out.res,
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf,
  LU               = c("cropland", "pasture", "otherland", "forest", "newforest"),
  limits           = c(-50,50)
)

print(p_LUC)

ggplot2::ggsave(
  filename = here("Output", cfg$country, paste0("LUC_", cfg$pathway, ".tiff")),
  plot = p_LUC,
  units = "in",
  height = 6.5, width = 3.5, dpi = 300)

# ------------------------------------------------------------
#  Downscaled CO2 emissions per land cover type over 2020-2050
# (Mt CO2/5 year)
# ------------------------------------------------------------

GHG <- fdr_plot_downscaled_GHG(
  out_res          = out.res,
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf,
  limits           = c(-8, 8)
)

print(GHG$plot)

ggplot2::ggsave(
  filename = here("Output", cfg$country, paste0("LUC_CO2_emissions_", cfg$pathway, ".tiff")), 
  plot = GHG$plot, 
  units = "in", 
  height = 3, width = 8.5, dpi = 300)

#total LUC CO2 emissions per 5 year period
TotalCO2_5year <- GHG$GHG %>% 
  group_by(times) %>% 
  summarize(GHG_biomass = sum(GHG_biomass))

#per transition CO2emissions per 5 year period
Transition_CO2_em <-  out.res %>% 
  group_by(times, lu.from, lu.to) %>% 
  summarise(CO2 = sum(GHG_biomass))

# average emission factor per land transition in 2020 using spatially explicit data fro ESA-CCI
average_EF <- out.res %>% 
  filter(times == 2020) %>% 
  select(-value) %>% 
  left_join(luc, by = join_by(ns, lu.to, lu.from)) %>% 
  group_by(times, lu.to, lu.from) %>% 
  summarise(EF_weighted = weighted.mean(ef_biomass * 3.6667, value, na.rm = T)) %>% 
  filter(times == 2020) %>% 
  filter(!is.na(EF_weighted))

#full grid emission factor per land transition (NASA 2010)
distribution_EF <- out.res %>% 
  mutate(ef_biomass = ef_biomass *3.6667) %>% 
  filter(times == 2020) %>% 
  select(ns, lu.to, lu.from, ef_biomass) %>% 
  tidyr::pivot_wider(names_from = c("lu.to", "lu.from"), values_from = ef_biomass)

# Standard set of stats (summary() sometimes adds "NA's" if there are missing values)
stat_names <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "NA's")

get_stats <- function(x) {
  if (is.numeric(x)) {
    s <- summary(x)
    full <- setNames(rep(NA_real_, length(stat_names)), stat_names)
    full[names(s)] <- as.numeric(s)
    full
  } else {
    setNames(rep(NA, length(stat_names)), stat_names)
  }
}

summary_mat <- sapply(distribution_EF, get_stats)          # variables as columns, stats as rows
summary_df  <- as.data.frame(summary_mat)
summary_df  <- cbind(Statistic = rownames(summary_df), summary_df)
rownames(summary_df) <- NULL

write_xlsx(
  list(
    TotalCO2_5year    = TotalCO2_5year,
    Transition_CO2_em = Transition_CO2_em,
    average_EF        = average_EF,
    StatDesc_EF       = summary_df
  ),
  path = here("Output", cfg$country, paste0("LUC_CO2_emissions_", cfg$pathway, ".xlsx"))
)


# ------------------------------------------------------------
#  Downscaled CO2 emissions per land cover type cumulative over 2020-2050
# (Mt CO2/5 year)
# ------------------------------------------------------------

GHG_cum <- fdr_plot_downscaled_GHG_cum(
  out_res          = out.res,
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf,
  year             = 2050,
  limits           = c(-25, 30)
  
)
print(GHG_cum$plot)

GHG_cum_newforest <- fdr_plot_downscaled_GHG_cum(
  out_res          = out.res,
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf,
  LU               = c("newforest"),
  year             = 2050
  
)
print(GHG_cum_newforest$plot)
GHG_cum_cropland <- fdr_plot_downscaled_GHG_cum(
  out_res          = out.res,
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf,
  LU               = c("cropland"),
  year             = 2050
  
)
print(GHG_cum_cropland$plot)
GHG_cum_otherland <- fdr_plot_downscaled_GHG_cum(
  out_res          = out.res,
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf,
  LU               = c("otherland"),
  year             = 2050
  
)
print(GHG_cum_otherland$plot)

View(GHG_cum_otherland$cum_GHG)

ggplot2::ggsave(
  filename = here("Output", cfg$country, paste0("Cummulative_LUC_CO2_emissions_", cfg$pathway, ".tiff")), 
  plot = GHG_cum$plot, 
  units = "in", 
  height = 5, width = 6, dpi = 300)

GHG_transition  <-  fdr_plot_downscaled_GHG_transition(
  out_res          = out.res,
  rasterized_layer = rasterized_layer,
  ns_map           = ns_map,
  border_sf        = border_sf,
  LU = c("cropland", "newforest", "otherland")
)
print(GHG_transition$plot)

