#### EF from mrorganic and downscled LUC for CHOICE D3.5
## Clara
## 20/04/2026

library(dplyr)
library(readr)
library(here)
library(readxl)

downscaled <- readRDS("C:/Users/Clara DOUZAL/Documents/Github Desktop/DownscalingFABLE/IND/260204_IND_CurrentTrends_HILDA_downscaled_LUC.rds")
ns_map     <- readRDS("C:/Users/Clara DOUZAL/Documents/Github Desktop/DownscalingFABLE/IND/260204_IND_CurrentTrends_HILDA_ns_map.rds")
EF         <- read_excel("C:/Users/Clara DOUZAL/Documents/Github Desktop/madrat/output/EF_Pools_transition_Ecoregion.xlsx")

EF_FABLE_cat <- EF %>% 
  mutate(from = ifelse(from == "grassland", "pasture",
                       ifelse(from %in% c("other", "othernatveg"), "otherland", from))) %>% 
  mutate(to = ifelse(to == "grassland", "pasture",
                       ifelse(to %in% c("other", "othernatveg"), "otherland", to))) %>% 
  mutate(to = ifelse(to == "forest", "newforest", to)) %>% 
  group_by(iso3, ECO_NAME, from, to) %>% 
  summarise(
    across(
      c(
        ef_soc,
        ef_biomass_above,
        ef_biomass_below,
        ef_total
      ),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  )

grid_table <- read_csv(here("Data", "global", "grid50_equal_area.csv")) %>% 
  left_join(ns_map) %>% 
  mutate(ns = id_c)


LUC <- downscaled %>% 
  left_join(grid_table) %>% 
  select(iso3, ECO_NAME, id_c, ns, times, lu.from, lu.to, value)

LUC_EF <- left_join(LUC, EF_FABLE_cat %>% 
                      rename(lu.from = from,
                             lu.to = to)) %>% 
  mutate(EM = value * (ef_biomass_above + ef_biomass_below))


### look at cropland abandonment
cropland_aband_EM <- LUC_EF %>% 
  filter(lu.from == "cropland",
         lu.to == "otherland") %>% 
  mutate(ef_biomass_CO2 = 3.67*(ef_biomass_above + ef_biomass_below)) %>% 
  group_by(ECO_NAME) %>% 
  mutate(average_EF = mean(ef_biomass_CO2)) %>% 
  ungroup() %>% 
  mutate(weighted_EF = weighted.mean(ef_biomass_CO2, value))

df <- cropland_aband_EM %>% 
  filter(times == 2025) %>% 
  select(ECO_NAME, lu.from, lu.to, average_EF, weighted_EF) %>% 
  unique()
summary(df)

df_transition <- cropland_aband_EM %>% 
  select(ECO_NAME, times, lu.from, lu.to, value, average_EF) %>%
  mutate(emission = value * average_EF) %>% 
  group_by(times) %>% 
  summarise(emission = sum(emission, na.rm = T))
