# ============================================================
# Compare cumulative spatial GHG emissions across pathways
#
# Rows = pathway
#
# Uses GHG_biomass already stored in downscaled_LUC.
# No EF calculation is repeated.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(here)
library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sf)
library(countrycode)
library(rnaturalearth)
library(FABLEDownscalR)



# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

country <- "IND"

pathways <- c(
  "CurrentTrends",
  "NationalCommitments"
)

start_map_source <- "HILDA"

start_year <- 2020
end_year   <- 2050


pathway_labels <- c(
  CurrentTrends        = "Current Trends",
  NationalCommitments  = "National Commitments",
  GlobalSustainability = "Global Sustainability"
)


# ------------------------------------------------------------
# 2. Locate country output folder
# ------------------------------------------------------------

country_dir_candidates <- c(
  here("Output", country),
  here(country)
)

country_dir <- country_dir_candidates[
  dir.exists(country_dir_candidates)
][1]

if (is.na(country_dir)) {
  stop(
    "Could not find country output folder.\n\nChecked:\n",
    paste(
      country_dir_candidates,
      collapse = "\n"
    )
  )
}


comparison_dir <- here(
  "Output",
  country,
  "pathway_comparison"
)

dir.create(
  comparison_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Helper: find latest run
# ------------------------------------------------------------

find_latest_run <- function(
    country,
    pathway,
    start_map_source,
    country_dir
) {
  
  expected_suffix <- paste0(
    "_",
    country,
    "_",
    pathway,
    "_",
    start_map_source,
    "_downscaled_LUC.rds"
  )
  
  candidate_files <- list.files(
    country_dir,
    pattern = paste0(
      expected_suffix,
      "$"
    ),
    full.names = TRUE
  )
  
  if (length(candidate_files) == 0) {
    stop(
      "No saved DownscalR result found for:\n",
      "  country          = ", country, "\n",
      "  pathway          = ", pathway, "\n",
      "  start_map_source = ", start_map_source
    )
  }
  
  candidate_info <- file.info(
    candidate_files
  )
  
  latest_file <- candidate_files[
    which.max(
      candidate_info$mtime
    )
  ]
  
  run_tag <- sub(
    "_downscaled_LUC\\.rds$",
    "",
    basename(
      latest_file
    )
  )
  
  message(
    "\nSelected run for ",
    pathway,
    ":\n  ",
    run_tag
  )
  
  list(
    run_tag = run_tag,
    luc_file = latest_file
  )
}


# ------------------------------------------------------------
# 4. Load one pathway
# ------------------------------------------------------------

load_pathway <- function(pathway) {
  
  run <- find_latest_run(
    country = country,
    pathway = pathway,
    start_map_source = start_map_source,
    country_dir = country_dir
  )
  
  run_tag <- run$run_tag
  
  
  # Downscaled result
  out_res <- readRDS(
    run$luc_file
  )
  
  
  if (!"GHG_biomass" %in% names(out_res)) {
    stop(
      "`GHG_biomass` is missing from:\n",
      run$luc_file,
      "\n\nAvailable columns:\n",
      paste(
        names(out_res),
        collapse = ", "
      )
    )
  }
  
  
  # ns map
  ns_map <- readRDS(
    file.path(
      country_dir,
      paste0(
        run_tag,
        "_ns_map.rds"
      )
    )
  )
  
  
  # Rasterized layer
  rasterized_layer <- terra::rast(
    file.path(
      country_dir,
      paste0(
        run_tag,
        "_rasterized_layer.tif"
      )
    )
  )
  
  
  list(
    pathway = pathway,
    out_res = out_res,
    ns_map = ns_map,
    rasterized_layer = rasterized_layer
  )
}


# ------------------------------------------------------------
# 5. Load pathways
# ------------------------------------------------------------

results <- purrr::map(
  pathways,
  load_pathway
)

names(results) <- pathways


# ------------------------------------------------------------
# 6. Reference raster
# ------------------------------------------------------------

reference_raster <- results[[1]]$rasterized_layer


# Check same geometry
if (length(results) > 1) {
  
  for (i in 2:length(results)) {
    
    if (
      !terra::compareGeom(
        reference_raster,
        results[[i]]$rasterized_layer,
        stopOnError = FALSE
      )
    ) {
      stop(
        "The pathways do not use the same spatial grid."
      )
    }
  }
}


# ------------------------------------------------------------
# 7. Prepare cumulative GHG for one pathway
# ------------------------------------------------------------
#
# This follows the same core logic as
# fdr_plot_downscaled_GHG_cum():
#
#   out_res already contains GHG_biomass
#   ns is converted to ns_int using ns_map
#   GHG is summed cumulatively over time
#
# ------------------------------------------------------------

prepare_cumulative_GHG <- function(result) {
  
  out_res <- result$out_res
  ns_map  <- result$ns_map
  
  
  # Convert model ns to raster ns_int
  out_int <- fdr_to_ns_int(
    out_res,
    ns_map
  )
  
  
  # Restrict to 2020-2050
  out_period <- out_int |>
    filter(
      times > start_year,
      times <= end_year
    )
  
  
  # Cumulative emissions across all transitions
  GHG_cumulative <- out_period |>
    group_by(
      ns
    ) |>
    summarise(
      value = sum(
        GHG_biomass,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    mutate(
      pathway = result$pathway
    )
  
  
  GHG_cumulative
}


# ------------------------------------------------------------
# 8. Combine pathways
# ------------------------------------------------------------

GHG_all <- purrr::map_dfr(
  results,
  prepare_cumulative_GHG
)


# ------------------------------------------------------------
# 9. Raster dataframe
# ------------------------------------------------------------
#
# Same plotting structure used by the package functions.
# ------------------------------------------------------------

df_pix <- terra::as.data.frame(
  reference_raster,
  xy = TRUE,
  na.rm = FALSE
)


names(df_pix)[3] <- "ns"


df_pix <- df_pix |>
  filter(
    !is.na(ns)
  )


# ------------------------------------------------------------
# 10. Join GHG results to raster pixels
# ------------------------------------------------------------

plot_df <- df_pix |>
  left_join(
    GHG_all,
    by = "ns"
  ) |>
  filter(
    !is.na(pathway)
  )


# ------------------------------------------------------------
# 11. Pathway labels
# ------------------------------------------------------------

plot_df <- plot_df |>
  mutate(
    pathway = factor(
      pathway,
      levels = pathways,
      labels = pathway_labels[
        pathways
      ]
    )
  )


# ------------------------------------------------------------
# 12. Country border
# ------------------------------------------------------------

country_name <- countrycode::countrycode(
  country,
  origin = "iso3c",
  destination = "country.name"
)


border_sf <- rnaturalearth::ne_countries(
  country = country_name,
  scale = "medium",
  returnclass = "sf"
)


raster_crs <- terra::crs(
  reference_raster
)


border_use <- sf::st_transform(
  border_sf,
  crs = raster_crs
)


# ------------------------------------------------------------
# 13. White mask outside country
# ------------------------------------------------------------

bbox_poly <- sf::st_as_sfc(
  sf::st_bbox(
    border_use
  )
)


outside_poly <- sf::st_difference(
  bbox_poly,
  sf::st_union(
    border_use
  )
)


# ------------------------------------------------------------
# 14. Common symmetric colour scale
# ------------------------------------------------------------

max_abs <- max(
  abs(
    plot_df$value
  ),
  na.rm = TRUE
)


limits <- c(
  -max_abs,
  max_abs
)


# ------------------------------------------------------------
# 15. Plot
# ------------------------------------------------------------

p_GHG_comparison <- ggplot2::ggplot(
  plot_df
) +
  
  ggplot2::geom_raster(
    ggplot2::aes(
      x = x,
      y = y,
      fill = value
    )
  ) +
  
  ggplot2::scale_fill_gradient2(
    low = "#00B300",
    mid = "white",
    high = "#FF0000",
    midpoint = 0,
    limits = limits,
    na.value = "grey90",
    name = expression(
      "Cumulative Mt CO"[2]
    )
  ) +
  
  ggplot2::coord_equal(
    expand = FALSE
  ) +
  
  ggplot2::facet_grid(
    pathway ~ .,
    switch = "y"
  ) +
  
  ggplot2::geom_sf(
    data = outside_poly,
    fill = "white",
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ggplot2::geom_sf(
    data = border_use,
    fill = NA,
    colour = "grey60",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  
  ggplot2::labs(
    title = paste0(
      "Cumulative land-use-change CO2 emissions in ",
      country_name
    ),
    subtitle = paste0(
      start_year,
      "-",
      end_year,
      " | Green = sequestration, red = emissions"
    )
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    strip.placement = "outside",
    
    strip.text.y.left =
      ggplot2::element_text(
        angle = 0,
        size = 10,
        face = "bold"
      ),
    
    legend.position = "bottom",
    
    plot.title =
      ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 9
      )
  )


print(
  p_GHG_comparison
)


# ------------------------------------------------------------
# 16. Save
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    comparison_dir,
    paste0(
      country,
      "_",
      start_year,
      "_",
      end_year,
      "_cumulative_GHG_pathway_comparison.tiff"
    )
  ),
  plot = p_GHG_comparison,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)


# ============================================================
# END
# ============================================================