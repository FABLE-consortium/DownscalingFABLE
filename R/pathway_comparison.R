# ============================================================
# Compare DownscalR pathways
#
# Cumulative land-use gains and losses over a selected period
#
# Rows    = pathway
# Columns = land-cover class
#
# Negative values = cumulative loss
# Positive values = cumulative gain
#
# Plotting follows the same logic as:
#   fdr_plot_downscaled_LUC()
#
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

#devtools::load_all("../FABLEDownscalR")


# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

country <- "IND"

pathways <- c(
  "CurrentTrends",
  "NationalCommitments"#,
  #"GlobalSustainability"
)

start_map_source <- "HILDA"

start_year <- 2020
end_year   <- 2050


# Land-cover classes to display
LU <- c(
  "cropland",
  "pasture",
  "otherland",
  "newforest"
)


# Optional fixed scale.
# Leave NULL to calculate a common symmetric scale
# from all selected pathways.
limits <- NULL


pathway_labels <- c(
  CurrentTrends        = "Current Trends",
  NationalCommitments  = "National Commitments",
  GlobalSustainability = "Global Sustainability"
)


# ------------------------------------------------------------
# 2. Output folders
# ------------------------------------------------------------

country_dir <- here(country)

if (!dir.exists(country_dir)) {
  stop(
    "Country results folder not found:\n",
    country_dir
  )
}


figure_directory <- here(
  "Output",
  country
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Country border
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


# ------------------------------------------------------------
# 4. Helper: find latest run for one pathway
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
      "No saved DownscalR results found for:\n",
      "  country          = ", country, "\n",
      "  pathway          = ", pathway, "\n",
      "  start_map_source = ", start_map_source, "\n",
      "in:\n",
      country_dir
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
    basename(latest_file)
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
# 5. Helper: load one pathway
# ------------------------------------------------------------

load_pathway_results <- function(pathway) {
  
  run <- find_latest_run(
    country = country,
    pathway = pathway,
    start_map_source = start_map_source,
    country_dir = country_dir
  )
  
  
  run_tag <- run$run_tag
  
  
  # Downscaled LUC
  out_res <- readRDS(
    run$luc_file
  )
  
  
  # Raster containing ns_int
  rasterized_layer <- terra::rast(
    file.path(
      country_dir,
      paste0(
        run_tag,
        "_rasterized_layer.tif"
      )
    )
  )
  
  
  # Mapping between original id_c / ns and ns_int
  ns_map <- readRDS(
    file.path(
      country_dir,
      paste0(
        run_tag,
        "_ns_map.rds"
      )
    )
  )
  
  
  list(
    pathway = pathway,
    run_tag = run_tag,
    out_res = out_res,
    rasterized_layer = rasterized_layer,
    ns_map = ns_map
  )
}


# ------------------------------------------------------------
# 6. Load selected pathways
# ------------------------------------------------------------

results <- purrr::map(
  pathways,
  load_pathway_results
)

names(results) <- pathways


# ------------------------------------------------------------
# 7. Check grids are identical
# ------------------------------------------------------------

reference_raster <- results[[1]]$rasterized_layer
reference_ns_map <- results[[1]]$ns_map


if (length(results) > 1) {
  
  for (i in 2:length(results)) {
    
    same_grid <- terra::compareGeom(
      reference_raster,
      results[[i]]$rasterized_layer,
      stopOnError = FALSE
    )
    
    
    if (!same_grid) {
      
      stop(
        "Spatial grids differ between pathways.\n",
        "Pathway comparison requires the same grid."
      )
    }
  }
}


# ------------------------------------------------------------
# 8. Raster base
# ------------------------------------------------------------
#
# This follows fdr_plot_downscaled_LUC().
#
# rasterized_layer contains ns_int.
# Within the plotting function it is renamed to `ns`.
# ------------------------------------------------------------

df_pix <- terra::as.data.frame(
  reference_raster,
  xy = TRUE,
  na.rm = FALSE
)


names(df_pix)[3] <- "ns"


df_pix <- df_pix |>
  dplyr::filter(
    !is.na(ns)
  )


# ------------------------------------------------------------
# 9. Prepare cumulative LUC for one pathway
# ------------------------------------------------------------

prepare_pathway_LUC <- function(result) {
  
  pathway <- result$pathway
  out_res <- result$out_res
  ns_map  <- result$ns_map
  
  
  # ----------------------------------------------------------
  # Convert original ns / id_c to integer ns used by raster
  #
  # Exactly the same approach as fdr_plot_downscaled_LUC()
  # ----------------------------------------------------------
  
  out_int <- fdr_to_ns_int(
    out_res,
    ns_map
  )
  
  
  # ----------------------------------------------------------
  # Restrict to 2020-2050
  #
  # times = 2025 represents 2020-2025, etc.
  # ----------------------------------------------------------
  
  out_period <- out_int |>
    dplyr::filter(
      times > start_year,
      times <= end_year,
      lu.from != lu.to
    )
  
  
  # ----------------------------------------------------------
  # Gains
  #
  # Same logic as fdr_plot_downscaled_LUC():
  # land entering a land-cover class is positive
  # ----------------------------------------------------------
  
  gains <- out_period |>
    dplyr::group_by(
      lu.to,
      ns
    ) |>
    dplyr::summarise(
      gain = sum(
        value,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      lu = lu.to
    )
  
  
  # ----------------------------------------------------------
  # Losses
  #
  # Same logic as fdr_plot_downscaled_LUC():
  # land leaving a land-cover class is negative
  # ----------------------------------------------------------
  
  losses <- out_period |>
    dplyr::group_by(
      lu.from,
      ns
    ) |>
    dplyr::summarise(
      loss = -sum(
        value,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      lu = lu.from
    )
  
  
  # ----------------------------------------------------------
  # Combine gain + loss
  #
  # This gives cumulative net LUC by:
  #   cell x land-cover class
  # ----------------------------------------------------------
  
  inputs <- gains |>
    dplyr::select(
      ns,
      lu,
      gain
    ) |>
    dplyr::full_join(
      losses |>
        dplyr::select(
          ns,
          lu,
          loss
        ),
      by = c(
        "ns",
        "lu"
      )
    ) |>
    dplyr::mutate(
      gain = tidyr::replace_na(
        gain,
        0
      ),
      loss = tidyr::replace_na(
        loss,
        0
      ),
      value = gain + loss,
      pathway = pathway
    ) |>
    dplyr::rename(
      lu.to = lu
    )
  
  
  inputs
}


# ------------------------------------------------------------
# 10. Prepare all pathways
# ------------------------------------------------------------

inputs_all <- purrr::map_dfr(
  results,
  prepare_pathway_LUC
)


# ------------------------------------------------------------
# 11. Select land-cover classes
# ------------------------------------------------------------

inputs_all <- inputs_all |>
  dplyr::filter(
    lu.to %in% LU
  )


# ------------------------------------------------------------
# 12. Join results to raster pixels
# ------------------------------------------------------------
#
# Same direction of join as fdr_plot_downscaled_LUC().
# ------------------------------------------------------------

plot_df <- df_pix |>
  dplyr::left_join(
    inputs_all,
    by = "ns"
  ) |>
  dplyr::filter(
    !is.na(lu.to),
    !is.na(pathway)
  )


# ------------------------------------------------------------
# 13. Order facets
# ------------------------------------------------------------

lu_order <- c(
  "cropland",
  "newforest",
  "otherland",
  "pasture",
  "forest",
  "urban"
)


plot_df$lu.to <- factor(
  plot_df$lu.to,
  levels = lu_order
)


plot_df$pathway <- factor(
  plot_df$pathway,
  levels = pathways,
  labels = pathway_labels[pathways]
)


lu_labels <- c(
  cropland  = "Cropland",
  forest    = "Forest",
  newforest = "New\nforest",
  otherland = "Other\nland",
  pasture   = "Pasture",
  urban     = "Urban"
)


# ------------------------------------------------------------
# 14. Common colour limits
# ------------------------------------------------------------

if (is.null(limits)) {
  
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
}


message(
  "\nColour limits: ",
  round(limits[1], 2),
  " to ",
  round(limits[2], 2),
  " (1000 ha)"
)


# ------------------------------------------------------------
# 15. Plot
# ------------------------------------------------------------
#
# Same style and colour logic as
# fdr_plot_downscaled_LUC().
# ------------------------------------------------------------

p_LUC_comparison <- ggplot2::ggplot(
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
    low = "#FF0000",
    mid = "white",
    high = "#00B300",
    midpoint = 0,
    limits = limits,
    na.value = "grey90",
    name = "1000 ha"
  ) +
  
  ggplot2::coord_equal(
    expand = FALSE
  ) +
  
  ggplot2::facet_grid(
    pathway ~ lu.to,
    labeller = ggplot2::labeller(
      lu.to = lu_labels
    ),
    switch = "y"
  )+
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    strip.text = ggplot2::element_text(
      size = 10,
      face = "bold"
    ),
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(
      size = 9
    ),
    legend.text = ggplot2::element_text(
      size = 8
    ),
    plot.title = ggplot2::element_text(
      size = 12,
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10
    ),
    panel.spacing = grid::unit(
      0.15,
      "lines"
    )
  )

# ------------------------------------------------------------
# 16. Add border and white mask
# ------------------------------------------------------------
#
# Same approach as fdr_plot_downscaled_LUC().
# ------------------------------------------------------------

raster_crs <- terra::crs(
  reference_raster
)


border_use <- sf::st_transform(
  border_sf,
  crs = raster_crs
)


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


p_LUC_comparison <- p_LUC_comparison +
  
  ggplot2::geom_sf(
    data = outside_poly,
    fill = "white",
    color = NA,
    linewidth = 0
  ) +
  
  ggplot2::geom_sf(
    data = border_use,
    fill = NA,
    color = "grey60",
    linewidth = 0.3
  )


# ------------------------------------------------------------
# 17. Add title
# ------------------------------------------------------------

p_LUC_comparison <- p_LUC_comparison +
  
  ggplot2::labs(
    title = paste0(
      "Cumulative land-use gains and losses in ",
      country_name, " ", start_year, "-", end_year
    )
  )


print(
  p_LUC_comparison
)


# ------------------------------------------------------------
# 18. Save figure
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = here(
    "Output",
    country,
    paste0(
      "LUC_pathway_comparison_",
      start_year,
      "_",
      end_year,
      ".tiff"
    )
  ),
  plot = p_LUC_comparison,
  units = "in",
  height = 5.5,
  width = 9,
  dpi = 300
)


# ------------------------------------------------------------
# 19. Save comparison data
# ------------------------------------------------------------

saveRDS(
  inputs_all,
  here(
    "Output",
    country,
    paste0(
      "LUC_pathway_comparison_",
      start_year,
      "_",
      end_year,
      ".rds"
    )
  )
)


# ============================================================
# END
# ============================================================

