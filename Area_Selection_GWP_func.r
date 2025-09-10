library(terra)
library(stringr)
library(sf)
library(smoothr)

#' Load and name GWP raster files
#' @param gwp_folder Folder containing GWP .tif files
#' @param indices Indices of files to select
#' @param crs Coordinate reference system string
#' @return Named list of raster objects
load_gwp_rasters <- function(gwp_folder, indices, crs) {
  list.files(gwp_folder, pattern = "\\.tif$", full.names = TRUE) |>
    (\(gwp_list) {
      gwp_list_filter <- gwp_list[indices]
      gwp_basename <- basename(gwp_list_filter)
      gwp_names <- paste0(str_extract(gwp_basename, "\\d{4}"), "_GWP")
      gwp_rasters <- lapply(gwp_list_filter, rast)
      names(gwp_rasters) <- gwp_names
      lapply(gwp_rasters, function(x) project(x, crs))
      gwp_rasters
    })()
}

#' Filter raster pixels by value range
#' @param r Raster object
#' @param lower Lower threshold
#' @param upper Upper threshold
#' @return Raster with filtered values (others set to NA)
filter_pixels <- function(r, lower = 50, upper = 150){
  r |> (\(x) ifel(x >= lower & x <= upper, x, NA))()
}

#' Clip rasters to a vector boundary, aligning CRS if needed
#' @param rasters List of raster objects
#' @param clip_vect Vector object to clip to
#' @return List of clipped rasters
clip_rasters <- function(rasters, clip_vect) {
  vect_crs <- crs(clip_vect)
  rasters |> lapply(\(x) {
    if (crs(x) != vect_crs) {
      x <- project(x, vect_crs)
    }
    crop(x, clip_vect)
  })
}

#' Calculate number of non-NA neighbors for each pixel
#' @param masks_clip List of masked rasters
#' @return List of rasters with neighbor counts
calculate_neighbors <- function(masks_clip) {
  w <- matrix(1, 5, 5)
  w[3, 3] <- 0
  masks_clip |> lapply(\(x) focal(!is.na(x), w = w, fun = "sum", na.rm = TRUE, na.policy = "omit", fillvalue = 0))
}

#' Mask out pixels with insufficient neighbors
#' @param masks_clip List of masked rasters
#' @param neighbors List of neighbor count rasters
#' @param threshold Minimum neighbor count
#' @return List of filtered rasters
apply_neighbor_filter <- function(masks_clip, neighbors, threshold = 15) {
  mapply(function(x, nbh) ifel(nbh > threshold, x, NA), masks_clip, neighbors, SIMPLIFY = FALSE)
}

#' Convert rasters to polygons
#' @param filtered_rasters List of filtered rasters
#' @return List of polygon vector objects
rasters_to_polygons <- function(filtered_rasters) {
  filtered_rasters |> lapply(as.polygons)
}

#' Merge, buffer, smooth polygons, and create convex hulls
#' @param polygons_list List of polygons
#' @param crs Coordinate reference system string
#' @param buffer_dist Buffer distance (default 500)
#' @return List containing smoothed vector object and convex hulls
merge_smooth_and_hull <- function(polygons_list, crs, buffer_dist = 500) {
  smoothed <- polygons_list |>
    vect() |>
    aggregate(dissolve = TRUE) |>
    fillHoles() |>
    project(crs) |>
    (\(x) buffer(x, buffer_dist))() |>
    aggregate(dissolve = TRUE) |>
    (\(x) smoothr::smooth(x, method = "ksmooth"))() |>
    project(crs) |>
    makeValid()
  
  convex_hulls <- polygons_list |> lapply(\(x) {
    hull <- convHull(x)
    hull_filled <- fillHoles(hull)
    makeValid(hull_filled)
  })
  
  list(smoothed = smoothed, convex_hulls = convex_hulls)
}

#' Write output shapefiles
#' @param smooth_vector_proj Smoothed AOI vector
#' @param poly_gwp_areas List of polygons
#' @param out_dir_vect Output directory for shapefiles
#' @param lower Lower pixel threshold
#' @param upper Upper pixel threshold
write_outputs <- function(smooth_vector_proj, poly_gwp_areas, out_dir_vect, lower = 50, upper = 150) {
  # Replace the smoothed vector with its convex hull
  smooth_vector_proj <- convHull(smooth_vector_proj)
  smooth_vector_proj <- makeValid(smooth_vector_proj)
  
  # Add lower and upper as attributes to the convex hull vector
  smooth_vector_proj$lower <- lower
  smooth_vector_proj$upper <- upper
  out_name <- paste0(out_dir_vect, "AOI_Smooth_", as.character(lower), "_", as.character(upper), ".shp")
  writeVector(smooth_vector_proj, out_name, overwrite = TRUE)
  
  for (i in seq_along(poly_gwp_areas)) { 
    out_name <- paste0(out_dir_vect, "PerYear/AOI_", names(poly_gwp_areas)[i], as.character(lower), "_", as.character(upper), ".shp")
    writeVector(poly_gwp_areas[[i]], out_name, overwrite = TRUE)
  }
}

#' Run the full GWP area selection workflow
#' @param gwp_folder Folder with GWP .tif files
#' @param aral_delta_shp Path to Aral delta shapefile
#' @param out_dir_vect Output directory for shapefiles
#' @param indices Indices of files to use (default 18:22)
#' @param crs CRS string (default "EPSG:4326")
#' @param lower Lower pixel threshold (default 50)
#' @param upper Upper pixel threshold (default 150)
#' @param neighbor_threshold Minimum neighbors (default 15)
#' @param buffer_dist Buffer distance (default 500)
run_area_selection_gwp <- function(
  gwp_folder,
  aral_delta_shp,
  out_dir_vect,
  indices = 18:22,
  crs = "EPSG:4326",
  lower = 60,
  upper = 199,
  neighbor_threshold = 20,
  buffer_dist = 500
) {
  dir.create(out_dir_vect, showWarnings = FALSE, recursive = TRUE)
  aral_delta <- vect(aral_delta_shp)
  gwp_rasters <- load_gwp_rasters(gwp_folder, indices, crs)
  gwp_masks <- gwp_rasters |> lapply(filter_pixels, lower = lower, upper = upper)
  gwp_masks_clip <- gwp_masks |> clip_rasters(aral_delta)
  gwp_neighbors <- gwp_masks_clip |> calculate_neighbors()
  gwp_filtered <- apply_neighbor_filter(gwp_masks_clip, gwp_neighbors, neighbor_threshold)
  poly_gwp_areas <- gwp_filtered |> rasters_to_polygons()
  poly_gwp_areas_clip <- poly_gwp_areas |> clip_rasters(aral_delta)
  merged_results <- merge_smooth_and_hull(poly_gwp_areas_clip, crs, buffer_dist)
  write_outputs(merged_results$smoothed, merged_results$convex_hulls, out_dir_vect, lower = lower, upper = upper)
}

# Usage:
run_area_selection_gwp(
  gwp_folder = "G:/Coen/Projects/Aral_Sea/Data/GWP",
  aral_delta_shp = "G:/Coen/Projects/Aral_Sea/Data/locations/Aral sea delta.shp",
  out_dir_vect = "G:/Coen/Projects/Aral_Sea/Data/locations/GWP/Output_func/",
)
