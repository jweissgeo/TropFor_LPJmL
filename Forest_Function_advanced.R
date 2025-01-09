### Funktion zur Schwellenwertberechnung ###
find_threshold <- function(target_percent_loss, forest_area, expansion_resampled, max_iterations = 100, tolerance = 0.02) {
  # Startwerte für die Iteration
  lower_bound <- 0
  upper_bound <- maxValue(expansion_resampled)
  iteration <- 0
  best_index <- NA
  achieved_percent_loss <- NA
  
  # Hilfsfunktion zur Flächenberechnung
  calc_total_area <- function(area_raster) {
    coords <- coordinates(area_raster)
    cell_areas <- calc_cellarea(coords[, 2], return_unit = "km2")
    total_area <- sum(values(area_raster) * cell_areas, na.rm = TRUE)
    return(total_area)
  }
  
  # Gesamtwaldfläche berechnen
  total_forest_area <- calc_total_area(forest_area)
  
  while (iteration < max_iterations) {
    current_index <- (lower_bound + upper_bound) / 2
    expansion_mask <- expansion_resampled > current_index
    
    # Waldflächen innerhalb der Expansion Mask berechnen
    forest_area_expansion <- forest_area * expansion_mask
    total_forest_expansion <- calc_total_area(forest_area_expansion)
    achieved_percent_loss <- (total_forest_expansion / total_forest_area) * 100
    
    # Debugging-Ausgaben
    cat("Iteration:", iteration, 
        "Lower Bound:", lower_bound, 
        "Upper Bound:", upper_bound, 
        "Current Index:", current_index, 
        "Achieved Percent Loss:", achieved_percent_loss, 
        "\n")
    
    # Überprüfung der Toleranz
    if (abs(achieved_percent_loss - target_percent_loss) < tolerance) {
      best_index <- current_index
      break
    }
    
    if (achieved_percent_loss < target_percent_loss) {
      upper_bound <- current_index
    } else {
      lower_bound <- current_index
    }
    
    iteration <- iteration + 1
  }
  
  if (!is.na(best_index)) {
    return(list(index = best_index, percent_loss = achieved_percent_loss))
  } else {
    stop("Kein Schwellenwert gefunden, maximale Iterationen erreicht.")
  }
}