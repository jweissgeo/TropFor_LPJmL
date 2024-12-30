library(raster) # Alternativ: library(terra)

# Beispiel: RasterStack mit den benötigten Variablen einlesen
#Raster mit readio einlesen
raster_stack <- stack("Klima Raster", "CFT Raster", "FPC Raster", "vegC Raster" "etc")

###Beispiel zum Schreiben einer Funktion. Achtung Werte und Konditionen sind noch NICHT geprüft

### Biomass = 2*VegC

# Funktion zur Klassifikation eines Pixels
classify_pixel <- function(values) {
  Urban <- values[1]
  Pop <- values[2]
  Pastoral <- values[3]
  Cropland <- values[4]
  Flooded <- values[5]
  Tropical_trees <- values[6]
  Temperate_trees <- values[7]
  Grass <- values[8]
  MeanTemp <- values[9]
  Precip <- values[10]
  PET <- values[11]
  Biomass <- values[12]
  
  # Entscheidungslogik
  if (Urban > 0.2 || Pop > 500) {
    return(1)  # Urban
  } else if (Pastoral > 0.5) {
    return(2)  # Pasture
  } else if (Cropland > 0.5) {
    return(3)  # Cropland
  } else if (Pastoral + Cropland > 0.5) {
    if (Pastoral > Cropland) {
      return(2)  # Pasture
    } else {
      return(3)  # Cropland
    }
  } else if (Flooded > 0.2) {
    return(4)  # Flooded
  } else if (Tropical_trees > 0.8) {
    if (MeanTemp < 24 && PET / Precip < 1) {
      if (Biomass > 40) {
        if (Precip > 4000 || PET / Precip < 0.5) {
          return(5)  # Closed Rainforest
        } else if (Precip > 2000 || PET / Precip < 1) {
          return(6)  # Closed Moist Forest
        } else {
          return(7)  # Open Moist Forest
        }
      } else {
        return(8)  # Open Dry Forest
      }
    } else {
      return(9)  # Shrubland
    }
  } else if (Tropical_trees > 0.4) {
    return(9)  # Shrubland
  } else if (Temperate_trees > 0.9) {
    return(10)  # Temperate Forest
  } else if (Temperate_trees > 0.1) {
    return(9)  # Shrubland
  } else if (Grass > 0.4) {
    return(11)  # Savannah/Grassland
  } else {
    return(12)  # Alpine or Desert
  }
}

# Anwendung der Klassifikation auf den RasterStack
classified_raster <- calc(raster_stack, classify_pixel)

# Optional: Klassifikation als Kategorien speichern
classified_raster <- ratify(classified_raster)
levels(classified_raster) <- data.frame(
  ID = 1:12,
  Class = c("Urban", "Pasture", "Cropland", "Flooded", 
            "Closed Rainforest", "Closed Moist Forest", 
            "Open Moist Forest", "Open Dry Forest", 
            "Shrubland", "Temperate Forest", 
            "Savannah/Grassland", "Alpine or Desert")
)

# Ergebnis speichern
writeRaster(classified_raster, "classified_biome.tif", overwrite = TRUE)