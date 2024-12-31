library(raster)
library(lpjmlkit)

### local path setzen, andere immer auskommentieren ###
#local_path <- "C:/Users/philipp/Documents/TropFor_LPJmL_LokalData/" #philipp
local_path <- "/Users/epigo/Documents/LPJmL_Lokal/" #Julius
#local_path <- "C:/Dokumente/Umweltsysteme/integrierte_modellierung/"  # Mareike

# Grid-Datei laden
lpjml_grid <- read_io(paste0(local_path, "CRUDATA_TS3_23/grid.bin"))

# MAT-Daten einlesen
clim_in_MAT <- read_io(
  paste0(local_path, "CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.clm"),
  name = "CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.clm",
  descr = "LPJCLIM",
  firstcell = 0,
  ncell = 67420,
  firstyear = 1901,
  nyear = 114,
  scalar = 1,
  nbands = 12,
  nstep = 1,
  timestep = 1,
  cellsize_lon = 0.5,
  cellsize_lat = 0.5,
  datatype = "short",
  order = "cellyear",
  endian = "little",
  version = 2,
  subset = list(year = 105)
)

# In Lon-Lat transformieren
climMAT_lonlat <- transform(clim_in_MAT, to = "lon_lat")

# RasterStack für MAT erstellen
raster_stackMAT <- stack()
for (band in 1:12) {
  band_data <- subset(climMAT_lonlat, band = band)
  band_raster <- as_raster(band_data)
  raster_stackMAT <- addLayer(raster_stackMAT, band_raster)
}

# RasterStack speichern
writeRaster(raster_stackMAT, filename = "modified_clim_stack_mat.tif", format = "GTiff", overwrite = TRUE)

# Plotten zur Überprüfung
spplot(raster_stackMAT[[12]], main = "Monthly Average Temperature Dezember") # Dezember
spplot(raster_stackMAT[[6]], main = "Monthly Average Temperature Juni") # Juni

# Durchschnittliche globale Temperatur berechnen und plotten
mean_global_temperature <- mean(raster_stackMAT)
spplot(mean_global_temperature, main = "Global MAT [°C] for 2005")

# Precip-Daten einlesen
clim_in_Precip <- read_io(
  paste0(local_path, "CRUDATA_TS3_23/gpcc_v7_cruts3_23_precip_1901_2013.clm"),
  name = "CRUDATA_TS3_23/gpcc_v7_cruts3_23_precip_1901_2013.clm",
  descr = "LPJCLIM",
  firstcell = 0,
  ncell = 67420,
  firstyear = 1901,
  nyear = 113,
  scalar = 1,
  nbands = 12,
  nstep = 1,
  timestep = 1,
  cellsize_lon = 0.5,
  cellsize_lat = 0.5,
  datatype = "short",
  order = "cellyear",
  endian = "little",
  version = 2,
  subset = list(year = 105)
)

# In Lon-Lat transformieren
climPrecip_lonlat <- transform(clim_in_Precip, to = "lon_lat")

# RasterStack für Precip erstellen
raster_stackPrecip <- stack()
for (band in 1:12) {
  band_data <- subset(climPrecip_lonlat, band = band)
  band_raster <- as_raster(band_data)
  raster_stackPrecip <- addLayer(raster_stackPrecip, band_raster)
}

# RasterStack speichern
writeRaster(raster_stackPrecip, filename = "modified_clim_stack_precip.tif", format = "GTiff", overwrite = TRUE)

# Plotten zur Überprüfung
spplot(raster_stackPrecip[[12]], main = "Monthly Precipitation Dezember") # Dezember
spplot(raster_stackPrecip[[6]], main = "Monthly Precipitation Juni") # Juni

# Durchschnittliche globale Niederschläge berechnen und plotten
global_precip <- sum(raster_stackPrecip)
spplot(global_precip, main = "Global Precipitation [mm] for 2005")


### Variable (PET)
read_meta(paste0(local_path, "Model_Output/FPC_test/mpet.bin.json"))                    # Metadaten anzeigen
clim_out_PET <- read_io(
  paste0(local_path, "Model_output/FPC_test/mpet.bin.json"),
  firstcell = 0,
  ncell = 67420,
  scalar = 1,
  nbands = 12,
  nstep = 1,
  timestep = 1,
  cellsize_lon = 0.5,
  cellsize_lat = 0.5,
)

PET_lonlat <- transform(clim_out_PET, to = "lon_lat")



raster_stackPET <- stack()
for (band in 1:12) {
  band_data <- subset(PET_lonlat, band = band)
  band_raster <- as_raster(band_data)
  raster_stackPET <- addLayer(raster_stackPET, band_raster)
}

print(raster_stackPET)

# RasterStack speichern
writeRaster(raster_stackPET, filename = "modified_stack_PET.tif", format = "GTiff", overwrite = TRUE)

# Plotten zur Überprüfung
spplot(raster_stackPET[[12]], main = " PET Dezember [mm]") # Dezember
spplot(raster_stackPET[[6]], main = " PET Juni[mm/]") # Juni

# Durchschnittliche globale Niederschläge berechnen und plotten
global_PET <- sum(raster_stackPET)
spplot(global_PET, main = "global PET 2005 [mm]")
# Margen reduzieren (z. B. von c(5, 4, 4, 2) auf kleinere Werte)

# Plot neu erstellen
plot(global_PET, main = "Global annual PET for 2005 [mm/Jahr]")


#Variable PET /precip

# Beispiel-Daten
Budyko <- global_PET / global_precip
# Wertebereiche und Farben (invertierte Reihenfolge, Pink für >3.0)
breaks <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.7, 0.8, 0.9, 1.0, 1.1, 1.4, 1.5, 1.75, 2.0, 3.0, Inf)
colors <- rev(c("pink", "#FF0000", "#FF4500", "#FF6347", "#FFA500", "#FFD700", 
                "#ADFF2F", "#7FFF00", "#32CD32", "#228B22", "#00FF7F", 
                "#00FA9A", "#00CED1", "#4682B4", "#5F9EA0", "#0000FF"))

# Plot erstellen
par(mar = c(5, 4, 4, 5))  # Margen anpassen (Platz für die Legende rechts)
plot(Budyko,
     col = colors,                  # Farbskala
     breaks = breaks,               # Wertebereiche
     main = "PET / P Plot (2005)",
     legend = FALSE)                # Automatische Legende deaktivieren

# Benutzerdefinierte Legende hinzufügen
legend("topright",                  # Position der Legende
       legend = rev(c("<= 0.05", "<= 0.1", "<= 0.2", "<= 0.3", "<= 0.4", 
                  "<= 0.7", "<= 0.8", "<= 0.9", "<= 1.0", "<= 1.1", 
                  "<= 1.4", "<= 1.5", "<= 1.75", "<= 2.0", "<= 3.0")), 
       fill = rev(colors[-1]),      # Farben ohne Pink für >3.0
       cex = 0.8,                   # Verkleinere die Schriftgröße der Legende
       title = "Aridity Index", 
       bty = "n")                   # Kein Rahmen um die Legende

### Berechnung Pastoral & Croplands 

### Daten vorbereiten ###
# 1. CFT-Raster vorbereiten
cft_in <- read_io(
  paste0(local_path, "gampe_baseline/cft1700_2005_irrigation_systems_64bands.bin"),
  name = "cft1700_2005_irrigation_systems_64bands.bin",
  descr = "LPJLUSE",
  firstcell = 0, ncell = 67420, 
  firstyear = 1700, nyear = 306,
  scalar = 0.001, nbands = 64, 
  datatype = 1, order = 1,
  cellsize_lat = 0.5, cellsize_lon = 0.5,
  nstep = 1, endian = "little",
  subset = list(year = 306)
)

cft_stack <- stack(lapply(1:64, function(band) {
  as_raster(subset(cft_in, band = band))
}))


### Pasture Fraction
# Definiere die Bänder für "pasture/managed grass (C3/C4 mixed)"
pasture_bands <- c(14, 30, 46, 62)


# Extrahiere und summiere die Bänder für Pasture
pasture_fraction <- calc(cft_stack[[pasture_bands]], sum, na.rm = FALSE)

# Plot der Pasture-Anteile
plot(pasture_fraction, main = "Fraction of Pasture in Each Grid Cell")

### Crop fraction

crop_fraction = sum(cft_stack) - pasture_fraction

# Plot der Crop Fraction
plot(crop_fraction, main = "Fraction of all Crops (non-Pasture")

#Variable  (VegC)
read_meta(paste0(local_path, "Model_Output/FPC_test/vegc.bin.json"))                    # Metadaten anzeigen
vegc_trop = read_io(paste0(local_path, "Model_Output/FPC_test/vegc.bin.json"))    # Datei einlesen
vegc_trop = transform(vegc_trop, to = "lon_lat")          # räumliches
vegc_trop = transform(vegc_trop, to = "year_month_day")   # zeitliches

# Konvertierung in ein RasterLayer-Objekt
vegc_trop_raster <- vegc_trop$as_raster()

Biomass <- vegc_trop_raster * 2

plot(Biomass, main = "Biomass (all PFTs) [gC/m^2]")

#Variable (FPC)
read_meta(paste0(local_path, "Model_Output/FPC_Test/fpc.bin.json"))                    # Metadaten anzeigen
vegFPC_trop = read_io(paste0(local_path, "Model_Output/FPC_Test/fpc.bin.json"))   # Datei einlesen
print(vegFPC_trop)
vegFPC_trop = transform(vegFPC_trop, to = "lon_lat")          # räumliches
vegFPC_trop = transform(vegFPC_trop, to = "year_month_day")   # zeitliches

### Fraction Für einzelne PFTs ausrechnen

free_vegetation <- 1 - calc(cft_stack, sum, na.rm = FALSE)
fpc_evergreen <- as_raster(subset(vegFPC_trop, band = "tropical broadleaved evergreen tree"))
fpc_raingreen <- as_raster(subset(vegFPC_trop, band = "tropical broadleaved raingreen tree"))

# Summieren der beiden Fraktionen zu tropical tree fraction
fpc_tropforest <- fpc_evergreen + fpc_raingreen

#Tropical tress
tropical_tree_fraction <- fpc_tropforest * free_vegetation


# Temperate Bäume (alle temperate PFTs)
fpc_temperate_needle <- as_raster(subset(vegFPC_trop, band = "temperate needleleaved evergreen tree"))
fpc_temperate_broadleaved_evergreen <- as_raster(subset(vegFPC_trop, band = "temperate broadleaved evergreen tree"))
fpc_temperate_broadleaved_summergreen <- as_raster(subset(vegFPC_trop, band = "temperate broadleaved summergreen tree"))
fpc_temperate_trees <- fpc_temperate_needle + fpc_temperate_broadleaved_evergreen + fpc_temperate_broadleaved_summergreen

# Boreale Bäume (alle boreale PFTs)
fpc_boreal_needleleaved_evergreen <- as_raster(subset(vegFPC_trop, band = "boreal needleleaved evergreen tree"))
fpc_boreal_broadleaved_summergreen <- as_raster(subset(vegFPC_trop, band = "boreal broadleaved summergreen tree"))
fpc_boreal_needleleaved_summergreen <- as_raster(subset(vegFPC_trop, band = "boreal needleleaved summergreen tree"))
fpc_boreal_trees <- fpc_boreal_needleleaved_evergreen + fpc_boreal_broadleaved_summergreen + fpc_boreal_needleleaved_summergreen

# Grasland (ohne Polar Grass)
fpc_tropical_c4_grass <- as_raster(subset(vegFPC_trop, band = "Tropical C4 grass"))
fpc_temperate_c3_grass <- as_raster(subset(vegFPC_trop, band = "Temperate C3 grass"))
fpc_grasslands <- fpc_tropical_c4_grass + fpc_temperate_c3_grass

# Polar Grass (extra)
fpc_polar_c3_grass <- as_raster(subset(vegFPC_trop, band = "Polar C3 grass"))

# Berechnung der Wüstenfraktion (free vegetation ohne PFTs)
total_fpc <- tropical_tree_fraction + fpc_temperate_trees + fpc_boreal_trees + fpc_grasslands + fpc_polar_c3_grass


# Tropische Baumfraktion
tropical_tree_fraction <- free_vegetation * fpc_forest

# Temperate Baumfraktion
temperate_tree_fraction <- free_vegetation * fpc_temperate_trees

# Boreale Baumfraktion
boreal_tree_fraction <- free_vegetation * fpc_boreal_trees

# Graslandfraktion (ohne Polar Grass)
grassland_fraction <- free_vegetation * fpc_grasslands

# Polar Grass Fraktion
polar_grass_fraction <- free_vegetation * fpc_polar_c3_grass

# Wüstenfraktion
desert_fraction <- free_vegetation - total_fpc * free_vegetation


# Tropische Baumfraktion
plot(tropical_tree_fraction, 
     main = "Tropical Tree Fraction")

# Temperate Baumfraktion
plot(temperate_tree_fraction, 
     main = "Temperate Tree Fraction")

# Boreale Baumfraktion
plot(boreal_tree_fraction, 
     main = "Boreal Tree Fraction")

# Graslandfraktion
plot(grassland_fraction, 
     main = "Grassland Fraction (except Polar Grass)")

# Polar Grass Fraktion
plot(polar_grass_fraction, 
     main = "Polar Grass Fraction")

# Wüstenfraktion
plot(desert_fraction, 
     main = "Desert Fraction")


# Klassifikationsfunktion
land_cover_classes <- c(
  "Pasture" = 1,
  "Cropland" = 2,
  "Montane Forest" = 3,
  "Closed Rainforest" = 4,
  "Closed Moist Forest" = 5,
  "Open Moist Forest" = 6,
  "Open Dry Forest" = 7,
  "Shrubland" = 8,
  "Temperate Forest" = 9,
  "Savannah/Grassland" = 10,
  "Alpine or Desert" = 11
)

classify_land_cover <- function(pasture, crop, tropical_trees, temperate_trees, grass, temp, Budyko, Biomass) {
  if (is.na(pasture) || is.na(crop) || is.na(tropical_trees) || 
      is.na(temperate_trees) || is.na(grass) || is.na(temp) || 
      is.na(Budyko) || is.na(Biomass)) {
    return(NA) # NA für Meerflächen oder fehlende Daten
  }
  
  if (pasture > 0.5) {
    return(1) # "Pasture"
  } else if (crop > 0.5) {
    return(2) # "Cropland"
  } else if ((pasture + crop) > 0.5) {
    if (pasture > crop) {
      return(1) # "Pasture"
    } else {
      return(2) # "Cropland"
    }
  } else if (tropical_trees > 0.8) {
    if (temp < 24 && Budyko < 1) {
      return(3) # "Montane Forest"
    } else {
      if (Biomass > 40000 && (Budyko < 0.5)) {
        return(4) # "Closed Rainforest"
      } else if (Biomass > 40 && (Budyko < 0.5)) {
        return(5) # "Closed Moist Forest"
      } else if (Biomass <= 40000 && Budyko < 1) {
        return(6) # "Open Moist Forest"
      } else {
        return(7) # "Open Dry Forest"
      }
    }
  } else if (tropical_trees > 0.4) {
    return(8) # "Shrubland"
  } else if (temperate_trees > 0.9) {
    return(9) # "Temperate Forest"
  } else if (temperate_trees > 0.1) {
    return(8) # "Shrubland"
  } else if (grass > 0.4) {
    return(10) # "Savannah/Grassland"
  } else {
    return(11) # "Alpine or Desert"
  }
}

# Rasterüberlagerung und Klassifikation
input_stack <- stack(pasture_fraction, crop_fraction, tropical_tree_fraction, temperate_tree_fraction, 
                     grassland_fraction, mean_global_temperature, Budyko, Biomass)

classified_raster <- calc(input_stack, fun = function(x) {
  classify_land_cover(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8])
})

# Definiere den geographischen Bereich der Tropen
tropen_extent <- extent(-179.75, 179.75, -23.25, 23.25)

# Aktualisierte Farbzuordnung mit hellgrauem Desert
land_cover_colors <- c(
  "#F5DEB3",  # Pasture (Weizenfarben)
  "#FFD700",  # Cropland (Gelb)
  "#228B22",  # Montane Forest (Dunkelgrün)
  "#006400",  # Closed Rainforest (Waldgrün)
  "#32CD32",  # Closed Moist Forest (Limettengrün)
  "#ADFF2F",  # Open Moist Forest (Hellgrün)
  "#808000",  # Open Dry Forest (Olivegrün, grünbraun)
  "#D2691E",  # Shrubland (Schokoladenbraun)
  "#6B8E23",  # Temperate Forest (Olivgrün)
  "#FFFF00",  # Savannah/Grassland (Gelb)
  "#D3D3D3"   # Alpine or Desert (Hellgrau)
)

land_cover_classes <- c(
  "Pasture",
  "Cropland",
  "Montane Forest",
  "Closed Rainforest",
  "Closed Moist Forest",
  "Open Moist Forest",
  "Open Dry Forest",
  "Shrubland",
  "Temperate Forest",
  "Savannah/Grassland",
  "Alpine or Desert"
)

# Plot der Karte ohne Skala (kein raster::plot-Farbbalken)
plot(classified_raster_tropen, 
     col = land_cover_colors, 
     main = "Classified Land Cover in the Tropics",
     legend = FALSE)  # Keine Skala anzeigen

# Legende horizontal anordnen: 6 Klassen in der ersten Reihe und 5 Klassen in der zweiten Reihe
legend("bottomright", 
       legend = land_cover_classes, 
       fill = land_cover_colors, 
       title = "Land Cover Classes", 
       bty = "n", 
       cex = 0.7,   # Kleinere Schriftgröße
       ncol = 2)    # Anzahl der Spalten (erste 6 in einer Reihe, letzte 5 in der anderen)
# Speichern des klassifizierten Rasters
writeRaster(classified_raster, "classified_land_cover.tif", format = "GTiff", overwrite = TRUE)