### Pakete laden ###
require(lpjmlkit)
require(raster)

### local path setzen, andere immer auskommentieren ###
#local_path <- "C:/Users/philipp/Documents/TropFor_LPJmL_LokalData/" #philipp
local_path <- "/Users/epigo/Documents/LPJmL_Lokal/" #Julius
#local_path <- "C:/Dokumente/Umweltsysteme/integrierte_modellierung/"  # Mareike


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

# 2. Tropenausdehnung definieren und Raster stacken
tropen_extent <- extent(-179.75, 179.75, -23.25, 23.25)
cft_stack <- stack(lapply(1:64, function(band) {
  crop(as_raster(subset(cft_in, band = band)), tropen_extent)
}))

# 3. Expansion Potential vorbereiten
expansion_potential <- crop(
  raster("Expansion_Potential/Figure_S3/integrated_expansion_potential.bil"),
  tropen_extent
)
expansion_resampled <- resample(expansion_potential, cft_stack, method = "bilinear")


# 4. Freie Vegetation und Waldfläche berechnen
free_vegetation <- 1 - calc(cft_stack, sum, na.rm = TRUE)
fpc_tropical_forest <- read_io(paste0(local_path, "Model_Output/FPC_Test/fpc.bin.json"))
fpc_forest <- crop(as_raster(subset(fpc_tropical_forest, band = "tropical broadleaved evergreen tree")), tropen_extent) +
  crop(as_raster(subset(fpc_tropical_forest, band = "tropical broadleaved raingreen tree")), tropen_extent)
forest_area <- free_vegetation * fpc_forest

### Funktion laden ###
source("forest_function_advanced.R")

### Schwellenwert berechnen ###
target_percent_loss <-15  # Ziel-Prozentverlust
result <- find_threshold(
  target_percent_loss = target_percent_loss,
  forest_area = forest_area,
  expansion_resampled = expansion_resampled
)

### Ergebnisse ausgeben ###
print(paste("Schwellenwert des Index für ", target_percent_loss, "% Verlust: ", round(result$index, 2)))

calculated_index <- result$index  # Der berechnete Schwellenwert
expansion_mask <- expansion_resampled > calculated_index  # Maske basierend auf dem Index

### Plot der Expansion Mask ###
plot(expansion_mask, main = paste("Expansion Mask mit Index", round(calculated_index, 2), ", bei Waldverlust", target_percent_loss, "%"))

# zum Vergleich
percent_forest <- forest_area * 100  # In Prozent umwandeln
plot(percent_forest, main = "Prozentuale Gesamtwaldfläche")

### Weiter geht es mit dem kNN Algorythmus, der bestimmt, welche Croptypes angebaut werden sollen in der Flöche. Siehe anderes Srkipt
### Wird dann hier eingefügt


### Letzter Schritt: Binärdatei für LPJmL so schreiben, dass dann wirklich auch nur Waldfläche in CFTs umgewandelt wird (keine Grasslands)
