### Funktion zur Threshold-Bestimmung ###

### Pakete laden
require(lpjmlkit)
require(raster)


### local path setzen, andere immer auskommentieren ###
#local_path <- "C:/Users/philipp/Documents/TropFor_LPJmL_LokalData/" #philipp
local_path <- "/Users/epigo/Documents/LPJmL_Lokal/" #Julius
#local_path <- "C:/Dokumente/Umweltsysteme/integrierte_modellierung/"  # Mareike


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
  
  # 2. Tropenausdehnung definieren
  tropen_extent <- extent(-179.75, 179.75, -23.25, 23.25)
  
  # 3. Alle Bänder in Raster konvertieren und auf Tropen beschränken
  cft_stack <- stack(lapply(1:64, function(band) {
    crop(as_raster(subset(cft_in, band = band)), tropen_extent)
  }))
  
  # 4. Berechnung der freien Vegetation pro Pixel
  free_vegetation <- 1 - calc(cft_stack, sum, na.rm = TRUE)
  
  # 5. FPC-Wert laden und auf Tropen beschränken
  fpc_tropical_forest <- read_io(paste0(local_path, "Model_Output/FPC_Test/fpc.bin.json"))
  fpc_rainforest <- crop(as_raster(subset(fpc_tropical_forest, band = "tropical broadleaved evergreen tree")), tropen_extent)
  fpc_dryforest <- crop(as_raster(subset(fpc_tropical_forest, band = "tropical broadleaved raingreen tree")), tropen_extent)
  
  # 6. Waldflächen pro Pixel berechnen
  rainforest_area <- free_vegetation * fpc_rainforest
  dryforest_area <- free_vegetation * fpc_dryforest
  forest_area <- rainforest_area + dryforest_area  # Gesamter tropischer Wald
  

  # Karten erstellen
  percent_forest <- forest_area * 100  # In Prozent umwandeln
  percent_rainforest <- rainforest_area * 100
  percent_dryforest <- dryforest_area * 100
  percent_cft <- calc(cft_stack, sum, na.rm = TRUE) * 100  # In Prozent umwandeln
  percent_free_veg <- free_vegetation * 100
  

  plot(percent_dryforest, main = "Prozentuale Trockenwaldfläche")
  plot(percent_forest, main = "Prozentuale Gesamtwaldfläche")
  plot(percent_rainforest, main = "Prozentuale Regenwaldfläche")
  plot(percent_cft, main = "Prozentuale Anbaufläche (CFT)")
  plot(percent_free_veg, main = "Prozentuale freie Vegetationsfläche")
  