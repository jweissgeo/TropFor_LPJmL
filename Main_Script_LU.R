### Pakete laden
require(lpjmlkit)
require(raster)
require(FNN)
library(leaflet)
#
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
# Main Raster Stack
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

## Funktion laden
source("forest_function_advanced.R")

# Ziel-Prozentverlust definieren
target_percent_loss <- 15

# Maske für Waldflächen mit Werten > 0.8 erstellen
forest_mask <- forest_area > 0.8
plot(forest_mask)
forest_area_filtered <- forest_area * forest_mask
plot(forest_area_filtered)

# Schwellenwert berechnen
result <- find_threshold(
  target_percent_loss = target_percent_loss,
  forest_area = forest_area_filtered,  # Gefilterte Waldflächen übergeben
  expansion_resampled = expansion_resampled
)

# Ergebnis ausgeben
print(result)

### Ergebnisse ausgeben ###
print(paste("Schwellenwert des Index für ", target_percent_loss, "% Verlust: ", round(result$index, 2)))

calculated_index <- result$index  # Der berechnete Schwellenwert
expansion_mask <- (expansion_resampled > calculated_index) & forest_mask  # Maske basierend auf dem Index

### Plot der Expansion Mask ###
plot(expansion_mask, main = paste("Expansion Mask mit Index", round(calculated_index, 2), ", bei Waldverlust", target_percent_loss, "%"))

# zum Vergleich
percent_forest <- forest_area * 100  # In Prozent umwandeln
plot(percent_forest, main = "Prozentuale Gesamtwaldfläche")

##calc area
calc_total_area <- function(area_raster) {
  coords <- coordinates(area_raster)
  cell_areas <- calc_cellarea(coords[, 2], return_unit = "km2")
  total_area <- sum(values(area_raster) * cell_areas, na.rm = TRUE)
  return(total_area)
}

# Gesamtfläche vor der Maskierung (in km²)
total_forest_area_km2 <- calc_total_area(forest_area * forest_mask)

# Gesamtfläche der entfernten Waldflächen (in km²)
removed_area_km2 <- calc_total_area(forest_area * forest_mask * expansion_mask)

# Prozentualer Verlust
percent_loss <- (removed_area_km2 / total_forest_area_km2) * 100

# Ergebnisse anzeigen
cat("Gesamtwaldfläche: ", total_forest_area_km2, "km²\n")
cat("Entfernte Waldfläche: ", removed_area_km2, "km²\n")
cat("Prozentualer Verlust: ", percent_loss, "%\n")

# Karte erstellen: Gesamtwaldfläche und rote Pixel für Abholzung

library(RColorBrewer)
# Basiswaldkarte
forest_map <- forest_area * forest_mask

# Abholzungskarte (nur relevante Pixel behalten)
deforestation_map <- forest_area * forest_mask * expansion_mask

# Hintergrund- und Null-Werte in der Abholzungskarte entfernen
deforestation_map[deforestation_map == 0] <- NA

# Anpassung der Farbpalette für die Waldkarte
# Basierend auf 'Greens' von RColorBrewer, aber erweitert
forest_colors <- c("lightgray", brewer.pal(9, "Greens"))

# Plot der Gesamtwaldfläche
plot(forest_map, 
     col = forest_colors, 
     main = "Gesamtwaldfläche mit Abholzung", 
     legend = TRUE)

# Hinzufügen der Abholzungspixel (nur relevante Pixel rot zeichnen)
plot(deforestation_map, 
     col = "red", 
     add = TRUE, 
     legend = FALSE)

### Trainingsdaten für k- nearest neighbour vorbereiten ###
# Indizes der auszuschließenden und eingeschlossenen CFTs
exclude_cfts <- c(13, 14, 15, 16, 29, 30, 31, 32, 45, 46, 47, 48, 61, 62, 63, 64)
include_cfts <- setdiff(1:64, exclude_cfts)

# Extrahieren der Werte aus den inkludierten Bändern
train_values <- getValues(cft_stack[[include_cfts]])
train_coords <- coordinates(cft_stack[[1]]) # holt sich die Koordinaten aus dem ersten Band

# Filtern der gültigen Pixel (keine fehlenden Werte in den inkludierten Bändern)
valid_train <- complete.cases(train_values)
train_values <- train_values[valid_train, ]
train_coords <- train_coords[valid_train, ]

# Definieren der dominierenden CFT-Klasse (mit Schwellenwert, um unrealistische Ergebnisse zu verhindern; nur jetzt schon hochskalierte Crops)
train_classes <- apply(train_values, 1, function(row) {
  max_val <- max(row, na.rm = TRUE)
  if (max_val >= 0.15) which.max(row) else NA #Schwellenwert 0.15
})

# Entfernen von Zellen ohne gültige dominierende Klasse
valid_train <- !is.na(train_classes)
train_values <- train_values[valid_train, ]
train_coords <- train_coords[valid_train, ]
train_classes <- train_classes[valid_train]

# Erstelle eine Liste der ursprünglichen Indizes der eingeschlossenen CFTs
original_cfts <- include_cfts  # Diese enthält die ursprünglichen Indizes der Bänder

# Mappe train_classes zurück zu den ursprünglichen CFT-Indices
train_classes_original <- original_cfts[train_classes]

"# Konvertiere in einen DataFrame mit den originalen CFT-Indices
result_table <- data.frame(
  train_coords,  # Koordinaten
  train_classes_original,  # Dominierende Klasse mit Originalindex
  train_values  # Werte der Bänder
)

# Tabelle anzeigen
View(result_table)"

### Zielzellen vorbereiten ###
# Extrahieren der Zielzellen basierend auf der Expansion Potential-Map
target_mask <- expansion_mask[] & !is.na(expansion_mask[])
target_coords <- coordinates(expansion_mask)[target_mask, ]

### kNN-Algorithmus ausführen ###
k <- 5  # Anzahl der Nachbarn die in die Berechnung einfließen
knn_result <- knn(
  train = train_coords,
  test = target_coords,
  cl = train_classes_original,
  k = k
)

### Erstellung des Rasterstacks mit aktualisierten Werten ###
# Erstelle einen leeren Raster-Stack für die neuen Werte
expanded_raster_stack <- stack()

# Iterieren über alle 64 CFTs und Einsetzen der Werte basierend auf den kNN-Ergebnissen
for (i in 1:64) {
  cft_raster <- raster(expansion_mask)
  values(cft_raster) <- NA
  values(cft_raster)[target_mask] <- as.numeric(knn_result == i)
  expanded_raster_stack <- addLayer(expanded_raster_stack, cft_raster)
}

### Ergebnis plotten und speichern ###
# Plot eines Beispiels
plot(expanded_raster_stack[[3]], main = "Expansion für CFT3 - Maize")

# Summe über alle Layer berechnen und plotten
#summed_raster <- calc(expanded_raster_stack, sum, na.rm = TRUE)
#plot(summed_raster, main = "Summierte Verteilung über alle CFTs")

# Balkendiagramm der Flächen um CFT-Expansiomn zu kontrollieren
area_by_cft <- sapply(1:nlayers(expanded_raster_stack), function(i) {
  sum(values(expanded_raster_stack[[i]]) == 1, na.rm = TRUE)
})
barplot(area_by_cft, names.arg = 1:64, 
        main = "Pixelzuwachs für jedes CFT", 
        xlab = "CFT Index", ylab = "Anzahl Zellen", col = "lightblue")

# Saven des neuen Raster-Stacks
writeRaster(expanded_raster_stack, 
            filename = "Expansion_Potential/expanded_cropland_stack_simple.tif",
            format = "GTiff", overwrite = TRUE) #Achtung!!!!!

### Karte zur Darstellung ###
# Kombinieren der Raster der Klassen 1-12 und Klasse 18
cft_expansion <- subset(expanded_raster_stack, c(1:12, 18))

# Zusammenfügen der Klassen (jede Zelle wird der höchsten Klasse zugeordnet)
combined_cft_raster <- which.max(cft_expansion)

### Definition der Legende ###
# Namen der CFT-Klassen (1-12 und 18)
cft_names <- c(
  "Temperate Cereals", "Rice", "Maize", "Tropical Cereals",
  "Pulses", "Temperate Roots", "Tropical Roots",
  "Oil Crops Sunflower", "Oil Crops Soybean", 
  "Oil Crops Groundnut", "Oil Crops Rapeseed", "Sugar Cane",
  "Rice (Surface Watering)"
)

# Farben für die CFT-Klassen (1-12 und 18)
cft_colors <- c(
  "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
  "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a",
  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99",
  "#41ab5d"  # Farbe für "Rice (Surface Watering)"
)

# Überprüfe die Werte im Raster
unique_values <- sort(unique(values(combined_cft_raster)))  # Eindeutige Werte im Raster (z. B. 1-13)


###Plot der CFT-Karte mit grauem hintergrund
# Plot der Waldkarte als Hintergrund (geht vlt)
plot(
  forest_map, 
  col = c("lightgray"), 
  legend = FALSE,
  main = "CFT Expansion der Klassen 1-12 und 18"
)

# Hinzufügen der kombinierten CFT-Karte
plot(
  combined_cft_raster, 
  col = cft_colors[unique_values], 
  legend = FALSE, 
  add = TRUE
)

# Hinzufügen der Legende (nur für vorhandene Klassen)
legend(
  "topright", 
  legend = cft_names[unique_values], 
  fill = cft_colors[unique_values], 
  title = "CFT Klassen", 
  cex = 0.6, 
  bty = "n"
)

#### Leaflet Karte zum überprüfen


### Konvertiere das kombinierte Raster in ein SpatialPointsDataFrame
spdf <- as(combined_cft_raster, "SpatialPointsDataFrame")

### Füge die Daten in ein DataFrame
df <- as.data.frame(spdf)
df$x <- coordinates(spdf)[, 1]  # x-Koordinaten (Longitude)
df$y <- coordinates(spdf)[, 2]  # y-Koordinaten (Latitude)

# Farben basierend auf der CFT-Klasse zuweisen
df$color <- cft_colors[df$layer]  # Farbe basierend auf der dominierenden Klasse

### Erstelle eine interaktive Leaflet-Karte
leaflet(data = df) %>%
  addTiles() %>%  # Basemap hinzufügen
  addCircleMarkers(
    ~x, ~y,  # Koordinaten
    radius = 2,  # Größe der Punkte
    color = ~color,  # Farbe basierend auf der Klasse
    stroke = FALSE,  # Keine Umrandung
    fillOpacity = 0.7,
    popup = ~paste(
      "<b>Koordinaten:</b>", x, y, "<br>",
      "<b>Dominierende Klasse:</b>", cft_names[df$layer]
    )
  )






##### Input für das Modell vorbereiten (.bin)
##### Angelehnt an Skript, das wir am 5.12. im Seminar gemacht haben

writeRaster(combined_cft_raster, filename = "cft_in_tropics.tif", format = "GTiff", overwrite = TRUE)
raster_data = raster("cft_in_tropics.tif")
grid_lpjml = read_io(paste0(local_path, "gampe_baseline/grid.bin.json"))

# Ein langer Vektor; die X und Y Werte sollen zusammen gebracht werden; Vektor kann dann verwendet werden um nachzuschauen wo welche Zelle ist
# Für grid und raster Koordinaten (grid und raster Auslesereihenfolge der Zellen ist unterschiedlich. Deswegen muss man über die Koordinaten vorgehen.)
grid_coords = paste(round(grid_lpjml$data[1,,],2), round(grid_lpjml$data[,,2],2), sep = "_")           # round(...,2) --> Zwei Nachkommastellen; Erst die ersten Werte, dann die zweiten Werte
ras_coords = paste(round(coordinates(raster_data)[,1],2), round(coordinates(raster_data)[,2],2), sep = "_")

# Vektor mit gleichen Sachen und richtiger Reihenfolge
cft_out = array(0,dim = (c(1,67420,64)))    # 1 Jahr, 67420 Zellen, 64 Bänder

# convert landuse grid (In R) to format needed in LPJmL
# R: cell - year - band
# LPJmL: year - cell - band
y = 1                   # man kann auch gleich 1 für y verwenden; die Zeile ist nur falls man mehrere Jahre machen möchte
for(i in 1:67420){
  for(b in 1:64){
    cft_out[y,i,b] = as.numeric(cft_in$data[i,y,b])    #cft_out[year, cell, band] = als Zahl (cft_in[cell, year, band])
    names(cft_out[y,i,b]) = names(cft_in$data[i,y,b])  # Name ändern
  }
}

# cft_out soll an Landnutzungdfile angehangen werden
for(i in 1:length(ras_coords)){
  if(is.na(match(ras_coords[i], grid_coords))){        # no data Abfrage
    next
  }else{
    cft_out[y, match(ras_coords[i], grid_coords), include_cfts] = raster_data[i]
  }
}

# land use file verketten
# use raster-year as new year (y=307)
# Jahr öffnen, nehmen, an f.out hängen (306 Jahre)
# es gibt 306 Jahre. Das 307te Jahr ist so wie das 306te Jahr mit dem Unterschied, dass unsere individuelle Veränderung (z.B. sugar cane in Tropen) dabei ist.
f.out <- file("cft_intmod24_deforestation.bin", "wb")     # file erstellen und öffnen, ist so lange geöffnet bis es wieder geschlossen wird
n_year_out = 307
for(y in 1:n_year_out){
  print(y)
  if(y<n_year_out){
    cft_in_tmp = read_io(paste0(local_path, "gampe_baseline/cft1700_2005_irrigation_systems_64bands.bin"),       # die folgenden Infos kennen wir aus dem header
                         name = "cft1700_2005_irrigation_systems_64bands.bin",
                         descr = "LPJLUSE",
                         firstcell = 0, ncell = 67420,                  # es gibt 67420 Zellen auf Land
                         firstyear = 1700, nyear = 306,
                         scalar = 0.001, nbands = 64,
                         datatype = 1, order = 1,                       # Datatype 1 ist integer
                         cellsize_lat = 0.5, cellsize_lon = 0.5,
                         nstep = 1, endian = "little",                  # endian Sortierung (kleine nach groß)
                         subset = list(year=306))                       # year 306 ist 2005 (es geht bei 1700 los)
    for(i in 1:67420){
      writeBin(as.integer(round(cft_in_tmp$data[i,1,]*1000)),f.out, size = 2, endian = "little", useBytes = TRUE)      # jedes Jahr wird hier hinten dran gehangen     | [iter Wert, ein Jahr, alle Bänder]
    }
  }else{
    for(i in 1:67420){
      writeBin(as.integer(round(cft_out[1,i,]*1000)), f.out, size = 2, endian ="little", useBytes = TRUE)
    }
  }
}

close(f.out)         # nach dem Ausführen dieser Zeile ist die Datei erst richtig da im Ordner

# write function to merge two binary files
merge_binary_files = function(file1, file2, output_file){
  f1 = readBin(file1, what = "raw", n = file.info(file1)$size)     # raw heißt es gibt keinen header, nur einlesen
  f2 = readBin(file2, what = "raw", n = file.info(file2)$size)
  out_f = c(f1, f2)       # f1 und f2 hintereinander gehangen
  writeBin(out_f, output_file)
}

# header schreiben; der soll am Ende vor die Daten geschrieben werden, die davor erstellt wurden
header_cft_out = create_header(
  name = "LPJLUSE",
  version = 2,
  order = 1,               # order 1 ist für land use files
  firstyear = 1700,
  nyear = 307,             # 307 jahre insgesamt
  firstcell = 0,
  ncell = 67420,
  nband = 64,
  cellsize_lon = 0.5,
  scalar =  0.001,         # für internes Umrechnen
  cellsize_lat = 0.5,
  datatype = 1,            # 1 ist integer
  nstep = 1,               # jähricher Zeitschritt
  endian = "little",
  verbose = TRUE
)

write_header(filename ="landuse_header.bin", header = header_cft_out, overwrite = TRUE)
merge_binary_files("landuse_header.bin", "cft_intmod24_deforestation.bin", "cft_intmod24_deforestation_final.bin")

header_cft_new = read_header("cft_intmod24_deforestation_final.bin")

cft_new = read_io("cft_intmod24_deforestation_final.bin",
                  name = "cft_intmod24_deforestation_final.bin",
                  descr = "LPJLUSE",
                  firstcell = 0, ncell = 67420,                  # es gibt 67420 Zellen auf Land
                  firstyear = 1700, nyear = 307,
                  scalar = 0.001, nbands = 64,
                  datatype = 1, order = 1,                       # Datatype 1 ist integer
                  cellsize_lat = 0.5, cellsize_lon = 0.5,
                  nstep = 1, endian = "little",                  # endian Sortierung (kleine nach groß)
                  subset = list(year=307))                       # year 307 ist 2006 (erstes Jahr mit modifikation)

plot(subset(cft_new))
