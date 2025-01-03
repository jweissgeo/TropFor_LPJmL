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

### Funktion laden ###
source("forest_function_advanced.R")

### Schwellenwert berechnen ###
target_percent_loss <- 15  # Ziel-Prozentverlust
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

# Plotten der kombinierten Karte
plot(
  combined_cft_raster, col = cft_colors[unique_values], legend = FALSE,
  main = "CFT Expansion der Klassen 1-12 und 18"
)

# Hinzufügen der Legende (nur für vorhandene Klassen)
legend(
  "topright", legend = cft_names[unique_values], fill = cft_colors[unique_values], 
  title = "CFT Klassen", cex = 0.6, bty = "n"
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


### Letzter Schritt: Binärdatei für LPJmL so schreiben, dass die CFTs anteikig umgewandelt werden (CFTs erhöhen so viel wie Free Vegetation da ist) 



#Dummy script
"
forest_area_expansion <- forest_area * expansion_mask
combined_raster_CFT_FORESTEXPANSION <- stack(forest_area_expansion, combined_cft_raster)
print(combined_raster_CFT_FORESTEXPANSION)

plot(forest_area_expansion + calc(cft_stack, sum, na.rm = TRUE) * expansion_mask + (free_vegetation - forest_area) * expansion_mask)

# 1. Konvertiere das kombinierte Raster in ein SpatialPointsDataFrame
combined_raster_CFT_FORESTEXPANSION <- stack(forest_area_expansion, combined_cft_raster)

# Konvertiere den Rasterstack in einen DataFrame
df <- as.data.frame(combined_raster_CFT_FORESTEXPANSION, xy = TRUE, na.rm = TRUE)

# Spaltennamen anpassen
colnames(df) <- c("Longitude", "Latitude", "ForestChange", "CFT_Class")

# Zeige die ersten Zeilen der Tabelle an
head(df)

# Optional: Gesamte Tabelle anzeigen lassen
#View(df)


library(scales)
library(viridis)

# 1. Werte filtern: Zeilen mit ForestChange < 0.01 entfernen
df <- df[df$ForestChange >= 0.8, ]
head(df)
View(df)
# 2. Farbpalette für CFT-Klassen (z. B. viridis)
cft_classes <- unique(df$CFT_Class)
num_classes <- length(cft_classes)
palette <- colorFactor(viridis(num_classes), domain = cft_classes)

# 3. Normalisiere ForestChange für Transparenz und Farbe
df$opacity <- rescale(df$ForestChange, to = c(0.3, 1))  # Opazität von 0.3 bis 1
df$color <- mapply(function(cft, op) {
  adjustcolor(palette(cft), alpha.f = op)
}, df$CFT_Class, df$opacity)

# 4. Leaflet-Karte erstellen
leaflet(data = df) %>%
  addTiles() %>%  # Basemap hinzufügen
  addCircleMarkers(
    ~Longitude, ~Latitude,  # Koordinaten
    radius = 5,  # Punktgröße
    color = ~color,  # Farbe mit Intensität
    fillOpacity = ~opacity,  # Transparenz basierend auf ForestChange
    stroke = FALSE,  # Keine Umrandung
    popup = ~paste(
      "<b>Koordinaten:</b>", Longitude, ", ", Latitude, "<br>",
      "<b>ForestChange:</b>", round(ForestChange, 2), "<br>",
      "<b>CFT_Class:</b>", CFT_Class
    )
  )

"
