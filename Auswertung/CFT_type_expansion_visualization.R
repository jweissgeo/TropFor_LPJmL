library(raster)
library(ggplot2)
library(dplyr)
library(colorspace)
library(ggspatial)
library(sf)
library(rnaturalearth)

# CSV-Datei einlesen
Crops_Classes <- read.csv("reclassified_CFTs_area.csv", sep = ",")

# Raster erstellen
Crop_class_expansion <- rasterFromXYZ(Crops_Classes[, c("x", "y", "layer")])
crs(r) <- "+proj=longlat +datum=WGS84"  # Koordinatensystem setzen

# Raster speichern
writeRaster(Crop_class_expansion, "CFT_Expansion_per_class", format="GTiff", overwrite=TRUE)

# Relevante CFT-Legende und zugehörige Farben (angepasst an Atlas-Farbgebung ohne Rot/Grün)
CFT_legend <- c("Rice (Rainfed)", "Maize (Rainfed)", "Tropical Cereals (Rainfed)", 
                "Pulses (Rainfed)", "Tropical Roots (Rainfed)", "Oil Crops Soybean (Rainfed)",
                "Sugar Cane (Rainfed)", "Rice (Surface Watering)")

CFT_layers <- c(2, 3, 4, 5, 7, 9, 12, 18)

CFT_colors <- c("#FFD700", "#FF8C00", "#D2691E", "#8B5A2B", "#A0522D", "#9370DB", "#5F9EA0", "#E6B800")

# Barplot-Daten vorbereiten
barplot_data <- Crops_Classes %>%
  filter(layer %in% CFT_layers) %>%
  group_by(layer) %>%
  summarise(total_areacrop = sum(areacrop, na.rm = TRUE)) %>%
  mutate(class = factor(layer, levels = CFT_layers, labels = CFT_legend)) %>%
  arrange(desc(total_areacrop))

# Barplot erstellen
barplot <- ggplot(barplot_data, aes(x = reorder(class, -total_areacrop), y = total_areacrop, fill = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CFT_colors, guide = FALSE) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Crop Expansion by CFT", x = "CFT Classes", y = "Area Increase (km²)")

# Weltkarte laden
world <- st_as_sf(ne_countries(scale = "medium", returnclass = "sf"))

# Bounding Box
lon_min <- -110
lon_max <- 55
lat_min <- -25
lat_max <- 25

# Raster mit Hintergrundkarte plotten
map_plot <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "black") +
  geom_tile(data = filter(Crops_Classes, layer %in% CFT_layers), aes(x = x, y = y, fill = factor(layer))) +
  scale_fill_manual(values = CFT_colors, labels = CFT_legend, name = "CFT Classes") +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Spatial Distribution of Crop Expansion")

# Einzelne Plots anzeigen
print(map_plot)
print(barplot)

