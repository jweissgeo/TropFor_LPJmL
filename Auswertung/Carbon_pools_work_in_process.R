require(lpjmlkit)
require(raster)
require(caTools)
require(trend)
require(ggplot2)
require(dplyr)
require(scales)

# Set working directory
setwd("C:/Users/philipp/Documents/TropForLPJmL_outputs/global_pools/")

# Load forest mask
forest_mask <- raster("forest_mask.tif")
if (is.null(forest_mask)) stop("Fehler: forest_mask konnte nicht geladen werden.")

# Load data
agb_newLU <- read_io("newLUconst2011/agb.bin.json")
vegC_newLU <- read_io("newLUconst2011/vegc.bin.json")
litC_newLU <- read_io("newLUconst2011/litc.bin.json")
soilC_newLU <- read_io("newLUconst2011/soilc.bin.json")

agb_oldLU <- read_io("oldLUconst2011/agb.bin.json")
vegC_oldLU <- read_io("oldLUconst2011/vegc.bin.json")
litC_oldLU <- read_io("oldLUconst2011/litc.bin.json")
soilC_oldLU <- read_io("oldLUconst2011/soilc.bin.json")

dataset_names <- c("AGB NewLU", "VegC NewLU", "LitC NewLU", "SoilC NewLU", "AGB OldLU", "VegC OldLU", "LitC OldLU", "SoilC OldLU")

# Funktion zur Transformation
transform_LPJmLData <- function(df, name) {
  print(paste("Transformiere", name))
  if (!inherits(df, "LPJmLData")) stop(paste("Fehler: Die Datei", name, "ist kein LPJmLData-Objekt."))
  df <- df$transform(to = "lon_lat")$transform(to = "year")
  return(df)
}

# Transformierte Daten
agb_newLU <- transform_LPJmLData(agb_newLU, "AGB NewLU")
vegC_newLU <- transform_LPJmLData(vegC_newLU, "VegC NewLU")
litC_newLU <- transform_LPJmLData(litC_newLU, "LitC NewLU")
soilC_newLU <- transform_LPJmLData(soilC_newLU, "SoilC NewLU")

agb_oldLU <- transform_LPJmLData(agb_oldLU, "AGB OldLU")
vegC_oldLU <- transform_LPJmLData(vegC_oldLU, "VegC OldLU")
litC_oldLU <- transform_LPJmLData(litC_oldLU, "LitC OldLU")
soilC_oldLU <- transform_LPJmLData(soilC_oldLU, "SoilC OldLU")

# Funktion zur Berechnung des Mittelwerts für 2001-2011 und Umrechnung in kg/m²
calculate_mean_2001_2011 <- function(df, name) {
  print(paste("Berechne Mittelwert für 2001-2011:", name))
  df_period <- df$subset(year = as.character(2001:2011))
  raster_stack <- stack(df_period$as_raster())
  df_mean <- calc(raster_stack, mean, na.rm = TRUE)
  
  # Umrechnung von gC/m² in kg/m²
  df_mean <- df_mean * 0.001
  
  return(df_mean)
}

datasets_newLU_LA <- list(
  calculate_mean_2001_2011(agb_newLU, "AGB NewLU") %>% crop(extent(-120, -30, -60, 30)),
  calculate_mean_2001_2011(vegC_newLU, "VegC NewLU") %>% crop(extent(-120, -30, -60, 30)),
  calculate_mean_2001_2011(litC_newLU, "LitC NewLU") %>% crop(extent(-120, -30, -60, 30)),
  calculate_mean_2001_2011(soilC_newLU, "SoilC NewLU") %>% crop(extent(-120, -30, -60, 30))
)


datasets_newLU_AF <- list(
  calculate_mean_2001_2011(agb_newLU, "AGB NewLU") %>% crop(extent(-30, 60, -40, 35)),
  calculate_mean_2001_2011(vegC_newLU, "VegC NewLU") %>% crop(extent(-30, 60, -40, 35)),
  calculate_mean_2001_2011(litC_newLU, "LitC NewLU") %>% crop(extent(-30, 60, -40, 35)),
  calculate_mean_2001_2011(soilC_newLU, "SoilC NewLU") %>% crop(extent(-30, 60, -40, 35))
)

datasets_oldLU_LA <- list(
  calculate_mean_2001_2011(agb_oldLU, "AGB OldLU") %>% crop(extent(-120, -30, -60, 30)),
  calculate_mean_2001_2011(vegC_oldLU, "VegC OldLU") %>% crop(extent(-120, -30, -60, 30)),
  calculate_mean_2001_2011(litC_oldLU, "LitC OldLU") %>% crop(extent(-120, -30, -60, 30)),
  calculate_mean_2001_2011(soilC_oldLU, "SoilC OldLU") %>% crop(extent(-120, -30, -60, 30))
)

datasets_oldLU_AF <- list(
  calculate_mean_2001_2011(agb_oldLU, "AGB OldLU") %>% crop(extent(-30, 60, -40, 35)),
  calculate_mean_2001_2011(vegC_oldLU, "VegC OldLU") %>% crop(extent(-30, 60, -40, 35)),
  calculate_mean_2001_2011(litC_oldLU, "LitC OldLU") %>% crop(extent(-30, 60, -40, 35)),
  calculate_mean_2001_2011(soilC_oldLU, "SoilC OldLU") %>% crop(extent(-30, 60, -40, 35))
)



# Maskierung mit der Waldmaske
apply_forest_mask <- function(df, mask, name) {
  print(paste("Wende Maske auf", name, "an."))
  if (!compareRaster(df, mask, extent = FALSE, stopiffalse = FALSE)) {
    print("Extent stimmt nicht überein - Resampling...")
    df <- resample(df, mask, method = "bilinear")
  }
  masked_data <- mask(df, mask)
  return(masked_data)
}

datasets_newLU <- Map(apply_forest_mask, datasets_newLU, mask = list(forest_mask), name = dataset_names[1:4])
datasets_oldLU <- Map(apply_forest_mask, datasets_oldLU, mask = list(forest_mask), name = dataset_names[5:8])

# Funktion zum Ausschneiden von Kontinenten (Asien entfernen)
subset_region <- function(raster_data, region, name) {
  extent_filter <- switch(region,
                          "Latin America" = extent(-120, -30, -60, 30),  
                          "Africa" = extent(-30, 60, -40, 35)
  )
  print(paste("Erstelle Subset für", name, "-", region))
  raster_subset <- crop(raster_data, extent_filter)
  return(raster_subset)
}

# Subsets für Lateinamerika und Afrika
atasets_newLU_LA <- Map(subset_region, datasets_newLU, 
                        region = rep("Latin America", length(datasets_newLU)), 
                        name = dataset_names[1:4])

datasets_newLU_AF <- Map(subset_region, datasets_newLU, 
                         region = rep("Africa", length(datasets_newLU)), 
                         name = dataset_names[1:4])

datasets_oldLU_LA <- Map(subset_region, datasets_oldLU, 
                         region = rep("Latin America", length(datasets_oldLU)), 
                         name = dataset_names[5:8])

datasets_oldLU_AF <- Map(subset_region, datasets_oldLU, 
                         region = rep("Africa", length(datasets_oldLU)), 
                         name = dataset_names[5:8])

# Funktion zur Berechnung von Below Ground Biomass (BGB)
calculate_bgb <- function(vegC, agb, litC, name) {
  print(paste("Berechne BGB für:", name))
  bgb <- vegC - agb - litC
  return(bgb)
}

# Berechnung von BGB für beide Szenarien und Regionen
bgb_newLU_LA <- calculate_bgb(datasets_newLU_LA[[2]], datasets_newLU_LA[[1]], datasets_newLU_LA[[3]], "BGB NewLU Latin America")
bgb_newLU_AF <- calculate_bgb(datasets_newLU_AF[[2]], datasets_newLU_AF[[1]], datasets_newLU_AF[[3]], "BGB NewLU Africa")

bgb_oldLU_LA <- calculate_bgb(datasets_oldLU_LA[[2]], datasets_oldLU_LA[[1]], datasets_oldLU_LA[[3]], "BGB OldLU Latin America")
bgb_oldLU_AF <- calculate_bgb(datasets_oldLU_AF[[2]], datasets_oldLU_AF[[1]], datasets_oldLU_AF[[3]], "BGB OldLU Africa")

# Funktion zur Berechnung des Mittelwerts und Umrechnung in kg/m²
safe_mean <- function(var) {
  if (!is.null(var) && length(values(var)) > 0) {
    return(mean(values(var), na.rm = TRUE))
  } else {
    return(0)
  }
}

# Datenframe für das Diagramm
df_bar <- data.frame(
  Region = rep(c("Latin America", "Latin America", "Africa", "Africa"), each = 4),
  LandUse = rep(c("Old LU", "New LU", "Old LU", "New LU"), each = 4),
  CarbonPool = rep(c("Above Ground Biomass", "Litter C", "Below Ground Biomass", "Soil C"), times = 4),
  Values = c(
    # Old LU Latin America
    safe_mean(datasets_oldLU_LA[[1]]), safe_mean(datasets_oldLU_LA[[3]]), 
    -abs(safe_mean(bgb_oldLU_LA)), -abs(safe_mean(datasets_oldLU_LA[[4]])),
    
    # New LU Latin America
    safe_mean(datasets_newLU_LA[[1]]), safe_mean(datasets_newLU_LA[[3]]), 
    -abs(safe_mean(bgb_newLU_LA)), -abs(safe_mean(datasets_newLU_LA[[4]])),
    
    # Old LU Africa
    safe_mean(datasets_oldLU_AF[[1]]), safe_mean(datasets_oldLU_AF[[3]]), 
    -abs(safe_mean(bgb_oldLU_AF)), -abs(safe_mean(datasets_oldLU_AF[[4]])),
    
    # New LU Africa
    safe_mean(datasets_newLU_AF[[1]]), safe_mean(datasets_newLU_AF[[3]]), 
    -abs(safe_mean(bgb_newLU_AF)), -abs(safe_mean(datasets_newLU_AF[[4]]))
  )
)

# Gesamtwerte berechnen
df_totals <- df_bar %>%
  group_by(Region, LandUse) %>%
  summarise(TotalCarbon = sum(abs(Values), na.rm = TRUE), .groups = "drop")

require(grid)  # Für Pfeile

# Korrekte Beschriftung der Karten mit Einheiten
par(mfrow = c(2, 4))  # 2 Zeilen, 4 Spalten für bessere Übersicht

# Funktion zum Plotten mit Einheitenbeschriftung
plot_with_units <- function(raster_data, title) {
  plot(raster_data, main = paste(title, "\n(kg/m²)"))
}

# Karten für Lateinamerika - New Land Use
for (i in 1:4) {
  plot_with_units(datasets_newLU_LA[[i]], 
                  paste("New Land Use - Latin America -", 
                        c("Above Ground Biomass", "Vegetation Carbon", "Litter Carbon", "Soil Carbon")[i]))
}

# Karten für Lateinamerika - Old Land Use
for (i in 1:4) {
  plot_with_units(datasets_oldLU_LA[[i]], 
                  paste("Old Land Use - Latin America -", 
                        c("Above Ground Biomass", "Vegetation Carbon", "Litter Carbon", "Soil Carbon")[i]))
}

# Karten für Afrika - New Land Use
for (i in 1:4) {
  plot_with_units(datasets_newLU_AF[[i]], 
                  paste("New Land Use - Africa -", 
                        c("Above Ground Biomass", "Vegetation Carbon", "Litter Carbon", "Soil Carbon")[i]))
}

# Karten für Afrika - Old Land Us
for (i in 1:4) {
  plot_with_units(datasets_oldLU_AF[[i]], 
                  paste("Old Land Use - Africa -", 
                        c("Above Ground Biomass", "Vegetation Carbon", "Litter Carbon", "Soil Carbon")[i]))
}


# Farbpalette für das Diagramm
color_palette <- c("Above Ground Biomass" = "#006400", "Below Ground Biomass" = "#8B4513", 
                   "Litter C" = "#E69F00", "Soil C" = "#A9A9A9")

# ggplot2-Plot ohne Pfeile
p <- ggplot(df_bar, aes(x = LandUse, y = Values, fill = CarbonPool)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_text(data = df_totals, aes(x = LandUse, y = -max(abs(df_bar$Values)) * 0.1, 
                                  label = paste0("Total: ", round(TotalCarbon, 1), " kg/m²")),
            vjust = 1.5, size = 3, fontface = "bold", inherit.aes = FALSE) +
  labs(title = "Carbon Pool Comparison: Above & Below Ground (2001-2011)", 
       y = "Carbon Storage (kg/m²)", 
       x = "") +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(expand = c(0, 0), labels = abs, breaks = pretty_breaks(n = 10)) +
  scale_x_discrete(limits = c("Old LU", "New LU")) +
  facet_wrap(~Region) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(0.6, "cm"),  # Kleinere Kästen in der Legende
        plot.title = element_text(size = 20, face = "bold", margin = margin(b = 30)),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, margin = margin(t = 10, b = 10)),
        legend.text = element_text(size = 14),
        plot.margin = margin(20, 20, 20, 20),
        #axis.line = element_line(size = 0.8),  # Achsen behalten, aber ohne Pfeile
        panel.border = element_blank())


# Diagramm anzeigen
print(p)

# Diagramm speichern
ggsave(filename = "carbon_pool_comparison_2001_kg_per_m2.png", plot = p, 
       width = 10, height = 6, dpi = 300, bg = "white")