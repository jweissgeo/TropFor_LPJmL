require(lpjmlkit)
require(raster)
library(leaflet)
library(dplyr)
require(sp)
require(RColorBrewer)
library(ggplot2)     # Für die Visualisierung
library(ggspatial)   # Für Hintergrundkarten
library(sf)          # Für Koordinatensysteme
library(rnaturalearth) # Weltkarten-Daten

# Lokalen Pfad setzen
local_path <- "/Users/epigo/Documents/LPJmL_Lokal/" # Julius

# Raster Daten laden
raster_data_change = stack("cft_in_tropics_no_australia.tif")

# CFT Raster vorbereiten
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

# Main Raster Stack
raster_data_nochange <- stack(lapply(1:64, function(band) {
  as_raster(subset(cft_in, band = band))
}))

# Schwelle für echte Änderungen definieren
threshold <- 1e-4

# Berechne absolute Differenz
raster_diff <- abs(raster_data_nochange - raster_data_change)

# Erstelle Change-Maske
change_mask <- calc(raster_diff, function(x) ifelse(any(x > threshold, na.rm = TRUE), 1, 0))

# Plot der Maske
plot(change_mask, main = "Veränderungsmaske über alle 64 Bänder")
change_mask[change_mask < 0.9] <- NA  # Setzt 0er auf NA

# Speichern der Maske
writeRaster(change_mask, filename = paste0(local_path, "change_mask.tif"), format = "GTiff", overwrite = TRUE)

# Weltkarte als Hintergrund
world <- ne_countries(scale = "medium", returnclass = "sf")

# Gemeinsame Bounding Box
lon_min <- -110
lon_max <- 55
lat_min <- -25
lat_max <- 25

# Runoff
o_R = read_io(paste0(local_path, "final_oldLU_cru2011_const/runoff.bin.json"))
o_R = transform(o_R, to = "lon_lat")
o_R = transform(o_R, to = "year_month_day")
o_R_raster <- as_raster(o_R, aggregate = list(month = sum, year = mean), na.rm = TRUE)

n_R = read_io(paste0(local_path, "final_newLU_cru2011_const/runoff.bin.json"))
n_R = transform(n_R, to = "lon_lat")
n_R = transform(n_R, to = "year_month_day")
n_R_raster <- as_raster(n_R, aggregate = list(month = sum, year = mean), na.rm = TRUE)

o_R_masked <- mask(o_R_raster, change_mask)
n_R_masked <- mask(n_R_raster, change_mask)

o_R_df <- as.data.frame(o_R_masked, xy = TRUE, na.rm = TRUE)
n_R_df <- as.data.frame(n_R_masked, xy = TRUE, na.rm = TRUE)

change_df_R <- merge(o_R_df, n_R_df, by = c("x", "y"), suffixes = c("_old", "_new"))
change_df_R$change <- change_df_R$runoff_new - change_df_R$runoff_old

# Transpiration
o_tr = read_io(paste0(local_path, "final_oldLU_cru2011_const/transp.bin.json"))
o_tr = transform(o_tr, to = "lon_lat")
o_tr = transform(o_tr, to = "year_month_day")
o_tr_raster <- as_raster(o_tr, aggregate = list(month = sum, year = mean), na.rm = TRUE)

n_tr = read_io(paste0(local_path, "final_newLU_cru2011_const/transp.bin.json"))
n_tr = transform(n_tr, to = "lon_lat")
n_tr = transform(n_tr, to = "year_month_day")
n_tr_raster <- as_raster(n_tr, aggregate = list(month = sum, year = mean), na.rm = TRUE)

o_tr_masked <- mask(o_tr_raster, change_mask)
n_tr_masked <- mask(n_tr_raster, change_mask)

o_tr_df <- as.data.frame(o_tr_masked, xy = TRUE, na.rm = TRUE)
n_tr_df <- as.data.frame(n_tr_masked, xy = TRUE, na.rm = TRUE)

change_df_tr <- merge(o_tr_df, n_tr_df, by = c("x", "y"), suffixes = c("_old", "_new"))
change_df_tr$change <- change_df_tr$transp_new - change_df_tr$transp_old

# Evaporation
o_ev = read_io(paste0(local_path, "final_oldLU_cru2011_const/evap.bin.json"))
o_ev = transform(o_ev, to = "lon_lat")
o_ev = transform(o_ev, to = "year_month_day")
o_ev_raster <- as_raster(o_ev, aggregate = list(month = sum, year = mean), na.rm = TRUE)

n_ev = read_io(paste0(local_path, "final_newLU_cru2011_const/evap.bin.json"))
n_ev = transform(n_ev, to = "lon_lat")
n_ev = transform(n_ev, to = "year_month_day")
n_ev_raster <- as_raster(n_ev, aggregate = list(month = sum, year = mean), na.rm = TRUE)

o_ev_masked <- mask(o_ev_raster, change_mask)
n_ev_masked <- mask(n_ev_raster, change_mask)

o_ev_df <- as.data.frame(o_ev_masked, xy = TRUE, na.rm = TRUE)
n_ev_df <- as.data.frame(n_ev_masked, xy = TRUE, na.rm = TRUE)

change_df_ev <- merge(o_ev_df, n_ev_df, by = c("x", "y"), suffixes = c("_old", "_new"))
change_df_ev$change <- change_df_ev$evap_new - change_df_ev$evap_old


### Plots

# Runoff Plot
ggplot() +
  geom_tile(data = change_df_R, aes(x = x, y = y, fill = pmin(change, 1000))) +
  scale_fill_viridis(option = "rocket", name = "Runoff Change", direction = -1, limits = c(min(change_df_R$change, na.rm = TRUE), 1000)) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  labs(
    title = "Änderung des Runoff (Neotropis & Afrotropis)",
    x = "Längengrad", y = "Breitengrad"
  ) +
  theme_minimal()


# Transpiration Plot
ggplot() +
  geom_tile(data = change_df_tr, aes(x = x, y = y, fill = change)) +
  scale_fill_viridis(option = "rocket", name = "Transp. Change", "direction = -1") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  labs(
    title = "Änderung der Transpiration (Neotropis & Afrotropis)",
    x = "Längengrad", y = "Breitengrad"
  ) +
  theme_minimal()

# Evaporation Plot
ggplot() +
  geom_tile(data = change_df_ev, aes(x = x, y = y, fill = change)) +
  scale_fill_viridis(option = "rocket", name = "Evaporation Change", direction = -1) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  labs(
    title = "Änderung der Evaporation (Neotropis & Afrotropis)",
    x = "Längengrad", y = "Breitengrad"
  ) +
  theme_minimal()
# Lade das benötigte Paket
library(ggplot2)
library(patchwork)

library(ggplot2)
library(viridis)
library(sf)
library(patchwork)



# Runoff Plot
p1 <- ggplot() +
  geom_tile(data = change_df_R, aes(x = x, y = y, fill = pmin(change, 1000))) +
  scale_fill_viridis(option = "rocket", name = "Mean Δ Surface Runoff (mm/year)", direction = -1, 
                     limits = c(min(change_df_R$change, na.rm = TRUE), 1000)) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  labs(title = "Change in Surface Runoff") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

# Transpiration Plot
p2 <- ggplot() +
  geom_tile(data = change_df_tr, aes(x = x, y = y, fill = change)) +
  scale_fill_viridis(option = "viridis", name = "Mean Δ Transpiration (mm/year)", direction = -1) +
  geom_sf(data = world, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  labs(title = "Change in Transpiration") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

# Evaporation Plot
p3 <- ggplot() +
  geom_tile(data = change_df_ev, aes(x = x, y = y, fill = change)) +
  scale_fill_viridis(option = "viridis", name = "Mean Δ Evaporation (mm/year)", "direction = -1") +
  geom_sf(data = world, fill = NA, color = "black", size = 0.3) +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  labs(title = "Change in Evaporation") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

# Plots übereinander anordnen
(p2 / p3 / p1)  # Transpiration oben, dann Evaporation, dann Runoff


# Kombiniere die Daten zu einem Dataframe
change_df <- data.frame(
  Variable = rep(c("Transpiration", "Evaporation", "Runoff"), 
                 times = c(nrow(change_df_tr), nrow(change_df_ev), nrow(change_df_R))),
  Change = c(change_df_tr$change, change_df_ev$change, change_df_R$change)
)
library(ggplot2)

# Kombiniere die Daten zu einem Dataframe
change_df <- data.frame(
  Variable = rep(c("Transpiration", "Evaporation", "Runoff"), 
                 times = c(nrow(change_df_tr), nrow(change_df_ev), nrow(change_df_R))),
  Change = c(change_df_tr$change, change_df_ev$change, change_df_R$change)
)
# Kombiniere die Daten zu einem Dataframe
change_df <- data.frame(
  Variable = rep(c("Transpiration", "Evaporation", "Runoff"), 
                 times = c(nrow(change_df_tr), nrow(change_df_ev), nrow(change_df_R))),
  Change = c(change_df_tr$change, change_df_ev$change, change_df_R$change)
)

# Statistische Tests (t-Test gegen 0)
test_results <- change_df %>%
  group_by(Variable) %>%
  summarise(
    Mean = mean(Change, na.rm = TRUE),
    SD = sd(Change, na.rm = TRUE),
    t_value = t.test(Change, mu = 0)$statistic,
    p_value = t.test(Change, mu = 0)$p.value
  )

# Statistische Tests (t-Test gegen 0)
test_results <- change_df %>%
  group_by(Variable) %>%
  summarise(
    Mean = mean(Change, na.rm = TRUE),
    SD = sd(Change, na.rm = TRUE),
    t_value = t.test(Change, mu = 0)$statistic,
    p_value = t.test(Change, mu = 0)$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Erstelle den Boxplot mit verschiedenen Farben
p <- ggplot(change_df, aes(x = Variable, y = Change, fill = Variable)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2) +  # Klassischer Boxplot mit Whiskern und Ausreißern
  scale_fill_manual(values = c("Transpiration" = "turquoise3", "Evaporation" = "gold", "Runoff" = "forestgreen")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Referenzlinie
  labs(y = "Δ (mm/year)", x = "", title = "Changes in Transpiration, Evaporation and Runoff") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",  # Keine separate Legende
    plot.title = element_text(hjust = 0.5)  # Überschrift mittig ausrichten
  )

# Anzeige des Plots
print(p)

# Ausgabe der Testergebnisse
print(test_results)
