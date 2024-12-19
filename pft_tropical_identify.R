### Versuch das PFT tropical forest anzuzeigen
### Öffnen und Anzeigen des Modell-Outputs

# Aus R-skrip_LPJmL_Auswertung.R, der Anfang

# Vorsicht lai_trop_forest und vegc_trop_forest sind Ergebnisse der Modellierung ab 2005. Nochmal ordentlich machen ab 1901 (mache ich (Mareike))

### Pakete laden
require(lpjmlkit)
require(raster)
require(caTools)
require(trend)

### local path setzen, andere immer auskommentieren ###
#local_path <- "C:/Users/philipp/Documents/TropFor_LPJmL_LokalData/" #philipp
#local_path <- "/Users/epigo/Documents/LPJmL_Lokal/" #Julius
local_path <- "C:/Dokumente/Umweltsysteme/integrierte_modellierung/"  # Mareike

## Entwurf 1  (LAI)
read_meta(paste0(local_path, "Model_Output/test_varpft/lai_trop_forest.bin.json"))                  # Metadaten anzeigen
lai_trop = read_io(paste0(local_path, "Model_Output/test_varpft/lai_trop_forest.bin.json"))         # Datei einlesen
lai_trop = transform(lai_trop, to = "lon_lat")          # räumliches
lai_trop = transform(lai_trop, to = "year_month_day")   # zeitliches

#lai_trop = subset(lai_trop, lat = as.character(seq(-23.25,23.25,0.5)))     # nur der Tropenbereich
lai_trop = subset(lai_trop, year = as.character(c(2006)), month = as.character(c(7)), band = as.character("tropical broadleaved evergreen tree"))    # nur 2006

plot(lai_trop, main = "LAI of 'tropical broadleaved evergreen tree' [m^2/m^2]")


## Entwurf 2  (VegC)
read_meta(paste0(local_path, "Model_Output/test_varpft/vegc_trop_forest.bin.json"))                    # Metadaten anzeigen
vegc_trop = read_io(paste0(local_path, "Model_Output/test_varpft/vegc_trop_forest.bin.json"))          # Datei einlesen
vegc_trop = transform(vegc_trop, to = "lon_lat")          # räumliches
vegc_trop = transform(vegc_trop, to = "year_month_day")   # zeitliches

#lai_trop = subset(vegc_trop, lat = as.character(seq(-23.25,23.25,0.5)))     # nur der Tropenbereich
vegc_trop = subset(vegc_trop, year = as.character(c(2006)), band = as.character("tropical broadleaved evergreen tree"))    # nur 2006

plot(vegc_trop, main = "VegC of 'tropical broadleaved evergreen tree' [gC/m^2]")


###################################################################################
# cell_area gibt an wie viele m^2 jede Zelle (0,5°x0,5°) hat (hier zwischen 25 u. 25,5 Längengrad, da dort überall Land ist in den Tropen)
one_meridian = subset(lai_trop,lat = as.character(seq(-23.25,23.25,0.5)), lon = as.character(25.25))
cell_area = calc_cellarea(one_meridian)
