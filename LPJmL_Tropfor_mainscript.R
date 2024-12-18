# R-Skript LPJmL Auswertung

#test#
require(lpjmlkit)
require(raster)
require(caTools)
require(trend)

local_path <- "/Users/epigo/Documents/LPJmL_Lokal/"
globalflux <- read.csv(paste0(local_path, "gampe_baseline/globalflux.csv"), sep = ",")

plot(globalflux$GPP,col="darkgreen",t="l")
long_term_gpp= runmean(globalflux$GPP,20,endrule = "mean",align="center")
lines(long_term_gpp,col="black",lwd=2)

## Mann-Kendall Trendbestimmung
trend_mk = mk.test(long_term_gpp,alternative = "greater")

mk_vec = c()
for(i in 10:(length(long_term_gpp)-9)){
  act_vec = long_term_gpp[c((i-9):(i+9))]   
  mk_vec[i] = mk.test(act_vec,alternative="greater")$p.value
}
trend.col=ifelse(mk_vec <=0.05,"grey90","black")
points(long_term_gpp,col=trend.col)


### spatial binary files
## "tutorial": 
vignette("lpjml-data")

## read meta data  
read_meta(paste0(local_path, "gampe_baseline/mgpp.bin.json"))

## read GPP grid
## braucht: *.bin.json, *.bin, grid.bin, grid.bin.json
gpp = read_io(paste0(local_path, "gampe_baseline/mgpp.bin.json"))
gpp = transform(gpp,to = "lon_lat")  
gpp = transform(gpp, to ="year_month_day")

gpp =subset(gpp, month = as.character(c(1,2,12)))  
gpp =subset(gpp,lat = as.character(seq(-0.25,-60.75,-0.5)))

plot(subset(gpp,year = as.character(c(2024))))

ref = subset(gpp,year = as.character(1901:1950))
fut = subset(gpp,year = as.character(2071:2100))


ref_mean = as_raster(ref,
                     aggregate = list(month=sum, year = mean,band=1),
                     na.rm=T
)

fut_mean = as_raster(fut,
                     aggregate = list(month=sum, year = mean,band=1),
                     na.rm=T
)

gpp_dif = fut_mean - ref_mean

col.bar = colorRampPalette(c("brown","white","darkgreen"))

spplot(gpp_dif,col.regions = col.bar(21),
       zlim =c(-850,850),scales = list(draw=T),
       xlab="Lon (°)",ylab="Lat (°)")

##### now as line plot...
gpp_vec = as_array(gpp,
                   aggregate = list(
                     lon=mean,lat=mean,month=sum
                   ),na.rm=T
) 

plot(gpp_vec,col="darkgreen",bty="n",
     xlab="Year",ylab="GPP (gC/m2/season)",t="l")

bliblablubb = lm(gpp_vec~c(1:200))
abline(bliblablubb,col="black")

## alternativ: direkt
plot(gpp,
     aggregate= list(
       lon = mean, lat = mean, month=sum
     )
)

gpp = read_io(paste0(local_path, "gampe_baseline/mgpp.bin.json"),
              subset = list(year = as.character(1971:1990))
)

vignette("lpjml-data")


### now land use input

cft_header_in = read_header(paste0(local_path, "gampe_baseline/cft1700_2005_irrigation_systems_64bands.bin")) # LPJmL lokal

cft_in = read_io(paste0(local_path, "gampe_baseline/cft1700_2005_irrigation_systems_64bands.bin"),
                name= "cft1700_2005_irrigation_systems_64bands.bin",
                descr = "LPJLUSE",
                firstcell = 0, ncell = 67420, 
                firstyear = 1700,nyear =306,
                scalar = 0.001, nbands= 64, 
                datatype = 1, order = 1,
                cellsize_lat = 0.5,cellsize_lon = 0.5,
                nstep = 1,endian = "little",
                subset=list(year=306)
)

# Crop Functional Types (CFTs):
# 1. temperate cereals
# 2. rice
# 3. maize
# 4. tropical cereals
# 5. pulses
# 6. temperate roots
# 7. tropical roots
# 8. oil crops sunflower
# 9. oil crops soybean
# 10. oil crops groundnut
# 11. oil crops rapeseed
# 12. sugar cane
# 13. other crops
# 14. pasture/managed grass (C3/C4 mixed)
# 15. bio-energy grass
# 16. bio-energy tree

act_band = 12
lu_in = transform(cft_in,to ="lon_lat")
act_y = subset(lu_in, band = act_band)
lpjml_grid = read_io(paste0(local_path, "gampe_baseline/grid.bin"))

act_ras =as_raster(act_y)
spplot(act_ras)
writeRaster(act_ras,filename = paste0(local_path, "gampe_baseline/cft_2005_band12.tif"),
            format ="GTiff", overwrite = TRUE) # Achtung Overwrite = True !!!!

lon.def = seq(-179.75, 179.75, 0.5)
lat.def = seq(-23.25, 23.25, 0.5)
act_coordbox = expand.grid(lon.def,lat.def)
act_cell = cellFromXY(act_ras,act_coordbox)

act_ras[act_cell][!is.na(act_ras[act_cell])] = 1

spplot(act_ras)
writeRaster(act_ras, paste0(local_path, "gampe_baseline/band12_modified2.tif"), overwrite = TRUE) # Achtung Overwrite = True !!!

##### break woche 5 -6

raster_data2 = raster(paste0(local_path, "gampe_baseline/band12_modified2.tif"))
grid_lpjml = read_io(paste0(local_path, "gampe_baseline/grid.bin.json"))

grid_coords = paste(round(grid_lpjml$data[,,1],2), round(grid_lpjml$data[,,2],2),sep="_")
ras_coords  = paste(round(coordinates(raster_data2)[,1],2),round(coordinates(raster_data2)[,2],2),sep="_")

cft_out = array(0,dim=(c(1,67420,64)))
y = 1

### Convert Landuse grid (in R) to format needed in LPJmL
### R: cell - year - band
### LPJmL: year - cell -band
for(i in 1:67420){
  for(b in 1:64){
    cft_out[y,i,b] = as.numeric(cft_in$data[i,y,b])
    # names(cft_out[y,,b]) = names(cft_in$data[i,,b])
  }
}

### overlay cft_out with new edited values from raster
### use matching coordinates to ensure correct placement
### act_band: band of land use class to be edited (=the raster)

### much faster: write raster data as vector first:
ras.vec = as.numeric(raster_data2[])

for(i in 1:length(ras_coords)){
  if(is.na(match(ras_coords[i],grid_coords))){
    next
  }else{
    cft_out[y,match(ras_coords[i],grid_coords),act_band] = ras.vec[i]
  }
}

### concatenate land use file
### use raster-year as new year (y = 307)
f.out <- file(paste0(local_path, "gampe_baseline/cft_intmod24_new1.bin","wb"))
n_year_out = 307
for(y in 1:n_year_out){
  print(y)
  if(y<n_year_out){
    cft_in_tmp =read_io(paste0(local_path, "gampe_baseline/cft1700_2005_irrigation_systems_64bands.bin"),
                        name= "cft1700_2005_irrigation_systems_64bands.bin",
                        descr = "LPJLUSE",
                        firstcell = 0, ncell = 67420, 
                        firstyear = 1700,nyear =306,
                        scalar = 0.001, nbands= 64, 
                        datatype = 1, order = 1,
                        cellsize_lat = 0.5,cellsize_lon = 0.5,
                        nstep = 1,endian = "little",
                        subset=list(year=y))
    for(i in 1:67420){
      writeBin(as.integer(round(cft_in_tmp$data[i,1,]*1000)),f.out,size = 2, endian="little",useBytes=TRUE)
    }  
  }else{
    ### year 307: edited year (from raster)
    for(i in 1:67420){
      writeBin(as.integer(round(cft_out[1,i,]*1000)),f.out,size = 2,endian = "little",useBytes=TRUE)
    }
  }
}
close(f.out)

### write function to merge 2 binary files
merge_binary_file =function(file1,file2,output_file){
  f1 = readBin(file1,what = "raw", n = file.info(file1)$size)
  f2 = readBin(file2,what = "raw", n = file.info(file2)$size)
  out_f = c(f1,f2)
  writeBin(out_f,output_file)
}


### create header for land use file
header_cft_out = create_header(
  name = "LPJLUSE",
  version = 2,
  order = 1,
  firstyear = 1700,
  nyear = 307,
  firstcell = 0,
  ncell = 67420,
  nbands = 64,
  cellsize_lon = 0.5,
  scalar = 0.001,
  cellsize_lat = 0.5,
  datatype = 1,
  nstep = 1,
  timestep = 1,
  endian = "little",
  verbose = TRUE
)

write_header(filename = paste0(local_path, "gampe_baseline/landuse_header4.bin"),header = header_cft_out, overwrite=TRUE)
merge_binary_file(paste0(local_path, "gampe_baseline/landuse_header4.bin"), paste0(local_path, "gampe_baseline/cft_intmod24_new1.bin"), paste0(local_path, "gampe_baseline/cft_intmod24_final2.bin") )

header_cft_new = read_header(paste0(local_path, "gampe_baseline/cft_intmod24_final2.bin"))

cft_new =read_io(paste0(local_path, "gampe_baseline/cft_intmod24_final2.bin"),
                 name= "cft_intmod24_final2.bin",
                 descr = "LPJLUSE",
                 firstcell = 0, ncell = 67420, 
                 firstyear = 1700,nyear =307,
                 scalar = 0.001, nbands= 64, 
                 datatype = 1, order = 1,
                 cellsize_lat = 0.5,cellsize_lon = 0.5,
                 nstep = 1,endian = "little",
                 subset=list(year=307))

plot(subset(cft_new,band=12))

#################################
### Cropland Expansion Raster ###

# Raster laden
Expansion_potential <- raster("Expansion_Potential/Figure_S3/integrated_expansion_potential.bil")

# Extent für die Tropen (zwischen 23.25° S und 23.25° N)
tropen_extent <- extent(min(lon.def), max(lon.def), min(lat.def), max(lat.def))

# Raster auf die Tropen beschneiden
trop_Expansion_potential <- crop(Expansion_potential, tropen_extent)

# Kontinuierliche Farbpalette von Dunkelgrün über Gelb nach Rot
col_palette <- colorRampPalette(c("darkgreen", "yellow", "red"))

# Plotten des zugeschnittenen Rasters
plot(trop_Expansion_potential,
     main = "Agricultural Expansion Potential (Tropen) (Zabel 2024)",
     col = col_palette(256),
     legend = TRUE)

