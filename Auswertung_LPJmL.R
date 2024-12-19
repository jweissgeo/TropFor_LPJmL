### R-Skript LPJmL Auswertung ###


### Pakete laden
require(lpjmlkit)
require(raster)
require(caTools)
require(trend)

### local path setzen, andere immer auskommentieren ###
local_path <- "C:/Users/philipp/Documents/TropFor_LPJmL_LokalData/" #philipp
#local_path <- "/Users/epigo/Documents/LPJmL_Lokal/" #Julius

### Start of the script; Einlesen der Daten
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

