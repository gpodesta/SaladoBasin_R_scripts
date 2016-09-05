# Leer netCDF con series sinteticas de variables climaticas

Sys.setenv(TZ = "UTC") 
options(warn = -1)

# -----------------------------------------------------------------------------#
# --- Erase all objects before starting ----

rm(list = ls()); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Load necessary R packages ----

# --- List of packages to download and load...

list.of.packages <- c("ADGofTest", "boot", "Cairo", "DBI",  
  "Evapotranspiration", "fitdistrplus",
  "futile.logger", "geosphere", "Hmisc",
  "insol", "lattice", "locfit", "lubridate", "MASS", "maps", "mgcv",
  "mice", "ncdf4", "nortest", "RColorBrewer", "raster", "rasterVis",
  "reshape2", "rgdal","rgl","RPostgreSQL", "rts", "sp", "spdep",
  "SPEI", "stringr", "xts", "yaml", "zoo", "dplyr")    

for (pack in list.of.packages) {
  if (!require(pack, character.only = TRUE)) {
    install.packages(pack, dependencies = TRUE)
    require(pack, character.only = TRUE)
  }
}

library(ncdf4)  # This library currently needs to be installed from zip file

rm(list.of.packages, pack); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Define CRS strings for geographic projections -----

# --- CRS for Gauss-Kruger projection used by Argentina

GK.string <- "+proj=tmerc +lat_0=-90 +lon_0=-60 +k=1
+x_0=5500000 +y_0=0 +ellps=intl
+twogs84=-148,136,90,0,0,0,0
+units=m +no_defs"

# Latitude of natural origin +lat_0=-90
# Longitude of natural origin	+lon_0=-60
# Scale factor at natural origin +k=1
# False Easting	+x_0=5500000
# False Northing +y_0=0
# ref.system : Gauss-Kruger, Zone 5, Campo Inchauspe projection : Gauss-Kruger
# datum : Campo Inchauspe delta WGS84 : -148 136 90 ellipsoid : International 1924
# major s-ax : 6378388 minor s-ax : 6356912 origin long : -60 origin lat : -90 
# origin X : 5500000 origin Y : 0 scale fac : 1.0 units : m parameters : 0

# --- CRS for lat-lon projection

ll.string <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Read outline of Salado Basin A (MIKE active cells)  ----

dir.basin <- 'D:/ENSO/Projects/CNH3/data/shapes_cuenca'
#dir.basin <- 'C:/Data/ENSO/Projects/CNH3/data/shapes_cuenca'

tt0 <- file.info(dir.basin)   # Info about directory where basin shapefile is located
if (!tt0$isdir)
  stop("ERROR: Specified directory for basin shapefile does not exist... check name...\n")              

# --- Read shapefile with Basin A boundary

cuenca.A.GK <- rgdal::readOGR(dir.basin,
  layer="Salado_A_MIKE_boundary_GK",
  verbose = TRUE)

# --- Plot boundary of Salado A Basin

plot(cuenca.A.GK, axes = TRUE,
  lwd = 2, border = '#2ca25f',
  col = '#99d8c9')

rm(dir.basin, tt0); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Read MIKE-SHE simulation grid for Salado Basin A ----

dir.grid <- 'D:/ENSO/Projects/CNH3/data/grillas_modelo'
#dir.grid <- 'C:/Data/ENSO/Projects/CNH3/data/grillas_modelo'

tt0 <- file.info(dir.grid)   # Informacion sobre el directorio especificado
if (!tt0$isdir)
  stop("ERROR: The specified grid directory does not exist... check dir name...\n")              

# --- Read simulation grid used for MIKE-SHE 

grilla.5000.GK <- rgdal::readOGR(dir.grid,
  layer="Salado_A_grilla_5000_GK",
  verbose = TRUE)

# --- Plot both GK and lat-lon grids

plot(grilla.5000.GK, pch = 16, col = "steelblue", axes = TRUE, cex = 0.2)
plot(cuenca.A.GK, add = TRUE, lwd = 2, border = "tomato")

rm(dir.grid, tt0); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Open netCDF file with synthetic weather series ----

dir.climate <- 'D:/ENSO/Projects/CNH3/data/WeatherSeries'
#dir.climate <- 'C:/Data/ENSO/Projects/CNH3/data/WeatherSeries'

tt0 <- file.info(dir.climate)   # Info about directory where synth weather is located
if (!tt0$isdir)
  stop("ERROR: Specified directory for synth weather does not exist... check name...\n")              

file.climate <- paste0(dir.climate, '/Weather-Master.nc')

if (!file.exists(file.climate)) {
  stop('Specified netCDF file with synth weather does NOT exist...\n')
}

# --- Open netCDF file with synthetic weather and get start date 

nc <- ncdf4::nc_open(filename = file.climate,
  write = FALSE, verbose = TRUE)

rm(dir.climate, tt0, file.climate);gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Extract information from netCDF file ----

# --- Get all global attributes

glob.attr <- ncdf4::ncatt_get(nc, varid = 0)

# --- Get values of x and y coordinates

names(nc$dim$projection_x_coordinate)

x_coord_vals <- nc$dim$projection_x_coordinate$vals   # Values of x coordinates (longitude)
n.x_coord_vals <- nc$dim$projection_x_coordinate$len  # Number of x coordinates

y_coord_vals <- nc$dim$projection_y_coordinate$vals   # Values of y coordinates (latitude)
n.y_coord_vals <- nc$dim$projection_y_coordinate$len  # Number of x coordinates

# --- Get values of time coordinate and convert to "date" object
# --- First we extract the origin of time units from "units" component
# --- of time coordinate.
# --- The following code assumes origin is expressed as "days since YYYY-MM-DD".

tt1 <- stringr::str_locate(nc$dim$time$units, "since ")
origin.date <- stringr::str_sub(nc$dim$time$units,
  start = tt1[, "end"] + 1, end = -1L)

if (origin.date != glob.attr$start_date) {
  warning('Start date specified as global attribute
    does not coincide with description of time variable')
}

time_coord_vals <- nc$dim$time$vals
n.time_coord_vals <- length(time_coord_vals)

date_coords <- as.Date(nc$dim$time$vals, origin = origin.date)

rm(tt1, origin.date)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Extract values of weather variables ----  

prcp.matrix <- array(dim = c(n.y_coord_vals, n.x_coord_vals,
  n.time_coord_vals))

tmax.matrix <- array(dim = c(n.y_coord_vals, n.x_coord_vals,
  n.time_coord_vals))

tmin.matrix <- array(dim = c(n.y_coord_vals, n.x_coord_vals,
  n.time_coord_vals))

for (i in 1:n.time_coord_vals) {

  prcp.layer <- t(
    ncdf4::ncvar_get(nc, varid = "prcp",
    start = c(1, 1, 1, i),
    count = c(1, n.x_coord_vals, n.y_coord_vals, 1),
    verbose = FALSE) )
  
  tmax.layer <- t(
    ncdf4::ncvar_get(nc, varid = "tmax",
      start = c(1, 1, 1, i),
      count = c(1, n.x_coord_vals, n.y_coord_vals, 1),
      verbose = FALSE) )
  
  tmin.layer <- t(
    ncdf4::ncvar_get(nc, varid = "tmin",
      start = c(1, 1, 1, i),
      count = c(1, n.x_coord_vals, n.y_coord_vals, 1),
      verbose = FALSE) )
  
  prcp.matrix[,, i] <- prcp.layer
  tmax.matrix[,, i] <- tmax.layer
  tmin.matrix[,, i] <- tmin.layer
  
}
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Build raster bricks for each weather variable ----

half_delta_x <- (x_coord_vals[2] - x_coord_vals[1]) / 2
half_delta_y <- (y_coord_vals[2] - y_coord_vals[1]) / 2

x.min <- min(x_coord_vals) - half_delta_x
x.max <- max(x_coord_vals) + half_delta_x
y.min <- min(y_coord_vals) - half_delta_y
y.max <- max(y_coord_vals) + half_delta_y

prcp.brick <- raster::brick(prcp.matrix,
  crs = sp::CRS(GK.string),
  xmn = x.min,
  xmx = x.max,
  ymn = y.min,
  ymx = y.max)

tmax.brick <- raster::brick(tmax.matrix,
  crs = sp::CRS(GK.string),
  xmn = x.min,
  xmx = x.max,
  ymn = y.min,
  ymx = y.max)

tmin.brick <- raster::brick(tmin.matrix,
  crs = sp::CRS(GK.string),
  xmn = x.min,
  xmx = x.max,
  ymn = y.min,
  ymx = y.max)

# --- Mask cells outside the Salado A Basin

prcp.brick2 <- raster::mask(prcp.brick, cuenca.A.GK)
tmax.brick2 <- raster::mask(tmax.brick, cuenca.A.GK)
tmin.brick2 <- raster::mask(tmin.brick, cuenca.A.GK)

# --- Set z-coordinate as dates for masked rasters

prcp.brick2 <- raster::setZ(prcp.brick2, date_coords, 'dates')
tmax.brick2 <- raster::setZ(tmax.brick2, date_coords, 'dates')
tmin.brick2 <- raster::setZ(tmin.brick2, date_coords, 'dates')
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Plot a few slices of bricks ----

p <- rasterVis::levelplot(prcp.brick2, layer = 1, margin = FALSE,
  main = paste('Prcp - ', as.character(raster::getZ(prcp.brick2)[1])))
p <- p + latticeExtra::layer(sp::sp.lines(cuenca.A.GK, lwd = 2,
  col = 'steelblue'))
plot(p)

p <- rasterVis::levelplot(tmax.brick2, layer = 1, margin = FALSE,
  main = paste('Tmax -', as.character(raster::getZ(tmax.brick2)[1])))
p <- p + latticeExtra::layer(sp::sp.lines(cuenca.A.GK, lwd = 2,
  col = 'steelblue'))
plot(p)

p <- rasterVis::levelplot(tmin.brick2, layer = 1, margin = FALSE,
  main = paste('Tmin -', as.character(raster::getZ(tmin.brick2)[1])))
p <- p + latticeExtra::layer(sp::sp.lines(cuenca.A.GK, lwd = 2,
  col = 'steelblue'))
plot(p)
# ------------------------------------------------------------------------------


# --- Calculate monthly mean of Tmax 

uu3 <- raster::zApply(tmax.brick2,
  by = as.yearmon,
  fun = mean,
  name = 'Monthly Tmax means')

# --- Calculate monthly mean of Tmin 

uu3b <- raster::zApply(tmin.brick2,
  by = as.yearmon,
  fun = mean,
  name = 'Monthly Tmin means')

# --- Calculate monthly total of prcp 

uu4 <- raster::zApply(prcp.brick2,
  by = as.yearmon,
  fun = sum,
  name = 'Monthly Prcp Totals')

uu4 <- raster::mask(uu4, cuenca.A.GK) # Mask basin again

ttt <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))

lll <- 18 # Monthly slice
limites.valores <- seq(10, 180, 10)

p <- rasterVis::levelplot(uu4, layer = lll, margin = FALSE,
  main = paste('Prcp Total -', as.character(raster::getZ(uu4)[lll])),
  at = c(0, 1, limites.valores),
  col.regions = c('grey20', ttt(length(limites.valores) + 1)))
p <- p + latticeExtra::layer(sp::sp.lines(cuenca.A.GK, lwd = 2,
  col = 'tomato'))
plot(p)




# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Build RTS objects ----

prcp.rts <- rts::rts(x = prcp.brick2, time = date_coords)
tmax.rts <- rts::rts(x = tmax.brick2, time = date_coords)
tmin.rts <- rts::rts(x = tmin.brick2, time = date_coords)

uu5 <- rts::apply.monthly(x = prcp.rts, FUN = sum, na.rm = TRUE)

plot(x = tmax.rts)

# ------------------------------------------------------------------------------












# --- Read water table depth from netCDF file into a RasterBrick

wtd.brick <- raster::brick(file.to.read,
  values = TRUE,
  crs = sp::CRS(GK.string),
  varname = 'wtd')

# --- Read times as Z coordinates (or layers)

times <- raster::getZ(wtd.brick)  # Get back z-coordinates (time)
wtd.dates <- as.Date(times, origin = as.Date(origin.date))
wtd.brick <- raster::setZ(wtd.brick, wtd.dates, 'wtd.dates')

p <- rasterVis::levelplot(wtd.brick, layer = 100, margin = FALSE,
  main = as.character(getZ(wtd.brick)[100]))
p <- p + layer(sp.lines(cuenca.A.GK, lwd = 2, col = 'steelblue'))
plot(p)


# --- Mask cells outside the Salado A Basin
# --- Apparently the masking operation deleted the Z values as dates

wtd.brick2 <- raster::mask(wtd.brick, cuenca.A.GK)
wtd.brick2 <- raster::setZ(wtd.brick2, wtd.dates, 'wtd.dates')

getZ(wtd.brick2)



p <- rasterVis::levelplot(wtd.brick2, layer = 120, margin = FALSE,
  main = as.character(getZ(wtd.brick2)[120]))
p <- p + layer(sp.lines(cuenca.A.GK, lwd = 2, col = 'steelblue'))
plot(p)


# --- Calculate monthly mean of WTD

uu3 <- raster::zApply(wtd.brick2,
  by = as.yearmon,
  fun = mean,
  name = 'Monthly WTD means')

# --- Calculate monthly median of WTD

uu3 <- raster::zApply(wtd.brick2,
  by = as.yearmon,
  fun = median, na.rm = TRUE,
  name = 'Monthly WTD medians')


ttt <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))

lll <- 18

p <- rasterVis::levelplot(uu3, layer = lll, margin = FALSE,
  main = as.character(getZ(uu3)[lll]),
  at = c(0.25, seq(0, -5, -0.1)),
  col.regions = c(rev(ttt(51)), 'tomato'))
p <- p + layer(sp.lines(cuenca.A.GK, lwd = 2, col = 'tomato'))
plot(p)




ncdf4::nc_close(nc)



