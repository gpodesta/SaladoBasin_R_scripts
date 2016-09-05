# Read netCDF with GPCC Global Precipitation Data

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
  "insol", "lattice", "latticeExtra",
  "locfit", "lubridate", "MASS", "maps", "mgcv",
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
# --- Read outline of CRC-SAS region  ----

dir.basin <- 'D:/ENSO/Projects/IDB/scripts/R/datos_adicionales/shapefiles'

tt0 <- file.info(dir.basin)   # Info about directory where basin shapefile is located
if (!tt0$isdir)
  stop("ERROR: Specified directory for basin shapefile does not exist... check name...\n")              

# --- Read shapefile with Basin A boundary

CRC.outline <- rgdal::readOGR(dir.basin,
  layer="borde_externo_CRC_SAS",
  verbose = TRUE)

# --- Plot boundary of Salado A Basin

plot(CRC.outline, axes = TRUE,
  lwd = 2, border = '#2ca25f',
  col = '#99d8c9')

rm(dir.basin, tt0); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Open netCDF file with global precipitation ----

dir.gpcc <- 'D:/ENSO/Projects/CNH3/data/gpcc_precipitation'

tt0 <- file.info(dir.gpcc)   # Info about directory where synth weather is located
if (!tt0$isdir)
  stop("ERROR: Specified directory for synth weather does not exist... check name...\n")              

file.gpcc <- paste0(dir.gpcc, '/precip.mon.combined.total.v7.nc')

if (!file.exists(file.gpcc)) {
  stop('Specified netCDF file with synth weather does NOT exist...\n')
}

# --- Open netCDF file with synthetic weather and get start date 

nc <- ncdf4::nc_open(filename = file.gpcc,
  write = FALSE, verbose = TRUE)

rm(dir.gpcc, tt0);gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Extract information from GPCC netCDF file ----

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
# --- Create a raster brick with the global GPCC data ----

gpcc.brick <- raster::brick(x = file.gpcc,
  values = TRUE,
  crs = sp::CRS(ll.string))

# --- Longitudes in original data have range [0-360] (degree east).
# --- Change them to range [-180, +180]

gpcc.brick2 <- raster::rotate(x = gpcc.brick)

rm(gpcc.brick); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Crop global GPCC data to CRC-SAS area ----

# Construir un objeto 'extent' para seleccionar area de CRC.
# El extent se construye a partir de las latitudes y longitudes
# minimas/maximas del area del CRC-SAS, redondeadas para dar valores 'lindos'.

tt1 <- bbox(CRC.borde)
tt2 <- range(pretty(tt1['x', ], 20))
tt3 <- range(pretty(tt1['y', ], 20))
CRC.extent <- raster::extent(list(x = tt2, y = tt3))

# Use the extent created above to crop the global data

gpcc.brick3 <- raster::crop(x = gpcc.brick2, y = CRC.extent)

# Finally, mask the data using the CRC-SAS external contour

gpcc.brick <- raster::mask(x = gpcc.brick3, mask = CRC.outline)

rm(gpcc.brick2, gpcc.brick3,tt1,tt2,tt3,CRC.extent)  
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Visualize slices of GPCC data for CRC area ----

colors <- RColorBrewer::brewer.pal(9, 'Blues')
n.colors <- length(colors)

my.at <- seq(from = 0, to = 400, length.out = n.colors)

myColorkey <- list(space = 'right',
  at = my.at, ## where the colors change
  col = colors,
  labels = list(
      labels = as.character(my.at),
      at = my.at ## where to print labels
  )
)

use.layer <- 1000  # Brick layer to plot

p <- rasterVis::levelplot(gpcc.brick,
  layer = use.layer,
  margin = FALSE,
  colorkey = myColorkey,
  col.regions = colors,
  contour = FALSE,
  main = paste('GPCC Prcp - ', as.character(raster::getZ(gpcc.brick)[use.layer])))

p <- p + latticeExtra::layer(sp::sp.lines(CRC.outline,
  lwd = 2,
  col = 'tomato'))

plot(p)
# ------------------------------------------------------------------------------


ncdf4::nc_close(nc)
