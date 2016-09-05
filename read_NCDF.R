
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
  "mice", "ncdf4", "nortest", "RColorBrewer", "raster", "rasterVis", "reshape2",
  "rgdal","rgl","RPostgreSQL", "sp", "spdep", "SPEI", "stringr", "xts",
  "yaml", "zoo", "dplyr")    

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

tt0 <- file.info(dir.basin)   # Info about directory where basin shapefile is located
if (!tt0$isdir)
  stop("ERROR: Specified directory for basin shapefile does not exist... check name...\n")              

# --- Read shapefile with Basin A boundary

cuenca.A.GK <- rgdal::readOGR(dir.basin, layer="Salado_A_MIKE_boundary_GK",
  verbose = TRUE)

setwd(proj.dir)

# --- Plot boundary of Salado A Basin

plot(cuenca.A.GK, axes = TRUE)

rm(dir.basin, tt0); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Read MIKE-SHE simulation grid for Salado Basin A ----

dir.grid <- 'D:/ENSO/Projects/CNH3/data/grillas_modelo'

tt0 <- file.info(dir.grid)   # Informacion sobre el directorio especificado
if (!tt0$isdir)
  stop("ERROR: The specified grid directory does not exist... check dir name...\n")              

# --- Read simulation grid used for MIKE-SHE 

grilla.5000.GK <- rgdal::readOGR(dir.grid, layer="Salado_A_grilla_5000_GK",
  verbose = TRUE)

# --- Plot both GK and lat-lon grids

plot(grilla.5000.GK, pch = 16, col = "steelblue", axes = TRUE, cex = 0.2)
plot(cuenca.A.GK, add = TRUE, lwd = 2, border = "tomato")

rm(dir.grid, tt0); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Open netCDF file ----

# --- Define netCDF file to be read 

file.to.read <- 'C:/Users/Guillermo Podesta/Desktop/HydroModelOutput-up-to-2003.nc'

if (!file.exists(file.to.read)) {
  stop('Specified netCDF file does NOT exist...\n')
}

# --- Open file

nc <- ncdf4::nc_open(filename = file.to.read,
  write = FALSE, verbose = TRUE)

rm(file.to.read)
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
# --- First we extract the origin of time units from "units" component of time coordinate.
# --- The following code assumes origin is expresed as "days since YYYY-MM-DD".

tt1 <- stringr::str_locate(nc$dim$time$units, "since ")
origin.date <- stringr::str_sub(nc$dim$time$units, start=tt1[, "end"]+1, end =-1L)

if (origin.date != glob.attr$start_date) {
  # Start date specified as global attribute
  # does not coincide with description of time variable
  warning('Start date specified as global attribute
    does not coincide with description of time variable')
}

time_coord_vals <- nc$dim$time$vals
n.time_coord_vals <- length(time_coord_vals)

date_coords <- as.Date(nc$dim$time$vals, origin = origin.date)

rm(tt1, origin.date)

# --- Expand grid of x- and y-coordinates and
# --- create spatial object with Gauss-Krueger coordinates
# --- (lats and lons are in these coordinates).

uu1 <- expand.grid(lon = x_coords, lat = y_coords)

uu2 <- expand.grid(time = time_coord, lat = y_coords, lon = x_coords)

locations.GK <- sp::SpatialPoints(uu1,
  proj4string = sp::CRS(GK.string))

# --- Convert Gauss-Kruger to decimal lat-lon coordinates

locations.ll <- sp::spTransform(locations.GK,
  sp::CRS(ll.string))

# --- Try extracting auxiliary coordinate variables

# --- Extract auxiliary coordinate variables

dec_lon <- ncdf4::ncvar_get(nc, varid = "dec_lon",
  start = c(1, 1),
  count = c(n.x_coord_vals, n.y_coord_vals),
  verbose = TRUE)

dec_lat <- ncdf4::ncvar_get(nc, varid = "dec_lat",
  start = c(1, 1),
  count = c(n.x_coord_vals, n.y_coord_vals),
  verbose = TRUE)

elevation <- ncdf4::ncvar_get(nc, varid = "elevation",
  start = c(1, 1),
  count = c(n.y_coord_vals, n.x_coord_vals),
  verbose = TRUE)

# --- Extract values of MIKE-SHE output variables

wtd <- ncdf4::ncvar_get(nc, varid = "wtd",
  start = c(1, 1, 1),
  count = c(n.x_coord_vals,  n.y_coord_vals, n.time_coord_vals), verbose = TRUE)

# --- Close original netCDF file

ncdf4::nc_close(nc)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Load MIKE-SHE netCDF output into raster objects ----

file.to.read <- 'C:/Users/Guillermo Podesta/Desktop/HydroModelOutput-up-to-2003.nc'

if (!file.exists(file.to.read)) {
  stop('Specified netCDF file does NOT exist...\n')
}

# --- Open netCDF file and get start date 

nc <- ncdf4::nc_open(filename = file.to.read,
  write = FALSE, verbose = TRUE)

glob.attr <- ncdf4::ncatt_get(nc, varid = 0)

tt1 <- stringr::str_locate(nc$dim$time$units, "since ")
origin.date <- stringr::str_sub(nc$dim$time$units, start=tt1[, "end"]+1, end =-1L)

if (origin.date != glob.attr$start_date) {
  # Start date specified as global attribute
  # does not coincide with description of time variable
  warning('Start date specified as global attribute
    does not coincide with description of time variable')
}

time_coord_vals <- nc$dim$time$vals
date_coords <- as.Date(nc$dim$time$vals, origin = origin.date)


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


rasterVis::plot3D(uu3[[18]], zfac = 2)




# -- Create a raster layer with terrain elevation

elev <- raster::raster(file.to.read,
  values = TRUE,
  crs = GK.string,
  varname = 'elevation')

elev <- raster::mask(elev, cuenca.A.GK)



# --- Plot elevation; add basin boundary

p <- rasterVis::levelplot(elev, margin = FALSE,
  col.regions  = terrain.colors(25))
p <- p + layer(sp.lines(cuenca.A.GK, lwd = 2, col = 'tomato'))
plot(p)

rasterVis::plot3D(elev, 
  zfac=1.5, drape=NULL, col=terrain.colors,
  useLegend=TRUE)







xyFromCell(elev, 1:75)
xyFromCell(elev, 1)





# ----------------------------

library(gdalUtils)

hdf_file <- "C:/Users/Guillermo Podesta/Downloads/DSI_halfdeg_2011337.hdf"

gdalUtils::gdalinfo(hdf_file, verbose = TRUE, mm = TRUE, stats = TRUE,
  proj4 = TRUE)

# --- Convert to GeoTIFF

outfile <- "testout"

gdalUtils::gdal_translate(hdf_file, outfile,
  of="GTiff", sds=TRUE, verbose=TRUE)

file.rename(outfile, paste(outfile,".tif",sep=""))

gdalUtils::gdalinfo('testout.tif', verbose=TRUE)

# --- Create a raster

uuu <- raster::raster(paste(outfile,".tif", sep=""))
uu2 <- raster::projectRaster(uuu, crs = ll.string)
NAvalue(uu2) <- 3e+35
e <- raster::extent(-77, -50, -45, -10)
rc <- raster::crop(uu2, e)

rasterVis::levelplot(rc, margin = FALSE)











  , method="bilinear", 
  alignOnly=FALSE, over=FALSE, filename="", ...) 


Use layout=c(1, 1). With your data: 
  
  levelplot(s, layout=c(1, 1)) 

If you need to create a movie use this code: 
  
  trellis.device(png, file='Rplot%02d.png') 
levelplot(s, layout=c(1, 1)) 
dev.off() 

## ffmpeg must be installed in your system 
movieCMD <- 'ffmpeg -i Rplot%02d.png output.mp4' 
system(movieCMD) 