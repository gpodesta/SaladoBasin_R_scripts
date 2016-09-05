
# --- Extraer informacion del atlas de suelos de INTA
# --- para cada punto en la grilla de simulacion del modelo MIKE-SHE
# --- en la subcuenca A del Rio Salado del Sur.

Sys.setenv(TZ = "UTC") 

# -----------------------------------------------------------------------------#
# --- Borrar todos los objetos antes de empezar ----

rm(list = ls()); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Cargar paquetes de R necesarios ----

options(repos=c(CRAN="http://lib.stat.cmu.edu/R/CRAN/"))
if (!require(lubridate)) {install.packages("lubridate"); library(lubridate)}
if (!require(stringr)) {install.packages("stringr"); library(stringr)}
if (!require(reshape2)) {install.packages("reshape2"); library(reshape2)}
if (!require(zoo)) {install.packages("zoo"); library(zoo)}
if (!require(xts)) {install.packages("xts"); library(xts)}
if (!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if (!require(Hmisc)) {install.packages("Hmisc"); library(Hmisc)}
if (!require(MASS)) {install.packages("MASS"); library(MASS)}
if (!require(lattice)) {install.packages("lattice"); library(lattice)}
if (!require(sp)) {install.packages("sp"); library(sp)}
if (!require(spdep)) {install.packages("spdep"); require(spdep)}
if (!require(maps)) {install.packages("maps"); library(maps)}
if (!require(mapdata)) {install.packages("mapdata"); library(mapdata)}
if (!require(maptools)) {install.packages("maptools"); library(maptools)}
if (!require(rgdal)) {install.packages("rgdal"); library(rgdal)}
if (!require(PBSmapping)) {install.packages("PBSmapping"); library(PBSmapping)}
if (!require(yaml)) {install.packages("yaml"); require(yaml)}
if (!require(Cairo)) {install.packages("Cairo"); require(Cairo)}
if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Define Coordinate Reference Systems for various projections ----

gk.crs.string <- '+proj=tmerc +lat_0=-90 +lon_0=-60 +k=1 +x_0=5500000 +y_0=0
+ellps=intl +twogs84=-148,136,90,0,0,0,0 +units=m +no_defs'

utm.crs.string <- '+proj=utm +zone=20H +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

ll.crs.string <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Define base directory depending on nodename (Macbook or PCs) ----

tt0 <- Sys.info()
tt1 <- which(names(tt0) == "nodename")
tt2 <- tt0[tt1]

if (tt2 == "GP-MACBOOK") {
  base.dir <- "C:/Data/"
} else {
  base.dir <- "D:/"
}

rm(tt0,tt1,tt2)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Leer contorno de subcuenca A  ----

dir.shape <- paste0(base.dir,"ENSO/Projects/CNH3/data/shapes_cuenca")

tt0 <- file.info(dir.shape)   # Informacion sobre el directorio especificado
if (!tt0$isdir)
  stop("ERROR: El directorio especificado para la cuenca no existe... verificar el nombre...\n")            	

# --- Leer contornos de la cuenca

cuenca.A.GK <- rgdal::readOGR(dir.shape, layer="Salado_A_GK", verbose = TRUE)
cuenca.A.ll <- rgdal::readOGR(dir.shape, layer="Salado_A_latlon", verbose = TRUE)

rm(dir.shape, tt0)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Leer los datos del Atlas de Suelos de INTA ----
# --- URL: http://geointa.inta.gov.ar/publico/INTA_SUELOS/suelos_500000_v9.zip.
# --- Este archivo tiene suelos agrupados por provincia.

# --- Especificamos el directorio donde estan los datos del atlas de suelos

dir.suelos <- paste0(base.dir, "ENSO/Projects/CNH3/data/Suelos Salado Atlas INTA")

tt0 <- file.info(dir.suelos) 	# Informacion sobre el directorio especificado
if (!tt0$isdir)
  stop("ERROR: El directorio especificado no existe... verificar el nombre...\n")            	

# --- Leemos los datos del atlas de suelos

suelos.atlas <- rgdal::readOGR(dir.suelos, layer="suelos_500000_v9", verbose = TRUE)

# --- Convertimos poligonos de suelos a coordenadas Gauss-Kruger

suelos.GK <- sp::spTransform(suelos.atlas, sp::CRS(gk.crs.string))

rm(dir.suelos, tt0, suelos.atlas)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Leer shapefiles de partidos/departamentos en Argentina ----

dir.depts <- paste0(base.dir, "ENSO-Data/Other Data/ARG_departments")

tt0 <- file.info(dir.depts)   # Informacion sobre el directorio especificado
if (!tt0$isdir)
  stop("ERROR: El directorio especificado no existe... verificar el nombre...\n")            	

dept.archivo <- paste(dir.depts, "deptscorrected.shp", sep="")
if (!file.exists(dept.archivo)) {
  stop("Specified SHP file for departments", dept.archivo, "does not exist...\n")
}

depts <- rgdal::readOGR(dir.depts, layer="deptscorrected", verbose = TRUE) # Read departments

# --- Convertir poligonos a coordenadas Gauss-Kruger (in meters)

depts.GK <- sp::spTransform(depts, CRS = sp::CRS(gk.crs.string)) 

rm(dir.depts, dept.archivo, tt0, depts); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Generar grilla densa de puntos ALREDEDOR de la cuenca A del Salado ----

point.spacing.x <- 1000  		# point separation along x-axis (in meters)
point.spacing.y <- 1000			# point separation along y-axis (in meters)

point.offset.x <- point.spacing.x / 2
point.offset.y <- point.spacing.y / 2

# --- Generate regular grid of points inside merged polygon
# --- The number of points in the grid is roughly equal to the area
# --- of the target region divided by the area encompassed by one grid point.

gg3 <- sp::bbox(cuenca.A.GK)		# bounding box of basin A
range.x <- range(pretty(range(gg3["x",]), n = 50, high.u.bias = 100))	# range of lon (x) values
range.y <- range(pretty(range(gg3["y",]), n = 50, high.u.bias = 100))	# range of lat (y) values

# --- The x- and y-values generated exceed the boundaries of basin A

x.values <- seq(from=(floor(range.x[1]) - (2 * point.spacing.x) + point.offset.x),
  to=(ceiling(range.x[2] + (2 * point.spacing.x))),
  by=point.spacing.x)

y.values <- seq(from=(floor(range.y[1]) - (2 * point.spacing.y) + point.offset.y),
  to=ceiling(range.y[2] + (2 * point.spacing.y)),
  by=point.spacing.y)

gg2 <- as.data.frame(expand.grid(lon = x.values, lat = y.values))

# --- Reorder the points so they are ordered along rows of constant latitude.
# --- Latitude is sorted in reverse so the first line is the northernmost lat.

gg2 <- dplyr::tbl_df(gg2) %>%
  dplyr::arrange(desc(lat), lon)
row.names(gg2) <- 1:nrow(gg2)	# matrix of coordinates for points in grid

# --- Create a spatial points object with what will be the centers
# --- of the tiles to be created later.

points.info <- data.frame(point.ID = row.names(gg2),
  x.coord = gg2[, "lon"], y.coord = gg2[, "lat"])

gg3 <- sp::SpatialPointsDataFrame(coords = as.data.frame(gg2),
  data = points.info,
  proj = sp::CRS(gk.crs.string), match.ID = TRUE)

# --- If plot.results == TRUE, plot grid of points

if (plot.results) {
  sp::plot(gg3, axes = TRUE, pch=".", col = "tomato", cex = 1)
  sp::plot(cuenca.A.GK, lwd = 2, border = "steelblue4", add = TRUE)
  title("Dense grid in Salado A Basin")
}

remove(gg1, gg2, range.x, range.y, x.values, y.values)
remove(point.offset.x, point.offset.y); gc()
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Usamos la grilla densa para extraer departamentos ----

pp5 <- sp::over(x = gg3, y = depts.GK)

pp6 <-dplyr::tbl_df(pp5) %>%
  dplyr::mutate(pointID = as.integer(rownames(pp5))) %>%
  dplyr::select(pointID, prov, dept, code) %>%
  dplyr::mutate(prov = as.character(prov), dept = as.character(dept))
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Usamos la grilla densa para extraer datos de suelos del Atlas ----

pp3 <- sp::over(x = gg3, y = suelos.GK)

# --- Construimos tabla con tipos de suelos suelos para grilla densa

tabla.suelos<- dplyr::tbl_df(pp3) %>%
  dplyr::select(SIMBC,IND_PROD,SGRUP_SUE1,PORC_SUE1,GGRUP_SUE1,
    SGRUP_SUE2,PORC_SUE2,GGRUP_SUE2,
    SGRUP_SUE3, PORC_SUE3,GGRUP_SUE3) %>%
  dplyr::transmute(soilUnitId = as.character(SIMBC),
    prodIndex = as.numeric(IND_PROD),
    nameSoil1 = as.character(SGRUP_SUE1),
    propSoil1 = as.numeric(PORC_SUE1),
    grpSoil1 = as.character(GGRUP_SUE1),
    nameSoil2 = as.character(SGRUP_SUE2),
    propSoil2 = as.numeric(PORC_SUE2),
    grpSoil2 = as.character(GGRUP_SUE2),
    nameSoil3 = as.character(SGRUP_SUE3),
    propSoil3 = as.numeric(PORC_SUE3),
    grpSoil3 = as.character(GGRUP_SUE3)) %>%
  dplyr::mutate(pointID = as.numeric(rownames(pp3))) %>%
  dplyr::arrange(soilUnitId)

tabla.suelos[tabla.suelos == '-'] <- NA
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Armar tabla con indice de productividad por unid de suelo ----

tt4 <- dplyr::tbl_df(tabla.suelos) %>%
  dplyr::select(pointID, soilUnitId, prodIndex) 

tt5 <- dplyr::left_join(pp6, tt4, by = "pointID")

partidos.ba.area5 <- c("GENERAL VILLEGAS",
  "TRENQUE LAUQUEN","RIVADAVIA","PEHUAJO","SALLIQUELO")

partidos.ba.area4 <- c("BRAGADO","LEANDRO N. ALEM","HIPOLITO YRIGOYEN",
  "CARLOS CASARES","BOLIVAR","25 de MAYO",
  "9 DE JULIO","GENERAL VIAMONTE","LINCOLN","GENERAL PINTO")

tt6 <- tt5 %>% dplyr::filter(dept %in% partidos.ba.area5)

tt7 <- tt5 %>% dplyr::filter(dept %in% partidos.ba.area4)


par(cex.axis = 0.5)
boxplot(split(tt6$prodIndex, tt6$dept),
  horizontal = TRUE, range = 0, las = 1,
  ylim = c(0, 100))

par(cex.axis = 0.5)
boxplot(split(tt7$prodIndex, tt7$dept),
  horizontal = TRUE, range = 0, las = 1,
  ylim = c(0, 100))


boxplot(tt6$prodIndex,
  horizontal = TRUE, range = 0, las = 1,
  ylim = c(0, 100))

quantile(tt5$prodIndex, probs = c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95, 1))
quantile(tt6$prodIndex, probs = c(0.10, 0.25, 0.75, 0.90))

range(tt5$prodIndex)

uuu <- tapply(X = tt7$prodIndex, INDEX = tt7$dept, FUN = mean)

tapply(X = tt6$prodIndex, INDEX = tt6$dept,
  FUN = range)















# -----------------------------------------------------------------------------#
# --- Construir tabla con suelos y proporciones para cada punto de grilla -----

tt1 <- dplyr::tbl_df(tabla.suelos) %>%
  dplyr::select(pointID, nameSoil1, grpSoil1, propSoil1, soilUnitId, pointID) %>%
  dplyr::rename(soilName = nameSoil1,
    soilGroup = grpSoil1,
    soilProp = propSoil1)

tt2 <- dplyr::tbl_df(tabla.suelos) %>%
  dplyr::select(pointID, nameSoil2, grpSoil2, propSoil2, soilUnitId, pointID) %>%
  dplyr::rename(soilName = nameSoil2,
    soilGroup = grpSoil2,
    soilProp = propSoil2)

tt3 <- dplyr::tbl_df(tabla.suelos) %>%
  dplyr::select(pointID, nameSoil3, grpSoil3, propSoil3, soilUnitId, pointID) %>%
  dplyr::rename(soilName = nameSoil3,
    soilGroup = grpSoil3,
    soilProp = propSoil3)

tt4 <- dplyr::bind_rows(tt1, tt2, tt3) %>%
  dplyr::select(pointID, soilUnitId, soilName, soilGroup, soilProp) %>%
  dplyr::mutate(pointID = as.integer(pointID)) %>%
  dplyr::arrange(pointID) %>%
  dplyr::mutate(soilName = stringr::str_replace_all(soilName,
    pattern = 'oles ', replacement = 'ol '))  %>%
  dplyr::mutate(soilName =stringr::str_replace_all(soilName,
    pattern = 'es ', replacement = ' '))

rm(tt1,tt2,tt3)
# ------------------------------------------------------------------------------

# -----------------------------------------------------------------------------#
# --- Unir tablas de suelos y partidos ----

pp7 <- dplyr::left_join(pp6, tt4, by = "pointID")
# ------------------------------------------------------------------------------



pp7 %>% dplyr::filter(dept == "GENERAL VILLEGAS")









# -----------------------------------------------------------------------------#

# ------------------------------------------------------------------------------











# -----------------------------------------------------------------------------#
# --- Construir tabla SOILUNIT ----

soilunit <- dplyr::distinct_(pp4) %>%
  dplyr::select(soilUnitId,
    soilUnitType,
    productivityIndex,
    limitation1,
    limitation2,
    surfaceTextureSoil1,
    subSurfaceTextureSoil1,
    drainageSoil1,
    positionSoil1,
    positionSoil2,
    positionSoil3) %>%
  as.data.frame(stringsAsFactors = FALSE)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# --- Construir tabla SOILUNIT_SOILTYPE ----

#  soilUnitId,
#  soilTypeId,
#  percentage

uu1 <- dplyr::tbl_df(pp4) %>%
  dplyr::select(soilUnitId,
    nameSoil1, propSoil1) %>%
  dplyr::rename(soilName = nameSoil1,
    proportion = propSoil1) %>%
  dplyr::group_by(soilUnitId, soilName, proportion) %>%
  dplyr::distinct()

uu2 <- dplyr::tbl_df(pp4) %>%
  dplyr::select(soilUnitId,
    nameSoil2, propSoil2) %>%
  dplyr::rename(soilName = nameSoil2,
    proportion = propSoil2) %>%
  dplyr::group_by(soilUnitId, soilName, proportion) %>%
  dplyr::distinct()

uu3 <- dplyr::tbl_df(pp4) %>%
  dplyr::select(soilUnitId,
    nameSoil3, propSoil3) %>%
  dplyr::rename(soilName = nameSoil3,
    proportion = propSoil3) %>%
  dplyr::group_by(soilUnitId, soilName, proportion) %>%
  dplyr::distinct()

uu4 <- dplyr::rbind_list(uu1,uu2,uu3)

uu5 <- dplyr::tbl_df(tt4) %>%
  dplyr::select(soilTypeId, soilName)

soilunit_soiltype <- dplyr::left_join(uu4, uu5) %>%
  dplyr::select(soilUnitId, soilTypeId, proportion) %>%
  dplyr::arrange(soilUnitId, soilTypeId) %>%
  as.data.frame(stringsAsFactors = FALSE)

rm(uu1,uu2,uu3,uu4,uu5,tt4)
# ------------------------------------------------------------------------------










# --- Llenamos variables para suelos en data.frame suelo.grilla.DF

suelo.grilla.DF[, "soilUnit"] <- as.character(pp3[ ,"SIMBC"])
suelo.grilla.DF[, "unitType"] <- as.character(pp3[ ,"TIPO_UC"])
suelo.grilla.DF[, "prodIndx"] <- as.numeric(pp3[ ,"IND_PROD"])
suelo.grilla.DF[, "prcntSl3"] <- as.numeric(pp3[ ,"PORC_SUE3"])
suelo.grilla.DF[, "prcntSl2"] <- as.numeric(pp3[ ,"PORC_SUE2"])
suelo.grilla.DF[, "prcntSl3"] <- as.numeric(pp3[ ,"PORC_SUE3"])
suelo.grilla.DF[, "grpSl3"] <- as.character(pp3[ ,"SGRUP_SUE3"])
suelo.grilla.DF[, "grpSl2"] <- as.character(pp3[ ,"SGRUP_SUE2"])
suelo.grilla.DF[, "grpSl3"] <- as.character(pp3[ ,"SGRUP_SUE3"])

# --- Verificamos que los porcentajes de los tres suelos principales
# --- en cada unidad cartografica de suelos no sumen mas de 300.

tt3 <- suelo.grilla.DF[, "prcntSl3"] +
  suelo.grilla.DF[, "prcntSl2"] +
  suelo.grilla.DF[, "prcntSl3"]

if (any(tt3 != 300)) {
  tt2 <- which(tt3 != 300)
  warning("Error en porcentajes de suelos ... Suma no da 300%\n")
}

bbb <- dplyr::filter(suelo.grilla.DF,
  prcntSl3 + prcntSl2 + prcntSl3 != 100)  %>%
  dplyr::select(., undCart,
    grpSl3, prcntSl3,
    grpSl2, prcntSl2,
    grpSl3, prcntSl3) %>%
  dplyr::group_by(., undCart)


rm(tt0, oldPath, dir.suelos,pp2,pp3); gc()
# ----------------------------------------------------------------------------------------




dplyr::select(SIMBC,TIPO_UC,LIMIT_PPAL,LIMIT_SECU,
  IND_PROD,
  PORC_SUE3,ORDEN_SUE3,GGRUP_SUE3,SGRUP_SUE3,
  PORC_SUE2,ORDEN_SUE2,GGRUP_SUE2,SGRUP_SUE2,
  PORC_SUE3,ORDEN_SUE3,GGRUP_SUE3,SGRUP_SUE3,
  TEXT_SUPS3, TEXT_BS3) %>%


# ---------------------------------------------------------------------------------------#
# --- Creamos un objeto para guardar datos de suelos ----
# --- Creamos un data.frame vacio donde guardaremos
# --- informacion sobre suelos a extraerse del Atlas de Suelos de INTA.

n.puntos <- nrow(grilla.5000.GK)

suelo.grilla.DF <- data.frame(puntoId = rep(NA, n.puntos),
  cncId = rep(NA, n.puntos),
  cncNom = rep(NA, n.puntos),
  subCncNom = rep(NA, n.puntos),
  GKlon = rep(NA, n.puntos),
  GKlat = rep(NA, n.puntos),
  deglon = rep(NA, n.puntos),
  deglat = rep(NA, n.puntos),
  UTMlon = rep(NA, n.puntos),
  UTMlat = rep(NA, n.puntos),
  provId = rep(NA, n.puntos),
  provNom = rep(NA, n.puntos),  		 
  deptId = rep(NA, n.puntos),
  deptNom = rep(NA, n.puntos),							  
  undCart = rep(NA, n.puntos),
  undTipo = rep(NA, n.puntos),
  prodInd = rep(NA, n.puntos),						 										 
  prcntSl3 = rep(NA, n.puntos),					 
  grpSl3 = rep(NA, n.puntos),							 
	prcntSl2 = rep(NA, n.puntos),					 
  grpSl2 = rep(NA, n.puntos),							 
	prcntSl3 = rep(NA, n.puntos),					 
  grpSl3 = rep(NA, n.puntos) )
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Llenamos algunas de las variables en objeto suelo.grilla.DF ----

suelo.grilla.DF[ ,"puntoId"] <- as.character(grilla.5000.GK@data[, "ID"])
suelo.grilla.DF[ ,"GKlon"] <- round(coordinates(grilla.5000.GK)[, "coords.x3"], 0)
suelo.grilla.DF[ ,"GKlat"] <- round(coordinates(grilla.5000.GK)[, "coords.x2"], 0)
suelo.grilla.DF[ ,"deglon"] <- round(coordinates(grilla.5000.ll)[, "coords.x3"], 2)
suelo.grilla.DF[ ,"deglat"] <- round(coordinates(grilla.5000.ll)[, "coords.x2"], 2)

# --- Extraer valores de cuenca y subcuenca

tt3 <- sp::over(x = grilla.5000.GK, y = cuenca.A.GK)

table(tt3$cuenca, useNA = "always")
table(tt3$subcuenca, useNA = "always")

suelo.grilla.DF[ ,"cncId"] <- round(as.numeric(tt3$ID), 3)
suelo.grilla.DF[ ,"cncNom"] <- tt3$cuenca
suelo.grilla.DF[ ,"subCncNom"] <- tt3$subcuenca

rm(tt3)
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Ahora extraer provincia (ID y nombre) para cada punto de la grilla ----
# --- Primero leemos los shape files CORREGIDOS de provincias argentinas

oldPath <- getwd()

prov.dir <- "D:/ENSO-Data/Other Data/ARG_provinces/"

tt0 <- file.info(prov.dir)    # Informacion sobre el directorio especificado
if (tt0$isdir) {							# Es un directorio?
	setwd(prov.dir)						  # Cambiar al directorio especificado
} else {
	stop("ERROR: El directorio especificado para suelos no existe... verificar el nombre...\n")
}          			
rm(tt0)

prov.file <- paste(prov.dir, "provs.corrected.shp", sep="")
if (!file.exists(prov.file))
  stop("ERROR: El archivo especificado ", prov.file, "no existe...\n")

pp3 <- rgdal::readOGR(".", layer="provs.corrected", verbose = TRUE);
setwd(oldPath)

# --- Convertimos poligonos de provincias a coordenadas Gauss-Kruger

pp2 <- spTransform(pp3,
	CRS("+proj=tmerc +lat_0=-90 +lon_0=-60 +k=3 +x_0=5500000 +y_0=0 +ellps=intl +units=m +no_defs"))

# --- Usamos los puntos de la grilla para extraer datos sobre la provincia.

pp3 <- sp::over(x = grilla.5000.GK, y = pp2)

# --- Llenar variables para provincia en data.frame suelo.grilla.DF

suelo.grilla.DF[ ,"provId"] <- floor(as.numeric(pp3[ ,"indec.code"]))
suelo.grilla.DF[ ,"provNom"] <- as.character(pp3[ ,"prov.name"])

remove(prov.dir, prov.file, oldPath, pp3, pp2, pp3); gc()
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Ahora extraemos datos para departamentos (ID y nombre) ----
# --- para cada punto de la grilla
# --- Primero leemos los shape files CORREGIDOS de departamentos/partidos argentinos

oldPath <- getwd()

dept.dir <- "D:/ENSO-Data/Other Data/Argentine departments/corrected/"
tt0 <- file.info(dept.dir)  	# Informacion sobre el directorio especificado
if (tt0$isdir) {							# Es un directorio?
	setwd(dept.dir)							# Cambiar al directorio especificado
} else {
	stop("ERROR: El directorio especificado para suelosno existe... verificar el nombre...\n")
}          			
rm(tt0)

dept.file <- paste(dept.dir, "deptscorrected.shp", sep="")
if (!file.exists(dept.file))
  stop("ERROR: Archivo SHP para departamentos", dept.file, "no existe...\n")

pp3 <- rgdal::readOGR(".", layer="deptscorrected", verbose = TRUE)
setwd(oldPath)

# --- Convertimos poligonos de departamentos/partidos a coordenadas Gauss-Kruger

pp2 <- sp::spTransform(pp3,
	CRS("+proj=tmerc +lat_0=-90 +lon_0=-60 +k=3 +x_0=5500000 +y_0=0 +ellps=intl +units=m +no_defs"))

# --- Usamos los puntos de la grilla recortada para extraer datos sobre departamentos

pp3 <- sp::over(x=grilla.5000.GK, y=pp2)

# --- Llenamos variables para departamentos en data.frame suelo.grilla.DF

suelo.grilla.DF[ ,"deptId"] <- floor(as.numeric(pp3[ ,"code"]))
suelo.grilla.DF[ ,"deptNom"] <- as.character(pp3[ ,"dept"])

remove(dept.dir, dept.file, oldPath,pp3, pp2, pp3); gc()
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Finalmente, agregamos datos de suelos para cada punto de la grilla
# ---------------------------------------------------------------------------------------#

# ---------------------------------------------------------------------------------------#
# --- Leemos los datos del Atlas de Suelos de INTA ----
# --- URL: http://geointa.inta.gov.ar/publico/INTA_SUELOS/suelos_500000_v9.zip.
# --- Este archivo tiene suelos agrupados por provincia.

oldPath <- getwd()

# --- Especificamos el directorio donde estan los datos del atlas de suelos

dir.suelos <- "D:/ENSO-Data/Other Data/ARG_soils_by_province/"	

tt0 <- file.info(dir.suelos) 	# Informacion sobre el directorio especificado
if (tt0$isdir) {							# Es un directorio?
	setwd(dir.suelos)						# Cambiar al directorio especificado
} else {
	stop("ERROR: El directorio especificado para suelos no existe... verificar nombre...\n")
}          			

# --- Leemos los datos del atlas de suelos

suelos.atlas <- rgdal::readOGR(".", layer="suelos_500000_v9", verbose = TRUE)
setwd(oldPath)	# Volver a directorio (carpeta) original

# --- Convertimos poligonos de suelos a coordenadas Gauss-Kruger

pp2 <- sp::spTransform(suelos.atlas,
	CRS("+proj=tmerc +lat_0=-90 +lon_0=-60 +k=3 +x_0=5500000 +y_0=0 +ellps=intl +units=m +no_defs"))

# --- Usamos los puntos de la grilla recortada para extraer datos sobre departamentos

pp3 <- over(x = grilla.5000.GK, y = pp2)

# --- Llenamos variables para suelos en data.frame suelo.grilla.DF

suelo.grilla.DF[, "undCart"] <- as.character(pp3[ ,"SIMBC"])
suelo.grilla.DF[, "undTipo"] <- as.character(pp3[ ,"TIPO_UC"])
suelo.grilla.DF[, "prodInd"] <- as.numeric(pp3[ ,"IND_PROD"])
suelo.grilla.DF[, "prcntSl3"] <- as.numeric(pp3[ ,"PORC_SUE3"])
suelo.grilla.DF[, "prcntSl2"] <- as.numeric(pp3[ ,"PORC_SUE2"])
suelo.grilla.DF[, "prcntSl3"] <- as.numeric(pp3[ ,"PORC_SUE3"])
suelo.grilla.DF[, "grpSl3"] <- as.character(pp3[ ,"SGRUP_SUE3"])
suelo.grilla.DF[, "grpSl2"] <- as.character(pp3[ ,"SGRUP_SUE2"])
suelo.grilla.DF[, "grpSl3"] <- as.character(pp3[ ,"SGRUP_SUE3"])

# --- Verificamos que los porcentajes de los tres suelos principales
# --- en cada unidad cartografica de suelos no sumen mas de 300.

tt3 <- suelo.grilla.DF[, "prcntSl3"] +
	suelo.grilla.DF[, "prcntSl2"] +
	suelo.grilla.DF[, "prcntSl3"]

if (any(tt3 != 300)) {
	tt2 <- which(tt3 != 300)
	warning("Error en porcentajes de suelos ... Suma no da 300%\n")
}

bbb <- dplyr::filter(suelo.grilla.DF,
              prcntSl3 + prcntSl2 + prcntSl3 != 300)  %>%
  dplyr::select(., undCart,
                grpSl3, prcntSl3,
                grpSl2, prcntSl2,
                grpSl3, prcntSl3) %>%
  dplyr::group_by(., undCart)


rm(tt0, oldPath, dir.suelos,pp2,pp3); gc()
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Creamos objeto de clase SpatialPointsDataFrame con la info de suelos ----

grilla.suelos.GK <- sp::SpatialPointsDataFrame(coordinates(grilla.5000.GK),
  data =	suelo.grilla.DF,
  match.ID = TRUE, proj4string = CRS("+proj=tmerc +lat_0=-90 +lon_0=-60 +k=3 +x_0=5500000 +y_0=0 +ellps=intl
  +twogs84=-348,336,90,0,0,0,0 +units=m +no_defs"))

# --- Creamos otros objetos SpatialPointsDataFrame con coordenadas lat-lon y UTM

grilla.suelos.ll <- spTransform(grilla.suelos.GK,
	CRS("+proj=longlat +datum=WGS84 +no_defs"))

grilla.suelos.UTM <- spTransform(grilla.suelos.GK,
	CRS("+proj=utm +zone=20H +south +datum=WGS84 +units=m +no_defs"))

# --- Escribir puntos de grilla con datos de suelos para cuenca A3 en coordenadas lat-lon y UTM 

suelos.dir <- "D:/ENSO/Projects/CNH-3 Project/Data/info_suelos/"

oldPath <- getwd()
setwd(suelos.dir)

writeOGR(grilla.suelos.GK, dsn=".", layer="Salado_A_suelos_5000_GK", driver="ESRI Shapefile",
	check_exists=TRUE, overwrite_layer=TRUE, verbose = TRUE)

writeOGR(grilla.suelos.ll, dsn=".", layer="Salado_A3_suelos_5000_latlon", driver="ESRI Shapefile",
	check_exists=TRUE, overwrite_layer=TRUE, verbose = TRUE)

writeOGR(grilla.suelos.UTM, dsn=".", layer="Salado_A3_suelos_5000_UTM", driver="ESRI Shapefile",
	check_exists=TRUE, overwrite_layer=TRUE, verbose = TRUE)

setwd(oldPath)
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Crear lista de tipos de suelo ----

uu3 <- dplyr::select(suelo.grilla.DF, -c(UTMlon, UTMlat, undCart, undTipo, prodInd))

# --- Suelos en toda la grilla de simulacion

uu2 <- reshape2::melt(uu3,
  id.vars = c("puntoId", "cncId"),     
  measure.vars = c("grpSl3", "grpSl2", "grpSl3"))

uu3 <- dplyr::arrange(uu2, puntoId)

uu3b <- sort(unique(uu3$value)) # Lista de suelos en TODA la grilla de simulacion

# --- Suelos dentro de la cuenca A

uu4 <- dplyr::filter(uu3, cncId == 3)

uu5 <- reshape2::melt(uu4,
  id.vars = c("puntoId", "cncId"),     
  measure.vars = c("grpSl3", "grpSl2", "grpSl3"))

uu6 <- dplyr::arrange(uu5, puntoId)

uu7 <- sort(unique(uu6$value))  # Lista de suelos DENTRO de la cuenca A

rm(uu3,uu2,uu3,uu4,uu5,uu6)
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Estimar proporcion de suelos ---- 
# --- Para la grilla, usamos el data.frame creado para estimar proporciones
# --- de departmentos, tipos de suelo, etc.

# --- Primero, proporcion de departamentos y provincias

tt3 <- suelo.grilla.DF[ , c("deptId","deptNom","provNom")]
tt2 <- table(tt3$deptId) 	# Tabla de numero de puntos de grilla por partido/depto
tt3 <- match(names(tt2), tt3$deptId)
tt4 <- tt3[tt3, ]
tt5 <- data.frame(tt4, npnts=as.vector(tt2))
tt6 <- tt5[order(tt5$npnts), ]
tt6$deptNom <- ordered(tt6$deptNom, levels=tt6$deptNom)

# --- Dotplot basico

dotplot(deptNom ~ npnts,
	data=tt6,
	xlab="Area",
	main="Area de Departamentos en Salado A3",
	pch=23, col="grey50", fill="tomato", cex=3.7)

# --- Dotplot algo mas sofisticado, usando funcion dotchart2()
# --- en paquete Hmisc. Los departamentos se agrupan por provincia.

dotchart2(tt6$npnts,
	labels=tt6$deptNom,
	groups=as.factor(tt6$provNom),
	auxdata=tt6$npnts,
	auxtitle="Nro\npuntos",
	xlab="Area (km2)",
	main="Area de Departamentos en Salado A3",
	dotsize=3.3,
	cex.labels=0.45,
	cex.group.labels=0.6,
	sort.=FALSE,
	pch=36, col="tomato")			

# --- Segundo, proporcion de suelos.
# --- Hacemos una larga columna con los tres suelos principales
# --- en cada unidad cartografica y sus proporciones respectivas

tt3a <- suelo.grilla.DF[, c("grpSl3","prcntSl3")]
tt3b <- suelo.grilla.DF[, c("grpSl2","prcntSl2")]
tt3c <- suelo.grilla.DF[, c("grpSl3","prcntSl3")]

colnames(tt3a) <- c("grp","prcnt")
colnames(tt3b) <- c("grp","prcnt")
colnames(tt3c) <- c("grp","prcnt")

tt3 <- rbind(tt3a, tt3b, tt3c)
tt3[,"prcnt"] <- (tt3[,"prcnt"] / 300)	# Convertir porcentaje (0-300) de cada suelo a proporcion (0-3)

tt2 <- round(tapply(tt3$prcnt, INDEX=tt3$grp, FUN=sum, simplify=TRUE),0)
tt3 <- tt2[order(tt2, decreasing=TRUE)]
tt4 <- data.frame(suelo=names(tt3), area=tt3)
rownames(tt4) <- NULL
tt4$suelo <- ordered(tt4$suelo, levels=tt4$suelo)

# --- Dotplot basico

dotplot(suelo ~ area,
	data=tt4,
	xlab="Area",
	main="Area de Suelos en Salado A3",
	pch=23, col="grey50", fill="tomato", cex=3.7)

dotchart2(tt4$area,
	labels=tt4$suelo,
	auxdata=tt4$area,
	auxtitle="Nro\npuntos",
	xlab="Area (km2)",
	main="Area de suelos en Salado A3",
	dotsize=3.25,
	cex.labels=0.45,
	cex.group.labels=0.6,
	sort.=FALSE,
	pch=36, col="tomato")			
# ----------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------#
# --- Seleccionar puntos en grilla de MIKE  DENTRO de cuenca ----

puntos.en.cuenca <- grilla.suelos.GK[!is.na(grilla.suelos.GK@data$cncId), ]

table(puntos.en.cuenca@data$prodInd, useNA= "always")
hist(puntos.en.cuenca@data$prodInd)

plot(cuenca.A.GK, axes = TRUE)
plot(puntos.en.cuenca[puntos.en.cuenca@data$prodInd < 3, ],
     add = TRUE)

# --- Categorias de indice de productividad de FAO
# --- TO DO: Referencia (Federico)

bb2 <- cut(puntos.en.cuenca@data$prodInd,
    breaks = c(0, 7, 39, 34, 64, 300),
    levels = c("Extr. pobre",
               "Pobre",
               "Medio",
               "Bueno",
               "Excelente"),
    include.lowest = TRUE, right = TRUE,
    ordered_result = TRUE) 

colores <- RColorBrewer::brewer.pal(5, "Spectral")

plot(cuenca.A.GK, axes = TRUE)
plot(puntos.en.cuenca, pch = 35, cex = 3.2,
     col = colores[bb2],
     add = TRUE)

# --- Definir usos posibles asociados a rangos de indice de productividad

bb0 <- puntos.en.cuenca@data$prodInd


# --- Categorias de indice de productividad de FAO
# --- TO DO: Referencia (Federico)
  
bb3 <- cut(puntos.en.cuenca@data$prodInd,
             breaks = c(0, 7, 39, 34, 64, 300),
             levels = c("Execrable",
                        "Solo cria",
                        "Cria o Invernada",
                        "Invernada o Agricultura",
                        "Solo Agricultura"),
             include.lowest = TRUE, right = TRUE,
             ordered_result = TRUE) 

colores <- RColorBrewer::brewer.pal(5, "RdYlGn")


Cairo::CairoPDF(file = "tipos_de_uso.pdf",
        onefile = TRUE,
        width = 4, height=6)

plot(puntos.en.cuenca, pch = 35, cex = 3.2,
     col = colores[bb3],
     axes = TRUE,
     xlim = c(5342500, 5457500),
     ylim = c(5922500, 6262500))
plot(cuenca.A.GK, add = TRUE, border = "steelblue", lwd = 2)
plot(pp2, add = TRUE, border = "grey 60")

legend(x = "bottom",
       legend = c("Execrable",
          "Solo cria",
          "Cria o Invernada",
          "Invernada o Agricultura",
          "Solo Agricultura"),
       fill = colores,
       title = "Tipos de uso",
       ncol = 3)
       
dev.off()


# --- Elegir y graficar puntos agricolas (ind prod > 40)

bb0 <-  puntos.en.cuenca@data[puntos.en.cuenca@data$prodInd >= 34, ]

puntos.agricolas <-  puntos.en.cuenca[puntos.en.cuenca@data$prodInd >= 34, ]

plot(puntos.agricolas, axes = TRUE, pch=36)
plot(cuenca.A.GK, add = TRUE, border = "steelblue", lwd = 2)
plot(pp2, add = TRUE, border = "grey 60")





# ----------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------#
# --- Crear lista de tipos de suelo en cuenca A ----

uu0 <-  puntos.en.cuenca@data[puntos.en.cuenca@data$prodInd >= 34, ]

uu3 <- dplyr::select(uu0, undCart,grpSl3,grpSl2,grpSl3)
                     
uu2 <- reshape2::melt(uu3,
    id.vars = "undCart",     
    measure.vars = c("grpSl3", "grpSl2", "grpSl3"))

uu3 <- dplyr::filter(uu2, !is.na(value)) %>%
  dplyr::select(., value) %>%
  dplyr::group_by(., value) %>%
  dplyr::summarise(., numero = n()) %>%
  dplyr::arrange(., desc(numero))




# ----------------------------------------------------------------------------------------




# ---------------------------------------------------------------------------------------#
# --- Clean up ----

rm(list = ls()); gc();
# ------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------#
# --- Read coastline and administrative units for plotting results

# --- Read GSHHS database containing shape files.
# --- Using version 2.3
# --- Available from: http://www.soest.hawaii.edu/pwessel/gshhs/index.html
# --- Reference: Wessel, P. and Smith, W.H.F., 3996. 
# --- A global, self-consistent, hierarchical, high-resolution shoreline database.
# --- J. Geophys. Res., 303(B4): 8743-8743.

gshhs.3 <- readShapeSpatial("D:/ENSO-Data/Other Data/gshhs_shapefiles/i/GSHHS_i_L3.shp",
      repair=TRUE, force_ring=TRUE, proj4string=CRS("+proj=longlat ellps=WGS84"));
gshhs.2 <- readShapeSpatial("D:/ENSO-Data/Other Data/gshhs_shapefiles/i/GSHHS_i_L2.shp",
    	repair=TRUE, force_ring=TRUE, proj4string=CRS("+proj=longlat ellps=WGS84"));

# --- Read WDBII files with political borders and rivers
# --- Available from: http://www.soest.hawaii.edu/pwessel/gshhs/index.html

wdb2.3 <- readShapeSpatial("D:/ENSO-Data/Other Data/WDBII_shapefiles/h/WDBII_border_h_L3.shp",
      repair=TRUE, force_ring=TRUE, proj4string=CRS("+proj=longlat ellps=WGS84"));
wdb2.2 <- readShapeSpatial("D:/ENSO-Data/Other Data/WDBII_shapefiles/f/WDBII_river_f_L3.shp",
    	repair=TRUE, force_ring=TRUE, proj4string=CRS("+proj=longlat ellps=WGS84"));
wdb2.3 <- readShapeSpatial("D:/ENSO-Data/Other Data/WDBII_shapefiles/f/WDBII_river_f_L2.shp",
    	repair=TRUE, force_ring=TRUE, proj4string=CRS("+proj=longlat ellps=WGS84"));

# --- Read shape files with boundaries for Argentine provinces

prov.dir <- "D:/ENSO-Data/Other Data/ARG_provinces/";
prov.file <- paste(prov.dir, "provs.corrected.shp", sep="");
if (!exists("prov.file")) {
  stop("ERROR: Specified SHP file for provinces", prov.file, "does not exist...\n") }

oldPath <- getwd() ;
setwd(prov.dir);
provs <- readOGR(".", layer="provs.corrected", verbose = TRUE);
setwd(oldPath);

rm(prov.dir, prov.file);

# --- Read CORRECTED shape files with Argentine departments

dept.dir <- "D:/ENSO-Data/Other Data/Argentine departments/corrected/"
dept.file <- paste(dept.dir, "deptscorrected.shp", sep="")
if (!exists("dept.file")) {
  stop("ERROR: Specified SHP file for departments", dept.file, "does not exist...\n") }

oldPath <- getwd() 
setwd(dept.dir)
depts <- readOGR(".", layer="deptscorrected", verbose = TRUE)
setwd(oldPath)

rm(dept.dir, dept.file);
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --- Plot base map

plot(gshhs.3, col=rgb(246/255, 232/255, 395/255),
    border="darkgrey", bg="lightblue3", axes=TRUE, xaxs="i", yaxs="i",
		ylim=c(-35.5,-33.5), xlim=c(-64.0, -60.0));
plot(gshhs.2, add=TRUE, border="darkgrey",  col="lightblue3");
box();

plot(wdb2.3, add=TRUE, col="grey40", lwd=2);				# Country borders

plot(provs, add=TRUE, border="grey60", lwd=2);			# Province boundaries

plot(depts, add=TRUE, border="grey80", lwd=3);			# Province boundaries


plot(cuenca.A3.ll, add=TRUE, border="tomato", lwd=2)

plot(grilla.A3.recortada.ll, add=TRUE, pch="36", col="slateblue", cex=0.07)


# ------------------------------------------------------------------------------


