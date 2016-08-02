# Load necessary packages
pacman::p_load(rgdal, raster, rgeos, maptools)

LM_NDWI <- brick("./GIS/NAIP14/LM_30m.tif")
LM_NDWI <- calc_NDWI(LM_NDWI)
#writeRaster(LM_NDWI, "./Output/LM_NDWI.tif", overwrite = TRUE)

# Dichotomize at NDWI cutoff
# Value changes from application to application; requires inspection
# Looking at the histogram of NDWI values can suggest a dividing value
#hist(LM_NDWI, breaks = 50)
LM_NDWI <- LM_NDWI >= 0.4
#plot(LM_NDWI)

# Calculate clumps and extract two lake halves
LM <- clump(LM_NDWI, directions = 4)
lake <- which(table(getValues(LM)) >2e6 / 900) # Finds clumps > 2000000 m^2
# Data frame to set non-lake clumps to NA
LM <- calc(LM, function(x) ifelse(x %in% lake, 1, NA))

# Polygonize the raster
LMp <- polygonizer(LM)
# Combine two polygons into single polygon
LMp <- raster::aggregate(LMp, by='DN') 
proj4string(LMp) <- CRS(proj4string(LM))

# Get rid of small holes typically representing individual or small stands of cypress
# Here we apply a cutoff is 250 m2
LMp <- remove_holes(LMp, min_area = 250, name = "Lake Mattamuskeet")
LMp <- disaggregate(LMp)
LMp@data$basin <- c("west", "east")
writeOGR(LMp, "./Output", "LM_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)

# Convert back to raster for use in bathymetry prediction
LM_100 <- raster(crs = crs(LM), ext = extent(LM), res = 100, vals = 1)
LM_100 <- trim(rasterize(LMp, LM_100), values = c(NA, NaN))
writeRaster(LM_100, "./Output/LM_100.tif", overwrite = TRUE)
