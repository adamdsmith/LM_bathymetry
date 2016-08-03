# Get list of NAIP 1 m resolution files
NAIP_files <- list.files("./GIS/NAIP14/1_m/", full.names = TRUE)

# Mosaic and aggregate NAIP14
NAIP <- lapply(NAIP_files, brick)
NAIP$fun <- mean
NAIP$na.rm <- TRUE

NAIP <- do.call(mosaic, NAIP)

# Aggregate to 3 m resolution
NAIP_3 <- aggregate(NAIP, fact = 3, expand = FALSE)
writeRaster(NAIP_3, "./GIS/NAIP14/LM_3m.tif")

# Aggregate to 30 m resolution
NAIP_30 <- aggregate(NAIP, fact = 30, expand = FALSE)
writeRaster(NAIP_30, "./GIS/NAIP14/LM_30m.tif")
