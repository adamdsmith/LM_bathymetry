pacman::p_load(ggplot2, ggmap, viridis)

# Reproject lake boundary and points to lat/lon
if (!exists("LM")) {
  pacman::p_load(rgdal, rgeos)
  LM <- readOGR("./Output", "LM_poly", verbose = FALSE)
  LM <- gBuffer(LM, byid = TRUE, width = 0) # To correct invalid geometry
}
LM <- spTransform(LM, CRS("+proj=longlat +datum=WGS84"))
LM_bb <- bbox(LM) + matrix(c(-0.01, 0.01, -0.01, 0.01), nrow=2, byrow=TRUE)
LM <- fortify(LM, regions = "basin")

# Create map
LM_map <- get_map(location = LM_bb, zoom = 11, maptype = "satellite", source = "google")
p <- ggmap(LM_map) +
  geom_path(data=LM, aes(x = long, y = lat, group = group), size=0.5,
            colour = viridis(1, b = 0.9)) + 
  #  geom_map(data = LM, map = LM,
#           aes(x = long, y = lat, map_id = id),
#           color = viridis(1, b = 1), fill = NA) +
  coord_fixed(xlim = LM_bb[1, ], ylim = LM_bb[2, ]) +
  xlab(NULL) + ylab(NULL)



