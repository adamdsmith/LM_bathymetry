pacman::p_load(leaflet, RColorBrewer)

# Reproject lake boundary and points to lat/lon
if (!exists("LM")) {
  pacman::p_load(rgdal, rgeos)
  LM <- readOGR("./Output", "LM_poly", verbose = FALSE)
  LM <- gBuffer(LM, byid = TRUE, width = 0) # To correct invalid geometry
}
LM_ll <- spTransform(LM, CRS("+proj=longlat +datum=WGS84"))

# Create map
p <- leaflet() %>%
  # Base map group
  addTiles("http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
           group = "Aerial") %>%
  addTiles("http://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}",
           group = "Terrain") %>%
  # Add Lake Mattamuskeet boundary
  addPolygons(data = LM_ll, fillOpacity = 0, smoothFactor = 0.5,
              color = "yellow", weight = 3) %>%
  addLayersControl(baseGroups = c("Aerial", "Terrain")) %>%
  addScaleBar(position = "bottomright") 
