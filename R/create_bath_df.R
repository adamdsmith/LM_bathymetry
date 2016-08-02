create_bath_df <- function(sav_bath = "./Data/depths_1315.csv",
                           use_sav_yrs = 2013:2015,
                           epicol_csv = "./Data/matta_bath.csv",
                           output_bath = TRUE) {
  
  if (!("pacman" %in% installed.packages())) install.packages("pacman", quiet = TRUE)
  pacman::p_load(dplyr, lubridate, sp)

  # Load depth measurements from SAV sample plots
  sav_bath <- read.csv(sav_bath, header = TRUE, stringsAsFactors = FALSE) %>%
    na.omit() %>%
    mutate(year = year(mdy(date)),
           date = as.character(mdy(date))) %>%
    dplyr::select(date, year, lat, lon, depth, basin)
  
  # Correct depth measurements by year since waterData won't allow query over 30? days
  tmp_bath <- lapply(use_sav_yrs, function(yr) {
    tmp <- sav_bath[sav_bath$year == yr, ]
    tmp <- navd_depth(tmp)
  })
  sav_bath <- do.call("rbind", tmp_bath) %>%
    dplyr::select(-date)
    
  
  # Load supplemental bathymetry data
  bath <- read.csv(epicol_csv, stringsAsFactors = FALSE) %>%
    dplyr::select(col_date, col_time, gps_loc_lat, gps_loc_lon, basin, lake_depth, gps_loc_acc) %>%
    # Remove points with 100" depth - deep scour channels at culverts that soap 
    # film smoother can't capture at current resolution
    filter(lake_depth != 100) %>%
    mutate(date = as.character(mdy(col_date)),
           year = year(mdy(col_date)),
           lat = gps_loc_lat,
           lon = gps_loc_lon,
           gps_acc = gps_loc_acc,
           basin = fill_var(basin),
           depth = lake_depth) %>%
    dplyr::select(date, year, lat, lon, depth, basin)
  # Use same day (average) gage data to correct bathymetry points
  bath <- navd_depth(bath) %>% dplyr::select(-date)
  
  # Combine with 2013-2015 points
  bath <- rbind(bath, sav_bath)
  
  # Extract UTM coordinates for sampling points
  bath_pts <- SpatialPointsDataFrame(coords = as.matrix(bath[, c("lon", "lat")]),
                                     data = dplyr::select(bath, -lon, -lat),
                                     proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  bath_pts <- spTransform(bath_pts, CRS("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  
  bath <- bath %>% 
    mutate(x = coordinates(bath_pts)[, 1],
           y = coordinates(bath_pts)[, 2])
  
  if (output_bath) {
    write.csv(bath, "./Output/final_bath_data.csv", row.names = FALSE)
  }
  return(bath)
  
}