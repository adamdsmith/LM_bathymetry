fit_LM_bath <- function(bath_data = NULL,
                        basins = c("west", "east"),
                        check_soap = TRUE,
                        knot_spacing = 900,
                        print_summary = TRUE,
                        gam_check = TRUE,
                        plot_smooth = TRUE,
                        pred_rast = TRUE,
                        seed = NULL,
                        export_gis = TRUE) {
  
  if (is.null(bath_data)) stop("You must specify the data frame containing bathymetric data.")
  
  if (is.null(seed)) r_seed <- sample(2^30, 1) else r_seed <- seed
  
  bath_data$depth <- -1 * bath_data$depth
  
  if (!("pacman" %in% installed.packages())) 
    install.packages("pacman", quiet = TRUE)
  pacman::p_load(raster, rgdal, rgeos, sp, dplyr, lubridate, mgcv, viridis)

  # Identify basins to model
  west <- "west" %in% basins
  east <- "east" %in% basins
  
  if (!any(west, east)) stop("The `basins` argument must be one of: 'west', 'east', or c('west', 'east').")
  
  # Load prediction area (100 m resolution grid)
  # West basin = 1, East basin = 2
  LM <- raster("./Output/LM_100.tif")
  LM_pts <- rasterToPoints(LM, spatial = TRUE)
  LM_pts@data <- LM_pts@data %>% 
    rename(basin = LM_100) %>%
    mutate(basin = ifelse(basin == 1, "west", "east"))
  
  # Load polygon shapefile of extracted shoreline
  LM_poly <- readOGR("./Output", "LM_poly", verbose = FALSE)
  
  # Crude means of getting valid geometry, as it doesn't seem to be
  #gIsValid(LM_poly)
  LM_poly <- gBuffer(LM_poly, byid = TRUE, width = 0)
  
  # Get lake elevations on date used to generate lake boundary
  elev <- LM_lake_elev("2014-10-05")

  # Set up soap film smoother
  LM_line <- as(LM_poly, "SpatialLinesDataFrame")
  
  if (west) {
    # Define basin boundary
    bnd_west <- soap_prep(LM_line[LM_line$basin == "west", ],
                          bound_val = elev[elev$basin == "west", ]$elev,
                          check = FALSE)
    # Create knot grid
    knots_west <- make_knots(LM_poly[LM_poly$basin == "west", ], 
                             spacing = knot_spacing, seed = r_seed,
                             use_data = FALSE, plot = FALSE,
                             data_coords = bath_data[bath_data$basin == "west", ])
    if (check_soap) {
      ok <- soap_check(bnd_west, 
                       knots = knots_west, 
                       data = bath_data[bath_data$basin == "west", ])
      if (!ok) stop("Soap check for west basin didn't pan out.")
      need_key("West basin soap film check.  Press [enter] to continue.")
    }

    # Fit soap-film smoother...
    w_soap <- gam(depth ~ s(x, y, bs = "so", xt = list(bnd = bnd_west)),
                  data = bath_data[bath_data$basin == "west", ], 
                  method = "REML", knots = knots_west)
    if (print_summary) {
      print(summary(w_soap))
      need_key("West soap film model summary above.  Press [enter] to continue.")
    }
    if (gam_check) {
      gam.check(w_soap)
      need_key("West soap film model check above.  Press [enter] to continue.")
    }
    if (plot_smooth) {
      plot(w_soap)
      need_key()
    }
    
    if (pred_rast) {
      cat("Generating predicted bathymetry raster for west basin.  Please wait...\n")
      w_soap_pred <- as.data.frame(coordinates(LM_pts[LM_pts$basin == "west", ]))
      w_soap_pred <- transform(w_soap_pred, 
                               depth = predict(w_soap, newdata = w_soap_pred))
      w_soap_pred <- transform(w_soap_pred, 
                               depth_se = predict(w_soap, newdata = w_soap_pred, 
                                                  se.fit = TRUE)$se.fit)
      w_soap_pred_rast <- raster(extent(LM), res = 100)
      w_soap_pred <- na.omit(w_soap_pred)
      w_soap_pred_rast <- rasterize(w_soap_pred[, c("x", "y")], 
                                    w_soap_pred_rast, 
                                    as.numeric(w_soap_pred[, "depth"]))
      w_soap_pred_se <- rasterize(w_soap_pred[, c("x", "y")], 
                                  w_soap_pred_rast, 
                                  w_soap_pred[, "depth_se"])
    }   
    
  }
  
  if (east) {
    # Define basin boundary
    bnd_east <- soap_prep(LM_line[LM_line$basin == "east", ],
                          bound_val = elev[elev$basin == "east", ]$elev,
                          check = FALSE)
    # Create knot grid
    knots_east <- make_knots(LM_poly[LM_poly$basin == "east", ], 
                             spacing = knot_spacing, seed = r_seed,
                             use_data = FALSE, plot = FALSE, 
                             data_coords = bath_data[bath_data$basin == "east", ])
    if (check_soap) {
      ok <- soap_check(bnd_east, 
                       knots = knots_east, 
                       data = bath_data[bath_data$basin == "east", ])
      if (!ok) stop("Soap check for east basin didn't pan out.")
      need_key("East basin soap film check.  Press [enter] to continue.")
    }   
    
    # Fit soap-film smoother...
    e_soap <- gam(depth ~ s(x, y, bs = "so", xt = list(bnd = bnd_east)),
                  data = bath_data[bath_data$basin == "east", ], 
                  method = "REML", knots = knots_east)
    if (print_summary) {
      print(summary(e_soap))
      need_key("East soap film model summary above.  Press [enter] to continue.")
    }
    if (gam_check) {
      gam.check(e_soap)
      need_key("East soap film model check above.  Press [enter] to continue.")
    }
    if (plot_smooth) {
      plot(e_soap)
      need_key()
    }
    
    if (pred_rast) {
      cat("Generating predicted bathymetry raster for east basin.  Please wait...\n")
      e_soap_pred <- as.data.frame(coordinates(LM_pts[LM_pts$basin == "east", ]))
      e_soap_pred <- transform(e_soap_pred, 
                               depth = predict(e_soap, newdata = e_soap_pred))
      e_soap_pred <- transform(e_soap_pred, 
                               depth_se = predict(e_soap, newdata = e_soap_pred, 
                                                  se.fit = TRUE)$se.fit)
      e_soap_pred_rast <- raster(extent(LM), res = 100)
      e_soap_pred <- na.omit(e_soap_pred)
      e_soap_pred_rast <- rasterize(e_soap_pred[, c("x", "y")], 
                                    e_soap_pred_rast, 
                                    as.numeric(e_soap_pred[, "depth"]))
      e_soap_pred_se <- rasterize(e_soap_pred[, c("x", "y")], 
                                  e_soap_pred_rast, 
                                  e_soap_pred[, "depth_se"])
    }

  }
  
  if (pred_rast && east && west) {
    lake_soap <- mosaic(e_soap_pred_rast, w_soap_pred_rast, fun = mean)
    lake_se <- mosaic(e_soap_pred_se, w_soap_pred_se, fun = mean)
  } else if (west) {
    lake_soap <- trim(w_soap_pred_rast, values = c(NA, NaN))
    lake_se <- trim(w_soap_pred_se, values = c(NA, NaN))
  } else {
    lake_soap <- trim(e_soap_pred_rast, values = c(NA, NaN))
    lake_se <- trim(e_soap_pred_se, values = c(NA, NaN))
  }
  
  if (pred_rast) {
    plot(lake_soap, col = viridis(256))
    need_key("Plotted predicted bathymetry raster.\nPress [enter] for standard error raster.\n")
    plot(lake_se, col = viridis(256))
    
    if (export_gis) {
      # Define coordinate system
      bath_rasts <- stack(lake_soap, lake_se)
      projection(bath_rasts) <- CRS("+proj=utm +zone=18 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
      writeRaster(bath_rasts[[1]], filename = "./Output/LM_bath.tif")
      writeRaster(bath_rasts[[2]], filename = "./Output/LM_bath_se.tif")
    }
    
  }
  
  out <- list()
  if (west) out$west_model = w_soap
  if (east) out$east_model = e_soap
  if (pred_rast) out$pred_rast = bath_rasts
  attr(out, "seed") <- r_seed
  return(out)
  
}
