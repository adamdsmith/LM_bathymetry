# Calculate NDWI
calc_NDWI <- function(cir_raster) {
  green <- cir_raster[[2]]
  nir <- cir_raster[[4]]
  NDWI <- (green - nir) / (green + nir)
  NDWI
}

# Fill variable down column until non-missing value detected
fill_var <- function(var, missing = "") {
  good <- var != missing
  fill_vals <- var[good]
  fill_vals[cumsum(good)]
}

need_key <- function(prompt = "Press [enter] to continue") invisible(readline(prompt=prompt))

convert_length <- function(data, from = c("in", "ft", "cm", "m"), to = c("in", "ft", "cm", "m")) {
  conv <- 1
  from <- match.arg(from)
  to <- match.arg(to)
  if (from == "ft") {
    from <- "in"
    data <- data * 12
  }
  if (from == "m") {
    from <- "cm"
    data <- data * 100
  }
  if (to == "ft") {
    to <- "in"
    conv <- 12
  }
  if (to == "m") {
    to <- "cm" 
    conv <- 100
  }
  out <- switch(from,
         'in' = grid::convertX(grid::unit(data, "in"), to), 
         'cm' = grid::convertX(grid::unit(data, "cm"), to))
  out <- as.numeric(out) / conv
  return(out)
}

# Function to use Python raster to polygon file from GDAL
# modified from https://gist.github.com/Pakillo/c6b076eceb0ef5a70e3b
# Forked for posterity at https://gist.github.com/adamdsmith/3f3d2245f7e7f2ededad
polygonizer <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile', 
                        pypath=NULL, readpoly=TRUE, quietish=TRUE) {
  # x: an R Raster layer, or the file path to a raster file recognised by GDAL
  # outshape: the path to the output shapefile (if NULL, a temporary file will be created)
  # gdalformat: the desired OGR vector format
  # pypath: the path to gdal_polygonize.py (if NULL, an attempt will be made to determine the location
  # readpoly: should the polygon shapefile be read back into R, and returned by this function? (logical)
  # quietish: should (some) messages be suppressed? (logical)
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- "C:\\OSGeo4W64\\bin\\gdal_polygonize"
  }
  ## The line below has been commented:
  # if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.") 
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)) 
      stop(sprintf('File already exists: %s', 
                   toString(paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.asc')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  
  ## Now 'python' has to be substituted by OSGeo4W
  #system2('C:\\Python27\\ArcGISx6410.3\\python.exe',
  system2('C:\\OSGeo4W64\\OSGeo4W.bat',
          args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"', 
                        pypath, rastpath, gdalformat, outshape)))
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quietish)
    return(shp) 
  }
  return(NULL)
}

# Get rid of small holes typically representing individual or small stands of cypress
# Default cutoff is 250 m2
remove_holes <- function(SPlyDF, min_area = 250, name = "poly") {
  library(sp)
  p <- slot(SPlyDF, "polygons") 
  holes <- lapply(p, function(x) sapply(slot(x, "Polygons"), slot, "hole") &
                    (sapply(slot(x, "Polygons"), slot, "area") < min_area)) 
  res <- lapply(1:length(p), function(i) slot(p[[i]], "Polygons")[!holes[[i]]]) 
  IDs <- row.names(SPlyDF) 
  out <- SpatialPolygonsDataFrame(
    SpatialPolygons(lapply(1:length(res), function(i) 
      Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(SPlyDF))),
    data = data.frame(ID = name))
}

# Modify waterData function to allow single day query
importDVs_1day <- function (staid, code = "00060", stat = "00003", sdate = "1851-01-01", 
                       edate = as.Date(Sys.Date(), format = "%Y-%m-%d")) 
{
  base_url <- "http://waterservices.usgs.gov/nwis/dv?"
  url <- paste(base_url, "site=", staid, "&parameterCd=", code, 
               "&statCd=", stat, sep = "")
  url <- paste(url, "&startDt=", sdate, "&endDt=", edate, sep = "")
  doc <- XML::xmlTreeParse(url, getDTD = FALSE, useInternalNodes = TRUE)
  r <- XML::xmlRoot(doc)
  val <- as.numeric(XML::xmlValue(r[[2]][[3]][[1]]))
  Attribute <- XML::xmlApply(r[[2]][[3]], XML::xmlAttrs)
  N <- length(val)
  NoDataValue <- XML::xmlValue(r[["timeSeries"]][["variable"]][["NoDataValue"]])
  NoDataValue <- as.integer(NoDataValue)
  dates <- vector(mode = "character", length = 1)
  dates <- as.Date(dates, "%Y-%m-%d")
  df <- data.frame(staid, val)
  df
}

# Generate potential knot locations
make_knots <- function(SPlyDF, spacing = 1500, 
                       random_start = TRUE, seed = NULL, 
                       data_coords = NULL,
                       knot_dist = max(spacing)/2,
                       clear_bnd = TRUE,
                       plot = TRUE,
                       use_data = TRUE, n_knots = 20) {
  # library(sp)
  # library(rgeos)
  # library(raster)
  # library(dplyr)
  ## library(rmapshaper)
  stopifnot(is.numeric(spacing) || is.integer(spacing))
  
  if (!use_data) {
    if (length(spacing) == 1) spacing <- rep(spacing, 2)
    if (random_start) {
      if (!is.null(seed)) set.seed(seed)
      shift <- runif(2, 0, spacing) 
    } else shift <- c(0, 0)
    
    # Make SpatialPointsDF based on extent of SPlyDF and spacing
    ext <- extent(SPlyDF)
    x <- seq(ext@xmin - spacing[1], ext@xmax + spacing[1], by = spacing[1]) - shift[1]
    y <- seq(ext@ymin - spacing[2], ext@ymax + spacing[2], by = spacing[2]) - shift[2]
    knots <- expand.grid(x = x, y = y)
    knotsSP <- SpatialPoints(knots, proj4string = CRS(proj4string(SPlyDF)))
    if (clear_bnd) {
      lines <- as(SPlyDF, "SpatialLines")
      too_close <- gWithinDistance(knotsSP, lines, dist = knot_dist/4, byid=TRUE)
      knots[too_close, ] <- NA
    }
    outside <- which(is.na((knotsSP %over% as(SPlyDF, "SpatialPolygons"))))
    knots[outside, ] <- NA
    
    if (!is.null(data_coords) && c("x", "y") %in% names(data_coords)) {
      prj <- proj4string(SPlyDF)
      dataSP <- SpatialPointsDataFrame(data_coords[, c("x", "y")], 
                                       dplyr::select(data_coords, -x, -y),
                                       proj4string = CRS(prj))
      knot_buff <- gBuffer(SpatialPoints(na.omit(knots), proj4string = CRS(prj)), 
                           byid = TRUE, width = knot_dist)
      knot_data <- knot_buff %over% dataSP
      keep_knot <- rep(NA, nrow(knots))
      keep_knot[-outside] <- 1
      knots <- knots * keep_knot
    }
  } else {
    if (is.null(data_coords)) stop("You have to specify input data to use it to generate knots.")
    library(spsurvey)
    prj <- proj4string(SPlyDF)
    crds <- unique(data_coords[, c("x", "y")])
    dataSP <- SpatialPointsDataFrame(crds, 
                                     data = data.frame(id = 1:nrow(crds)),
                                     proj4string = CRS(prj))
    dsgn <- list(None=list(panel=c(PanelOne=n_knots), seltype="Equal"))
    sites <- grts(design=dsgn,
                  DesignID="knots",
                  type.frame="finite",
                  src.frame="sp.object",
                  sp.object=dataSP,
                  shapefile = FALSE)
    knots <- data.frame(x = sites@data$xcoord, y = sites@data$ycoord)
  }
  
  if (plot) {
    plot(SPlyDF)
    if (!is.null(data_coords)) with(data_coords, points(x, y))
    points(knots$x, knots$y, pch = "+", col = "red")
  }
  
  return(na.omit(knots))
  
}

# Check soap film surface
# Modified slightly from David Miller
# https://github.com/dill/soap_checker
soap_check <- function(bnd, knots=NULL, data=NULL, plot=TRUE,
                       tol=sqrt(.Machine$double.eps)){
  
  # library(rgeos)
  ## check that the boundary makes sense
  # check that boundary is a list
  stopifnot(is.list(bnd))
  
  # check that the boundary (or boundary part) have x and y elements
  lapply(bnd, function(x) stopifnot(c("x","y") %in% names(x)))
  # each boundary part must have at least 4 elements!
  lapply(bnd, function(x) stopifnot(length(x$x)>3, length(x$y)>3))
  
  # check that the boundary loops are actually loops
  check_ends <- function(x, tol){
    all.equal(c(x$x[1], x$y[1]),
              c(x$x[length(x$y)], x$y[length(x$y)]),tolerance=tol)
  }
  end_check <- unlist(lapply(bnd, check_ends, tol=tol))
  end_check_logical <- is.character(end_check)
  if(any(end_check_logical)){
    stop(paste("Boundary loop(s)",which(end_check_logical),
               "don't have identical start & end points",collapse=" "))
  }
  
  # Retain only x-y columns for the test
  bnd <- lapply(bnd, function(x) x[c("x", "y")] )
  
  islands <- FALSE
  # check for intersections
  if(length(bnd)>1){
    inds <- combn(1:length(bnd),2)
    
    ## make the bnds into polys here
    make_bnd_poly <- function(bnd){
      
      bnd$x[length(bnd$x)] <-  bnd$x[1]
      bnd$y[length(bnd$y)] <-  bnd$y[1]
      SpatialPolygons(list(Polygons(list(Polygon(bnd)),ID=1)))
    }
    bnd_poly <- lapply(bnd, make_bnd_poly)
    
    # function to see if two polygons intersect
    intersects <- function(this.ind, bnd){
      poly1 <- bnd[[this.ind[1]]]
      poly2 <- bnd[[this.ind[2]]]
      
      gIntersects(poly1, poly2)
    }
    # apply over all the combinations
    inter <- apply(inds, 2, intersects, bnd=bnd_poly)
    
    if(any(inter)){
      # get the index for the prospective "outer" loop
      outer_ind <- which.max(unlist(lapply(bnd_poly, gArea)))
      outer_bnd <- bnd_poly[[outer_ind]]
      
      other_bnd <- bnd_poly
      other_bnd[[outer_ind]] <- NULL
      
      # is everything else inside that?
      islands <- unlist(lapply(other_bnd, gWithin, spgeom2=outer_bnd))
      
      if(!all(islands)){
        stop(paste("Polygon parts",
                   paste0(apply(inds[,inter, drop=FALSE], 2, paste0,
                                collapse=" and "),
                          collapse=", "), "intersect"))
      }
      
      islands <- all(islands)
    }
  }
  
  ## plot what the boundary is
  # highlighting the area to be modelled
  if(plot){
    # colourblind-safe colours from colorbrewer2 "qualitative" map
    red <- "#d95f02"
    # if the boundary is only 1 part, plotting is rather easier
    if(!islands){
      plot(bnd[[1]], type="l", xlab = "x", ylab = "y", 
           main="Red indicates soap film surface", asp=1)
      lapply(bnd, polygon, col=red)
    }else{
      outer_bnd <- bnd[[outer_ind]]
      other_bnd <- bnd
      other_bnd[[outer_ind]] <- NULL
      plot(outer_bnd, type="n", xlab = "x", ylab = "y", 
           main="Red indicates soap film surface", asp=1)
      # plot the outer loop
      polygon(outer_bnd, col=red)
      # plot the other polygons on top in white
      lapply(other_bnd, polygon, col="white")
    }
  }
  
  # function to check if points are inside the boundary
  point_check <- function(bnd, x, y, type){
    
    if(length(bnd)>1){
      # inSide doesn't deal with edge points very well
      # but does handle multiple rings better
      inout <- inSide(bnd, x, y)
    }else{
      # use sp::point.in.polygon
      # see ?point.in.polygon for returned codes, 1 is inside
      pip <- function(bnd, x, y){
        point.in.polygon(x, y, bnd$x, bnd$y)==1
      }
      # apply over the parts of the polygon
      inout <- pip(bnd[[1]], x, y)
    }
    if(!all(inout)){
      stop(paste(type, paste(which(!inout),collapse=", "),
                 "are outside the boundary."))
    }
  }
  
  ## check the knots
  if(!is.null(knots)){
    # check that the points have x and y elements
    stopifnot(c("x","y") %in% names(knots))
    point_check(bnd, knots$x, knots$y, "Knots")
    if(plot) points(knots, col="#1b9e77", pch=19)
  }
  
  ## check the data
  if(!is.null(data)){
    # check that the points have x and y elements
    stopifnot(c("x","y") %in% names(data))
    point_check(bnd, data$x, data$y, "Data points")
    if(plot) points(data$x, data$y, col="#7570b3", pch=19)
  }
  
  return(TRUE)
}

soap_prep <- function(SpLinDF, bound_val = 0, check = TRUE) {
  crds <- coordinates(SpLinDF)
  bnd <- lapply(crds[[1]], function(ply) {
    tmp <- list(x = ply[, 1], y = ply[, 2], f = rep(bound_val, nrow(ply)))
  })
  
  if (check) {
    soap_check(bnd)
  }
  
  return(bnd)
}

navd_depth <- function(depth_df, in_units =  c("in", "ft", "m", "cm"), 
                       out_units = c("ft", "m", "cm", "in")) {
  
  if (!("pacman" %in% installed.packages())) install.packages("pacman", quiet = TRUE)
  pacman::p_load(dplyr, lubridate, waterData)
  
  if (!any(c("date", "depth", "basin") %in% names(depth_df))) 
    stop("Input data frame must contain columns 'date', 'depth', and 'basin'")
  in_units <- match.arg(in_units)
  out_units <- match.arg(out_units)
  
  # Check if date in acceptable format
  if (is.character(depth_df$date)) {
    if (!any(sapply(gregexpr("-|/", depth_df$date), function(x) identical(as.vector(x), c(5L,8L)))))
      stop("At least one date didn't parse correctly.  Check that they are formatted YYYY-MM-DD.")
  } else {
    if (!any(sapply(gregexpr("-|/", as.character(depth_df$date)), function(x) identical(as.vector(x), c(5L,8L)))))
      stop("At least one date didn't parse correctly.  Check that they are formatted YYYY-MM-DD.")
  }
  
  stations <- c("0208458892", "0208458893"); basins <- c("west", "east")
  
  # Surface elevation available only in feet, so convert
  # Get date range
  start <- as.character(min(ymd(depth_df$date)))
  end <- as.character(max(ymd(depth_df$date)))
  # Gages are -2.0 feet above NAVD88
  gage <- lapply(stations, function(stn) {
    b <- basins[which(stations == stn)]
    if (start == end) {
      out <- importDVs_1day(stn, code = "62615", stat = "00003",
                       sdate = start, edate = end) %>%
        mutate(date = start)
    } else {
      out <- importDVs(stn, code = "62615", stat = "00003",
                       sdate = start, edate = end) %>%
        mutate(date = as.character(dates))
    }
    out <- out %>%
      mutate(basin = b,
             # Default units: feet
             wl = convert_length(val, from = "ft", to = out_units))
  })
  
  gage <- as.data.frame(do.call(rbind, gage)) %>% dplyr::select(date, basin, wl) 
  
  out <- depth_df %>% 
    mutate(depth = convert_length(depth, in_units, out_units) - 
             left_join(depth_df, gage, by = c("date", "basin"))$wl)

  return(out)

}

LM_lake_elev <- function(start_date = "2014-10-05",
                         end_date = start_date) {
  
  if (lubridate::ymd(end_date) < lubridate::ymd(start_date)) 
    stop("`end_date` must be the same or after `start_date`.")
  
  # Lake elevation in NAVD88 (feet)
  # Used to set the NADV88 elevation for the lake boundary
  stations <- c("0208458892", "0208458893")
  basins <- c("west", "east")
  
  elev <- lapply(stations, function(stn) {
    b <- basins[which(stations == stn)]
    if (identical(start_date, end_date)) {
      out <- importDVs_1day(stn, code = "62615", stat = "00003", 
                            sdate = start_date, edate = end_date)
      out$dates <- start_date
    } else {
      out <- waterData::importDVs(stn, code = "62615", stat = "00003", 
                                  sdate = start_date, edate = end_date)
    }
    
    out <- out %>%
      mutate(basin = b,
             elev = val,
             date = as.character(dates)) %>% # convert to meters
      dplyr::select(date, basin, elev)
  })
  
  elev <- as.data.frame(do.call(rbind, elev))
  
  return(elev)
  
}

points_to_line <- function(data, lon, lat, id_field = NULL, sort_field = NULL) {
  
  # Convert to SpatialPointsDataFrame
  coordinates(data) <- c(lon, lat)
  
  # If there is a sort field...
  if (!is.null(sort_field)) {
    if (!is.null(id_field)) {
      data <- data[order(data[[id_field]], data[[sort_field]]), ]
    } else {
      data <- data[order(data[[sort_field]]), ]
    }
  }
  
  # If there is only one path...
  if (is.null(id_field)) {
    
    lines <- SpatialLines(list(Lines(list(Line(data)), "id")))
    
    return(lines)
    
    # Now, if we have multiple lines...
  } else if (!is.null(id_field)) {  
    
    # Split into a list by ID field
    paths <- sp::split(data, data[[id_field]])
    
    sp_lines <- SpatialLines(list(Lines(list(Line(paths[[1]])), "line1")))
    
    # I like for loops, what can I say...
    for (p in 2:length(paths)) {
      id <- paste0("line", as.character(p))
      l <- SpatialLines(list(Lines(list(Line(paths[[p]])), id)))
      sp_lines <- spRbind(sp_lines, l)
    }
    
    return(sp_lines)
  }
}

CB_classify <- function(sav_vals) {
  pacman::p_load(dplyr)
  # Calculate estimate of total area covered by SAV
  # i.e., estimated proportion per cell * cell area
  # cells are 1 ha, so simple sum = total ha
  est_ha <- sum(sav_vals)
  CB_bins <- c(0, 0.1, 0.4, 0.7, 1.01)
  CB_labels = data.frame(class = c("1", "2", "3", "4"),
                         label = c("Very sparse (< 10%)", 
                                   "Sparse (10 - 40%)", 
                                   "Moderate (40 - 70%)", 
                                   "Dense (> 70%)"),
                         stringsAsFactors = FALSE)
  sav_mat <- cut(sav_vals, breaks = CB_bins, labels = c("1", "2", "3", "4"), right = FALSE)
  sav_df <- data.frame(table(sav_mat)) %>% 
    mutate(class = as.character(sav_mat),
           area_ha = Freq) %>%
    dplyr::select(class, area_ha)
  out <- left_join(CB_labels, sav_df, by = "class") %>%
    mutate(area_ha = ifelse(is.na(area_ha), 0, area_ha),
           LM_SAV_ha = est_ha)
  return(out)

}

every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}