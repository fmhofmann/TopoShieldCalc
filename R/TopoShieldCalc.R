################################################
### Compute topographic shielding with a DEM ###
################################################

# 1. Define the functions #

## Make sure that the terra and readODS packages are installed ##

# 1.1 Create a function to load the required geodata interactively #

#' Load geodata for topographic shielding factor calculation
#'
#' This function allows for loading geodata (digital elevation model (DEM) and shapefile) for topographic shielding factor calculation.
#' @param radius Numeric. The radius (in metres) around the sampling site relevant for shielding factor calculation.  
#' @details
#' The function adds two objects to the global environment: "dem" and "points", i.e., the digital elevation model (DEM) and the shapefile with the sampling site. 
#' For this step, the [terra::rast()] and [terra::vect()] commands are used, respectively.
#' As soon as the "Select the shapefile" message appears, the shapefile (.shp) can be interactively selected. 
#' The function then performs a quality check of the input-shapefile with [terra::crs()]. 
#' If the unit of the coordinate reference system (CRS) is NOT in metres, the function returns an error message and the function execution stops.
#' The function then checks whether the attribute table of the shapefile contains all necessary data for shielding factor calculation.
#' Note that the attribute table of the input file must contain the following columns: Name, Strike_deg, Dip_deg, and BouldHt.
#' If this is not the case, an error message is returned and the function execution stops.
#' 
#' The function will then prompt the user to interactively load the digital elevation model (DEM). 
#' The function then performs a quality check of the input-DEM with [terra::crs()].
#' If the unit of the coordinate reference system (CRS) is NOT in metres, the function returns an error message and the function execution stops.
#' The function then assesses with [terra::crs()] whether the input-DEM and the input-shapefile have the same CRS. 
#' If this is not the case, an error message is returned and the function execution stops.
#' The function subsequently assesses with [terra::xmin()], [terra::xmax()], [terra::ymin()], and [terra::ymax()] whether the DEM covers the whole area relevant for shielding factor calculation.
#' If this is not the case, an error message is returned and the function execution stops.
#' If the input-shapefile and the input-raster meet the requirements for shielding factor calculation, the function returns the following message: "The input-shapefile and the input-DEM have the same coordinate reference system (CRS) and the unit of the CRS is metres. 
#' They are suitable for topographic shielding factor calculation".
#' 
#' The function adds two objects to the global environment: "dem" and "points". 
#' These are of class "SpatRast" and "SpatVect", respectively. 
#' For this step, the [terra::rast()] and [terra::vect()] commands are used, respectively.
#' @examples 
#' load_geodata(radius = 10000)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
load_geodata = function(radius){ 
  message("Select the shapefile")
  points = terra::vect(file.choose()) # Load the shapefile
  if (terra::linearUnits(points) == 0){
    stop("The linear unit of the coordinate reference system of the input-shapefile needs to be metres")
  } else {}
  colnames = colnames(terra::as.data.frame(points)) # Retrieve the column names
  if ("Name" %in% colnames & "Strike_deg" %in% colnames & "Dip_deg" %in% colnames & "BouldHt" %in% colnames){ # Check if all columns exist
    message("The attribute table of the input-shapefile contains all necessary columns for shielding factor calculation")
  } else {
    stop("The input-shapefile does not contain all columns needed for shielding factor calculation: Name, Strike_deg, Dip_deg, and BouldHt")
  }
  message("Select the DEM")
  dem = terra::rast(file.choose()) # Load the dem
  if(terra::linearUnits(dem) == 0){
    stop("The linear unit of the input-DEM must be metres")
  } else {}
  if (terra::crs(points, proj = TRUE) != terra::crs(dem, proj = TRUE)){
    stop("The input-DEM and the input-shapefile do not have the same coordinate reference system")
  } else {}
  if(terra::xmin(points) - radius - 500 < terra::xmin(dem) || terra::xmax(points) + radius + 500 > terra::xmax(dem) || terra::ymin(points) - radius - 500 < terra::ymin(dem) || terra::ymax(points) + radius + 500 > terra::ymax(dem)){ # Checks if dem covers the whole area relevant for shielding factor calculation
    stop("The DEM does not cover the whole area relevant for shielding factor calculation")
  } else {}  
  message("The input-shapefile and the input-DEM fulfil all requirements")
  list = list(dem = dem, points = points)
  return(list2env(list, envir = globalenv()))
}

# 1.2 Create a function to crop the dem #

#' Crop the input-DEM for shielding factor calculation
#'
#' This function allows for cropping the selected DEM for shielding factor calculation.
#' @param dem The input-DEM (an object of type "SpatRaster", created with [terra::rast] or [TopoShieldCalc::load_geodata]).
#' @param point The input-shapefile (object of type "SpatVect", created with [terra::vect] or [TopoShieldCalc::load_geodata]).
#' @param plot Logical. Should the results be plotted?
#' @param radius Numeric. The radius around the sampling site(s) relevant for shielding factor calculation in metres.
#' @details
#' This function extracts the portion of the DEM which is relevant for shielding factor calculation for a single sampling site. 
#' This depends on the selected radius around the sampling site. 
#' For this step, [terra::crop] is used. 
#' Please note that 500 metres are added to radius to make sure that all all raster cells relevant for shielding factor calculation are extracted. 
#' If plot = TRUE, a map is created with [terra::plot]. 
#' This map shows both the DEM and the input-point.
#' The cropped DEM (dem) is added to the global environment.
#' @examples 
#' crop_dem(dem, point, radius = 10000, plot = TRUE)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
crop_dem = function(dem, # SpatRaster object
                    point, # SpatVector object
                    radius, # Numeric: maximum distance around the sampling sites
                    plot = FALSE){ # Logical. If plot = TRUE, plot the cropped DEM
  mask = terra::ext(terra::xmin(point) - radius - 500, # Create a mask to crop the DEM
                    terra::xmax(point) + radius + 500,
                    terra::ymin(point) - radius - 500,
                    terra::ymax(point) + radius + 500)
  dem_2 = terra::crop(dem, # Crop the DEM with the mask
                      mask,
                      touches = TRUE)
  if (plot == TRUE){
    terra::plot(dem_2) # Plot the result
    terra::points(point, col = "red")
  } else {}
  list = list(dem_2 = dem_2, point = point)
  return(list2env(list, envir = globalenv())) # Return dem_2
}

# 1.4 Function to extract the x-, y-, and z-coordinates of the input-points #

#' Get the coordinates of the sampling site
#' 
#' This function allows for obtaining the \emph{x}-, \emph{y}-, and \emph{z}-coordinates of a sampling site.  
#' @param dem The input-DEM (an object of type "SpatRaster", created with [terra::rast] or [TopoShieldCalc::load_geodata]).
#' @param point The input-shapefile (object of type "SpatVect", created with [terra::vect] or [TopoShieldCalc::load_geodata]).
#' @param boulder_height Numeric. The boulder height in metres.
#' @details
#' For the first step, the \emph{x}- and \emph{y}-coordinates of the sampling site are extracted with the aid of [terra::geom].
#' For the second step, the \emph{z}-coordinate is retrieved from dem with [terra::extract] (method = "bilinear").
#' The boulder height in metres is added to the \emph{z}-coordinate.
#' The \emph{x}-, \emph{y}-, and \emph{z}-coordinates of the sampling site (vectors) are finally added to the global environment (point_x, point_y, and point_z, respectively).
#' @examples 
#' point_xyz(dem, point, boulder_height = 0.3)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
point_xyz = function(dem, # SpatRaster. The cropped DEM
                     point, # SpatVect. Sampling sites
                     boulder_height){ # Numeric. Boulder height in m
  point_x = terra::geom(point, df = TRUE)
  point_x = point_x$x
  point_y = terra::geom(point, df = TRUE)
  point_y = point_y$y
  point_z = terra::extract(dem, y = point, method = "bilinear")
  colnames(point_z) = c("ID","z") # Rename the columns in points_z
  point_z = as.numeric(point_z$z)
  point_z = point_z + boulder_height
  list = list(point_x = point_x, point_y = point_y, point_z = point_z)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

# 1.5 Function for deriving a DEM relative to the z-coordinate of an input-point #

#' Elevation with respect to the sampling site 
#' 
#' This function allows for creating a raster with the elevation in metres above sea-level with respect to the \emph{z}-coordinate of the sampling site. 
#' @param dem The input-DEM (an object of type "SpatRaster", created with [terra::rast] or [TopoShieldCalc::load_geodata]).
#' @param point_z Numeric. The \emph{z}-coordinate, i.e., the elevation of the sampling site. 
#' @param plot Logical. Create a map with the raster using the [terra::plot] function? Default to FALSE.
#' @details
#' To make this calculation as fast as possible, the cropped digital elevation model should preferably used.
#' For this step, the [TopoShieldCalc::crop_dem] function can be used. 
#' point_z can be obtained with [TopoShieldCalc::point_xyz].
#' For the first step, the function converts the input-DEM into a dataframe with the aid of [terra::as.data.frame].
#' For the second step, the difference between the elevation of the sampling site and each cell of the DEM is calculated.
#' Note that if z < 0, the value is replaced by point_z.
#' If cells in the input-DEM contain "NA" values, an error message is returned.
#' An object of type "SpatRaster" (dem_rel) is added to the environment. 
#' This raster is populated with the difference between the elevation of the sampling site and each cell of the DEM through the use of [terra::init]. 
#' The raster (dem_rel; object of type "SpatRaster") is finally added to the global environment.
#' If plot = TRUE, a map is created with [terra::plot]. 
#' @examples 
#' dem_relative(dem = dem, point_z = 1356, plot = FALSE)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
dem_relative = function(dem, # SpatRaster. The cropped DEM
                        point_z, # Numeric. The z-coordinate of the input-point
                        plot = FALSE){ # Numeric. z-coordinates of the input-point
  dem_dataframe = terra::as.data.frame(dem, xy = TRUE, na.rm = FALSE) # Convert the DEM into a dataframe
  colnames(dem_dataframe) = c("x","y","z") # Rename the columns of the dataframe
  z = dem_dataframe$z
  if(NA %in% z){ # If z contains 'NA'...
    stop(paste("Site #",point$Name,": The portion of the DEM relevant for shielding factor calculation contains NA values. Select a DEM that entirely covers the area relevant for shielding factor calculation."), call. = FALSE)
  } else {}
  z = replace(z,
              z < point_z, # If dem_dataframe$z is smaller than point_z
              point_z) # Replace the value with point_z
  z = z - point_z
  empty_raster = terra::rast(res = terra::res(dem), # Create an empty raster
                             nlyrs = 1,
                             extent = terra::ext(dem),
                             crs = terra::crs(dem))
  dem_rel = terra::init(empty_raster, fun = z) # Populate the raster with z
  if(plot == TRUE){
    terra::plot(dem_rel)
  } else {}
  list = list(dem_rel = dem_rel)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

# 1.6 Function for the creation of a raster with the distance from the sampling locality to the midpoints of the raster cells #

#' Distance from the the sampling site 
#' 
#' This function allows for generating a raster showing the distance between the sampling site and the cells of the input-DEM.
#' The distance is needed for the calculation of the elevation angle for each raster cell.
#' [TopoShieldCalc::raster_elevation] performs this task.
#' @param dem The input-DEM (an object of type "SpatRaster", created with [terra::rast] or [TopoShieldCalc::load_geodata]).
#' @param point_x Numeric. The \emph{x}-coordinate of the sampling site. 
#' @param point_y Numeric. The \emph{y}-coordinate of the sampling site. 
#' @details
#' To perform this calculation as fast as possible, the cropped digital elevation model should preferably be used.
#' For this step, the [TopoShieldCalc::crop_dem] function can be used.
#' point_z can be obtained with [TopoShieldCalc::point_xyz].
#' The raster (raster_distance; object of type "SpatRaster") is finally added to the global environment.
#' @examples 
#' raster_dist(dem = dem, point_x = 414416.5, point_y = 5316445.8)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
raster_dist = function(dem, # SpatRaster. The cropped DEM.
                       point_x, # Numeric. x-coordinate of the sampling locality
                       point_y){ # Numeric. y-coordinate of the sampling locality
  dem_dataframe = terra::as.data.frame(dem, xy = TRUE, na.rm = FALSE) # Convert the DEM into a dataframe
  colnames(dem_dataframe) = c("x","y","z") # Rename the columns of the dataframe
  x = dem_dataframe$x
  y = dem_dataframe$y
  distance = sqrt((x - point_x)^2 + (y - point_y)^2)
  empty_raster = terra::rast(res = terra::res(dem), # Create an empty raster
                             nlyrs = 1,
                             extent = terra::ext(dem),
                             crs = terra::crs(dem))
  raster_distance = terra::init(empty_raster, fun = distance) # Populate the raster with the values in distance
  list = list(raster_distance = raster_distance)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

# 1.7 Function to derive a raster with elevation angle relative to the sampling site #

#' Create a raster showing the elevation angle with respect to the sampling site 
#' 
#' This function allows for generating a raster showing the elevation angle (degrees) with respect to the sampling site.
#' 
#' @param dem_rel The raster (an object of type "SpatRaster", created with [TopoShieldCalc::dem_relative]) showing the elevation angle with respect to the sampling site.
#' @param raster_distance The raster (an object of type "SpatRaster", created with [TopoShieldCalc::raster_dist]) showing the distance between the raster cells and the sampling site.
#' @param plot Logical. Should the raster with the elevation angle be plotted? Default to FALSE.
#' @details
#' For the first step, the function calculates the difference in elevation for each cell of the input-DEM.
#' The slope of the with respect to the sampling site is then calculated for each raster cell.
#' Elevation angles are calculated as follows: atan(slope) * 180 / pi.
#' The function adds an object (raster_elev) of class "SpatRaster" to the environment.  
#' @examples 
#' raster_elevation(dem_rel = dem_relative, raster_distance = raster_distance, plot = TRUE)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
raster_elevation = function(dem_rel,
                            raster_distance,
                            plot = FALSE){
  dataframe_dem_relative = terra::as.data.frame(dem_rel)
  colnames(dataframe_dem_relative) = "z_rel"
  elevation_relative = dataframe_dem_relative$z_rel
  dataframe_distance = terra::as.data.frame(raster_distance)
  colnames(dataframe_distance) = "dist"
  distance = dataframe_distance$dist
  slope = elevation_relative / distance 
  elevation = atan(slope) * 180 / pi
  empty_raster = terra::rast(res = terra::res(raster_distance), # Create an empty raster
                             nlyrs = 1,
                             extent = terra::ext(raster_distance),
                             crs = terra::crs(raster_distance))
  raster_elev = terra::init(empty_raster, fun = elevation) # Populate the raster with the values in elevation
  if(plot == TRUE){
    terra::plot(raster_elev)
  } else {}
  list = list(raster_elev = raster_elev)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

# 1.8 Define a function to compute pairs of elevation and azimuth angles #

#' Calculate elevation angles for each azimuth  
#' 
#' This function allows for calculating elevation angles for each azimuth (1:360). 
#' @param radius Numeric. The radius around the sampling site (in metres) relevant for shielding factor calculation.
#' @param raster_elev The raster (an object of type "SpatRaster", created with [TopoShieldCalc::dem_relative]) showing the elevation angle with respect to the sampling site.
#' @param point_x Numeric. \emph{x}-coordinate of the sampling site. Can be obtained with [TopoShieldCalc::point_xyz].
#' @param point_y Numeric. \emph{y}-coordinate of the sampling site. Can be obtained with [TopoShieldCalc::point_xyz].
#' @param plot Logical. Should pairs of azimuth and elevation angles be plotted? Default to FALSE.
#' @details
#' This function creates points along polylines that extend from the sampling site.
#' The length of the polylines is equal to radius.
#' Polylines are created for each azimuth (1:360)
#' The spacing of the points corresponds to the x--y resolution of raster_elev.
#' The points are then converted into an object of type "SpatVect" with the aid of [terra::vect].
#' The elevation at the points is extracted from raster_elev with the aid of [terra::extract] (method = "bilinear").
#' The x--y coordinates of the points with the maximum elevation at each polyline are obtained to write n ESRI shapefile to the working directory with the aid of [terra::writeVector] (filetype = "ESRI Shapefile").
#' If plot = TRUE, raster_elev and the skyline are plotted with [terra::plot].
#' @examples 
#' TEXT.
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
azimuth_elevation = function(radius, # Numeric. The maximum distance from the sampling locality in metres 
                             raster_elev, # SpatRaster. The elevation with respect to the sampling locality
                             point_x, # Numeric. The x-coordinate of the sampling site
                             point_y, # Numeric. The x-coordinate of the sampling site
                             plot = FALSE){ # Logical. Enable plotting. Default to FALSE.
  resolution_dem = terra::res(raster_elev)[1] # Determine the xy-resolution of the DEM
  number_points = floor(radius / resolution_dem) # Determine the number of points (distance between the points: resolution_dem)
  radius_2 = number_points * resolution_dem # Maximum distance from the sampling site (should be roughly equal to radius) 
  azimuth = 1 # Start with azimuth = 1
  for (i in 1:360){
    if (azimuth >= 90){ # Convert the azimuth
      azimuth_2 = 360 - abs(azimuth - 90)
    } else {
      azimuth_2 = 90 - azimuth
    }
    pointmaxdist_x = radius_2 * cos(azimuth_2 * pi / 180) + point_x # Determine the x-coordinate of the point situated at the maximum distance apart from the observer point 
    pointmaxdist_y = radius_2 * sin(azimuth_2 * pi / 180) + point_y # Determine the y-coordinate of the point situated at the maximum distance apart from the observer point
    coord_point_x = seq.int(from = point_x, to = pointmaxdist_x, length.out = number_points) # Create x-coordinates of points between the observer point and the point at the maximum distance around the observer point)
    coord_point_y = seq.int(from = point_y, to = pointmaxdist_y, length.out = number_points) # Create y-coordinates of points between the observer point and the point at the maximum distance around the observer point)
    coord = data.frame(coord_point_x, coord_point_y) # Put the xy-coordinates of the points in a dataframe
    extract = terra::extract(raster_elev,
                             coord,
                             method = "bilinear",
                             na.rm = FALSE,
                             xy = TRUE,
                             ID = FALSE)
    max_elevation = max(extract[1]) # Extract the maximum elevation at the points
    index = which(extract$lyr.1 == max_elevation) # Retrieve the index number of the point with the maximum elevation angle
    if (length(index) > 1){ # If multiple cells contain max_elevation (unlikely but theoretically possible)
      index = index[1]
    } else {}
    max_elevation_x = extract$x[as.numeric(index)] # Determine the x-coordinate of the cell with the maximum elevation angle
    max_elevation_y = extract$y[as.numeric(index)] # Determine the y-coordinate of the cell with the maximum elevation angle
    if (i == 1){ # During the first iteration...
      elevation = max_elevation
      skyline_x = extract$x[as.numeric(index)] # Determine the x-coordinate of the cell with the maximum elevation angle
      skyline_y = extract$y[as.numeric(index)] # Determine the y-coordinate of the cell with the maximum elevation angle
      azimuth = azimuth + 1 # Move to the next azimuth
    } else { # During all subsequent iterations
      elevation = append(elevation, max_elevation)
      skyline_x = append(skyline_x, extract$x[as.numeric(index)])
      skyline_y = append(skyline_y, extract$y[as.numeric(index)])
      if (i < 360){ # If this is not the last iteration...
        azimuth = azimuth + 1 # Move to the next azimuth
      } else { # If i == 360...
        skyline_coordinates = cbind(id = 1, part = 1, skyline_x, skyline_y)
        skyline = terra::vect(skyline_coordinates, # Create a SpatVector with the skyline (polygon)
                              type = "polygons", 
                              crs = terra::crs(raster_elev))
        terra::writeVector(skyline, # Export the skyline as ESRI shapefile
                           filename = paste(point$Name,"_skyline.shp",sep=""),
                           filetype = "ESRI Shapefile",
                           overwrite = TRUE)
        if(plot == TRUE){
          terra::plot(raster_elev)
          terra::polys(skyline, col = NA, border = "red") # Visualise the results
        } else {}
      }
    }
  }
  list = list(elevation = elevation)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

# 1.9 Create a function to determine self-shielding for all azimuths #  

#' Calculate elevation angles for each azimuth (self-shielding)  
#' 
#' This function calculates for every azimuth (0--360 degrees) the corresponding elevation angle in degrees.
#' @param strike Numeric. The strike of the sampling surface (in degrees).
#' @param dip Numeric. The dip of the sampling surface (in degrees)
#' @param plot Logical. Should the results be plotted? Default to FALSE.
#' @details
#' This function is largely based on the MATLAB code of Greg Balco's topographic shielding calculator, available at \url{http://stoneage.ice-d.org/math/skyline/skyline.m} (last access: 16 October 2024).
#' It adds two vectors to the global environment: self_shielding_azimuth and self_shielding_elevation, containing the azimuth in degrees from north and the corresponding elevation angle in degrees, respectively.
#' @examples 
#' self_shield(strike = 20, dip = 40, plot = TRUE)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
self_shield = function(strike, # Numeric. Strike of the sampling surface
                       dip, # Numeric. Dip of the sampling surface
                       plot = FALSE){ # Logical. Enable plotting. Default to FALSE
  if (dip == 0){
    self_shielding_elevation = rep.int(0, 360)
    self_shielding_azimuth = 1:360
  } else {
    azimuth_radians = seq(from = pi / 180,
                          to = 2 * pi,
                          by = pi / 180)
    dip_radians = (dip / 360) * (2 * pi) # Convert dip to radians
    strike_radians = (strike / 360) * (2 * pi) # Convert strike to radians
    a = azimuth_radians - (strike_radians - (pi / 2)) # updip direction = strike direction - pi/2 by convention
    elevation_radians = atan(tan(dip_radians) * cos(a))
    elevation_radians = replace(elevation_radians, # Replace all values < 0 with 0
                                elevation_radians < 0,
                                0)
    self_shielding_elevation =  elevation_radians * (180 / pi)
    self_shielding_azimuth = 1:360
  }
  if(plot == TRUE){
    plot(self_shielding_azimuth,
         self_shielding_elevation,
         type = "l",
         col = "red",
         xlab = "Azimuth (degrees from North)",
         ylab = "Elevation (degrees)")
  } else {}
  list = list(self_shielding_azimuth = self_shielding_azimuth,
              self_shielding_elevation = self_shielding_elevation)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

# 1.10 Create a function to calculate topographic shielding factors #

#' Calculate the topographic shielding factor
#' 
#' This function calculates the unitless topographic shielding factor for a sampling site.
#' @param self_shielding_elevation A vector with elevation angles for each azimuth (1:360). Can be generated with [TopoShieldCalc::self_shield]
#' @param elevation 
#' @details 
#' This functions compares the elevation angles obtained with [TopoShieldCalc::self_shield] (dipping surface) and [TopoShieldCalc::azimuth_elevation] (horizon).
#' Following Greg Balco's topographic shielding calculator (MATLAB code available at \url{http://stoneage.ice-d.org/math/skyline/skyline.m}, last access: 16 October 2024), the larger elevation angle is chosen for shielding factor calculation.
#' The unitless shielding factor is calculated in the same way as in Greg Balco's topographic shielding calculator.
#' @examples 
#' shielding_factor(self_shielding_elevation = self_shielding_elevation, elevation  = elevation)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
shielding_factor = function(self_shielding_elevation,
                            elevation){
  elevation_final = pmax(elevation, # Select the larger elevation value for each azimuth
                         self_shielding_elevation)
  pdf(file = paste(point$Name, ".pdf", sep = ""),
      paper = "a4",
      pointsize = 6)
  plot(1:360, # Plot self_shielding_elevation
       self_shielding_elevation,
       type = "l",
       col = "red",
       xlab = "Azimuth (degrees from N)",
       ylab = "Elevation (degrees)",
       ylim = c(0,90))
  lines(1:360, # Plot the elevation derived from the DEM
        elevation,
        col = "blue")
  lines(1:360,
        elevation_final,
        lty = 2,
        lwd = 2,
        col = "black")
  legend("topright",
         legend = c("Elevation (self-shielding)",
                    "Elevation (horizon)",
                    "Elevation (final)"),
         lty = c(1,1,2),
         lwd = c(1,1,2),
         col = c("red","blue","black"))
  dev.off()
  azimuth = (1 / 360) * (2 * pi) # One in radians 
  elevation_final_radians = elevation_final * pi / 180
  shielding = (azimuth / (2 * pi)) * (sin(elevation_final_radians)^3.3) # Integration formula gives the shielding in %
  shielding_fact = 1 - sum(shielding)
  list = list(shielding_fact = shielding_fact)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

# 1.11 Create a function that combines all previous functions #

#' Calculate topographic shielding factors
#' 
#' This function allows for calculating topographic shielding factors for one or multiple sampling sites.
#' @param radius Numeric. The radius around the sampling site (in metres) relevant for shielding factor calculation. 
#' @details
#' This function requires a shapefile (.shp) with the sampling sites and a digital elevation model (DEM; .tif)
#' The following columns need to be added to the attribute table and populated with data: "Name", "Strike_deg", "Dip_deg", and "BouldHt", i.e., the name of the sampling site, the strike of the sampling surface (degrees), the dip of the sampling surface (degrees), and the boulder height (metres), respectively.
#' The function interactively loads the required geodata and performs quality checks. See [TopoShieldCalc::load_geodata] for further details.
#' Subsequently, the digital elevation model (DEM) is cropped.
#' See [TopoShieldCalc::crop_dem] for further details.
#' The \emph{x}-, \emph{y}-, and \emph{z}-coordinates of the sampling sites are then retrieved.
#' See [TopoShieldCalc::point_xyz] for further details.
#' The function then creates a raster with the elevation in metres above sea-level with respect to the \emph{z}-coordinates of the sampling sites. 
#' See [TopoShieldCalc::dem_relative] for further details.
#' A raster showing the distance between the sampling site and the cells of the input-DEM is subsequently generated.
#' See the documentation of [TopoShieldCalc::raster_dist] for further details.
#' For the next step, a raster showing the elevation angle (degrees) with respect to the sampling site is created.
#' See the documentation of [TopoShieldCalc::raster_elevation] for further details.
#' Elevation angles for each azimuth (1:360) are then calculated. 
#' See the documentation of [TopoShieldCalc::azimuth_elevation] for further details.
#' Elevation angles for each azimuth (self-shielding) are then calculated.
#' See the documentation of [TopoShieldCalc::self_shield] for further details.
#' Finally, thhe unitless topographic shielding factor for the sampling sites is determined.
#' See the documentation of [TopoShieldCalc::shielding_factor] for further details.
#' This function adds plots (azimuth versus elevation) to the working directory (.pdf files).
#' It also adds the skylines (ESRI shapefiles, i.e., .shp and associated files) to the working directory via [terra::writeVector].
#' The function writes an open document spreadsheet (.ods) entitled "shielding.ods" to the working directory via [readODS::write_ods].
#' This spreadsheet contains to columns entitled "Sampling site" and "Topographic shielding factor". 
#' @examples 
#' TopoShieldFact(radius = 10000)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
TopoShieldFact = function(radius){ # Numeric. Radius around the point.
  load_geodata(radius = radius) # Load the DEM and the shapefile
  for (i in 1:length(points$Name)){ # Length: number of points
    point = points[i,] # Subset the point
    crop_dem(dem = dem,
             point = point,
             radius = radius) # Crop the DEM
    point_xyz(dem = dem_2,
              point = point,
              boulder_height = point$BouldHt) # Determine x-, y-, and z-coordinates
    dem_relative(dem = dem_2,
                 point_z = point_z) # Create a DEM relative to point_z
    raster_dist(dem = dem_2, # Create a raster showing the distance between the sampling site and the midpoints of the cells of the DEM
                point_x = point_x,
                point_y = point_y)
    raster_elevation(dem_rel = dem_rel, # Compute the elevation with respect to the sampling site for all cells in the DEM
                     raster_distance = raster_distance)
    azimuth_elevation(radius = radius, # Calculate the maximum elevation for each azimuth
                      raster_elev = raster_elev,
                      point_x = point_x,
                      point_y = point_y)
    self_shield(strike = point$Strike_deg,
                dip = point$Dip_deg)
    shielding_factor(self_shielding_elevation = self_shielding_elevation,
                     elevation = elevation)
    if(i == 1){ # During the first iteration
      ShieldFact = shielding_fact
      Name = point$Name
      message(paste("Shielding factor for site #",point$Name," calculated", sep = ""))
      if (length(points$Name) == 1){ # If the shapefile only contains one point
        ShieldFact = shielding_fact
        Name = point$Name
        ShieldFact_all = data.frame(Name,ShieldFact)
        colnames(ShieldFact_all) = c("Sampling site","Topographic shielding factor")
        readODS::write_ods(ShieldFact_all,
                           path = "Shielding.ods")
        message(paste("Shielding factors for ",i," sampling sites successfully calculated", sep = ""))
      } else {}
    } else { # During the subsequent iterations
      if(i < length(points$Name)){
        ShieldFact = append(ShieldFact,shielding_fact)
        Name = append(Name,point$Name)
        message(paste("Shielding factor for site #",point$Name," calculated", sep = ""))
      } else { # During the last iteration
        ShieldFact = append(ShieldFact,shielding_fact)
        Name = append(Name,point$Name)
        ShieldFact_all = data.frame(Name,ShieldFact)
        colnames(ShieldFact_all) = c("Sampling site","Topographic shielding factor")
        readODS::write_ods(ShieldFact_all,
                           path = "shielding.ods")
        message(paste("Shielding factor for site #",point$Name," calculated", sep = ""))
        message(paste("Shielding factors for ",i," sampling site(s) successfully calculated", sep = ""))
      }
    }
  }
}
