################################################
### Compute topographic shielding with a DEM ###
################################################

## Make sure that the terra and readODS packages are installed ##

#' Load geodata for topographic shielding factor calculation
#'
#' This function allows for loading geodata (digital elevation model (DEM) and shapefile) for topographic shielding factor calculation.
#' @param radius Numeric. The radius (in metres) around the sampling site relevant for shielding factor calculation. Default to 10000 metres.  
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
load_geodata = function(radius = 10000){ 
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
                     point, # SpatVect. Sampling site
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

#' Obtain 360 pairs of azimuth and corresponding elevation angles for a sampling site.
#' 
#' This function allows for obtaining 360 pairs of azimuth and the corresponding elevation angles with the aid of a shapefile and a digital elevation model (loaded with [TopoShieldCalc::load_geodata]).
#' @param dem The DEM (object of class "SpatRaster", see [terra::SpatRaster] for further details).
#' @param point Shapefile of the sampling site (object of class "SpatVect, see [terra::SpatVector]).
#' @param point_x Numeric. The \emph{x}-coordinate of the sampling site.
#' @param point_y Numeric. The \emph{y}-coordinate of the sampling site.
#' @param point_z Numeric. The \emph{z}-coordinate of the sampling site.
#' @param radius Numeric. The radius around the sampling site relevant for topographic shielding calculation. Default to 10000 metres.
#' @param plot Logical. Create a map with the skyline? 
#' @details
#' For the first step, this function determines the \emph{xy}-resolution of the input-DEM with [terra::res].
#' This function then creates point along polylines from the sampling site to the point at the maximum distance (defined by radius) for each azimuth.
#' The function subsequently creates points along the polylines. 
#' The distance between the points corresponds to the \emph{xy}-resolution of the input-DEM.
#' The elevation at the points is then retrieved with [terra::extract] via bilinear interpolation.
#' \emph{xy}-coordinates of the points are also returned to construct the skyline around the sampling site.
#' This function then determines the elevation angle with respect to the sampling site for all points.
#' The maximum elevation angle for each azimuth is subsequently retrieved. 
#' The skyline is constructed with the \emph{xy}-coordinates of the points with the maximum elevation angles.
#' The skyline (created with the aid of [terra::vect]) is then exported to the working directory with [terra::writeVector]. 
#' A vector entitled "elevation" is added to the global environment. 
#' This vector contains the elevation angles for all azimuths (in 1:360).
#' Note that the exported shapefile is an "ESRI shapefile" which can be opened with geographic information system (GIS) software, such as QGIS.
#' @examples 
#' point_xyz(dem, point_x = 414416.5, point_y = 5316445.8, point_z = 279, radius = 10000, plot = TRUE)
#' @author Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@@gmail.com})
#' @export
azimuth_elevation_horizon = function(dem,
                                     point,
                                     point_x,
                                     point_y,
                                     point_z,
                                     radius = 10000,
                                     plot = FALSE){
  geom = terra::geom(point) # Get the x- and y-coordinate of the sampling site
  mask = terra::ext(geom[3] - radius - 200, 
                    geom[3] + radius + 200,
                    geom[4] - radius - 200,
                    geom[4] + radius + 200) # Create a mask 
  
  dem_2 = terra::crop(dem,
                      mask,
                      touches = TRUE) # Crop dem
  resolution_raster = terra::res(dem_2)[1] # Determine the xy-resolution of the DEM
  number_points = floor(radius / resolution_raster) + 1 # Determine the number of points (distance between the points: resolution_raster)
  radius_2 = number_points * resolution_raster - 1 # Maximum distance from the sampling site (should be roughly equal to radius)
  azimuth = 1 # Begin with azimuth = 1
  for (i in 1:360){
    if (azimuth >= 90){ # Convert the azimuth
      azimuth_2 = 360 - abs(azimuth - 90)
    } else {
      azimuth_2 = 90 - azimuth
    }
    pointmaxdist_x = radius_2 * cos(azimuth_2 * pi / 180) + point_x # Determine the x-coordinate of the point situated at the maximum distance apart from the observer point 
    pointmaxdist_y = radius_2 * sin(azimuth_2 * pi / 180) + point_y # Determine the y-coordinate of the point situated at the maximum distance apart from the observer point
    coord_point_x = seq.int(from = point_x, 
                            to = pointmaxdist_x, 
                            length.out = number_points) # Create x-coordinates of points between the observer point and the point at the maximum distance around the observer point)
    coord_point_y = seq.int(from = point_y, 
                            to = pointmaxdist_y, 
                            length.out = number_points) # Create y-coordinates of points between the observer point and the point at the maximum distance around the observer point)
    extract = terra::extract(dem_2, 
                             data.frame(coord_point_x, coord_point_y), 
                             method = "bilinear", 
                             xy = TRUE, 
                             ID = FALSE)
    x = extract$x
    y = extract$y
    z = extract[,1]
    slope = (z - point_z) / sqrt((x - point_x)^2 + (y - point_y)^2) # Slope is defined as: change in elevation / distance
    elevation_angle = atan(slope) * 180 / pi # Compute the elevation angle
    max_elevation = max(elevation_angle) # Extract the maximum elevation at the points
    index = which(elevation_angle == max_elevation) # Retrieve the index number of the point with the maximum elevation angle
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
                              crs = terra::crs(dem_2))
        terra::writeVector(skyline, # Export the skyline as ESRI shapefile
                           filename = paste(point$Name,"_skyline.shp",sep=""),
                           filetype = "ESRI Shapefile",
                           overwrite = TRUE)
        if(plot == TRUE){
          terra::plot(dem_2)
          terra::polys(skyline, col = NA, border = "red") # Visualise the results
        } else {}
      }
    }
  }
  elevation = replace(elevation, # If present, replace elevation < 0 by 0
                      elevation < 0,
                      0)
  list = list(elevation = elevation)
  return(list2env(list, envir = globalenv())) # Return the vectors to the global environment
}

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
    a = azimuth_radians - (strike_radians - (pi / 2)) # updip direction equals strike direction - pi/2 by convention
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

#' Calculate the topographic shielding factor
#' 
#' This function calculates the unitless topographic shielding factor for a sampling site.
#' @param point Shapefile of the sampling site (object of class "SpatVect, see [terra::SpatVector]).
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
shielding_factor = function(point, 
                            self_shielding_elevation,
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

#' Calculate topographic shielding factors
#' 
#' This function allows for calculating topographic shielding factors for one or multiple sampling sites.
#' @param radius Numeric. The radius around the sampling site (in metres) relevant for shielding factor calculation. Default to 10000 m. 
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
TopoShieldFact = function(radius = 10000){ # Numeric. Radius around the point.
  load_geodata(radius = radius) # Load the DEM and the shapefile
  for (i in 1:length(points$Name)){ # Length: number of points
    point = points[i,] # Subset the point
    point_xyz(dem = dem_2,
              point = point,
              boulder_height = point$BouldHt) # Determine x-, y-, and z-coordinates
    azimuth_elevation_horizon(dem = dem,
                              point = point,
                              point_x = point_x,
                              point_y = point_y,
                              point_z = point_z,
                              radius = radius)
    self_shield(strike = point$Strike_deg,
                dip = point$Dip_deg)
    shielding_factor(point = point,
                     self_shielding_elevation = self_shielding_elevation,
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
