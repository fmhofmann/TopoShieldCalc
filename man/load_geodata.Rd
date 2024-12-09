% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TopoShieldCalc.R
\name{load_geodata}
\alias{load_geodata}
\title{Load geodata for topographic shielding factor calculation}
\usage{
load_geodata(radius = 10000)
}
\arguments{
\item{radius}{Numeric. The radius (in metres) around the sampling site relevant for shielding factor calculation. Default to 10000 metres.}
}
\description{
This function allows for loading geodata (digital elevation model (DEM) and shapefile) for topographic shielding factor calculation.
}
\details{
The function adds two objects to the global environment: "dem" and "points", i.e., the digital elevation model (DEM) and the shapefile with the sampling site.
For this step, the \code{\link[terra:rast]{terra::rast()}} and \code{\link[terra:vect]{terra::vect()}} commands are used, respectively.
As soon as the "Select the shapefile" message appears, the shapefile (.shp) can be interactively selected.
The function then performs a quality check of the input-shapefile with \code{\link[terra:crs]{terra::crs()}}.
If the unit of the coordinate reference system (CRS) is NOT in metres, the function returns an error message and the function execution stops.
The function then checks whether the attribute table of the shapefile contains all necessary data for shielding factor calculation.
Note that the attribute table of the input file must contain the following columns: Name, Strike_deg, Dip_deg, and BouldHt.
If this is not the case, an error message is returned and the function execution stops.

The function will then prompt the user to interactively load the digital elevation model (DEM).
The function then performs a quality check of the input-DEM with \code{\link[terra:crs]{terra::crs()}}.
If the unit of the coordinate reference system (CRS) is NOT in metres, the function returns an error message and the function execution stops.
The function then assesses with \code{\link[terra:crs]{terra::crs()}} whether the input-DEM and the input-shapefile have the same CRS.
If this is not the case, an error message is returned and the function execution stops.
The function subsequently assesses with \code{\link[terra:xmin]{terra::xmin()}}, \code{\link[terra:xmin]{terra::xmax()}}, \code{\link[terra:xmin]{terra::ymin()}}, and \code{\link[terra:xmin]{terra::ymax()}} whether the DEM covers the whole area relevant for shielding factor calculation.
If this is not the case, an error message is returned and the function execution stops.
If the input-shapefile and the input-raster meet the requirements for shielding factor calculation, the function returns the following message: "The input-shapefile and the input-DEM have the same coordinate reference system (CRS) and the unit of the CRS is metres.
They are suitable for topographic shielding factor calculation".

The function adds two objects to the global environment: "dem" and "points".
These are of class "SpatRast" and "SpatVect", respectively.
For this step, the \code{\link[terra:rast]{terra::rast()}} and \code{\link[terra:vect]{terra::vect()}} commands are used, respectively.
}
\examples{
load_geodata(radius = 10000)
}
\author{
Felix Martin Hofmann, University of Freiburg, Germany (\email{fmhofmann9892@gmail.com} and Stefan Hergarten, University of Freiburg, Germany)
}
