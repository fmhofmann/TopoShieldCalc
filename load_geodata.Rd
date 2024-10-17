\name{load_geodata}
\alias{load_geodata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{~~ Dummy title ~~}
\description{~~ Dummy description ~~}
\usage{
load_geodata()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function () 
{
    message("Select the shapefile")
    points = terra::vect(file.choose())
    if (terra::linearUnits(points) == 0) {
        stop("The linear unit of the coordinate reference system of the input-shapefile needs to be metres")
    }
    else {
    }
    message("Select the DEM")
    dem = terra::rast(file.choose())
    if (terra::linearUnits(dem) == 0) {
        stop("The linear unit of the input-DEM needs to be metres")
    }
    else {
        if (terra::crs(points, proj = TRUE) != terra::crs(dem, 
            proj = TRUE)) {
            stop("The input-DEM and the input-shapefile do not have the same coordinate reference system")
        }
        else {
            colnames = colnames(terra::as.data.frame(points))
            if ("Name" \%in\% colnames & "Strike_deg" \%in\% colnames & 
                "Dip_deg" \%in\% colnames & "BouldHt" \%in\% colnames) {
                message("The attribute table of the input-shapefile contains all necessary data for topographical shielding factor calculation")
            }
            else {
                stop("The input-shapefile does not contain all columns needed for shielding factor calculation: Name, Strike_deg, Dip_deg, and BouldHt")
            }
        }
    }
    message("The input-shapefile and the input-DEM have the same coordinate reference system (CRS) and the unit of the CRS is metres. They are suitable for topographic shielding factor calculation")
    list = list(dem = dem, points = points)
    return(list2env(list, envir = globalenv()))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
