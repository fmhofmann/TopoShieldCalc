# TopoShieldCalc 

<img align="right" width="250" src="man/Logo_TopoShieldCalc.jpg">

<p align="right"; style="font-size:11px">logo by Felix Martin Hofmann</p>

`TopoShieldCalc` is an *R* package for topographical shielding factor calculation for discrete sampling sites with the aid of a digital elevation model (DEM)
## Installation

`TopoShieldCalc` is available from GitHub, so you can use the following code to get the current *released version*:

1. Install `devtools` on your machine

```
install.packages("devtools")
```

2. Load the devtools package

```
library(devtools)
```

3. Install the package with the `install_github()`command
 
```
install_github("fmhofmann/TopoShieldCalc")
```
Make sure that the `terra` and `readODS`packages are installed. 

4. Put the package into the library
 
```
library(TopoShieldCalc)
```
Now you should be able to use the package. Have fun and enjoy!

## Preparation of the input-ESRI shapefile 

Note that the unit of the coordinate reference system (CRS) of the input-shapefile has to be metres. If this is not the case, the load_geodata() function will return an error message.

The attribute table of the ESRI shapefile has to be structured as follows:

Name                      | Strike_deg                       | Dip_deg                         | BouldHt                      |
--------------------------|----------------------------------|---------------------------------|------------------------------|
Name of the sampling site | The strike in degrees from north | The dip of the sampling surface | The boulder height in metres |

Example:

Name   | Strike_deg | Dip_deg | BouldHt |
-------|------------|---------|---------|
FS-01a | 210        | 5       | 1       |
