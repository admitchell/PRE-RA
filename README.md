PRE-RA: Probabilistic Runout Estimator – Rock Avalanche 
Copyright (C) 2019  Andrew Mitchell

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.

Developer contact: Andrew Mitchell <amitchell@eoas.ubc.ca>

OVERVIEW:

Definitions of all terms, a description of the data used in the statistical analysis, and details of the analyses implemented in this code are described in Mitchell et al. (2019). This code was implemented using the R programming language.

The user has two options for analysis, path or profile. A path analysis calculates the probability of runout exceedance probability ranges along the entire digitized path for a given volume, and a point analysis provides an estimate the volume that corresponds to the user-specified target probability of runout exceedance. The code uses a user-input topography grid to create a hillshade that the user digitizes a potential rock avalanche runout path on. For the path analysis the program will produce an image with probability ranges for the runout distance and mean path width on the supplied topography. With the point analysis the program supplies an image with the digitized path, and a numerical output of the calculated volume corresponding the the target probability of runout exceedance.

Before running the program R packages 'sp' and 'raster' are required for the code.

USER INPUTS:

Line 25: Working directory where the rock avalanche case history data and the topography for the analysis are saved.

Line 29: Confinement state, TRUE if the analysis follows a laterally confined path, FALSE if the path is frontally confined or unconfined.

Line 32: Analysis type, "path" to analyse the probability of runout exceedance probability ranges along the entire digitized path for a given volume, "point" to estimate the volume that corresponds to the user-specified target probability of exceedance.

Line 35: Volume to be used in the path analysis, with units of millions meters cubed.

Line 38: Target probability of runout exceedance used in the point analysis.

Line 46: Topography data to be used for the analysis, in an ascii grid format. The analysis is set up to use UTM coordinates. The Coordinate Reference System (CRS) information for the grid must be supplied, see https://www.rdocumentation.org/packages/sp/versions/1.3-1/topics/CRS-class for more details on specifying the CRS to be used.

OPTIONS:
Optional segments are currently commented out in the code.

Lines 109 - 112: Option to read in a csv file with X and Y (easting and northing) path coordinates. This will overwrite any digitized points.

Lines 558 - 559: Option to output the X, Y and Z path coordinates from the analysis.

REFERENCE:
Mitchell, A., McDougall, S., Nolde, N., Brideau, M.-A., Whittall, J., Aaron, J. (2019) Landslides. https://doi.org/10.1007/s10346-019-01331-3
