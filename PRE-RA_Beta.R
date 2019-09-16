# Probabilistic Runout Estimator - Rock Avalanche: PRE-RA
# Version 0.1 - Beta 
# Copyright (C) 2019  Andrew Mitchell
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# For a copy of the GNU General Public License see <https://www.gnu.org/licenses/>.
# 
# Developer contact: Andrew Mitchell <amitchell@eoas.ubc.ca>
#
######################################################################################################################

# Clear the workspace
rm(list = ls())

# set the wd and read in the table for complete cases
setwd("C:/Users/amitchell/Documents/Rock Avalanche Paper/Predictive")
all.data = read.csv("Can_data.csv")

# Select confined (TRUE) or unconfined (FALSE)
confinement = TRUE 

# select the analysis to be run, "path" or "point"
analysis = "path"

# if path is selected, define the event volume in M m^3
volume = 3

# if point is selected, define the target probability of exceedance
PE = 0.25

# Read in a topo file and draw a profile
# Libraries 'sp' and 'raster' must be installed
library(sp)
library(raster)

# create a spatialGridDataFrame from an ASCII grid of the data points
points = read.asciigrid("Joffre_clip.asc", proj4string = "+proj=utm +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

# create a raster from the points
DEM = raster(points)

# create a slope and aspect raster (must be in radians to make hillshade)
slope.data = terrain(DEM, opt = "slope", unit = "radians")
aspect.data = terrain(DEM, opt = "Aspect", unit = "radians")

# create a hillshade - use for plotting
hs = hillShade(slope.data, aspect.data, angle=45, direction=0, normalize=TRUE)

# store the min and max northing and easting extents for digitizing a path
northing.min = DEM@extent@ymin
northing.max = DEM@extent@ymax
easting.min = DEM@extent@xmin
easting.max = DEM@extent@xmax

# plot the hillshade and overlay terrain colours
#start a window to display the figure
SF=0.9
windows(width = SF*11, height = SF*8.5)
par(mar=c(4, 4, 0, 1))
par(oma = c(0, 0, 0, 2))
plot(hs, col = gray(seq(1, 0, -0.1)), alpha = 0.5, legend = FALSE)
mtext("Elevation (m)", side = 4, line = 4.5, outer = FALSE, adj = 0.5)
plot(DEM, col=terrain.colors(35), alpha = 0.5, add = TRUE)
contour(DEM, add = TRUE, col = grey(0.5))
title(main = "", xlab = "Easting (m)", ylab = "Northing (m)")

# turn off the locator bell
options(locatorBell = FALSE)

### next lines must be run all together ###
# digitize the runout path
legend("topright", legend = "Digitize runout path from high to low elevation", bg = "white")

TopData = locator(n = 100, type = "o")
# break with esc button

# digitize the max and min y axis (northing) values 
legend("topright", legend = "Digitize min then max northing value...........", bg = "white")

Yrange = c(northing.min, northing.max) # small to large, order of clicks matters
Ydig = locator(n = 2, type = "o")

# digitize the max/min x axis (easting) values
legend("topright", legend = "Digitize min then max easting value............", bg = "white")

Xrange = c(easting.min, easting.max) 
Xdig = locator(n = 2, type = "o")

#then fit a linear model to transform data from arbitary to scaled values
Ymod = lm(Yrange ~ Ydig$y)
Xmod = lm(Xrange ~ Xdig$x)

#and translate the digitized values into topographic data
Y = (Ymod$coef[1] + TopData$y * Ymod$coef[2]) # intercept + interpolated value of y*slope
X = (Xmod$coef[1] + TopData$x * Xmod$coef[2])

### end path digitization ###

# alternate to read in pre-digiitzed points
# digitized = read.csv("May16_path-analysis.csv", header = TRUE)
# 
# X = digitized$X
# Y = digitized$Y

### get H and L values to use in analysis ###
# create a dataframe of path coordinates
path.coords = data.frame(X, Y)

# extract the elevation profile from the digitized path
path.z = extract(DEM, path.coords, "simple")

# initialize vectors to create/store profile
L = numeric(length(X))
delta.x = numeric(length(L))
delta.y = numeric(length(L))

# use 3D coordinates to generate 2D path
for(i in 2:length(L)) {
  delta.x[i] = path.coords[i - 1, 1] - path.coords[i, 1]
  delta.y[i] = path.coords[i - 1, 2] - path.coords[i, 2]
  
  delta.L = sqrt(delta.x[i]^2 + delta.y[i]^2)
  
  L[i] = L[i - 1] + delta.L
}

# calculate H at all points
H = max(path.z) - path.z
# set H and L at first value to 1 to allow for interpolation 
H[1] = 1
L[1] = 1

# generate linear regression for runout length and width
logV = log10(all.data$Volume)
logH = log10(all.data$H)
logL = log10(all.data$L)
con = all.data$Confinement == "Lateral"

loglm.L = lm(logL ~ logV + logH + con)

# pull out the coefficients from the regression
coeff.L = coefficients(loglm.L)

# calculate the standard deviation of the residual distribution
stdev.L = sd(loglm.L$residuals)

# calculate the mean width for the runout path for all cases in the dataset
W = all.data$Area/all.data$L
logW = log10(W)

loglm.W = lm(logW ~ logV)

# pull out the coefficients from the regression
coeff.W = coefficients(loglm.W)

# calculate the standard deviation of the residual distribution
stdev.W = sd(loglm.W$residuals)

# function to find the probability of exceedance for a given volume, fall height and distance from source
# supply regression coefficients and standard deviations
# consider the effect of lateral confinement (TRUE) and frontal or unconfined (FALSE)
# Currently only implemented for Canadian dataset
####################################################################################################

runout.probabilities = function(Vol, H, L, coeff, stdev, conf) {
  
  # calculate the probability of exceedance
  z.score = (log10(L) - coeff[1] - coeff[2]*log10(Vol) - coeff[3]*log10(H) - coeff[4]*conf)/stdev 
  
  p.exceed = 1 - pnorm(z.score)
  
  return(p.exceed)  
}
###################################################################################################

# function to find the best-fit W value for a given volume and probability of exceedance
# Currently only implemented for Canadian dataset
####################################################################################################

width.fit = function(Vol, coeff, stdev, PE) {
  
  # starting guess for width
  w = 100
  
  # calculate the probability of exceedance for the initial guess
  z.score = (log10(w) - coeff[1] - coeff[2]*log10(Vol))/stdev 
  
  p.exceed = 1 - pnorm(z.score)
  
  # calculate the initial error, and define the upper and lower bounds for width
  error = PE - p.exceed
  w.upper = 10000
  w.lower = 1
  
  # iterate until the target PE is reached
  while(abs(error) > 0.001){
    if(error < 0){
      w.new = (w + w.upper)/2
      w.lower = w
    } else {
      w.new = (w + w.lower)/2
      w.upper = w
    }
    
    z.score = (log10(w.new) - coeff[1] - coeff[2]*log10(Vol))/stdev 
    p.exceed = 1 - pnorm(z.score)
    
    # calculate the initial error, and define the upper and lower bounds for width
    error = PE - p.exceed
    w = w.new
  }
  
  return(w)  
}
###################################################################################################

# function to do point analysis, find the volume that corresponds to the target probability of exceedance
# Currently only implemented for Canadian dataset
####################################################################################################

volume.probabilities = function(PE, H, L, coeff, stdev, conf) {
  
  # starting guess for volume
  v = 10
  
  # calculate the probability of exceedance for the initial guess
  z.score = (log10(L) - coeff[1] - coeff[2]*log10(v) - coeff[3]*log10(H) - coeff[4]*conf)/stdev 
  
  p.exceed = 1 - pnorm(z.score)
  
  # calculate the initial error, and define the upper and lower bounds for width
  error = PE - p.exceed
  v.upper = 1000
  v.lower = 0.1
  
  # iterate until the target PE is reached
  while(abs(error) > 0.001){
    if(error < 0){
      v.new = (v + v.lower)/2
      v.upper = v
    } else {
      v.new = (v + v.upper)/2
      v.lower = v
    }
    
    z.score = (log10(L) - coeff[1] - coeff[2]*log10(v.new) - coeff[3]*log10(H) - coeff[4]*conf)/stdev 
    p.exceed = 1 - pnorm(z.score)
    
    # calculate the initial error, and define the upper and lower bounds for width
    error = PE - p.exceed
    v = v.new
  }
  
  return(v)  
}
###################################################################################################

# Path analysis 
if(analysis == "path"){
  
  legend("topright", legend = "Calculating runout exceedance probabilities....", bg = "white")
  
  # define the probability of exceedance values
  target.PE = c(0.95, 0.75, 0.5, 0.25, 0.05, 0.01)
  
  # find the probability of exceedance for each point on the profile, create a vector to store data
  point.pred = numeric(length = length(L))
  
  # fill the vector with probability estimates
  for(i in 1:length(point.pred)) {
    point.pred[i] = runout.probabilities(volume, H[i], L[i], coeff.L, stdev.L, confinement)
  }
  
  # check if the calculated probabilities span the range of target.PE
  check.range = which(target.PE > min(point.pred))
  target.L = target.PE[check.range]
  
  # find the values of W corresponding to different P(E) levels, create a vector to store data
  W.pred = numeric(length = length(target.PE))
  
  # loop through values of target.PE and calculate predicted W for all values of PE
  for(i in 1:length(target.PE)) {
    W.pred[i] = width.fit(volume, coeff.W, stdev.W, target.PE[i])
  }
  
  # interpolate between digitized points to find the location corresponding to the target values  
  intersect.coords = matrix(nrow = length(target.L), ncol = 2) 
  # create a vector to store the index for the values (used later for plotting on topography)
  index.PE = numeric(length = length(target.L))
  
  # check that the runout path is long enough to do the next calculation
  error.flag = 0
  if(min(point.pred > 0.95)) {
    error.flag = 1
  } else {
    for(i in 1:length(target.L)){
      # find the index first value of L from the profile that is greater than the predicted value of L at this PE
      upper.index = min(which(point.pred < target.L[i]))
      # lower index is one less
      lower.index = upper.index - 1
      
      # interpolate to find H and L values of the profile where it crosses the prediction curve
      L.upper = L[upper.index]
      L.lower = L[lower.index]
      H.upper = H[upper.index]
      H.lower = H[lower.index]
      
      # calculate slopes and intercepts for linear interpolations
      slope.prof = (H.upper - H.lower)/(L.upper - L.lower)
      b.prof = H.lower - slope.prof*L.lower
      
      # calculate an initial error halfway between the digitized points
      H.test = (H.upper + H.lower)/2
      L.test = (L.upper + L.lower)/2
      
      p.exceed = runout.probabilities(volume, H.test, L.test, coeff.L, stdev.L, confinement)
      
      error = target.L[i] - p.exceed
      H.max = H.upper
      H.min = H.lower
      L.max = L.upper
      L.min = L.lower
      
      # iterate until the target PE is reached
      while(abs(error) > 0.001){
        if(error < 0){
          H.new = H.max - (H.test - H.min)/2
          L.new = L.max - (L.test - L.min)/2
          H.lower = H.test
          L.lower = L.test
        } else {
          H.new = H.min + (H.test - H.min)/2
          L.new = L.min + (L.test - L.min)/2
          H.upper = H.test
          L.upper = L.test
        }
        
        H.test = (H.upper + H.lower)/2
        L.test = (L.upper + L.lower)/2
        
        p.exceed = runout.probabilities(volume, H.test, L.test, coeff.L, stdev.L, confinement)
        
        # calculate the initial error, and define the upper and lower bounds for width
        error = target.L[i] - p.exceed
      }
      
      # store the result
      intersect.coords[i, ] = c(L.test, H.test)
      index.PE[i] = lower.index
    }
    
    # Convert profile values back to the northing easting coordinates
    # interpolate to find the point on the profile between them
    X.PE = numeric(length = length(target.L))
    Y.PE = numeric(length = length(target.L))
    
    for(i in 1:length(target.L)) {
      X.PE[i] = X[index.PE[i]] + (X[index.PE[i]+1] - X[index.PE[i]])*((L[index.PE[i]] - intersect.coords[i, 1])/(L[index.PE[i]] - L[index.PE[i]+1]))
      Y.PE[i] = Y[index.PE[i]] + (Y[index.PE[i]+1] - Y[index.PE[i]])*((L[index.PE[i]] - intersect.coords[i, 1])/(L[index.PE[i]] - L[index.PE[i]+1]))
    }
  }
  
  # make points offset from the path line segments to plot average width
  X.mid = numeric(length = length(X) - 1)
  Y.mid = numeric(length = length(Y) - 1)
  m = numeric(length = length(X.mid))
   
  for(i in 1:length(X.mid)){
    X.mid[i] = (X[i + 1] - X[i])/2 + X[i]
    Y.mid[i] = (Y[i + 1] - Y[i])/2 + Y[i]

    m[i] = (Y[i + 1] - Y[i])/(X[i + 1] - X[i])
  }
  
  # make a vertical plane containing each line segment 
  # calculate the normal vector to that plane for each line segment (used for plotting width offsets)
  norm.vec = matrix(nrow = length(X) - 1, ncol = 2)
  
  for(i in 2:length(X)){
    delta.z = H[i] - H[i - 1]
    norm.vec[i - 1, 1] = delta.z * delta.y[i]
    norm.vec[i - 1, 2] = -1 * delta.z * delta.x[i]
  }
  
  # convert to unit normal vectors
  unit.norm = matrix(nrow = length(X) - 1, ncol = 2)
  
  for(i in 2:length(X)){
    length.n = sqrt((norm.vec[i - 1, 1]^2) + (norm.vec[i - 1, 2])^2)
    unit.norm[i - 1, 1] = norm.vec[i - 1, 1]/length.n
    unit.norm[i - 1, 2] = norm.vec[i - 1, 2]/length.n
  }
  
  # make a series of offset points for each value of W at each midpoint
  upperoffset.X = matrix(nrow = length(X.mid), ncol = length(W.pred))
  loweroffset.X = matrix(nrow = length(X.mid), ncol = length(W.pred))
  upperoffset.Y = matrix(nrow = length(X.mid), ncol = length(W.pred))
  loweroffset.Y = matrix(nrow = length(X.mid), ncol = length(W.pred))
  
  # Fill in the offset points along each line segment
  for(i in 1:length(W.pred)) { 
    for(j in 1:length(X.mid)) {
      upperoffset.X[j, i] = unit.norm[j, 1] * W.pred[i]/2 + X.mid[j]
      upperoffset.Y[j, i] = unit.norm[j, 2] * W.pred[i]/2 + Y.mid[j]
      loweroffset.X[j, i] = -1 * unit.norm[j, 1] * W.pred[i]/2 + X.mid[j]
      loweroffset.Y[j, i] = -1 * unit.norm[j, 2] * W.pred[i]/2 + Y.mid[j]
    }
  }
  
  # check that the lines drawn from the offset points don't cross
  for(i in 2:length(upperoffset.X[ , 1])) {
    # calculate the slope of the upper and lower offset lines
    m.upper = (upperoffset.Y[i, 1] - upperoffset.Y[i - 1, 1])/(upperoffset.X[i, 1] - upperoffset.X[i - 1, 1])
    m.lower = (loweroffset.Y[i, 1] - loweroffset.Y[i - 1, 1])/(loweroffset.X[i, 1] - loweroffset.X[i - 1, 1])
    
    # check for cases with very small x components (very large slope values)
    if(abs(m.upper) > 100) m.upper = 100
    if(abs(m.lower) > 100) m.lower = 100
    
    # check that the two slopes are within 15% (allowing for round off error)
    m.check = (m.upper - m.lower)/m[i]
    
    if(abs(m.check) > 0.15) {
      temp.x = upperoffset.X[i, ]
      temp.y = upperoffset.Y[i, ]
      
      upperoffset.X[i, ] = loweroffset.X[i, ]
      loweroffset.X[i, ] = temp.x
      
      upperoffset.Y[i, ] = loweroffset.Y[i, ]
      loweroffset.Y[i, ] = temp.y
    }
  }

  # replot topography for final figure
  # SF=0.9
  windows(width = SF*11, height = SF*8.5)
  par(mar=c(4, 4, 0, 1))
  par(oma = c(0, 0, 0, 2))
  plot(hs, col = gray(seq(1, 0, -0.1)), alpha = 0.5, legend = FALSE)
  mtext("Elevation (m)", side = 4, line = 4.5, outer = FALSE, adj = 0.5)
  plot(DEM, col=terrain.colors(35), alpha = 0.5, add = TRUE)
  contour(DEM, add = TRUE, col = grey(0.5))
  title(main = "", xlab = "Easting (m)", ylab = "Northing (m)")
  
  #plot coloured lines over the runout path to indicate the probabilities of exceedance

  if(error.flag == 1) {
    lines(X[1:length(X)], Y[1:length(Y)], col = "black", lwd = 3)
  } else {
    if(length(X.PE) == 6) {
      lines(c(X.PE[6],X[1 + index.PE[6]:length(X)]), c(Y.PE[6],Y[1  +index.PE[6]:length(Y)]), col = "green", lwd = 3)
      lines(c(X.PE[5],X[1 + index.PE[5]:index.PE[6] - 1], X.PE[6]), c(Y.PE[5],Y[1 + index.PE[5]:index.PE[6] - 1], Y.PE[6]), col = "blue", lwd = 3)
      lines(c(X.PE[4],X[1 + index.PE[4]:index.PE[5] - 1], X.PE[5]), c(Y.PE[4],Y[1 + index.PE[4]:index.PE[5] - 1], Y.PE[5]), col = "yellow", lwd = 3)
      lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[4]), col = "orange", lwd = 3)
      lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
      lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
      lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
    }
    if(length(X.PE) == 5) {
      lines(c(X.PE[5],X[1 + index.PE[5]:length(X)]), c(Y.PE[5],Y[1  +index.PE[5]:length(Y)]), col = "blue", lwd = 3)
      lines(c(X.PE[4],X[1 + index.PE[4]:index.PE[5] - 1], X.PE[5]), c(Y.PE[4],Y[1 + index.PE[4]:index.PE[5] - 1], Y.PE[5]), col = "yellow", lwd = 3)
      lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[4]), col = "orange", lwd = 3)
      lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
      lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
      lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
    }
    if(length(X.PE) == 4) {
      lines(c(X.PE[4],X[1 + index.PE[4]:length(X)]), c(Y.PE[4],Y[1  +index.PE[4]:length(Y)]), col = "yellow", lwd = 3)
      lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[4]), col = "orange", lwd = 3)
      lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
      lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
      lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
    } 
    if(length(X.PE) == 3) {
      lines(c(X.PE[3],X[1 + index.PE[3]:length(X)]), c(Y.PE[3],Y[1  +index.PE[3]:length(Y)]), col = "orange", lwd = 3)
      lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
      lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
      lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
    }
    if(length(X.PE) == 2) {
      lines(c(X.PE[2],X[1 + index.PE[2]:length(X)]), c(Y.PE[2],Y[1  +index.PE[2]:length(Y)]), col = "red", lwd = 3)
      lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
      lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
    }
    if(length(X.PE) == 1) {
      lines(c(X.PE[1],X[1 + index.PE[1]:length(X)]), c(Y.PE[1],Y[1  +index.PE[1]:length(Y)]), col = "darkred", lwd = 3)
      lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
    }
    points(X.PE, Y.PE, pch = 4, lwd = 3)
  }
  
  legend("bottomright", title = "Probability of Exceedance", lty = c(1, 1, 1, 1, 1, 1, 1), lwd = c(3, 3, 3, 3, 3, 3, 3),
         legend = (c("> 0.95", "0.95 - 0.75", "0.75 - 0.5", "0.5 - 0.25", "0.25 - 0.05", "0.05 - 0.01")),
         col = c("black", "darkred", "red", "orange", "yellow", "blue", "green"), bg = "white")
  
  # add in lines for the width exceedance contours
  # points(upperoffset.X[ , 2], upperoffset.Y[ , 2], pch = 18, col = "darkred")
  # points(loweroffset.X[ , 2], loweroffset.Y[ , 2], pch = 18, col = "darkred")
  # 
  # points(upperoffset.X[ , 3], upperoffset.Y[ , 3], pch = 18, col = "red")
  # points(loweroffset.X[ , 3], loweroffset.Y[ , 3], pch = 18, col = "red")
  # 
  # points(upperoffset.X[ , 4], upperoffset.Y[ , 4], pch = 18, col = "yellow")
  # points(loweroffset.X[ , 4], loweroffset.Y[ , 4], pch = 18, col = "yellow")

  lines(upperoffset.X[ , 2], upperoffset.Y[ , 2], lty = 2, col = "darkred")
  lines(loweroffset.X[ , 2], loweroffset.Y[ , 2], lty = 2, col = "darkred")

  lines(upperoffset.X[ , 3], upperoffset.Y[ , 3], lty = 2, col = "red")
  lines(loweroffset.X[ , 3], loweroffset.Y[ , 3], lty = 2, col = "red")

  lines(upperoffset.X[ , 4], upperoffset.Y[ , 4], lty = 2, col = "yellow")
  lines(loweroffset.X[ , 4], loweroffset.Y[ , 4], lty = 2, col = "yellow")
}  
 
### Point analysis  ###
if(analysis == "point"){
  legend("topright", legend = "Calculating target volume..........................", bg = "white")

  # use the volume.probabilities function to calculate the minimum volume for the defined probability of exceedance
  min.volume = volume.probabilities(PE, H[length(H)], L[length(L)], coeff.L, stdev.L, confinement)

  vol.plot = as.character(signif(min.volume, 2))
  PE.plot = as.character(PE)
  output.text = paste(vol.plot, "M m^3 = volume for P(L>l) = ", PE.plot)

#   # replot topography for final figure
  SF=0.9
  windows(width = SF*11, height = SF*8.5)
  par(mar=c(4, 4, 0, 1))
  par(oma = c(0, 0, 0, 2))
  plot(hs, col = gray(seq(1, 0, -0.1)), alpha = 0.5, legend = FALSE)
  mtext("Elevation (m)", side = 4, line = 4.5, outer = FALSE, adj = 0.5)
  plot(DEM, col=terrain.colors(35), alpha = 0.5, add = TRUE)
  contour(DEM, add = TRUE, col = grey(0.5))
  title(main = "", xlab = "Easting (m)", ylab = "Northing (m)")

  legend("topright", legend = output.text, bg = "white")
  
  lines(X, Y)
  points(X[1], Y[1], pch = 4, lwd = 3)
  text(X[1], Y[1], "Source crest", pos = 4)
  points(X[length(X)], Y[length(Y)], pch = 4, lwd = 3)
  text(X[length(X)], Y[length(Y)], "Point of interest", pos = 3)
}

# # output X, Y and Z values to re-run an analysis if needed
# output = data.frame(X, Y, path.z)
# write.csv(output, file = "xyz_coords.csv")
