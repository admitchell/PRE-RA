# Stochastic Runout Estimation - Rock Avalanche: SRE-RA
# Version 0.1 - Beta
# 
# Read in user supplied topo
# Digitize a profile on the topo
# Find the probability of exceedance for points along the path using the multilinear regression
# Plot the probability of runout exceedance intervals on the path
#
# 18 January 2019 - AM multiple linear regression implemented
# 22 January 2019 - AM addition of script to predict width of flows
# 06 February 2019 - AM update for indicator variable and probability calculation
######################################################################################################################

# Clear the workspace and set library path
rm(list = ls())
.libPaths("D:/R Library/win-library/3.5")

# set the wd and read in the table for complete cases
setwd("D:/Rock Avalanche Database/Predictive")
all.data = read.csv("Can_data.csv")

# select the dataset to use in the statical analysis
# dataset = "Canadian"

# Select confined (TRUE) or unconfined (FALSE)
confinement = TRUE 

# select the analysis to be run, "path" or "point"
analysis = "path"

# if path was selected, define the event volume in Mm^3
volume = 3

# if point was selected, define the target probability of exceedance
PE = 0.25

# Read in a topo file and draw a profile
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
title(main = "", xlab = "Easting (m)", ylab = "Northing (m)")

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
digitized = read.csv("May16_path-analysis.csv", header = TRUE)

X = digitized$X
Y = digitized$Y

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
  
  ### NEW CODE FOR MULTILINEAR REGRESSION ###
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
  target.PE = target.PE[check.range]
  
  # find the values of W corresponding to different P(E) levels, create a vector to store data
  W.pred = numeric(length = length(target.PE))
  
  # loop through values of target.PE and calculate predicted W for all values of PE
  for(i in 1:length(target.PE)) {
    W.pred[i] = width.fit(volume, coeff.W, stdev.W, target.PE[i])
  }
  
  # interpolate between digitized points to find the location corresponding to the target values  
  intersect.coords = matrix(nrow = length(target.PE), ncol = 2) 
  # create a vector to store the index for the values (used later for plotting on topography)
  index.PE = numeric(length = length(target.PE))
  
  for(i in 1:length(target.PE)){
    # find the index first value of L from the profile that is greater than the predicted value of L at this PE
    upper.index = min(which(point.pred < target.PE[i]))
    # lower index is one less
    lower.index = upper.index - 1
    
    # interpolate to find H and L values of the profile where it crosses the prediction curve
    # pred.upper = point.pred[upper.index, i]
    # pred.lower = point.pred[lower.index, i]
    L.upper = L[upper.index]
    L.lower = L[lower.index]
    H.upper = H[upper.index]
    H.lower = H[lower.index]
    
    # calculate slopes and intercepts for linear interpolations
    # slope.pred = (H.upper - H.lower)/(pred.upper - pred.lower)
    slope.prof = (H.upper - H.lower)/(L.upper - L.lower)
    # b.pred = H.lower - slope.pred*pred.lower
    b.prof = H.lower - slope.prof*L.lower
    
    # calculate an initial error halfway between the digitized points
    H.test = (H.upper + H.lower)/2
    L.test = (L.upper + L.lower)/2
    
    p.exceed = runout.probabilities(volume, H.test, L.test, coeff.L, stdev.L, confinement)
    
    error = target.PE[i] - p.exceed
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
      error = target.PE[i] - p.exceed
    }
    
    # store the result
    intersect.coords[i, ] = c(L.test, H.test)
    index.PE[i] = lower.index
  }
 
  # Convert profile values back to the northing easting coordinates
  # interpolate to find the point on the profile between them
  X.PE = numeric(length = length(target.PE))
  Y.PE = numeric(length = length(target.PE))

  for(i in 1:length(target.PE)) {
    X.PE[i] = X[index.PE[i]] + (X[index.PE[i]+1] - X[index.PE[i]])*((L[index.PE[i]] - intersect.coords[i, 1])/(L[index.PE[i]] - L[index.PE[i]+1]))
    Y.PE[i] = Y[index.PE[i]] + (Y[index.PE[i]+1] - Y[index.PE[i]])*((L[index.PE[i]] - intersect.coords[i, 1])/(L[index.PE[i]] - L[index.PE[i]+1]))
  }
  
  # make points offset from the path line segments to plot average width
  X.mid = numeric(length = length(X) - 1)
  Y.mid = numeric(length = length(Y) - 1)
  # m = numeric(length = length(X.mid))
  # 
  for(i in 1:length(X.mid)){
    X.mid[i] = (X[i + 1] - X[i])/2 + X[i]
    Y.mid[i] = (Y[i + 1] - Y[i])/2 + Y[i]

  #  m[i] = (Y[i + 1] - Y[i])/(X[i + 1] - X[i])
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
  

  # replot topography for final figure
  SF=0.9
  windows(width = SF*11, height = SF*8.5)
  par(mar=c(4, 4, 0, 1))
  par(oma = c(0, 0, 0, 2))
  plot(hs, col = gray(seq(1, 0, -0.1)), alpha = 0.5, legend = FALSE)
  mtext("Elevation (m)", side = 4, line = 4.5, outer = FALSE, adj = 0.5)
  plot(DEM, col=terrain.colors(35), alpha = 0.5, add = TRUE)
  contour(DEM, add = TRUE, col = grey(0.5))
  title(main = "", xlab = "Easting (m)", ylab = "Northing (m)")

  #plot coloured lines over the runout path to indicate the probabilities of exceedance
  # green line for probability of exceedance less than 1%, blue line for probability of exceedance less than 5%, greater than 1%, check for NA's
  if(is.na(X.PE[6]) & !is.na(X.PE[5])) {
    lines(c(X.PE[5],X[1 + index.PE[5]:length(X)]), c(Y.PE[5],Y[1  +index.PE[5]:length(Y)]), col = "blue", lwd = 3)
    lines(c(X.PE[4],X[1 + index.PE[4]:index.PE[5] - 1], X.PE[5]), c(Y.PE[4],Y[1 + index.PE[4]:index.PE[5] - 1], Y.PE[5]), col = "yellow", lwd = 3)
    lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[4]), col = "orange", lwd = 3)
    lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
    lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
    lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
  }
  if(is.na(X.PE[5]) & !is.na(X.PE[4])) {
    lines(c(X.PE[4],X[1 + index.PE[4]:length(X)]), c(Y.PE[4],Y[1  +index.PE[4]:length(Y)]), col = "yellow", lwd = 3)
    lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[4]), col = "orange", lwd = 3)
    lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
    lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
    lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
  } else {
    lines(c(X.PE[6],X[1 + index.PE[6]:length(X)]), c(Y.PE[6],Y[1  +index.PE[6]:length(Y)]), col = "green", lwd = 3)
    lines(c(X.PE[5],X[1 + index.PE[5]:index.PE[6] - 1], X.PE[6]), c(Y.PE[5],Y[1 + index.PE[5]:index.PE[6] - 1], Y.PE[6]), col = "blue", lwd = 3)
    lines(c(X.PE[4],X[1 + index.PE[4]:index.PE[5] - 1], X.PE[5]), c(Y.PE[4],Y[1 + index.PE[4]:index.PE[5] - 1], Y.PE[5]), col = "yellow", lwd = 3)
    lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[4]), col = "orange", lwd = 3)
    lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
    lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
    lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)
  }
  # # yellow line for probability of exceedance less than 34 %, greater than 5, NA check
  # if(is.na(X.PE[5])) {
  #   lines(c(X.PE[4],X[1 + index.PE[4]:length(X)]), c(Y.PE[4],Y[1  +index.PE[4]:length(Y)]), col = "yellow", lwd = 3)
  # } else {
  #   lines(c(X.PE[5],X[1 + index.PE[5]:length(X)]), c(Y.PE[5],Y[1  +index.PE[5]:length(Y)]), col = "blue", lwd = 3)
  #   lines(c(X.PE[4],X[1 + index.PE[4]:index.PE[5] - 1], X.PE[5]), c(Y.PE[4],Y[1 + index.PE[4]:index.PE[5] - 1], Y.PE[5]), col = "yellow", lwd = 3)
  # }
  # # orange line for probability of exceedance less than 50 %, greater than 34%, NA check
  # if(is.na(X.PE[4])) {
  #   lines(c(X.PE[3],X[1 + index.PE[3]:length(X)]), c(Y.PE[3],Y[1  +index.PE[3]:length(Y)]), col = "yellow", lwd = 3)
  # } else {
  #   lines(c(X.PE[4],X[1 + index.PE[4]:length(X)]), c(Y.PE[4],Y[1  +index.PE[4]:length(Y)]), col = "blue", lwd = 3)
  #   lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[3]), col = "yellow", lwd = 3)
  # }
  # 
  # 
  # 
  # # yellow line for probability of exceedance less than 34 %, greater than 5%
  # lines(c(X.PE[4],X[1 + index.PE[4]:index.PE[5] - 1], X.PE[5]), c(Y.PE[4],Y[1 + index.PE[4]:index.PE[5] - 1], Y.PE[5]), col = "yellow", lwd = 3)
  # # orange line for probability of exceedance less than 50 %, greater than 34%
  # lines(c(X.PE[3],X[1 + index.PE[3]:index.PE[4] - 1], X.PE[4]), c(Y.PE[3],Y[1 + index.PE[3]:index.PE[4] - 1], Y.PE[4]), col = "orange", lwd = 3)
  # # red line for probability of exceedance less than 68 %, greater than 50%
  # lines(c(X.PE[2],X[1 + index.PE[2]:index.PE[3] - 1], X.PE[3]), c(Y.PE[2],Y[1 + index.PE[2]:index.PE[3] - 1], Y.PE[3]), col = "red", lwd = 3)
  # # darkred line for probability of exceedance less than 95 %, greater than 68%
  # lines(c(X.PE[1],X[1 + index.PE[1]:index.PE[2] - 1], X.PE[2]), c(Y.PE[1],Y[1 + index.PE[1]:index.PE[2] - 1], Y.PE[2]), col = "darkred", lwd = 3)
  # # black line for probability of exceedance greater than 95%
  # lines(c(X[1:index.PE[1]], X.PE[1]), c(Y[1:index.PE[1]], Y.PE[1]), col = "black", lwd = 3)

  legend("bottomright", title = "Probability of Exceedance", lty = c(1, 1, 1, 1, 1, 1, 1), lwd = c(3, 3, 3, 3, 3, 3, 3),
         legend = (c("> 0.95", "0.95 - 0.75", "0.75 - 0.5", "0.5 - 0.25", "0.25 - 0.05", "0.05 - 0.01")),
         col = c("black", "darkred", "red", "orange", "yellow", "blue", "green"), bg = "white")

  points(X.PE, Y.PE, pch = 4, lwd = 3)
  
  # add in lines for the width exceedance contours
  points(upperoffset.X[ , 2], upperoffset.Y[ , 2], pch = 18, col = "darkred")
  points(loweroffset.X[ , 2], loweroffset.Y[ , 2], pch = 18, col = "darkred")

  points(upperoffset.X[ , 3], upperoffset.Y[ , 3], pch = 18, col = "red")
  points(loweroffset.X[ , 3], loweroffset.Y[ , 3], pch = 18, col = "red")

  points(upperoffset.X[ , 4], upperoffset.Y[ , 4], pch = 18, col = "yellow")
  points(loweroffset.X[ , 4], loweroffset.Y[ , 4], pch = 18, col = "yellow")

  # lines(upperoffset.X[ , 2], upperoffset.Y[ , 2], lty = 2, col = "darkred")
  # lines(loweroffset.X[ , 2], loweroffset.Y[ , 2], lty = 2, col = "darkred")
  # 
  # lines(upperoffset.X[ , 3], upperoffset.Y[ , 3], lty = 2, col = "red")
  # lines(loweroffset.X[ , 3], loweroffset.Y[ , 3], lty = 2, col = "red")
  # 
  # lines(upperoffset.X[ , 4], upperoffset.Y[ , 4], lty = 2, col = "yellow")
  # lines(loweroffset.X[ , 4], loweroffset.Y[ , 4], lty = 2, col = "yellow")
}  
 
# # Point analysis 
if(analysis == "point"){
  legend("topright", legend = "Calculating target volume..........................", bg = "white")

  # # take the HL value for the last point on the path
  # target.HL = HL[length(HL)]
  # # check if there is a point with a lower H/L value elsewhere on the path
  # if(target.HL > min(HL)) target.HL = min(HL)

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
  text(X[length(X)], Y[length(Y)], "Point of interest", pos = 4)
}

png(filename = "May16Profile.png", width = 2244, height = 1628, res = 300)

plot(L[1:12],path.z[1:12], type = "l", main = "May 16 Event Profile", xlab = "Runout Distance (m)", ylab = "Elevation (m)")

dev.off()

# 
# ############################################################################################
# # create a plot showing the breakpoints along the profile
# 
# # find the path distance for each P(E) breakpoint
# # find the values of X.PE and Y.PE that are valid
# X.valid = which(is.finite(X.PE))
# Y.valid = X.valid
# 
# # calculate the distance between the previous digitized point on the profile and the P(E) breakpoint
# path.increment = sqrt((X[index.PE[X.valid]]-X.PE[X.valid])^2 + (Y[index.PE[Y.valid]]-Y.PE[Y.valid])^2)
# 
# # calculate the change in elevation 
# H.increment = numeric(length = length(path.increment))
# index = 1
# 
# for(i in min(X.valid):max(X.valid)) {
#   H.increment[index] = (path.z[index.PE[i]])+(path.z[index.PE[i]+1] - path.z[index.PE[i]])*((HL[index.PE[i]] - HL.optim[i])/(HL[index.PE[i]] - HL[index.PE[i]+1]))
#   index = index + 1
# } 
# 
# H.increment = max(path.z) - H.increment
# 
# # find the P(E) for the first line segment
# PE.valid = target.PE[X.valid]
# 
# # initialize segments for plotting
# L100to95 = 0
# L96to68 = 0
# L68to50 = 0
# L50to34 = 0
# L34to5 = 0
# L5to1 = 0
# Lless1 = 0
# 
# H100to95 = 0
# H96to68 = 0
# H68to50 = 0
# H50to34 = 0
# H34to5 = 0
# H5to1 = 0
# Hless1 = 0
# 
# 
# L100to95 = c(path.dist[1],path.dist[index.PE[1]]+path.increment[1])
# H100to95 = c(H[1], H.increment[1])
#   
# if(index.PE[1]==index.PE[2]) {
#   L95to68 = c(path.dist[index.PE[1]]+path.increment[1], path.dist[index.PE[2]]+path.increment[2])
#   H95to68 = c(H.increment[1], H.increment[2])
# } else {
#   L95to68 = c(path.dist[index.PE[1]]+path.increment[1], path.dist[(index.PE[1]+1):index.PE[2]], path.dist[index.PE[2]]+path.increment[2])
#   H95to68 = c(H.increment[1], H[(index.PE[1]+1):index.PE[2]], H.increment[2])
# }
#   
# if(index.PE[2]==index.PE[3]) {
#   L68to50 = c(path.dist[index.PE[2]]+path.increment[2], path.dist[index.PE[3]]+path.increment[3])
#   H68to50 = c(H.increment[2], H.increment[3])
# } else {
#   L68to50 = c(path.dist[index.PE[2]]+path.increment[2], path.dist[(index.PE[2]+1):index.PE[3]], path.dist[index.PE[3]]+path.increment[3])
#   H68to50 = c(H.increment[2], H[(index.PE[2]+1):index.PE[3]], H.increment[3])
# }
#   
# if(index.PE[3]==index.PE[4]) {
#   L50to34 = c(path.dist[index.PE[3]]+path.increment[3], path.dist[index.PE[4]]+path.increment[4])
#   H50to34 = c(H.increment[3], H.increment[4])
# } else {
#   L50to34 = c(path.dist[index.PE[3]]+path.increment[3], path.dist[(index.PE[3]+1):index.PE[4]], path.dist[index.PE[4]]+path.increment[4])
#   H50to34 = c(H.increment[3], H[(index.PE[3]+1):index.PE[4]], H.increment[4])
# }
#   
# if(index.PE[4]==index.PE[5]) {
#   L34to5 = c(path.dist[index.PE[4]]+path.increment[4], path.dist[index.PE[5]]+path.increment[5])
#   H34to5 = c(H.increment[4], H.increment[5])
# } else {
#   L34to5 = c(path.dist[index.PE[4]]+path.increment[4], path.dist[(index.PE[4]+1):index.PE[5]], path.dist[index.PE[5]]+path.increment[5])
#   H34to5 = c(H.increment[4], H[(index.PE[4]+1):index.PE[5]], H.increment[5])
# }
#   
# if(index.PE[5]==index.PE[6]) {
#   L5to1 = c(path.dist[index.PE[5]]+path.increment[5], path.dist[index.PE[6]])
#   H5to1 = c(H.increment[5], H[index.PE[6]])
# } else {
#   L5to1 = c(path.dist[index.PE[5]]+path.increment[5], path.dist[(index.PE[5]+1):index.PE[6]])
#   H5to1 = c(H.increment[5], H[(index.PE[5]+1):index.PE[6]])
# }
# 
# 
# 
# SF=0.9
# windows(width = SF*11, height = SF*5)
# par(mar=c(4, 4, 0, 1))
# par(oma = c(0, 0, 0, 2))
# 
# plot(path.dist, H, xlab = "L (m)", ylab = "H (m)", ylim = c(max(H), min(H)))
# lines(L100to95, H100to95, col = "black", lwd = 3)
# lines(L95to68, H95to68, col = "darkred", lwd = 3)
# lines(L68to50, H68to50, col = "red", lwd = 3)
# lines(L50to34, H50to34, col = "orange", lwd = 3)
# lines(L34to5, H34to5, col = "yellow", lwd = 3)
# lines(L5to1, H5to1, col = "blue", lwd = 3)
# 
# legend("topright", title = "Probability of Exceedance", lty = c(1, 1, 1, 1, 1, 1, 1), lwd = c(3, 3, 3, 3, 3, 3, 3),
#        legend = (c("> 0.95", "0.95 - 0.68", "0.68 - 0.5", "0.5 - 0.34", "0.34 - 0.05", "0.05 - 0.01")),
#        col = c("black", "darkred", "red", "orange", "yellow", "blue", "green"), bg = "white")
# 
# SF=0.9
# windows(width = SF*11, height = SF*5)
# par(mar=c(4, 4, 0, 1))
# par(oma = c(0, 0, 0, 2))
# 
# plot(path.dist, c(1,HL[2:length(HL)]), type = "l", xlab = "L (m)", ylab = "H/L", ylim = c(0.1, 1))
# lines(c(-500,max(L100to95, na.rm = TRUE)), c(HL.optim[1], HL.optim[1]), lty = 2)
# lines(c(-500,max(L95to68, na.rm = TRUE)), c(HL.optim[2], HL.optim[2]), lty = 2)
# lines(c(-500,max(L68to50, na.rm = TRUE)), c(HL.optim[3], HL.optim[3]), lty = 2)
# lines(c(-500,max(L50to34, na.rm = TRUE)), c(HL.optim[4], HL.optim[4]), lty = 2)
# lines(c(-500,max(L34to5, na.rm = TRUE)), c(HL.optim[5], HL.optim[5]), lty = 2)
# lines(c(-500,max(L5to1, na.rm = TRUE)), c(HL.optim[6], HL.optim[6]), lty = 2)
# 
# # output X, Y and Z values to re-run an analysis if needed
# output = data.frame(X, Y, path.z)
# write.csv(output, file = "xyz_coords.csv")