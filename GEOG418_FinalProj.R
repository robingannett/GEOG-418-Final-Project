#Libraries
install.packages("spgwr")
install.packages("spatstat")
install.packages("tmap")
install.packages("gstat")
install.packages("sf")
install.packages("raster")
install.packages("rgdal")
install.packages("e1071")
install.packages("spdep")
install.packages("gtable")
install.packages("gridExtra")
install.packages("grid")
install.packages("ggplot2")
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)
library(gtable)
library(gridExtra)
library(grid)
library(ggplot2)
library(rgeos)

#Set working directory
dir <- ""
setwd(dir)


# Data prep
#Reading in particulate matter dataset
#Read in PM2.5 data:
pm2.5 <- readOGR("Pm25Sample.shp") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))


#Reading in dissemination tract and income data
#Read in census income data:
income <- read.csv("Income.csv")  
#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income") 
#Read in dissemination tract shapefile:
census.tracts <- readOGR("BC_DA.shp") 
#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Determine the number of columns in the dataframe:
nrow(income.tracts)
#Remove NA values:
income.tracts <- income.tracts[!is.na(income.tracts$Income),]
#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))

#Create choropleth map of income:
map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income",
              style = "jenks",
              palette = "YlGnBu", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM")) + 
  tm_shape(pm2.5) +
  tm_dots(col = "PM25", title = "Mean Annual PM2.5", style = "jenks", palette = "YlOrRd", size = 0.2) +
  tm_layout(main.title.size = 1, main.title = "Median Income and Mean Annual PM2.5 in Vancouver (2016)", main.title.position = "center", bg.color = "grey73") +
  tm_compass(type = "arrow", position = c(0.9, 0.9), size = 3) +
  tm_scale_bar(breaks = c(0, 5, 10), text.size = 0.6, position = "LEFT", "TOP") 

map_Income



# Testing for SAC of Income using Global and Local Moran's I

#Defining Neighbourhood
income.nb <- poly2nb(income.tracts) #Defines the neighbours using Queen's case
income.net <- nb2lines(income.nb, coords=coordinates(income.tracts)) #Converts the neighbours list into a network/line graph
crs(income.net) <- crs(income.tracts)

#Creating weights matrices
income.lw <- nb2listw(income.nb, zero.policy = TRUE, style = "W") 
print.listw(income.lw, zero.policy = TRUE)

#runs a Global Moran's I test on the data using the income variable and the list of neighbour weights
mi <- moran.test(income.tracts$Income, income.lw, zero.policy = TRUE) 
mi

#Calculating the Z score
mI <- mi$estimate[[1]] 
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- ((mI-eI)/sqrt(var))


I <- signif(mI, digits=3)
EI <- signif(eI, digits=3)
Variance <- signif(var, digits=3)
Z <- signif(z, digits=3)

data.for.moran = data.frame(I, EI, Variance, Z)

tableMoran <- tableGrob(data.for.moran, rows = c(""), cols = c("I", "E(I)", "Variance", "Z")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Global Moran's I results using census data\nin Vancouver (2016) to test for\nspatial autocorrelation of income", gp = gpar(fontsize = 08))
padding <- unit(5, "mm")

tableMoran <- gtable_add_rows(tableMoran, 
                              heights = grobHeight(t1Caption) + padding, 
                              pos = 0)

tableMoran <- gtable_add_grob(tableMoran,
                              t1Caption, t = 1, l = 2, r = ncol(data.for.moran) + 1)


grid.arrange(tableMoran, newpage = TRUE)

#Calculating the Local Moran's I statistic using the Income variable and the list of neighbour weights
lisa.test <- localmoran(income.tracts$Income, income.lw, zero.policy = TRUE)

#Adding the Local Moran's I, E(Ii), variance, Z, and p values to the Census dataset
income.tracts$Ii <- lisa.test[,1]
income.tracts$E.Ii<- lisa.test[,2]
income.tracts$Var.Ii<- lisa.test[,3]
income.tracts$Z.Ii<- lisa.test[,4]
income.tracts$P<- lisa.test[,5]

#Mapping the Local Moran's I results across the study area using the Z value
tmap_mode("view")

map_LISA <- tm_shape(income.tracts, name = "Vancouver") + 
  tm_layout(main.title.size = 1, main.title = "Spatial Autocorrelation of Income in Vancouver (2016)", main.title.position = "center") +
  tm_polygons(col = "Z.Ii", 
              title = "Z Value", 
              style = "fixed", 
              auto.palette.mapping=FALSE,
              palette = "PuOr", breaks = c(-Inf, -1.96, 1.96, Inf )) #manually sets the breaks to desired intervals

map_LISA


# Creating a raster surface of air pollution by spatial interpolation of PM2.5 values


#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(pm2.5, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
# Create SpatialPixel object:
gridded(grd)     <- TRUE  
# Create SpatialGrid object:
fullgrid(grd)    <- TRUE  
#Reproject the grid:
proj4string(grd) <- proj4string(income.tracts)



##Spatial Interpolation with Kriging

# Define the trend model
f.0 <- as.formula(PM25 ~ 1) 

#Create variogram
var.smpl <- variogram(f.0, pm2.5, cloud = FALSE, cutoff = 50000) 
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=0.032, model="Sph", range=30500, nugget=0))
plot(var.smpl, dat.fit, xlab = "Distance (m)", ylab = "Semivariance")

# Perform the kriging interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg <- krige(f.0, pm2.5, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
rm.kriging <- mask(r, income.tracts)

tmap_mode("plot")

# Plot the map
tm_Krig <- tm_shape(rm.kriging) + 
  tm_raster(n=6, palette="OrRd",  
            title="Predicted PM2.5 \n(in ug/m3)") +
  tm_shape(pm2.5) + tm_dots(size=0.1) +
  tm_legend(legend.outside=TRUE) +
  tm_layout(main.title.size = 1, main.title = "Interpolated Surface of PM2.5 using Ordinary Kriging, Vancouver (2016)", main.title.position = "left")

tm_Krig

#### Looks at variance
r   <- raster(dat.krg, layer="var1.var")
r.m <- mask(r, income.tracts)

tm_var <- tm_shape(r.m) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n(in squared ppm)") +tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
tm_var

#### Converts variance into a confidence interval
r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m <- mask(r, income.tracts)

tm_KrigCI <- tm_shape(r.m) + 
  tm_raster(n=6, palette ="Reds",
            title="95% CI map \n(in ppm)") +tm_shape(pm2.5) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE) + 
  tm_layout(main.title.size = 1, main.title = "Interpolated surface of PM2.5 using Ordinary Kriging, Vancouver (2016)", main.title.position = "center")

tm_KrigCI


#Extract average pm2.5 for each polygon
income.tracts$Pm2.5 <- round(extract(rm.kriging, income.tracts, fun = mean)[,1], 5)

#Removing negative PM2.5 values
income.tracts <- income.tracts[!is.na(income.tracts@data$Pm2.5),]
for(i in 1:nrow(income.tracts)){
  if(income.tracts@data[i,30] < 0){
    income.tracts@data[i,30] <- 0
  }
}

View(income.tracts@data)


######Linear Regression##########

plot(income.tracts$Income~income.tracts$Pm2.5)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
income.tracts.no0 <-  income.tracts[which(income.tracts$Pm2.5 > 0), ]

#Now plot the data again
plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
#Add the regression model to the plot you created
plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5, xlab = "PM2.5 (ug/m3)", ylab = "Income")
abline(lm.model, col = "red")
#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
income.tracts.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
income.tracts.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(income.tracts.no0)

tmaptools::palette_explorer()

tmap_mode("view")
#Now, create choropleth map of residuals
map_resid <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "residuals",
              title = "Residuals",
              style = "jenks",
              palette = "PRGn", n = 6, auto.palette.mapping = TRUE) +
  tm_layout(main.title.size = 0.8, main.title = "Regression Analysis Residuals of Mean Annual PM2.5 and Median Income in Vancouver (2016)", main.title.position = "center")

map_resid


# Run Global Moran's I test for SAC of residuals

#Defining Neighbourhood
residuals.nb <- poly2nb(income.tracts.no0) #Defines the neighbours using Queen's case
residuals.net <- nb2lines(residuals.nb, coords=coordinates(income.tracts.no0)) #Converts the neighbours list into a network/line graph
crs(residuals.net) <- crs(income.tracts.no0)


#Creating weights matrices
residuals.lw <- nb2listw(residuals.nb, zero.policy = TRUE, style = "W") 
print.listw(residuals.lw, zero.policy = TRUE)

#runs a Global Moran's I test on the data using the income variable and the list of neighbour weights
miRes <- moran.test(income.tracts.no0$residuals, residuals.lw, zero.policy = TRUE) 
miRes


#Calculating the Z score
mIi <- miRes$estimate[[1]] 
eIi <- miRes$estimate[[2]]
vari <- miRes$estimate[[3]]

zi <- ((mIi-eIi)/sqrt(vari))


I <- signif(mIi, digits=3)
EI <- signif(eIi, digits=3)
Variance <- signif(vari, digits=3)
Z <- signif(zi, digits=3)

data.for.moran = data.frame(I, EI, Variance, Z)

tableMoran <- tableGrob(data.for.moran, rows = c(""), cols = c("I", "E(I)", "Variance", "Z")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Global Moran's I results for\nSpatial Autocorrelation test of Residuals from Regression Analysis of Mean Annual PM2.5 and Median Income,\nVancouver (2016)", gp = gpar(fontsize = 08))
padding <- unit(5, "mm")

tableMoran <- gtable_add_rows(tableMoran, 
                              heights = grobHeight(t1Caption) + padding, 
                              pos = 0)

tableMoran <- gtable_add_grob(tableMoran,
                              t1Caption, t = 1, l = 2, r = ncol(data.for.moran) + 1)


grid.arrange(tableMoran, newpage = TRUE)


####Geographically Weighted Regression
#Let's say you are continuing with 
#your data from the regression analysis. 
#The first thing you need to do is to add the 
#polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the 
#"coordinates" function from the sp library
income.tracts.no0.coords <- sp::coordinates(income.tracts.no0)
#Observe the result:
head(income.tracts.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
income.tracts.no0$X <- income.tracts.no0.coords[,1]
income.tracts.no0$Y <- income.tracts.no0.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                        data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

gwr#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
income.tracts.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style = "jenks",
              palette = "PuBuGn", n = 6) +
  tm_layout(main.title.size = 0.8, main.title = "R2 Values for GWR Analysis of Mean Annual PM2.5 and Median Income in Vancouver (2016)", main.title.position = "center")
map_r2

#Time for more magic. Let's map the coefficients
income.tracts.no0$coeff <- results$income.tracts.no0.Pm2.5

map_coef <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              auto.palette.mapping=FALSE,
              palette = "RdBu", n= 6)
map_coef

range(income.tracts.no0$coeff)


### Point pattern analysis to check if sample points of pm2.5 are randomly located or not


pm2.5.n <- as.numeric(pm2.5$PM25)

pm2.5.coords <- coordinates(pm2.5)

pm2.5.n$x <- pm2.5.coords[,1]
pm2.5.n$y <- pm2.5.coords[,2]

win.ext <- as.matrix(extent(income.tracts))
window <- as.owin(list(xrange = win.ext[1,], yrange = win.ext[2,]))
pm2.5.n.ppp <-ppp(x=pm2.5.n$x, y=pm2.5.n$y,window=window)



plot(pm2.5.n.ppp)

nearestNeighbour <- nndist(pm2.5.n.ppp)

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"

pm2.5.n.ppp$n
nrow(pm2.5)

N <- nrow(pm2.5)

##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
nnd = (sum(nearestNeighbour$Distance))/(N)

#mean nearest neighbour for random spatial distribution
VanArea.m <- gArea(income.tracts)

studyArea <- VanArea.m
studyArea/1000000
pointDensity <- N/studyArea

r.nnd = 1/(2*sqrt(pointDensity))

d.nnd = 1.07453/(sqrt(pointDensity))

R = nnd/r.nnd

SE.NND <- 0.26136/(sqrt(N*pointDensity))

z = (nnd-r.nnd)/(SE.NND)

NND <- round(nnd, 3)
NNDr <- round(r.nnd, 3)
NNDd <- round(d.nnd, 3)
Zn <- round(z, 3)

data.for.neigh = data.frame(NND, NNDr, NNDd, Zn)

tablen <- tableGrob(data.for.neigh, rows = c("")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 2: Nearest Neighbour Analysis Results for\nPM2.5 Sample Points, Vancouver (2016)", gp = gpar(fontsize = 08))
padding <- unit(5, "mm")

tablen <- gtable_add_rows(tablen, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

tablen <- gtable_add_grob(tablen,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.neigh) + 1)

grid.arrange(tablen, newpage = TRUE)

##K-FUNCTION 
#basic k-function
k.fun <- Kest(pm2.5.n.ppp, correction = "Ripley")
plot(k.fun, main = "K-Function Plot for PM2.5 Sample Points")


#use simulation to test the point pattern against CSR
k.fun.e <- envelope(pm2.5.n.ppp, Kest, nsim = 99, correction = "Ripley")
plot(k.fun.e, main = "K-Function with Confidence Intervals for PM2.5 Sample Points")


##### Descriptive statistics for income and air pollution
#Mean
meanIncome <- mean(income.tracts$Income, na.rm = TRUE) #This is going to produce a wrong value (NA) due to a single NA value in data
meanpm2.5 <- mean(pm2.5$PM25, na.rm = TRUE) #Use na.rm = TRUE to ignore NA values in calculation


#Standard Deviation
sdIncome <- sd(income.tracts$Income, na.rm = TRUE) #Calculate the SD, ignoring NA values
#sdSummer <- sd(subset(df_2020, IGN_Day >= 183 & IGN_Day <= 244)$CURRENT_SI) #Calculate the SD, ignoring NA values only for the summer months
sdpm2.5 <- sd(pm2.5$PM25, na.rm = TRUE)

#Median
medIncome <- median(income.tracts$Income, na.rm = TRUE)
medpm2.5 <- median(pm2.5$PM25, na.rm = TRUE)

#Min
minIncome <- min(income.tracts$Income)
minPm2.5 <- min(pm2.5$PM25)

#Max
maxIncome <- max(income.tracts$Income)
maxPm2.5 <- max(pm2.5$PM25)

#####
#Create a table of descriptive stats

Variable = c("Income", "PM2.5") #Create an object for the labels
Min = c(minIncome, minPm2.5)
Max = c(maxIncome, maxPm2.5)
means = c(meanIncome, meanpm2.5) #Create an object for the means
sd = c(sdIncome, sdpm2.5) #Create an object for the standard deviations
median = c(medIncome, medpm2.5) #Create an object for the medians

Mean <- round(means, 3)
SD <- round(sd, 3)
Median <- round(median, 3)

data.for.table1 = data.frame(Variable, Min, Max, Mean, Median, SD)
data.for.table2 = data.frame(samples, skewness, kurtosis, CoV, normality)

#Make table 1
table1 <- tableGrob(data.for.table1, rows = c("","")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Descriptive Statistics of Median Income\nand Mean Annual PM2.5 in Vancouver, 2016", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)

grid.arrange(table1, newpage = TRUE)

####THE END!!!!####





