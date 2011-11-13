##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process spatial data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-11-11
##################################################
## This script is run as part of `./lib/load.R`
if (FALSE) {        # do not run automatically
  rm(list=ls())     # house-keeping
  setwd('./ SECC/') # project directory
  getwd()           # check current wd
}

##================================================
## CHECK DATA
##================================================
## Should be 2 data frames:
## - Plot_xy: xy coordinates, in UTM mE & mN
## - Plot_vectors : distances (m) & directions (degrees) from plots with known coordinates (GPS)
## I need to use the data in PlotVectors to compute missing coordinates in SECC.spatial
if (FALSE) {
  str(Plot.xy)
  str(Plot.vectors)                     # column names may have to be cleaned (smybols may cause errors)
}

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.
names(Plot.xy)[names(Plot.xy)=="time"] <- "Time"
Plot.xy <- within(Plot.xy, {
                  GPS <- as.logical(GPS)
                  Plot <- as.character(Plot)
                  mx <- as.double(mx)
                  my <- as.double(my)
})
Plot.vectors <- within(Plot.vectors, {
                  comments <- as.character(comments)
                  From <- as.character(From)
                  To   <- as.character(To)
})

## Standardize ID column names & values
## SECC.xy <- checkSECCdata(Plot.xy, "SECC.spatial")




##################################################
## CALCULATIONS
##################################################
## still buggy.  Some plots overlapping.  Check the From & To vectors?
PlotGPS <- Plot.xy                     # save a copy of original GPS coordinates, just in case
CalcXY <- rep(NA, nrow(Plot.xy))       # to track which plots have had destinations calculated?

GPSXY <- which(Plot.xy$GPS == TRUE)
KnownXY <- which(!is.na(Plot.xy$GPS))  # first set of starting points
while (length(KnownXY) > 0) {
  for (i in KnownXY) {
    Destinations <- which(Plot.vectors$From == Plot.xy$Plot[i])
    Startxy <- Plot.xy[i, ]
    for (d in Destinations) {
      Pvector <- Plot.vectors[d, ]       # get vector data
      EndPlot <- Pvector$To              # which plot are we going to?
      Endi <- which(Plot.xy$Plot == EndPlot) # convert label to index (easier)
      Endxy <- Plot.xy[Endi, ]           # current data
      ## calculate destination coordinates
      Pdist   <- Pvector$distance
      Pangle  <- Pvector$direction * (pi/180)
      dx      <- sin(Pangle) * Pdist
      dy      <- cos(Pangle) * Pdist
      EndX    <- Startxy$mx + dx 
      EndY    <- Startxy$my + dy 
      ## average subsequent coordinates, or use only first set?
      ## ideally, I'd like to weight some vectors over others 
      ## i.e. vectors from GPS points take priority over calculated ones, 
      ## or closer vectors take priority over vectors farther away...
      ## store each calculated position separately, with fields for:
      ## - distance (from source)
      ## - was source a GPS position?
      ## computed weights, based on above information
      ##       if (is.na(Endxy$mx)) 
      Endxy$mx <- mean(c(EndX, Endxy$mx) , na.rm = TRUE)
      ##       if (is.na(Endxy$my)) 
      Endxy$my <- mean(c(EndY, Endxy$my) , na.rm = TRUE)
      if (is.na(Endxy$GPS)) Endxy$GPS <- FALSE
      if (Endxy$GPS != TRUE) Plot.xy[Endi, ] <- Endxy # update values (if not a GPS value)
      if (is.na(CalcXY[Endi]))  CalcXY[Endi] <- FALSE # flag for use as future Start point
    }
    CalcXY[i] <- TRUE
  }
  KnownXY <- which(CalcXY == FALSE & Plot.xy$GPS == FALSE) # new set of starting points
}


## calculate patch-level coordinates?
SECC.xy <- Plot.xy

##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.xy, "SECC columns") <- c('mx', 'my')
attr(SECC.xy, "labels") <- list("mx" = "UTM East",
                                "my" = "UTM North"
                                )
attr(SECC.xy, "units")  <- list("mx" = quote("m"),
                                "my" = quote("m")
                                )

##################################################
## CHECK DATA
##################################################
if (FALSE) {
  plot(SECC.xy$mx, SECC.xy$my, asp=1, type="n",
       xlab = attr(SECC.xy, "labels")$mx,
       ylab = attr(SECC.xy, "labels")$my,
       main = "Plot positions"
       )
  text(SECC.xy$mx, SECC.xy$my, labels=SECC.xy$Plot, cex=0.1)
}




##################################################
## SAVE DATA
##################################################

##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('Plot.xy', 'Plot.vectors', 'SECC.xy')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.xy' )

