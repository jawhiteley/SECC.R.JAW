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
  str(Plot.vectors)                    # column names may have to be cleaned (smybols may cause errors)
  Raw.xy <- Plot.xy                    # save a copy
  Raw.vectors <- Plot.vectors          # save a copy
  Plot.xy <- Raw.xy                    # restore
  Plot.vectors <- Raw.vectors          # restore
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
                  skip <- as.logical(skip)
                  skip[is.na(skip)] <- FALSE
})

## Standardize ID column names & values
## SECC.xy <- checkSECCdata(Plot.xy, "SECC.spatial")




##################################################
## CALCULATIONS
##################################################
ReverseVector <- function (Vector.df, Vector.i) {
      NewVector <- Vector.df[Vector.i, ]
      OldFrom   <- NewVector$From
      OldAngle  <- NewVector$direction
      NewVector$From      <- NewVector$To
      NewVector$To        <- OldFrom
      NewVector$direction <- OldAngle + 180
      Vector.df[Vector.i, ] <- NewVector
      Vector.df
}

Calc.xy <- function (Start.xy, Vector.xy) {
      ## calculate destination coordinates
      Pdist   <- Vector.xy$distance
      Pangle  <- Vector.xy$direction * (pi/180)
      dx      <- sin(Pangle) * Pdist
      dy      <- cos(Pangle) * Pdist
      EndX    <- Startxy$mx + dx 
      EndY    <- Startxy$my + dy 
      c(EndX, EndY)
}

PlotGPS <- Plot.xy                     # save a copy of original GPS coordinates, just in case
CalcXY <- rep(NA, nrow(Plot.xy))       # to track which plots have had destinations calculated?

GPSXY <- which(Plot.xy$GPS == TRUE)
Plot.xy[-GPSXY, c("mx", "my", "GPS")] <- NA # wipe non-GPS coordinates
## GPS coordinates have error associated with them as well.  
## Fortunately, I also have vector information between them.
## Perhaps I should start by adjusting these initial GPS coordinates first?
## Selected points have been "de-selected" by removing the GPS flag: their positions are open for re-calculation based on vector data
GPSfrom <- Plot.vectors$From %in% Plot.xy$Plot[GPSXY]
GPSto   <- Plot.vectors$To   %in% Plot.xy$Plot[GPSXY]
GPSvectors <- GPSfrom & GPSto
Plot.vectors[GPSvectors, ]
Bvectors <- grep(".+B", Plot.vectors$From)
Plot.vectors <- ReverseVector(Plot.vectors, Bvectors) # set B plots to 'destinations' only

KnownXY <- GPSXY # which(!is.na(Plot.xy$GPS))  # first set of starting points
while (length(KnownXY) > 0) {
  ## Make new, blank data frame to store new calculated positions?
  ## if positions are added to the same data frame throughout, would this upset the 'KnownXY' indices?  Unless new positions are added at the *END*
  for (i in KnownXY) {
    Startxy <- Plot.xy[i, ]            # average of current positions?
    ## Collect all plots whose positions depend on currently known positions
    Destinations <- which(Plot.vectors$From == Plot.xy$Plot[i])
    ## - include "From" AND "To" vectors?  Reverse "To" vectors?
    ReVectors    <- which(Plot.vectors$To   == Plot.xy$Plot[i])
    for (r in ReVectors) {
      NewVector <- Plot.vectors[r, ]
      OldFrom   <- NewVector$From
      OldAngle  <- NewVector$direction
      NewVector$From      <- NewVector$To
      NewVector$To        <- OldFrom
      NewVector$direction <- OldAngle + 180
      Plot.vectors[r, ]   <- NewVector
    }
    Destinations <- c(Destinations, ReVectors)
    KeepVectors  <- which(Plot.vectors[Destinations, "skip"] != TRUE)
    Destinations <- Destinations[KeepVectors] # drop vectors marked to skip
    for (d in Destinations) {
      Pvector <- Plot.vectors[d, ]       # get vector data
      EndPlot <- Pvector$To              # which plot are we going to?
      Endi    <- which(Plot.xy$Plot == EndPlot) # convert label to index (easier)
      Endxy   <- Plot.xy[Endi, ]           # current data
      ## calculate destination coordinates
      Pdist   <- Pvector$distance
      Pangle  <- Pvector$direction * (pi/180)
      dx      <- sin(Pangle) * Pdist
      dy      <- cos(Pangle) * Pdist
      EndX    <- Startxy$mx + dx 
      EndY    <- Startxy$my + dy 
      ## Reverse vectors to prioritize C over B chambers?
      ## - switch "From" & To" if "To" is a GPS coord, or "From" is a B plot.
      ## Adjust *both* Start & End positions if neither are GPS coords?
      ## Average subsequent coordinates, or use only first set?
      ## ideally, I'd like to weight some vectors over others 
      ## i.e. vectors from GPS points take priority over calculated ones, 
      ## or closer vectors take priority over vectors farther away...
      ## store each calculated position separately, with fields for:
      ## - distance (from source)
      ## - was source a GPS position?
      ## - computed weights, based on above information (=2/distance; GPS = 1)
      ## aggregate data frame using weighted means.
      ##       Endxy$mx <- mean(c(EndX, Endxy$mx) , na.rm = TRUE)
      ##       Endxy$my <- mean(c(EndY, Endxy$my) , na.rm = TRUE)
      ##       browser()
      if (is.na(Endxy$GPS) | is.na(Endxy$mx) | is.na(Endxy$my)) {
        ## & (!is.na(EndX) & !is.na(EndY))
        Endxy$mx  <- EndX
        Endxy$my  <- EndY
        Endxy$GPS <- FALSE
      }
      if (Endxy$GPS != TRUE) Plot.xy[Endi, ] <- Endxy # update values (if not a GPS value)
      if (is.na(CalcXY[Endi]))  CalcXY[Endi] <- FALSE # flag for use as future Start point
    }
    CalcXY[i] <- TRUE
  }
  ## merge new calculated positions to 'official' data frame
  KnownXY <- which(CalcXY == FALSE) # new set of starting points
}
## Note: comparing coordinates produced by averaging vs. first vector only
##       highlights possible data entry errors (plots that jump a lot because of errors in distance or direction)
## Plots with problematic coordinates / vectors
## Distances between GPS locations mismatch measured distances: 
## 14A-24A (10.7; 12), 34A-44A (7 ; 11.4), 44A-42A (8.4 ; 4.2), 54A-51A (6.5 ; 16.2), 64A-62A (13.5 ; 8.1), 84A-81A (7.5, 20.6)
## Change wildly if multiple vectors are averaged vs. only initial: 32C*, 33C, 83C*, 74C~ (71C~)
## Abnormally close to each other:  42B, 42C, 22B / 33A / 33B
## Have vectors with missing distance or direction (may end up with NA coords, depending on algorithm):  22C, 41B


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
row.names(SECC.xy) <- SECC.xy$Plot




##################################################
## CHECK DATA
##################################################
if (FALSE) {
  distxy <- dist(SECC.xy[, c("mx", "my")])
  dist.mat <- as.matrix(distxy)        # coerce to matrix for easier sub-setting
  dist.mat[GPSXY, GPSXY]
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

