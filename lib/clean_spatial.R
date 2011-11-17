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

cat('  - Processing Spatial Data.\n')
library(plyr)                          # applying functions to subsets of data
Objects.ls <- ls()                     # List of current objects in memory (for clean-up after)

##================================================
## CHECK DATA
##================================================
## Should be 2 data frames:
## - Plot_xy: xy coordinates, in UTM mE & mN
## - Plot_vectors : distances (m) & directions (degrees) from plots with known coordinates (GPS)
Raw.xy      <- Plot.xy                 # save a copy
Raw.vectors <- Plot.vectors            # save a copy
if (FALSE) {
  str(Plot.xy)
  str(Plot.vectors)                    # column names may have to be cleaned (smybols may cause errors)
  Plot.xy      <- Raw.xy               # restore
  Plot.vectors <- Raw.vectors          # restore
}

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
cat('    - Preparing Spatial Data.\n')
distBase <- 2                          # this distance, in m, has a weight of '1' for weighted averaging of calculated positions.
names(Plot.xy)[names(Plot.xy)=="time"] <- "Time"
Plot.xy <- within(Plot.xy, {
                  GPS <- as.logical(GPS)
                  Plot <- as.character(Plot)
                  mx <- as.double(mx)
                  my <- as.double(my)
                  distance <- 0        # distance from origin for calculated positions
                  GPSstart <- FALSE    # Was the origin a GPS point?
                  wt <- 1              # weights for averaging calculated positions
})
Plot.vectors <- within(Plot.vectors, {
                  comments <- as.character(comments)
                  From <- as.character(From)
                  To   <- as.character(To)
                  skip <- as.logical(skip)
                  skip[is.na(skip)] <- FALSE
                  calc <- 0            # How many times has this vector been calculated?
                  calcFrom <- FALSE    # Has it beed used to calculate position of "From" plot?
                  calcTo   <- FALSE    # Has it beed used to calculate position of "To" plot?
})




##################################################
## CALCULATIONS
##################################################
## compute destination positions for each vector twice:
## - once in each direction (From->To, To->From)
## Stop once all vectors have been used in each direction once
## - collect list of Known Positions (available starting points)
## - calculate destinations for remaining available vectors
##   - including GPS coords? (these have some error, too...)
## - Collect all calculated positions & average together for final position
## Ideally, I'd like to weight some vectors over others 
## i.e. vectors from GPS points take priority over calculated ones, 
## or closer vectors take priority over vectors farther away...
## Store each calculated position separately, with fields for:
## - distance (from source)
## - was source a GPS position?
## - computed weights, based on above information (=2/distance; GPS = 1)
## Aggregate data frame using weighted means.
cat('    - Calculating Plot Positions.  This could take a minute or two.\n')

## Reduce effects of GPS error by averaging across GPS readings
## Use 1 GPS coordinate position at a time (loop through all, 1 @ a time)
## Compute all possible dependent (relative) positions using vector data
## - There may be multiple calculated positions for each GPS starting point,
##   but I only want to use positions calculated from the same GPS start point at a time.
## Store all calculated positions, grouped by starting point
## Average positions derived from each GPS starting point

ReverseVector <- function (Vector.df) {
  NewVector <- Vector.df               # only pass in the rows you want to reverse
  NewVector$From      <- Vector.df$To
  NewVector$To        <- Vector.df$From
  NewVector$calcFrom  <- Vector.df$calcTo
  NewVector$calcTo    <- Vector.df$calcFrom
  NewVector$direction <- Vector.df$direction + 180
  NewVector
}

CalcCoords <- function (StartXY, VectorXY) {
  EndXY   <- CalcXYrow(Plot=VectorXY$To, GPS=FALSE, 
                       distance=VectorXY$distance, GPSstart=StartXY$GPS, 
                       wt=1) # Blank row for end pt.
  EndXY   <- PlotRowFactors(EndXY)
  ## calculate destination coordinates
  Pdist   <- VectorXY$distance
  Pangle  <- VectorXY$direction * (pi/180)
  dx      <- sin(Pangle) * Pdist
  dy      <- cos(Pangle) * Pdist
  EndX    <- StartXY$mx + dx 
  EndY    <- StartXY$my + dy 
  if (FALSE) {
    c(EndX, EndY)
  } else {
    EndXY$mx <- EndX
    EndXY$my <- EndY
    ## compute weights for weighted averaging.  
    ##   NB: ifelse() works with a vector of logical values (results element-wise)!
    EndXY$wt <- max(as.numeric(EndXY$GPSstart), distBase/EndXY$distance)
    EndXY
  }
}

CalcXYrow <- function (PlotLabel=NA, mx=NA, my=NA, comments="", InitPlot=NA,
                       GPS=NA, GPSstart=NA, distance=NA, wt=1) {
  ## create a (blank) row for new calculated Plot coordinates
  data.frame(Block=NA, Time=NA, Chamber=NA, Plot=PlotLabel, 
             mx=mx, my=my, GPS=GPS, GPSstart=GPSstart, distance=distance, wt=wt, 
             InitPlot=InitPlot, comments="", stringsAsFactors=FALSE)
}

PlotRowFactors <- function (PlotXY) {
  ## extract factor levels from Plot labels
  within(PlotXY, {
         Block   <- NA
         Time    <- NA
         Chamber <- NA
         Block   <- substr(Plot, 1, 1)
         Time    <- substr(Plot, 2, 2)
         Chamber <- substr(Plot, 3, 3)
         Block   <- factor(Block)
         Time    <- factor(Time)
         Chamber <- factor(Chamber)
         })
}

## if positions are added to the same data frame throughout, would this upset the 'KnownXY' indices?  Unless new positions are added at the *END*
## However, it might be best to track based on Plot label, rather than index of Plot.xy:
##   The starting location should include the average of all (current) positions for a given Plot (not a given position in the Plot.xy data frame)
## Or should I be tracking each vector instead of each plot?
##   calculate each vector once?  or twice (once in each direction)?
CalcXY <- rep(NA, nrow(Plot.xy))       # to track which plots have had destinations calculated?
CalcXY <- data.frame(Plot=Plot.xy$Plot, Calc=FALSE, stringsAsFactors=F) # track by Plot label, rather than index

GPSXY <- which(Plot.xy$GPS == TRUE)    # Plots flagged to use as GPS starting points
GPSXY <- which(!is.na(Plot.xy$mx) & !is.na(Plot.xy$my))     # Plots with existing (GPS) coordinates
GPSplots <- Plot.xy$Plot[GPSXY]        # Plot labels with GPS-based coordinates
PlotsIni <- which(!is.na(Plot.xy$mx) & !is.na(Plot.xy$my))
IniPlots <- Plot.xy$Plot[PlotsIni] # Plots with initial GPS coordinates (whether used or not)
## Plot.xy[-GPSXY, c("mx", "my", "GPS")] <- NA # wipe existing 'non-GPS' coordinates
## Plot.xy <- Plot.xy[PlotsIni, ]            # remove unknown coordinates? Plots with no calculated positions will disappear!
## GPS coordinates have error associated with them as well.  
## Fortunately, I also have vector information between them.
## Selected points have been "de-selected" by removing the GPS flag: their positions will be re-calculated using vector data
GPSfrom <- Plot.vectors$From %in% IniPlots
GPSto   <- Plot.vectors$To   %in% IniPlots
GPSvectors <- GPSfrom & GPSto
Plot.vectors[GPSvectors, ]             # vectors between GPS locations
## Reverse vectors to prioritize C over B chambers as starting points?
Bvectors <- grep(".+B", Plot.vectors$From)
Plot.vectors[Bvectors, ] <- ReverseVector(Plot.vectors[Bvectors, ])
PlotVectors.ini <- Plot.vectors        # store a copy of processed, unused vectors.

Calc.xy <- CalcXYrow()[0, ]            # empty data frame to store Calculated coordinates
for (initPlot in IniPlots) {
  Plot.vectors <- PlotVectors.ini      # Reset vector tracking
  InitXY <- Plot.xy[Plot.xy$Plot == initPlot, ] # get original GPS starting point
  ## initialize container for calculated positions, with starting point given appropriate attributes
  Calc.ini <- CalcXYrow(Plot = initPlot, mx = InitXY$mx, my = InitXY$my,
                         GPS = TRUE, wt = 1, InitPlot = initPlot)

  Destinations <- 1:nrow(Plot.vectors) # indices of vectors to process (start with all candidates to begin the loop)
  Prev.Destinations <- 1               # Fake start value (just has to be != Destinations)
  while ((length(Destinations) > 0) & !identical(Destinations, Prev.Destinations)) {
    ## stop when all destinations have been calculated in both directions
    ## OR New Loop has the exact same Destinations as the last one (no change)
    Prev.Destinations <- Destinations  # save values from last loop
    ## Make new, blank data frame to store new calculated positions
    New.xy  <- Calc.xy[0, ]            # 0-row data frame with same columns as Calc.xy
    ## Make a list of starting points, and associated vectors to calculate target locations
    KnownXY <- !is.na(Calc.ini[, "mx"]) & !is.na(Calc.ini[, "my"])
    KnownXY <- Calc.ini$Plot[KnownXY]  # Known starting points (may be duplicates)
    KnownXY <- unique(KnownXY)
    VectorsToKeep <- Plot.vectors$skip == FALSE 
    VectorsToCalc <- Plot.vectors$calcFrom==FALSE | Plot.vectors$calcTo==FALSE
    ## Keep only Vectors based on From / To with known locations
    VectorsFrom   <- Plot.vectors$From %in% KnownXY
    VectorsTo     <- Plot.vectors$To   %in% KnownXY
    Destinations  <- which(VectorsToKeep & VectorsToCalc & (VectorsFrom | VectorsTo))  # indices of vectors to process
    ## Remove starting points that have had all related vectors processed already
    ##   Starts <- (KnownXY %in% Plot.vectors[Destinations, "From"]) |
    ##   (KnownXY %in% Plot.vectors[Destinations, "To"])
    ##   KnownXY <- KnownXY[Starts]

    for (d in Destinations) {
      Pvector <- Plot.vectors[d, ]       # get vector data
      ## check origin & destination
      ## - if necessary, reverse vector to match known starting location
      ## It is perhaps safer to reverse vectors individually here, 
      ## to keep "From" & "To" fields consistent,
      ## and allow tracking based on each direction
      ## But how do I know which direction is the one that should be calculated?
      FromKnown <- Pvector$From %in% KnownXY 
      ToKnown <- Pvector$To %in% KnownXY 
      ReVector  <- FALSE
      if (ToKnown & FromKnown) {
        ## if both destinations known, use direction that has not been calculated
        if (Pvector$calcTo & Pvector$calcFrom) {
          ## if both destinations known, and both calculated, then don't bother; skip
          next
        } else if (!Pvector$calcFrom) {
          ReVector <- FALSE
        } else if (!Pvector$calcTo  ) {
          ReVector <- TRUE
        }
      } else {
        ## If one end is unknown, start with direction to unknown location
        if (!ToKnown & !FromKnown) {
          ## neither location is known: what are w doing here?
          next
        } else if (!ToKnown  ) {
          ReVector <- FALSE
        } else if (!FromKnown) {
          ReVector <- TRUE
        }
      }
      if (ReVector) {
        Dvector <- ReverseVector(Pvector)
      } else {
        Dvector <- Pvector
      }
      if (Dvector$calcFrom) next # Already done this.

      ## get starting coordinates
      i <- which(Calc.ini$Plot == Dvector$From)
      Startxy <- Calc.ini[i, ]            # get currently known positions
      ## average of current positions 
      ## this only works because there is only 1 subset: the better way to do this using the tools in plyr (as at the end of all these loops)
      Startxy <- aggregate(Startxy[, c("mx", "my")], 
                           by=list(Plot=Startxy$Plot), 
                           FUN=weighted.mean, w = Startxy$wt, 
                           na.rm=TRUE)                            
      Startxy$GPS <- Startxy$Plot %in% GPSplots

      ## calculate destination coordinates
      Endxy <- CalcCoords(Startxy, Dvector)
      Endxy$InitPlot <- initPlot       # Assign current GPS starting point
      ## store new coordinates
      New.xy <- rbind(New.xy, Endxy)
      ## Track progress: Update vector
      if (Startxy$Plot == Pvector$From) {
        Pvector$calcFrom <- TRUE
      } else if (Startxy$Plot == Pvector$To) {
        Pvector$calcTo   <- TRUE
      }
      Pvector$calc <- Pvector$calc +1
      Plot.vectors[d, ] <- Pvector
      CalcXY[CalcXY$Plot == Dvector$From, "Calc"] <- TRUE # plot has been used to calculate dependent locations
    }
    ## add new calculated positions to main data frame (for this starting position)
    Calc.ini <- rbind(Calc.ini, New.xy)
  }
  ## computed weighted averages for this run before adding to main data frame?
  ## I have to do this eventually, but it might be more efficient to wait and do it all at once.
  ## add new calculated positions to main data frame
  Calc.xy <- rbind(Calc.xy, Calc.ini)
}

## aggregate all calculated positions using a weighted average
## http://stackoverflow.com/questions/3367190/aggregate-and-weighted-mean-in-r
## weighted means within each InitPlot run
PlotCalc <- ddply(Calc.xy, c("Plot", "InitPlot"), summarise, 
                  mx = weighted.mean(mx, w = wt, na.rm = TRUE), 
                  my = weighted.mean(my, w = wt, na.rm = TRUE) )
## Means for each Plot, over all InitPlot runs
PlotCalc <- ddply(PlotCalc, "Plot", summarise, 
                  mx = mean(mx, na.rm = TRUE), 
                  my = mean(my, na.rm = TRUE) )
PlotCalc[is.nan(PlotCalc$mx), "mx"] <- NA
PlotCalc[is.nan(PlotCalc$my), "my"] <- NA
## merge new with Raw (or CalcXY), to ensure a row for every Plot, even if NA?
PlotC.xy <- merge(PlotCalc, Raw.xy[, c("Plot", "Block")], all=TRUE) # could also include "comments" from Raw.xy, if I really want to keep those.
## Clean-up Results
PlotC.xy <- PlotRowFactors(PlotC.xy)     # ensure all factors are properly specified
PlotC.xy <- PlotC.xy[, c("Block", "Time", "Chamber", "Plot", "mx", "my")] # re-order columns



##================================================
## Calculate Patch-level Coordinates
##================================================
cat('    - Calculating Sample Patch Positions.  This could take a few seconds.\n')
Cardinals <- c("N", "S", "E", "W")
PatchDirs <- c("NW", "NE", "SE", "SW")
NE <- pi*1/4                           # Direction in radians
SE <- pi*3/4
SW <- pi*5/4
NW <- pi*7/4
Pdiam <- 12.5/100                      # Patch diameter, in m
Pdist <- 10/100                        # inter-patch distance, in m
## distance to centre of edge plots from plot centre, along x OR y
EdgeDist <- Pdist/2 + Pdiam + Pdist + Pdiam/2   
## distance from centre of plot to centre of each Frag treatment
FragC <- Pdist/2 + Pdiam + Pdist/2  # to centre of Frag trts
Fdist <- sqrt( (FragC^2) + (FragC^2) )
## distance from centre of each Frag treatment to centre of Patches
FPlotC <- Pdist/2 + Pdiam/2
FPdist <- sqrt( (FPlotC^2) + (FPlotC^2) ) # Fdist/2 ;)
FCdist <- Pdiam/2                      # distance to centres of patches in Continuous Frag trts
FCPlot <- sqrt((FCdist^2)/2)           # distance along x & y for Continuous Frag trts

## Build list of all patches, and direction from Frag/Plot centre
Patch.d <- within(SECC.base, 
                  {
                    Xdir  <- NA        # X direction: W/E
                    Ydir  <- NA        # Y direction: N/S
                    FragD <- NA        # direction from Plot to Frag trt centre
                    FPD   <- NA        # direction from Frag to Patch centre
                    dx    <- NA        # adjustment to x from Plot centre to patch, in m
                    dy    <- NA        # adjustment to y from Plot centre to patch, in m
                  }
)
Patch.d <- ddply(Patch.d, c("Block", "Time", "Chamber", "Frag"), 
                 function(FragGroup) {
                   WhichCardinals  <- which(Cardinals %in% FragGroup$Pos)
                   WhichCardinalY  <- which(Cardinals[1:2] %in% FragGroup$Pos)
                   WhichCardinalX  <- which(Cardinals[3:4] %in% FragGroup$Pos)
                   FragGroup$FragD <- paste(Cardinals[WhichCardinals], collapse="")
                   FragGroup$Xdir  <- paste(Cardinals[WhichCardinalX+2], collapse="")
                   FragGroup$Ydir  <- paste(Cardinals[WhichCardinalY  ], collapse="")
                   FragGroup
                 }
)
## loop over every row.  There's probably a more efficient way to do this, I just haven't figured it out yet.
Patch.d <- ddply(Patch.d, "SampleID", 
                 function (Prow) {
                   CardX <- Cardinals %in% Prow$Xdir
                   CardY <- Cardinals %in% Prow$Ydir
                   Cards <- c(CardY[1:2], CardX[3:4])
                   if (Prow$Pos == "O") {
                     Prow$FragD
                   } else if (Prow$Pos == "I") {
                     Cards <- !Cards
                   } else if ( any(Prow$Pos %in% Cardinals[1:2]) ) { # N,S
                     Cards <- c( Cards[1:2], !Cards[3:4])
                   } else if ( any(Prow$Pos %in% Cardinals[3:4]) ) { # E,W
                     Cards <- c(!Cards[1:2],  Cards[3:4])
                   }
                   Prow$FPD <- paste(Cardinals[Cards], collapse="")
                   Prow
                 }
)

## Calculate distances for each Patch from Plot centre
for (i in 1:nrow(Patch.d)) {
  dx <- 0
  dy <- 0
  ## calculate relative position of Frag trt centre
  xdir <- Patch.d[i, "Xdir"]
  ydir <- Patch.d[i, "Ydir"]
  dx <- if (xdir == Cardinals[4]) -FragC else FragC
  dy <- if (ydir == Cardinals[2]) -FragC else FragC
  ## add relative position of Patch centre
  PatchD <- if (Patch.d[i, "Frag"] == 1) FCPlot else FPlotC
  xdir <- substr(Patch.d[i, "FPD"], 2, 2)
  ydir <- substr(Patch.d[i, "FPD"], 1, 1)
  dx <- dx + if (xdir == Cardinals[4]) -PatchD else PatchD
  dy <- dy + if (ydir == Cardinals[2]) -PatchD else PatchD
  ## save results
  Patch.d[i, "dx"] <- dx
  Patch.d[i, "dy"] <- dy
}

## Prepare final data frame
## SECC.xy <- PlotC.xy
SECC.xy <- within(SECC.base, 
                  {
                    my  <- NA        # Y, metres UTM North
                    mx  <- NA        # X, metres UTM East
                  }
)

## Add Patch distances to Plot positions
for (i in 1:nrow(SECC.xy)) {
  PatchID    <- SECC.xy[i, "SampleID"]
  PlotString <- substr(PatchID, 1, 3)
  PlotCoords <- PlotC.xy[PlotC.xy$Plot==PlotString, ]
  PatchRow   <- which(Patch.d$SampleID == SECC.xy[i, "SampleID"])
  PatchCoord <- Patch.d[PatchRow, ]
  SECC.xy[i, "mx"] <- PlotCoords$mx + PatchCoord$dx
  SECC.xy[i, "my"] <- PlotCoords$my + PatchCoord$dy
}




##################################################
## Assign Attributes
##################################################
cat('    - Assigning Attributes to Spatial Data.\n')
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.xy, "SECC columns") <- c('mx', 'my')
attr(SECC.xy, "labels") <- list("mx" = "UTM East",
                                "my" = "UTM North"
                                )
attr(SECC.xy, "units")  <- list("mx" = quote("m"),
                                "my" = quote("m")
                                )
row.names(SECC.xy) <- SECC.xy$SampleID




##################################################
## CHECK DATA
##################################################
## Standardize ID column names & values
SECC.xy <- checkSECCdata(SECC.xy, "SECC.xy")

if (FALSE) {
  row.names(PlotC.xy) <- PlotC.xy$Plot   # makes subsetting calculated distance matrix SO MUCH easier
  distxy <- dist(PlotC.xy[, c("mx", "my")])
  dist.mat <- as.matrix(distxy)        # coerce to matrix for easier sub-setting
  dist.mat[GPSXY, GPSXY]
  dist.mat[IniPlots, IniPlots]
  min(dist.mat[dist.mat>0], na.rm=TRUE) # should be >1 (51C -- 52A)
  max(dist.mat[dist.mat>0], na.rm=TRUE) # ~200m?
  Plot.vectors[GPSvectors, ]           # vectors between GPS locations
  ShortVectors <- which(Plot.vectors$distance < 2)
  Plot.vectors[ShortVectors, ]
  ClosePlots <- c(Plot.vectors[ShortVectors, "From"], Plot.vectors[ShortVectors, "To"])
  ClosePlots <- unique(ClosePlots)
  CloseDist <- which(rownames(dist.mat) %in% ClosePlots)
  print( dist.mat[CloseDist, CloseDist], digits=2 )

  plot(PlotC.xy$mx, PlotC.xy$my, asp=1, type="n",
       xlab = attr(SECC.xy, "labels")$mx,
       ylab = attr(SECC.xy, "labels")$my,
       main = "Plot positions"
       )
  text(PlotC.xy$mx, PlotC.xy$my, labels=PlotC.xy$Plot, cex=0.1)

  plot(SECC.xy$mx, SECC.xy$my, asp=1, type="n",
       xlab = attr(SECC.xy, "labels")$mx,
       ylab = attr(SECC.xy, "labels")$my,
       main = "Sample Patch positions"
       )
  text(SECC.xy$mx, SECC.xy$my, labels=SECC.xy$SampleID, cex=0.005)
}




##################################################
## SAVE DATA
##################################################
cat('    - Saving Spatial Data to file.\n')
write.csv(format(SECC.xy, digits = 3, nsmall=3), 
          file = paste("./save/", "SECC_xy", ".csv", sep=""), row.names=FALSE)


##================================================
## Housekeeping
##================================================
cat('    - Cleaning up after Spatial Calculations.\n')
## Remove old objects from memory
New.Objects <- setdiff(ls(), Objects.ls)
rm.objects <- setdiff(New.Objects, 'SECC.xy')
rm.objects <- c(rm.objects, 'Plot.xy', 'Plot.vectors')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.xy' )

