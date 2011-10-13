##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# standard functions used in most analyses
# Jonathan Whiteley      R v2.12      2011-03-28
##################################################
## access these functions in another file by using: 
## 	 source("path/to/this/file.R")

##################################################
## ANOVA: standardized analysis of individual response variables
##################################################
# Standard univariate Analysis of Variance on NESTED treatments: 
# No fancy mixed effects or error distribution families (GLMM) yet.
# The same experimental design applies for each response variable, 
# therefore a standard analytical procedure applies.
# Major differences include details such as:
#   - which Transformation is appropriate
#   - Graphs of significant Main Effects

##================================================
## DATA PROCESSING
##================================================
# See `/lib/fun.R` for data processing functions
# used to produce a data frame with standard columns
# for standardized analysis steps.

SECC_sampleID <- function( factors=NULL ) {
  ## Auto-concatenate factor codes into standard SampleIDs.
  ## the supplied data frame had BETTER contain all necessary columns.
  factors.names <- names(factors)
  factors.present <- c("Chamber", "Time", "Frag", "Pos", "Position") %in% names(factors)
  factors.present[4] <- any( factors.present[c(4,5)] )
  factors.present <- factors.present[1:4]
  if ( all( factors.present ) == FALSE ) {
    stop( paste("ID columns missing:", 
                factors.names[which(factors.present==FALSE)]
    ))
  }
  
  pos.sep <- "." # default
  if( "Pos" %in% names(factors) ) {
    pos.col <- factors$Pos
  } else if( "Position" %in% names(factors) ) {
    pos.col <- factors$Position  # would I ever really do this?
  } else {
    pos.sep <- ""
    pos.col <- ""
  }

  SampleIDs <-  with(factors,
                     paste(Block, Time, Chamber,
                           "-", Frag, pos.sep, pos.col,
                           sep = ""
                           )
                     )

  return(SampleIDs)
}

SECC_aggregate <- function( dataf = NULL , trt = "Frag", drop.trt = FALSE, ... ) {
  ## aggregate full SECC data to means for each level of specified treatment
  ## - default: Frag (for Meta-Community scale data)
  ## optional: specify a drop argument to aggregate across all other levels, 
  ##		   except for the dropped level(s)
  
  if (is.null(drop.trt)) {
	dropt.trt <- FALSE
  } else {
	if (drop.trt==TRUE) {
	  drop.trt <- trt                     # 'trt' will be dropped
	} else if (!is.logical(drop.trt)) {
	  drop.trt -> trt                     # use 'drop.trt' for 'trt'
	  drop.trt <- TRUE
	}
  }
  data.cols <- names(dataf)
  col.types <- lapply(dataf, class)
  agg.cols  <- names(dataf)[(col.types=="numeric" | col.types=="integer")]  # not sure if integer is necessary, but ...
  agg.by    <- with(dataf, list(Block    = Block,
								Time     = Time,
								Chamber  = Chamber,
								Frag     = Frag,
								Position = Position
								)
                    )
  if (trt %in% names(agg.by) == FALSE) {
    stop( paste("Treatment \'", trt, "\' is not a valid column.  Choose one of:\n",
                names(agg.by),
                sep = "" )
         )
  }
  agg.lvl <- match(TRUE, names(agg.by) %in% trt ) # == or %in%? %in% allows vector of trts
  if (drop.trt==FALSE) {
	agg.by <- agg.by[1:agg.lvl[1]]        # use only the first value, if more than one (avoids Warning messages).
  } else {
	agg.by <- agg.by[-agg.lvl]
  }
  
  data.agg <- with(dataf,
                  aggregate(dataf[, agg.cols], by=agg.by, mean, na.rm = TRUE, ... )
                  )
  ## note that the mean of only NAs is NaN.
  return(data.agg)
}


SECCclean <- function(data=NULL,
                      Time.lvls     = 1:3,
                      Chamber.lvls  = c("A", "B", "C"),
                      Frag.lvls     = 1:4,
                      Position.lvls = c("Inner", "other", "Outer")
                      )
{
  ## clean empty and unused data from SECC data frame

  ## strip empty rows (rows with only NAs)
  SECC.lvls <- strip_empty_dims( data, dim = 1, col.class = "numeric" )  

  ## Filter data by treatment levels.
  SECC.lvls <- SECC.lvls[SECC.lvls$Time     %in% Time.use     &
                         SECC.lvls$Chamber  %in% Chamber.use  &
                         SECC.lvls$Frag     %in% Frag.use     &
                         SECC.lvls$Position %in% Position.use 
                         , ]
  ## drop unused factor levels (for plotting)
  SECC.lvls <- within( SECC.lvls, {
                      Time     <- factor(Time,     levels = Time.lvls)
                      Chamber  <- factor(Chamber,  levels = Chamber.lvls)
                      Frag     <- factor(Frag,     levels = Frag.lvls)
                      Position <- factor(Position, levels = Position.lvls)
                      })
  return(SECC.lvls)
}


##================================================
## DATA EXPLORATION
##================================================



##================================================
## DEFINE MODEL FORMULA
##================================================




##================================================
## CHECK ASSUMPTIONS: residuals, standard diagnostic plots
##================================================




##================================================
## ANALYSIS: GET RESULTS
##================================================



##================================================
## PLOTS & GRAPHS
##================================================

plotMap <- function (factor = c("Chamber", "Frag", "Position"), 
                     labels = c() ) 
{
  ## allow labels to be a subset of factor levels, and auto-drop labels here?
  ## strings should *NOT* be stored as factors, or they will not be recognized as strings, and cause problems when trying to use them in graphing functions.  Alternatively, I could probably store these as matrices or tables, rather than data frames, but data frames make more sense to me for this purpose.
  
  factor <- match.arg(factor)

  if (factor == "Chamber") {
    if (is.null(labels)) labels <- c("A", "B", "C")
    PlotMap <- 
      data.frame(label = labels, 
                 col = c("#000000","#000099","#990000"), 
                 bg  = c("#FFFFFF","#FFFFFF","#FFFFFF"), 
                 pch = c(21, 23, 18), lty = c(3, 2, 1),
                 stringsAsFactors = FALSE
                 )
    ## A) Ambient = black, open circles with dotted line ; 
    ## B) Partial = blue, open diamonds with dashed line ; 
    ## C) Full    = red, solid diamond with solid line.
  }

  if (factor == "Frag") {
    if (is.null(labels)) labels <- c("1", "2", "3", "4")
    PlotMap <- 
      data.frame(label = labels, 
                 col = c("#000000", "#666666", "#000099", "#990000"), 
                 bg  = c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF"), 
                 pch = c(19, 15, 22, 21),
                 lty = c(1, 2, 3, 3),
                 stringsAsFactors = FALSE
                 )
    ## 1) Continuous         = black, filled circles with solid line ; 
    ## 2) Full Corridors     =  grey, filled squares with dashed line ; 
    ## 3) Pseudo-Corridors   =  blue, open squares with dotted line.
    ## 4) Isolated           =   red, open circles with dotted line.
  }

  if (factor == "Position") {
    if (is.null(labels)) labels <- c("I", "*", "O")
    PlotMap <- 
      data.frame(label = labels, 
                 col = c("#000000", "#000099", "#990000"), 
                 bg  = c("#FFFFFF", "#FFFFFF", "#FFFFFF"), 
                 pch = c(19, 8, 21), lty = c(2, 3, 1),
                 stringsAsFactors = FALSE
                 )
    ## 1) Inner   = black filled circles with dotted line ; 
    ## 2) other   = blue           stars with dashed line ; 
    ## 3) Outer   = red   open   circles with solid  line.
  }  

  return(PlotMap)
}

ggPts.SECC <- function (ptMap = plotMap("Chamber"), name = "Chamber Treatment") {
  ## Add Chamber shapes & colours to a ggplot scatterplot
  require(ggplot2)
  plotTemplate <- list(scale_colour_manual(name = name,
                                           values = ptMap$col, 
                                           breaks = levels(ptMap$label)
                                           ),
                       scale_shape_manual(name = name,
                                          values = ptMap$pch, 
                                          breaks = levels(ptMap$label)
                                          )
                       )
  ## should just be breaks = Chamber.map$label, but that produces right-aligned text :(
}

