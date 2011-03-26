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
  pos.sep <- "." # default
  if( "Pos" %in% names(factors) ) {
    pos.sep <- "."
    pos.col <- factors$Pos
  } else if( "Position" %in% names(factors) ) {
    pos.sep <- "."
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

SECC_aggregate <- function( data = NULL , trt = "Frag", ... ) {
  ## aggregate full SECC data to means for each level of specified treatment
  ## - default: Frag (Meta-Community)
  
  data.cols <- names(data)
  col.types <- lapply(data, class)
  agg.cols  <- names(data)[(col.types=="numeric" | col.types=="integer")]  # not sure if integer is necessary, but ...
  agg.by    <- with( data, list(Block    = Block,
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
  agg.lvl <- match(TRUE, names(agg.by) == trt )
  agg.by  <- agg.by[1:agg.lvl]
  
  data.mc <- with(data,
                  aggregate(data[, agg.cols], by=agg.by, mean, na.rm = TRUE, ... )
                  )
  ## note that the mean of only NAs is NaN.
  return(data.mc)
}

strip_empty_dims  <- function( data = NULL, dim = c(1, 2), 
							  rows = NULL, cols = NULL, col.class = NULL ) {
  ## strip whole rows or columns if they are entirely empty (NA or NaN)
  ## `data` is the data frame to strip empty values from.
  ## `rows` is a vector of indices or logical values of 
  ##        which rows to check for empty values (default all)
  ## `cols` see `rows` above (but applied to columns)
  ## `dim`  the dimensions to strip if all empty
  ## `col.class` is a vector of classes: columns that match will be included
  ##             This applies to the `cols` argument.

  if (!is.data.frame(data) && !is.matrix(data)) stop(
	  "The first argument must be a \"data.frame\" or \"matrix\"."
     )

  dim = match.args(dim)
  cols <- as.vector(cols)
  if ( is.null(cols) ) {
	cols <-  1:ncol(data)
  }
  if ( !is.null(col.class) ) {
	check.cols <- which( lapply(data[, cols], class) %in% as.vector(col.class) )
	check.cols <- cols[check.cols]
  } else {
	check.cols <- cols
  }
  rows <- as.vector(rows)
  if ( is.null(rows) ) {
	rows <-  1:nrow(data)
  }  
  
  empty.rows <- which( apply( data[, check.cols], 1, function(x) all(is.na(x)) ) )
  empty.cols <- which( apply( data[check.rows, ], 2, function(x) all(is.na(x)) ) )

  if (all(dim == 1)) {
	data <- data[-empty.rows, ]
  } else if (all(dim == 1))  {
	data <- data[, -empty.cols]
  } else if (all(dim %in% c(1, 2))) {
	data <- data[-empty.rows, -empty.cols]
  }
  
  return(data)
  
  if (FALSE) {  # Testing
	test.df <-  data.frame( )
	df.proc <- strip_empty_dims(test.df)
	df.proc
  }
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
chamberMap <- function ( labels=c("A", "B", "C") ) {
  Chamber.map <- data.frame( label=labels, 
    col = c("#000000","#000099","#990000"), 
    bg  = c("#FFFFFF","#FFFFFF","#FFFFFF"), 
    pch = c(21,23,18), lty = c(3,2,1) 
  )
    # A) Ambient = black, open circles with dotted line ; 
    # B) Partial = blue, open diamonds with dashed line ; 
    # C) Full    = red, solid diamond with solid line.
}


