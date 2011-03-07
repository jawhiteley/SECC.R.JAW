##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process data files stored in "./data/"
# Jonathan Whiteley		R v2.12		2011-01-26
##################################################

# setwd('./ SECC/')   # project directory
getwd()             # check current wd

## LOAD LIBRARIES
source("./lib/fun.R")   # define functions
library(car)		# load external package 'car', for recode()

##################################################
## LOAD DATA FILES
##################################################
load('./save/SECC_factors.R')	# includes SECC.base data.frame, 
  # and other vectors of standard column names and levels.
SECC.raw <- read.csv("./save/SECC_base.csv") # just to check how it is imported.  
  # Should be the same values as SECC.base, but raw column types (only character columns are factors).

Data_files <- dir('./data/')
Data_objects <- NULL          # empty object to hold vector of object names.
for (File_name in Data_files) {
  File_path <- paste("./data/", File_name, sep="")
  temp <- read.csv(File_path)
  Object_name <- cleanVarName(File_name)
  # remove file extension from the object name  
  Object_name <- gsub("\\.csv\\b", "", Object_name, perl=TRUE ) 
  # assign data to object with similar name as the file.
  assign(Object_name, temp)
  # collect a list of object names for further processing.
  Data_objects <- c(Data_objects, Object_name)
}
rm( list=c('File_name', 'File_path','temp','Object_name') )

##================================================
## CHECK DATA
##================================================
str(SECC.raw)
head(SECC.raw)

str(SECC.base)
head(SECC.base)


##################################################
## CHECK & CLEAN DATA
##################################################
# If files need Special Attention (individual cleaning), 
# that goes here; or in a separate script that is source()'d here.

Data_Check <- TRUE	# default value: set to FALSE if a check fails.

## Standardize ID column names and types (based on template)
ColNames_std <- colnames(SECC.base)	# could use colnames(SECC.base) or Trt_nest_order
if (  min( colnames(SECC.base) %in% ColNames_std ) == 0 ) {
  # min(boolean_vector) == 0 means there was at least on FALSE result
  # Check for most likely non-standard names
  # Attempt to rename based on likely matches.
  # If column names still do not match, throw an Error message: this file needs Special Attention.
}

## Standardize ID column types & values (levels)
  # Check that column types & values (levels) match template
  # Attempt to convert column types if necessary
  # Check for unexpected Warnings: 
  # Unexpected Warnings means this file needs Special Attention.


## Standardize factor levels
if ( min(levels(SECC.base$Pos) %in% Pos_all_lvls) == 0 ) {
  # min(boolean_vector) == 0 means there was at least on FALSE result
  
}

# Check that all IDs match EXACTLY with Template (or that all rows in data are represented in the template)
  # Check number of rows - data files may only be a subset
  # If fewer rows, are they a subset: can missing rows be replaced with NA's?
# If IDs still do not match, throw an Error message: this file needs Special Attention.

# If all checks pass, continue.  
# Otherwise, throw an Error message: this file needs Special Attention.

if Data_Check == FALSE stop(
  paste("ERROR: There is an unknown problem with data frame: ", DataFrame)
)

##################################################
## PROCESS DATA
##################################################
##================================================
## Merge all main response variables into a single data frame: SECC
##================================================


##================================================
## Generate factor columns, re-order factor levels, etc.
##================================================
SECC <- within( SECC.raw, {
  ## Rename columns /convert to standard informative names for the rest of this script
  Block   <- factor(Block)
  Time    <- factor(Time)
  Chamber <- factor(Warming)
  Frag    <- factor(Frag)
  Pos     <- factor(Pos)
  ## rename and reorder factor levels the easy way - maintains empty values if empty factor specified, otherwise converts to 'NA'.  Requires package 'car'
  Chamber <- recode( Warming, 
    "'A'='Ambient'; 'B'='Partial Chamber'; 'C'='Full Chamber'; else=''", 
    levels=c( "Ambient", "Partial Chamber", "Full Chamber" ),
    as.factor.result=TRUE
  )
	levels(Frag) <- list( 'Continuous'=1, 'Full Corridors'=2, 'Pseudo-Corridors'=3, 'Isolated'=4 )
	Position <- factor(Pos, levels=c('1', 'S', 'W', 'E', 'N', '0'))	# safely reorder factor levels
	levels(Position) <- list( 'I'='1', 'S'='S', 'W'='W', 'E'='E', 'N'='N', 'O'='0' )	# rename some factor levels (omitted levels are dropped and replaced with empty strings).
	# New factor with simplified recoded values for Patch Position
	Pos <- recode( Position, 
		"'I'='Inner'; 'O'='Outer'; else='other'", 
		levels=c( "Inner", "other", "Outer" ),
		as.factor.result=TRUE
	)
})


## Summarize data by means across different (or all) Positions to prevent unbalanced effects?
SECCr <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, TimePt=TimePt, Chamber=Chamber, Frag=Frag), mean ) )	# for regional-level analyses (ignoring Position)
SECC <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, TimePt=TimePt, Chamber=Chamber, Frag=Frag, Pos1=Pos1), mean ) )	# using cbind() on the response variables allows multiple columns to be summarized, and also preserves column names.


##################################################
## CHECK DATA
##################################################
## Patch-level data
head(SECC)		# have a peek at the first 6 rows & columns: is this what you expected?
str(SECC)		# check structure: are the appropriate variables factors, numeric, etc.?

## Regional data
head(SECCr)		# have a peek at the first 6 rows & columns: is this what you expected?
str(SECCr)		# check structure: are the appropriate variables factors, numeric, etc.?

# Save data to native R file in "./save/", or leave in memory.
