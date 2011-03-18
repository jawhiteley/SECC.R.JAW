##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process data files stored in "./data/"
# Jonathan Whiteley		R v2.12		2011-01-26
##================================================
## All Data files in './data' are loaded into memory.
## Load scripts for manual processing are source()'d
##   as necessary: cleaning, processing, calculations, etc.
##   This may also generate a new object suitable for analysis
##   and mergining, by combining several others.
##   - redundant objects should be removed from memory.
## Data objects whose names start with 'SECC.'
##   are checked and merged onto 'SECC.base'
##   -> SECC : Main data frame of response variables.
## Data objects starting with 'SECC_env'
##   are checked and merged onto 'SECC.base' -> SECC.env
##   -> SECC.env : Environmental data.
## Faunal data is loaded into 'SECC.fauna'
## Data objects starting with 'TRH.'
##   are checked and merged as a Time Series.
##   -> SECC.TRH : Temperature & Relative Humidity Time Series.
##################################################

# setwd('./ SECC/')   # project directory
getwd()             # check current wd

## LOAD LIBRARIES
source("./lib/fun.R")   # define functions
require(car)		# load external package 'car', for recode()
require(reshape)	# sort_df (sort data frame) wrapper for order

##################################################
## LOAD DATA FILES
##################################################
load('./save/SECC_factors.R')	# includes SECC.base data.frame, 
  # and other vectors of standard column names and levels.
SECC.raw <- read.csv("./save/SECC_base.csv") # just to check how it is imported.  
  # Should be the same values as SECC.base, but raw column types (only character columns are factors).

Data_files <- dir('./data/')  # get list of file names.
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
rm( list=c('File_name', 'File_path','temp','Object_name') ) # clean-up

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
##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
# If files need Special Attention (individual cleaning), 
# that goes here; or in a separate script that is source()'d here.
#   Specify which columns are to be merged with attribute "SECC columns"

source("./lib/load_ARA.R", echo=FALSE)  # Process ARA N-fixation data ; hide output?

##================================================
## AUTOMATICALLY CHECK & CLEAN DATA
##================================================
# Check & merge only objects beginning with "SECC."
merge.SECC <- Data_objects[substr(Data_objects, 1, 5)=="SECC."]  

for (DataObject in merge.SECC)
{
  assign( DataObject, checkSECCdata( get(DataObject), DataObject ) )
    # Returns an updated clean version of data object.
    # It will fail with an error message if the object is unsuitable for merging.
}
# What does this warning mean exactly?
#   In min(levels(data[[ColName]]) %in% Col_lvls) :
#   no non-missing arguments to min; returning Inf
# It might result from a column of NAs being compared?

##################################################
## PROCESS DATA
##################################################
##================================================
## Merge all main response variables into a single data frame: SECC
##================================================
SECC <- SECC.base  # Base template into which relevant frames will be merged.

attr(SECC, "labels") <- list("Time" = "Sample Time",
                             "Chamber" = "Chamber Treatment",
                             "Frag" = "Fragmentation Treatment",
                             "Pos" = "Patch Position"
                             )

SECC.by <- names(SECC.base)
for (ObjectName in merge.SECC)
{
  DataObject <- get(ObjectName)
  # Keep only columns for respons variables, denoted by 'SECC' prefix in column name
  KeepCols <- attr(DataObject, "SECC columns") # extract list of column names from attributes (set manually during load).
  KeepCols <- c(SECC.by, KeepCols)
  DataObject <- DataObject[,KeepCols]
  ## Merge data frame into SECC, keeping all rows of SECC, but NOT rows in target that do not match SECC.
  SECC <- merge( SECC, DataObject, by=SECC.by, all.x=TRUE, all.y=FALSE, sort=TRUE )
  # merge attributes: labels, units.
  attr(SECC, "labels") <- c( attr(SECC, "labels"), attr(DataObject, "labels") )
  attr(SECC, "units" ) <- c( attr(SECC, "units" ), attr(DataObject, "units" ) )
}
# Proper sort order
SECC <- within( SECC, {
  # levels in sort order.  The Position column will be replaced with recoded values later anyway.
  Position <- factor(Pos, levels=Pos_all_sort)
})
SECC <- sort_df( SECC, vars=c(Trt_sort_order, "Position") )

##================================================
## Calculate columns across source data frames in new object
##================================================
ARA.Nfix.ratio <- 1/3  # ratio of N-fixation : ARA.

SECC <- within( SECC, {
  # ARA per gram dry weight of sample.
# ARAg <- ARA / ARA.dwt
# Nfix <- ARA * ARA.Nfix.ratio
})


##================================================
## Generate factor columns, re-order factor levels, etc.
##================================================
SECC <- recodeSECC( SECC )  # standard function in `./lib/fun.r`  

### Summarize data by means across different (or all) Positions to prevent unbalanced effects?
# SECCr <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, Time=Time, Chamber=Chamber, Frag=Frag), mean ) )	# for regional-level analyses (ignoring Position)
# SECC <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, Time=Time, Chamber=Chamber, Frag=Frag, Position=Position), mean ) )	# using cbind() on the response variables allows multiple columns to be summarized, and also preserves column names.


##################################################
## CHECK DATA
##################################################
## Patch-level data
head(SECC)		# have a peek at the first 6 rows & columns: is this what you expected?
str(SECC)		# check structure: are the appropriate variables factors, numeric, etc.?

# ## Regional data
# head(SECCr)		# have a peek at the first 6 rows & columns: is this what you expected?
# str(SECCr)		# check structure: are the appropriate variables factors, numeric, etc.?

##################################################
## SAVE DATA
##################################################
# Save data to native R file "./save/SECC_data.R", or leave in memory:
# + SECC       - data from experimental patches.  **Main table for analysis**
# + SECC.env   - Environmental data corresponding to SECC.
# + SECC.fauna - Microarthropod community data corresponding to SECC.
# + SECC.TRH   - Temperature & Relative Humidity (time-series) data.
# + [Other]
load_export <- c( 'SECC', 'SECC.env', 'SECC.fauna', 'SECC.TRH' )
save( list=load_export, file="./save/SECC_data.R" )

# Export data frames to csv, just in case.
for (DataFrame in load_export)
{
  filename <- gsub('.', '_', DataFrame, perl = TRUE)  # replace '.' with '_' 
  write.csv( get(DataFrame), 
    file = paste("./save/", filename, ".csv", sep=""), 
    row.names=FALSE
  ) 
}
