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


##================================================
## AUTOMATICALLY CHECK & CLEAN DATA
##================================================

for (DataObject in Data_objects)
{
  assign( DataObject, checkSECCdata( get(DataObject), DataObject ) )
    # Returns an updated clean version of data object.
    # It will fail with an error message if the object is unsuitable for merging.
}

##################################################
## PROCESS DATA
##################################################
##================================================
## Merge all main response variables into a single data frame: SECC
##================================================
SECC <- SECC.base  # Base template into which relevant frames will be merged.

merge.SECC <- Data_objects[substr(Data_objects, 1, 5)=="SECC."]  # merge objects beginning with "SECC."
SECC.by <- colnames(SECC.base)
for (DataObject in Data_objects)
{
  SECC <- merge( SECC, DataObject, by=SECC.by, all=TRUE, sort=FALSE)
}

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
