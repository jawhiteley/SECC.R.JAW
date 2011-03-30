##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process data files stored in "./data/"
# Jonathan Whiteley		R v2.12		2011-03-23
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
# rm(list=ls())       # house-keeping
# setwd('./ SECC/')   # project directory
getwd()             # check current wd

## LOAD LIBRARIES
source("./lib/fun.R")   # define functions
require(car)            # load external package 'car', for recode()
require(reshape)        # sort_df (sort data frame) wrapper for order

##################################################
## LOAD DATA FILES
##################################################
cat('- Loading data files.\n')
load('./save/SECC_factors.R')	# includes SECC.base data.frame, 
  # and other vectors of standard column names and levels.
SECC.raw <- read.csv("./save/SECC_base.csv") # just to check how it is imported.  
  # Should be the same values as SECC.base, but raw column types (only character columns are factors).

Data_files <- dir('./data/')  # get list of file names.
Data_objects <- NULL          # empty object to hold vector of object names.
for (File_name in Data_files) {
  File_path <- paste("./data/", File_name, sep="")
  temp <- read.csv(File_path)  # na.strings = c("NA", "-", "", ".")
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
# str(SECC.raw)
# head(SECC.raw)

# str(SECC.base)
# head(SECC.base)


##################################################
## CHECK & CLEAN DATA
##################################################
cat('- Processing individual data objects.\n')
##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
# If files need Special Attention (individual cleaning), 
# that goes here; or in a separate script that is source()'d here.
#   Specify which columns are to be merged with attribute "SECC columns"

source("./lib/clean_cyanobacteria.R", echo=FALSE)  # Process Cyanobacteria data.
source("./lib/clean_ARA.R", echo=FALSE)            # Process ARA N-fixation data.

source("./lib/clean_H2O.R", echo=FALSE)            # Process Moisture data.


##================================================
## AUTOMATICALLY CHECK & CLEAN DATA
##================================================
cat('- Auto-checking data.\n')
# Check & merge only objects beginning with "SECC."
merge.SECC <- Data_objects[substr(Data_objects, 1, 5)=="SECC."]  
for (DataObject in merge.SECC) {
  assign( DataObject, checkSECCdata( get(DataObject), DataObject ) )
    # Returns an updated clean version of data object.
    # It will fail with an error message if the object is unsuitable for merging.
}


##################################################
## PROCESS DATA
##################################################
cat('- Merging data.\n')
## Prepare base data frame to merge others into
SECC <- SECC.base  # Base template into which relevant frames will be merged.

## Create containers for attributes.
 # attributes are lost after every merge,
 # so they are accumulated in a separate object to be added later.
SECC.labels <- list("Time" = "Sample Time",
                    "Chamber" = "Chamber Treatment",
                    "Frag" = "Fragmentation Treatment",
                    "Pos" = "Patch Position"
                    )
SECC.units  <- list()

##================================================
## Merge all main response variables into a single data frame: SECC
##================================================

SECC.by <- names(SECC.base)
  # [-1] - drop first row (SampleID is not re-coded or checked?)
  # I need SampleID to store extra information about samples:
  # especially duplicate IDs for controls.
  # In theory, this could be the *only* column I need to use to merge,
  # but I figure it's safer to use all the main ID columns.
for (ObjectName in merge.SECC) {
  DataObject  <- get(ObjectName)
  # merge attributes: labels, units.
  SECC.labels <- c( SECC.labels, attr(DataObject, "labels") )
  SECC.units  <- c( SECC.units , attr(DataObject, "units" ) )
  # Keep only columns for respons variables, denoted by 'SECC' prefix in column name
  KeepCols <- attr(DataObject, "SECC columns") # extract list of column names from attributes (set manually during load).
  KeepCols <- c(SECC.by, KeepCols)
  DataObject <- DataObject[,KeepCols]  # attributes lost?
  ## Merge data frame into SECC, keeping all rows of SECC, but NOT rows in target that do not match SECC.
  SECC <- merge( SECC, DataObject, by=SECC.by, all.x=TRUE, all.y=FALSE, sort=TRUE )
}
rm(list=c('SECC.by', 'ObjectName', 'KeepCols', 'DataObject')) # house-keeping

## Preferred sort order
SECC <- within( SECC, {
  # Levels in sort order: The Position column will be replaced with recoded values later anyway.
  Position <- factor(Pos, levels=Pos_all_sort)
})
# Sort, but also drops attributes :-(
SECC.sorted <- sort_df( SECC, vars=c(Trt_sort_order, "Position") )

## Re-build SECC data frame properties
SECC <- SECC.sorted
row.names(SECC) <- SECC$SampleID  # Change row.names to at least something meaningful.
attr(SECC, "labels") <- SECC.labels
attr(SECC, "units" ) <- SECC.units


##================================================
## Calculate columns across source data frames in new object
##================================================
cat('- Final calculations.\n')
Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.
sampleA  <- 6	# sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
ARA.m2	 <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * (12.5^2)      # patch area
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2

SECC <- within( SECC, {
  # ARA per gram dry weight of sample.
  # Dry weights are only available for samples with cyanobacteria data :(
  ARA.g <- ARA.ml / ARA.dwt  
  Nfix  <- ARA.m * Nfix.ARA.ratio
})

attr(SECC, "labels")[["Nfix"]] <- "N-fixation rate"
attr(SECC, "units" )[["Nfix"]] <- quote(mu*"mol" %.% m^-2 %.% d^-1)
## attr(SECC, "units" )[["Nfix"]] <- quote(mu*"mol" %.% g^-1 %.% d^-1)

##================================================
## Generate factor columns, re-order factor levels, etc.
##================================================
SECC.coded <- SECC          # keep a copy, just in case
SECC <- recodeSECC( SECC )  # standard function in `./lib/fun.r`  

### Summarize data by means across different (or all) Positions to prevent unbalanced effects?
# SECCr <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, Time=Time, Chamber=Chamber, Frag=Frag), mean ) )	# for regional-level analyses (ignoring Position)
# SECC <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, Time=Time, Chamber=Chamber, Frag=Frag, Position=Position), mean ) )	# using cbind() on the response variables allows multiple columns to be summarized, and also preserves column names.


##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d
  head(SECC)		# have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC)		# check structure: are the appropriate variables factors, numeric, etc.?
  SECCstr(SECC)
  ## For a more accurate summary, excluding unsampled Time Points:
  SECCstr(SECC[!is.na(SECC$Time), ])
  ## For a more accurate summary, excluding t3, B's & Other Pos's
  ## ALL Column counts here should be equal to the 'Total' counts
  SECCstr(SECC[SECC.coded$Time!=3 &
               SECC.coded$Chamber %in% c("A", "C") &
               SECC.coded$Pos     %in% c("I", "O"),
               ])

  invisible(edit(SECC.coded))

# ## Regional data
# head(SECCr)		# have a peek at the first 6 rows & columns: is this what you expected?
# str(SECCr)		# check structure: are the appropriate variables factors, numeric, etc.?
}

##################################################
## SAVE DATA
##################################################
# Save data to native R file "./save/SECC_data.R", or leave in memory:
# + SECC       - data from experimental patches.  **Main table for analysis**
# + SECC.env   - Environmental data corresponding to SECC.
# + SECC.fauna - Microarthropod community data corresponding to SECC.
# + SECC.TRH   - Temperature & Relative Humidity (time-series) data.
# + [Other]
cat('- Saving data & cleaning up.\n')

Load.export <- c( 'SECC', 'SECC.coded')  # , 'SECC.env', 'SECC.fauna', 'SECC.TRH' )
save( list=Load.export, file="./save/SECC_data.R" )

# Export data frames to csv, just in case.
for (DataFrame in Load.export)
{
  Filename <- gsub('\\.', '_', DataFrame, perl = TRUE)  # replace '.' with '_' 
  write.csv( get(DataFrame), 
    file = paste("./save/", Filename, ".csv", sep=""), 
    row.names=FALSE
  ) 
}

## House-keeping
rm(list=c('SECC.sorted', 'DataFrame', 'Filename', 'rm.objects'))
try( detach('package:car',     character.only = TRUE) )
try( detach('package:reshape', character.only = TRUE) )

cat('= FINISHED Loading data. =\n')
