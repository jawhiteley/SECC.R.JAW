##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process data files stored in "./data/"
# Jonathan Whiteley		R v2.12		2011-08-21
##================================================
## Faunal data is loaded into 'SECC.fauna'
## Fauna data is kept separate mostly because 
## of the number of individual variables (species),
## and to keep the data matrix separate for specific
## analyses of community data (e.g. ordination).
## Aggregate response variables (group totals) 
## may be included in 'SECC' data frame 
## for univariate analyses.
##################################################
if (FALSE) {        # do not run automatically
  rm(list=ls())     # house-keeping
  setwd('./ SECC/') # project directory
  getwd()           # check current wd

  ## LOAD LIBRARIES
  source("./lib/fun.R")   # define functions
  require(car)            # load external package 'car', for recode()
  require(reshape)        # sort_df (sort data frame) wrapper for order
}

##################################################
## LOAD DATA FILES
##################################################
cat('- Loading fauna data.\n')

SECC.fauna.raw <- read.csv(".data/fauna_t4_raw.csv", na.strings = c("NA", "", ".") )

##################################################
## CLEAN INPUT
##################################################
## strip rows with empty species names

##================================================
## Extract meta-data
##================================================
## clone first metadata columns to 'SECC.fauna.sp': species records, aliases, taxanomic categories, etc.

## strip species metadata columns

## strip empty columns (NA or 0s)

## strip 0-only rows

##################################################
## PROCESS DATA
##################################################
## Transpose input data table

##################################################
## CHECK & CLEAN DATA
##################################################
## Check SECC structure


