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

SECC.fauna.raw <- read.csv("./data/fauna_t4_raw.csv", na.strings = c("NA", "", "."), as.is=TRUE)

if (FALSE) {
  str(SECC.fauna.raw)
  head(SECC.fauna.raw)
}

## SECC.fauna <- SECC.fauna.raw    # make a copy to clean & process
if (FALSE) {
  str(SECC.fauna)
  head(SECC.fauna)
}

##################################################
## CLEAN INPUT
##################################################
## strip rows with empty species names
SECC.fauna.raw <- SECC.fauna.raw[!is.na(SECC.fauna.raw$ID), ]

##================================================
## Extract meta-data
##================================================
cat('- Extracting fauna species metadata.\n')
## clone first metadata columns to 'SECC.fauna.sp': species records, aliases, taxanomic categories, etc.
SECC.fauna.meta <- SECC.fauna.raw[, -grep("^X\\d", colnames(SECC.fauna.raw))]
                                        # keep only columns that do not start with "X" and a number (these are sample IDs).
## strip empty rows from meta-data
SECC.fauna.meta <- strip_empty_dims(SECC.fauna.meta, dim=1, cols = 2:ncol(SECC.fauna.meta))
Mdata.cols <- intersect( colnames(SECC.fauna.raw), colnames(SECC.fauna.meta))[-1]
Mdata.cols <- which( colnames(SECC.fauna.raw) %in% colnames(SECC.fauna.meta))[-1]

## strip species metadata columns
SECC.fauna <- SECC.fauna.raw[, -Mdata.cols]

## strip empty columns (NA or 0s)
SECC.fauna <- strip_empty_dims(SECC.fauna, dim=2)
sp.rows <- which(SECC.fauna$ID %in% SECC.fauna.meta$ID)
cols.0  <- which( apply( SECC.fauna[sp.rows, ], 2, function(x) all(x==0) ) )
SECC.fauna <- SECC.fauna[, -cols.0]

## strip 0-only rows
sample.cols <- 2:ncol(SECC.fauna)
rows.0 <- which( apply( SECC.fauna[, sample.cols], 1, function(x) all(x==0) ) )
SECC.fauna <- SECC.fauna[-rows.0, ]

##################################################
## PROCESS DATA
##################################################
## Transpose input data table

##################################################
## CHECK & CLEAN DATA
##################################################
## Check SECC structure


