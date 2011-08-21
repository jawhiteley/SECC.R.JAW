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
cat('  - Loading fauna data files.\n')

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
cat('  - Cleaning raw input.\n')
## strip rows with empty species names
SECC.fauna.raw <- SECC.fauna.raw[!is.na(SECC.fauna.raw$ID), ]

## clean species names
IDs.formal <- SECC.fauna.raw$ID
IDs.formal <- sub( "\\s\\(imm\\)$", "" , IDs.formal)
IDs.formal <- sub( "\\s\\(beetle\\)", "" , IDs.formal)
IDs.formal <- sub( "\\s\\(Diptera\\?\\)", "" , IDs.formal)
IDs.formal <- sub( "beetle", "Coleoptera spp." , IDs.formal)
SECC.fauna.raw$ID <- IDs.formal

##================================================
## Extract meta-data
##================================================
cat('  - Extracting fauna species metadata.\n')
## clone first metadata columns to 'SECC.fauna.sp': species records, aliases, taxanomic categories, etc.
SECC.fauna.meta <- SECC.fauna.raw[, -grep("^X\\d", colnames(SECC.fauna.raw))]
                                        # keep only columns that do not start with "X" and a number (these are sample IDs).
## strip empty rows from meta-data
SECC.fauna.meta <- strip_empty_dims(SECC.fauna.meta, dim=1, cols = 2:ncol(SECC.fauna.meta))
## Mdata.cols <- intersect( colnames(SECC.fauna.raw), colnames(SECC.fauna.meta))[-1]
Mdata.cols <- which( colnames(SECC.fauna.raw) %in% colnames(SECC.fauna.meta))[-1]

## strip species metadata columns
SECC.fauna <- SECC.fauna.raw[, -Mdata.cols]

##################################################
## PROCESS DATA
##################################################
cat('  - Processing fauna data.\n')
## strip empty columns (NA or 0s)
SECC.fauna <- strip_empty_dims(SECC.fauna, dim=2)
sp.rows <- which(SECC.fauna$ID %in% SECC.fauna.meta$ID)
cols.0  <- which( apply( SECC.fauna[sp.rows, ], 2, function(x) all(x==0) ) )
SECC.fauna <- SECC.fauna[, -cols.0]

## strip 0-only rows
sample.cols <- 2:ncol(SECC.fauna)
rows.0 <- which( apply( SECC.fauna[, sample.cols], 1, function(x) all(x==0) ) )
SECC.fauna <- SECC.fauna[-rows.0, ]

##================================================
## Transpose input data table
##================================================
cat('    - Transposing fauna data.\n')
## rownames(SECC.fauna) <- cleanVarName(SECC.fauna$ID)
sp.names   <- make.names(SECC.fauna$ID, unique = TRUE)
sample.IDs <- colnames(SECC.fauna)[-1]
rownames(SECC.fauna) <- sp.names
## SECC.fauna.mat <- as.matrix(SECC.fauna)
SECC.fauna.t <- t(SECC.fauna[, -1])         # omit ID column (replaced by rownames)
SECC.fauna.t <- cbind(SampleID = sample.IDs, SECC.fauna.t)
rownames(SECC.fauna.t) <- NULL  # remove annoying rownames (stored in sample.IDs)
## colnames(SECC.fauna.t)[1] <- "SampleID"
SECC.fauna.df <- as.data.frame(SECC.fauna.t, stringsAsFactors = FALSE)

## Convert character columns to appropriate data type
if (FALSE) {
  ## This does not work as expected
  ## I can either coerce everything, replacing character columns with NULL (omit else)
  ## or it returns the original object with the else statement.
  ## grr.
  ## ANSWER: Can not mix data types in a matrix (as returned by apply).  dur.  That's what data frames are for!
  SECC.fauna.df <- apply(SECC.fauna.df, 2, function(x) {
                         if (any(is.na(as.numeric(x))) == FALSE) as.numeric(x) else x
    ## convert to numeric if ALL values are numeric (no NAs introduced by coercion). 
  })
  str(SECC.fauna.df)
  head(SECC.fauna.df)

  cols.num <- which( apply(SECC.fauna.t, 2, function(x) any(!is.na(as.numeric(x)))) )
  SECC.fauna.num <- apply(SECC.fauna.t[, cols.num], 2, function(x) as.numeric(x) )
  SECC.fauna.df <- cbind(SECC.fauna.t[, -cols.num], SECC.fauna.num)  # merge numeric & non-numeric columns
  SECC.fauna.df <- as.data.frame(SECC.fauna.df)
}

for (i in 1:ncol(SECC.fauna.df)) {
    coli <- SECC.fauna.df[[i]]
    if (suppressWarnings( any(is.na(as.numeric(coli))) ) == FALSE) {
      coli <- as.numeric(coli)
    } else if (suppressWarnings( any(is.na(as.logical(coli))) ) == FALSE) {
      coli <- as.logical(coli)
    } else if (suppressWarnings( any(is.na((coli))) ) == FALSE) {
      coli <- as.factor(coli)
    }
    SECC.fauna.df[[i]] <- coli
}

## Same for Meta-data
for (i in 2:ncol(SECC.fauna.meta)) {    # skip first column (leave as character data) - all rows should be unique.
    coli <- SECC.fauna.meta[[i]]
    if (suppressWarnings( any(is.na(as.numeric(coli))) ) == FALSE) {
      coli <- as.numeric(coli)
    } else if (suppressWarnings( any(is.na(as.logical(coli))) ) == FALSE) {
      coli <- as.logical(coli)
    } else if (suppressWarnings( any(is.na((coli))) ) == FALSE) {
      coli <- as.factor(coli)
    }
    SECC.fauna.meta[[i]] <- coli
}


##################################################
## CHECK & CLEAN DATA
##################################################
cat('  - Checking & Cleaning fauna data structure.\n')
## Check SECC structure
SECC.fauna <- checkSECCdata(SECC.fauna.df, "SECC.fauna.df")
## re-generate Sample IDs
SECC.fauna <- within(SECC.fauna, {
                     SampleID <- paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep="")
  })
rownames(SECC.fauna) <- SECC.fauna$SampleID
## recode factors with meaningful labels
## SECC.fauna <- recodeSECC( SECC.fauna )  # standard function in `./lib/fun.r`  

## SECC.fauna still includes many extra columns with sample meta-data not needed for ordination analyses.
## stripping these should be fairly easy: keep only the first column (SampleIDs - or rely on these as rownames), and all numeric columns (species counts).


## House-keeping
rm(list=c('coli', 'rows.0', 'cols.0', 'i', 'IDs.formal', 'Mdata.cols', 'sample.IDs', 'SECC.fauna.t', 'SECC.fauna.df', 'SECC.fauna.raw', 'sp.rows'))
