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
  source("./lib/init.R")   # load basic objects (needed for environmental data)
  rm(list=ls())     # house-keeping
  setwd('./ SECC/') # project directory
  getwd()           # check current wd

  ## LOAD LIBRARIES
  source("./lib/fun.R")   # define functions
}

##################################################
## LOAD DATA FILES
##################################################
cat('  - Loading fauna data files.\n')

SECC.fauna.raw <- read.csv("./data/fauna_t4_raw.csv", na.strings = c("NA", "", "."), as.is=TRUE)

if (FALSE) {
  str(SECC.fauna.raw)
  head(SECC.fauna.raw)

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

## Make species labels R-friendly (prepare for transposing to column names)
## rownames(SECC.fauna) <- cleanVarName(SECC.fauna$ID)
sp.names   <- make.names(SECC.fauna.raw$ID, unique = TRUE)
rownames(SECC.fauna.raw) <- sp.names
SECC.fauna.raw$ID <- sp.names

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
sample.IDs <- colnames(SECC.fauna)[-1]
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

##================================================
## Check & clean data structure
##================================================
cat('    - Checking & Cleaning fauna data structure.\n')
## Check SECC structure & overwrite un-transposed data frame if successful
SECC.fauna <- checkSECCdata(SECC.fauna.df, "SECC.fauna.df")
SECC.fauna <- within(SECC.fauna, {
    Comments[is.na(Comments)] <- ""
    ## re-generate Sample IDs
    SampleID <- SECC_sampleID(SECC.fauna)
  })
## SECC.fauna <- within(SECC.fauna, {
##                      SampleID <- paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep="")
##   })

## remove missing sp IDs from meta-data
SECC.fauna.meta <- SECC.fauna.meta[which( SECC.fauna.meta$ID %in% intersect(colnames(SECC.fauna), SECC.fauna.meta$ID) ), ]
SECC.fauna.meta$sp_alias <- make.names(SECC.fauna.meta$sp_alias, unique = FALSE)

##################################################
## CALCULATIONS
##################################################
cat('  - Calculating fauna data.\n')
SECC.fauna.c <- SECC.fauna  # save count data

##================================================
## Convert species counts to # / g dwt
##================================================
## Find suitable environmental data
## (depending on where in the load script this script is called, 
## it might be stored in a different object)
Env.df <- c('SECC.H2O', 'SECC.env', 'SECC')
Env.df <- Env.df[which( sapply(Env.df, exists) )]
Env.df <- get( Env.df[1] ) # only the first one that exists
## merge in Patch.dwt from SECC(.env)
SECC.fauna.dwt <- merge(SECC.fauna, Env.df[, c('SampleID', 'Patch.dwt')], 
                        all.x = TRUE, all.y = FALSE 
)

## divide all species count columns by Patch.dwt
SECC.fauna.g <- SECC.fauna.dwt
for (i in which( sapply(SECC.fauna.g, is.numeric) ) ) {
  if (colnames(SECC.fauna.g)[i] != "Patch.dwt")
    SECC.fauna.g[[i]] <- SECC.fauna.dwt[[i]] / SECC.fauna.dwt$Patch.dwt
}

SECC.fauna <- SECC.fauna.g


##================================================
## Summary variables for univariate analyses
## & merging into main SECC dataframe
##================================================
## Could have been done partly using aggregate() if the species were rows instead of columns :P
SECC.fauna.sum <- SECC.fauna[, c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos')]
Taxa.groups <- unique(SECC.fauna.meta$Taxonomic.Group)

SECC.fauna.sum <- within(SECC.fauna.sum, {
### Mesostigs
### Collembola
### Prostigs
### Other
    for (taxa in Taxa.groups) {
      assign( taxa, apply(SECC.fauna[, SECC.fauna.meta$ID[which(SECC.fauna.meta$Taxonomic.Group == taxa)] ], 1, sum) )
    }
    rm(taxa)

### Uropodina (non-predatory Mesostigs)
    Uropodina <- apply(SECC.fauna[, SECC.fauna.meta$ID[which(SECC.fauna.meta$Major.Taxa == "  Uropodina")] ], 1, sum)
### Mesostig.preds = Mesostigs - Uropodina
    Mesostig.preds <- Mesostigmata - Uropodina
### Predators = Mesostigs.preds + Prostigs
    Predators <- Mesostig.preds + Prostigmata
### Grazers = Collembola + Uropodina (+ Oribatids)
    Grazers <- Collembola + Uropodina
### fauna.jaw = Mesostigs + Collembola (+ Prostigs?)
    fauna.jaw <- Mesostigmata + Collembola
### fauna = all (-Other) * including ZL data?
    fauna <- apply(SECC.fauna[, SECC.fauna.meta$ID], 1, sum)
})
SECC.fauna.sum <- SECC.fauna.sum[, c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos',
                              'Mesostigmata', 'Collembola', 'Prostigmata', 'Other',
                              'Uropodina', 'Mesostig.preds', 'Predators', 'Grazers', 
                              'fauna.jaw', 'fauna')]  # manually reorder columns 


##================================================
## Species Richness & Diversity metrics
##================================================



##################################################
## FINAL CLEANUP 
##################################################
rownames(SECC.fauna) <- SECC.fauna$SampleID
## recode factors with meaningful labels
## SECC.fauna <- recodeSECC( SECC.fauna )  # standard function in `./lib/fun.r`  

## SECC.fauna still includes many extra columns with sample meta-data not needed for ordination analyses.
## stripping these should be fairly easy: keep only the first column (SampleIDs - or rely on these as rownames), and all numeric columns (species counts).
## - remove Patch.dwt (numeric) column.


## House-keeping
rm(list=c('coli', 'rows.0', 'cols.0', 'i', 'IDs.formal', 'Mdata.cols', 'sample.IDs', 'SECC.fauna.t', 'SECC.fauna.df', 'SECC.fauna.raw', 'SECC.fauna.dwt', 'SECC.fauna.g', 'sp.rows'))
