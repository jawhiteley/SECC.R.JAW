##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process data files stored in "./data/"
# Jonathan Whiteley		R v2.12		2012-03-18
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
if (FALSE) {       d# do not run automatically
  source("./lib/init.R")   # load basic objects (needed for environmental data)
  rm(list=ls())     # house-keeping
  setwd('./ SECC/') # project directory
  setwd('../')      # relative to this file (\rd in Vim-R)
  getwd()           # check current wd

  ## LOAD LIBRARIES
  source("./lib/fun.R")   # define functions
}

cat('- Loading fauna data.\n')
library(vegan) # for diversity metrics

cleanSpNames <- function(sp.names = NULL, ...)
{
  new.names <- gsub("\\s", "_", sp.names)
  new.names <- make.names(new.names, ...)
  new.names
}


##################################################
## LOAD DATA FILES
##################################################
cat('  - Loading fauna data files.\n')

SECC.fauna.raw <- read.csv("./data/fauna_t4_raw.csv", na.strings = c("NA", "", "."), as.is=TRUE)
Fauna.JAW.raw <- read.csv("./data/fauna_JAW_raw.csv", na.strings = c("NA", "", "."), as.is=TRUE)
Fauna.ZL.raw  <- read.csv("./data/fauna_ZL_raw.csv",  na.strings = c("NA", "", "."), as.is=TRUE)
Fauna.species <- read.csv("./data/fauna_JAW-ZL_species.csv",  
                          na.strings = c("NA", "", "-"), as.is=TRUE)

Fauna.dataframes <- c("Fauna.JAW.raw", "Fauna.ZL.raw")

for (dat in Fauna.dataframes)
{
  SECC.fauna.raw <- get(dat)

##################################################
### CLEAN INPUT
##################################################
  cat('  -', dat, ': Cleaning raw input.\n')
  ## strip rows with empty species names
  SECC.fauna.raw <- SECC.fauna.raw[!is.na(SECC.fauna.raw$ID), ]

  ## Drop extra 'meta-data' columns (this is all in 'Fauna.species' now)
  SECC.fauna <- SECC.fauna.raw[, c(1, grep("^X\\d", colnames(SECC.fauna.raw)) )]

  if (FALSE) 
  {                                    # 'Meta-data' is now in 'Fauna.species'
    ##================================================
    ## Extract meta-data
    ##================================================
    cat('  - Extracting fauna species metadata.\n')
    ## clone first metadata columns to 'SECC.fauna.meta': species records, aliases, taxanomic categories, etc.
    SECC.fauna.meta <- SECC.fauna.raw[, -grep("^X\\d", colnames(SECC.fauna.raw))]
                                        # keep only columns that do not start with "X" and a number (these are sample IDs).
    ## strip empty rows from meta-data
    SECC.fauna.meta <- strip_empty_dims(SECC.fauna.meta, dim=1, cols = 2:ncol(SECC.fauna.meta))

    ## Mdata.cols <- intersect( colnames(SECC.fauna.raw), colnames(SECC.fauna.meta))[-1]
    Mdata.cols <- which( colnames(SECC.fauna.raw) %in% colnames(SECC.fauna.meta))[-1]
    ## strip species metadata columns
    SECC.fauna <- SECC.fauna.raw[, -Mdata.cols]
  }


##################################################
### PROCESS DATA
##################################################
  cat('  -', dat, ': Processing fauna data.\n')
  ## strip empty (NAs or 0s) columns (samples) 
  SECC.fauna <- strip_empty_dims(SECC.fauna, dim=2)
  sp.rows <- (grep("comments", SECC.fauna$ID, ignore.case = TRUE) +1):NROW(SECC.fauna)
  cols.0  <- which( apply( SECC.fauna[sp.rows, ], 2, function(x) all(x==0) ) )
  if (length(cols.0) > 0) SECC.fauna <- SECC.fauna[, -cols.0]

  ## strip 0-only rows (species)
  sample.cols <- 2:ncol(SECC.fauna)
  rows.0 <- which( apply( SECC.fauna[, sample.cols], 1, function(x) all(x==0) ) )
  if (length(rows.0) > 0) SECC.fauna <- SECC.fauna[-rows.0, ]

  ## Make species labels R-friendly (prepare for transposing to column names)
  ## rownames(SECC.fauna) <- cleanVarName(SECC.fauna$ID)
  sp.names   <- cleanSpNames(SECC.fauna$ID, unique = TRUE)
  rownames(SECC.fauna) <- sp.names
  SECC.fauna$ID <- sp.names


  ##================================================
  ## Transpose input data table
  ##================================================
  cat('    -', dat, ': Transposing fauna data.\n')
  sample.IDs <- colnames(SECC.fauna)[-1]
  ## SECC.fauna.mat <- as.matrix(SECC.fauna)
  SECC.fauna.t <- t(SECC.fauna[, -1])         # omit ID column (replaced by rownames)
  SECC.fauna.t <- cbind(SampleID = sample.IDs, SECC.fauna.t)
  rownames(SECC.fauna.t) <- NULL  # remove annoying rownames (stored in sample.IDs)
  ##   colnames(SECC.fauna.t)[-1] <- sp.names
  SECC.fauna.df <- as.data.frame(SECC.fauna.t, stringsAsFactors = FALSE)

  ## Convert character columns to appropriate data type
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

  if(F)
  {                                      ## Same for Meta-data
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

  ##================================================
  ## Check & clean data structure
  ##================================================
  cat('    -', dat, ': Checking & Cleaning fauna data structure.\n')
  ## Check SECC structure & overwrite un-transposed data frame if successful
  SECC.fauna <- checkSECCdata(SECC.fauna.df, "SECC.fauna.df")
  SECC.fauna <- within(SECC.fauna, 
                       {
                         Comments[is.na(Comments)] <- ""
                         ## re-generate Sample IDs
                         SampleID <- SECC_sampleID(SECC.fauna)
                       })

  if (F)
  {
    ## remove missing sp IDs from meta-data
    SECC.fauna.meta <- SECC.fauna.meta[which( SECC.fauna.meta$ID %in% intersect(colnames(SECC.fauna), SECC.fauna.meta$ID) ), ]
    SECC.fauna.meta$sp_alias <- cleanSpNames(SECC.fauna.meta$sp_alias, unique = FALSE)
  }

  assign(gsub("\\.raw", "", dat), SECC.fauna) # result to new object name without ".raw" @ end
}                                      # END Processing individual data frames...

##################################################
## MERGE JAW & ZL data
##################################################
## JAW: Mesostigs & Collembola from Frag 1, 2, & 4; Blocks 1, 3, 5, 7; t4
##  ZL: All microarthropods    from a subset of Blocks 1, 3, 5, 7; t1, t4
cat('  - Merging data sets.\n')

## STANDARDIZE species names using metadata
fauna_sp_alias  <- function (SECC.fauna = NULL, raw.sp.col = "JAW.sp.ID")
{
  col.class <- lapply(SECC.fauna, function(x) class(x) )
  sp.cols   <- which(col.class == "numeric")
  ## convert species names into same friendly format used for column names
  Species.labels <- within( Fauna.species,
                        {
                          sp_alias <- cleanSpNames(sp_alias)
                        })
  Species.labels[[raw.sp.col]] <- cleanSpNames( Species.labels[[raw.sp.col]] )
  Species.raw <- colnames(SECC.fauna)[sp.cols]
  Species.match <- match(Species.raw, Species.labels[[raw.sp.col]])
  Species.new <- Species.labels$sp_alias[ Species.match ]
  colnames(SECC.fauna)[sp.cols] <- Species.new
  ## there are likely now duplicate column names: these need to be aggregated
  SECC.fauna.new <- SECC.fauna[, -sp.cols]
  fcols <- seq_along( colnames(SECC.fauna.new) )
  Species.labels.new <- na.omit( unique(Species.new) ) # species not in the master list are dropped
  for (sp in Species.labels.new)
  {
    sp.col <- which( colnames(SECC.fauna) == sp )
    if (length(sp.col) > 1)
    {
      sp.new <- apply(SECC.fauna[, sp.col], 1, sum)
    } else {
      sp.new <- SECC.fauna[, sp.col]
    }
    SECC.fauna.new <- cbind(SECC.fauna.new, sp.new)
  }
  colnames(SECC.fauna.new)[-fcols] <- Species.labels.new

  SECC.fauna.new
}

Fauna.JAW.std <- fauna_sp_alias(Fauna.JAW, raw = "JAW.sp.ID")
Fauna.ZL.std  <- fauna_sp_alias(Fauna.ZL , raw = "ZL.sp.ID" )

## Prepare combined sample x species matrix
col.class <- lapply(Fauna.JAW.std, function(x) class(x) )
sp.cols   <- which(col.class == "numeric")
fcols     <- seq_along(col.class)[-sp.cols]
xp.cols   <- c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos')
All.spp <- union(colnames(Fauna.JAW.std)[-fcols], colnames(Fauna.ZL.std)[-fcols])
All.spp <- All.spp[order( match(All.spp, cleanSpNames(Fauna.species$sp_alias) ) )]
All.samples <- union(Fauna.JAW.std$SampleID, Fauna.ZL.std$SampleID)

SECC.fauna <- merge( Fauna.JAW.std[, xp.cols], Fauna.ZL.std[, xp.cols], all = TRUE)
SECC.fauna$Time <- factor(SECC.fauna$Time, levels = c(1, 2, 4) )
SECC.fauna <- SECC.fauna[order(SECC.fauna$Time, SECC.fauna$Block, SECC.fauna$Chamber, SECC.fauna$Frag, SECC.fauna$Pos), ]
## SECC.fauna.mat <- matrix(NA, nrow = NROW(SECC.fauna), ncol = length(All.spp))
## SECC.fauna.mat <- as.data.frame(SECC.fauna.mat)
## colnames(SECC.fauna.mat) <- All.spp
## SECC.fauna <- cbind(SECC.fauna, SECC.fauna.mat)

## Keep Max. count in each cell
dframes <- c('Fauna.JAW.std', 'Fauna.ZL.std')
for (sp in All.spp)
{
  sp.counts <- rep(NA, length.out = NROW(SECC.fauna) )
  for (IDsample in SECC.fauna$SampleID)
  {
    Counts <- 0                        # default
    for (datn in dframes)
    {
      dat   <- get(datn)
      IDrow <- which(dat$SampleID == IDsample)
      Counts <- c(Counts, dat[IDrow, sp])
    }
    ##     browser()
    ##     Count <- if (length(Counts) > 0) max(Counts, na.rm = TRUE) else 0
    IDrow <- which(SECC.fauna$SampleID == IDsample)
    sp.counts[IDrow] <- max(Counts, na.rm = TRUE)
  }
  SECC.fauna <- cbind(SECC.fauna, sp.counts)
}
colnames(SECC.fauna) <- c(xp.cols, All.spp)


if(F)
{
  ## clean species names
  IDs.formal <- SECC.fauna$ID
  IDs.formal <- sub( "\\s\\(imm\\)$", "" , IDs.formal)
  IDs.formal <- sub( "\\s\\(beetle\\)", "" , IDs.formal)
  IDs.formal <- sub( "\\s\\(Diptera\\?\\)", "" , IDs.formal)
  IDs.formal <- sub( "beetle", "Coleoptera spp." , IDs.formal)
  SECC.fauna$ID <- IDs.formal
}

##================================================
## Extract meta-data
##================================================
cat('  - Extracting fauna species metadata.\n')
Species.labels <- unique(Fauna.species$sp_alias)
Meta.sp.rows   <- match(Species.labels, Fauna.species$sp_alias)
SECC.fauna.meta <- Fauna.species[Meta.sp.rows, c("sp_alias", "Major.Taxa", "Taxonomic.Group")]
colnames(SECC.fauna.meta) <- gsub( "sp_alias", "Lowest.Level.ID", colnames(SECC.fauna.meta) )
SECC.fauna.meta$ID <- cleanSpNames(SECC.fauna.meta$Lowest.Level.ID) # column name friendly format for lookups with variable names in other data frames
rownames(SECC.fauna.meta) <- SECC.fauna.meta$ID
SECC.fauna.meta <- SECC.fauna.meta[, c("ID", "Lowest.Level.ID", "Major.Taxa", "Taxonomic.Group")]
## Document Trophic Groups
SECC.fauna.meta <- within(SECC.fauna.meta,
                          {
                            Trophic.Group <- NA
                            Trophic.Group <- ifelse(Taxonomic.Group %in% c("Mesostigmata", "Prostigmata"), 
                                                    "Predator", Trophic.Group)
                            Trophic.Group <- ifelse(Taxonomic.Group == "Collembola", "Grazer", Trophic.Group)
                            Trophic.Group <- ifelse(Major.Taxa == "Uropodina", "Grazer", Trophic.Group)
                            Trophic.Group <- ifelse(Taxonomic.Group == "Other", "Other", Trophic.Group)
                          })




##################################################
## CALCULATIONS
##################################################
cat('  - Calculating fauna data.\n')
SECC.fauna.counts <- SECC.fauna  # save count data

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
SECC.fauna <- SECC.fauna[order(SECC.fauna$Time, SECC.fauna$Block, SECC.fauna$Chamber, SECC.fauna$Frag, SECC.fauna$Pos), ]


##================================================
## Summary variables for univariate analyses
## suitable for merging into main SECC dataframe
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
      taxa.cols <- SECC.fauna.meta$ID[which(SECC.fauna.meta$Taxonomic.Group == taxa)]
      taxa.cols <- intersect(colnames(SECC.fauna), taxa.cols)
      assign( taxa, apply(SECC.fauna[, taxa.cols], 1, sum) )
    }

### Uropodina (non-predatory Mesostigs)
    Uropodina <- apply(SECC.fauna[, SECC.fauna.meta$ID[which(SECC.fauna.meta$Major.Taxa == "Uropodina")] ], 1, sum)
### Mesostig.preds = Mesostigs - Uropodina
    Mesostig.preds <- Mesostigmata - Uropodina
### Predators = Mesostigs.preds + Prostigs
    Predators <- Mesostig.preds + Prostigmata
### Grazers = Collembola + Uropodina (+ Oribatids)
    Grazers <- Collembola + Uropodina
### fauna.jaw = Mesostigs + Collembola (+ Prostigs?)
    fauna.jaw <- Mesostigmata + Collembola
### fauna = all (-Other) * including ZL data?
    taxa.cols <- intersect(colnames(SECC.fauna), SECC.fauna.meta$ID)
    taxa.cols <- setdiff(taxa.cols, SECC.fauna.meta$ID[which(SECC.fauna.meta$Taxonomic.Group == "Other")])
    fauna <- apply(SECC.fauna[, taxa.cols], 1, sum)
    rm(taxa, taxa.cols)
})
SECC.fauna.sum <- SECC.fauna.sum[, c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos',
                              'Mesostigmata', 'Collembola', 'Prostigmata', 'Other',
                              'Uropodina', 'Mesostig.preds', 'Predators', 'Grazers', 
                              'fauna.jaw', 'fauna')]  # manually reorder columns 



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
rm(list=c('coli', 'rows.0', 'cols.0', 'i', 'sample.IDs', 'SECC.fauna.t', 'SECC.fauna.df', 'SECC.fauna.raw', 'SECC.fauna.dwt', 'SECC.fauna.g', 'sp', 'sp.names', 'sp.rows', 'Species.labels', 'fcols', 'xp.cols', 'sp.cols', 'sp.counts', 'col.class', 'Counts', 'dat', 'datn', 'dframes', 'Meta.sp.rows', 'Taxa.groups')) # , 'IDs.formal', 'Mdata.cols'
