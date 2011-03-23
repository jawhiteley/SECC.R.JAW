##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process ARA data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-03-15
##################################################
## This script is run as part of `./lib/load.R`
require(car)		# load external package 'car', for recode()

##================================================
## CHECK DATA
##================================================
# str(SECC.cyanobacteria)  # Should already be in memory.  str() produces output on source().

# Standardize ID column names (but not values)
SECC.cyanobacteria <- SECCcolumnNames(SECC.cyanobacteria)
# Strip empty rows (before checking values).
SECC.cyanobacteria <- SECC.cyanobacteria[!is.na(SECC.cyanobacteria$Time),]
# Standardize column types (needed for aggregate & other processing based on column types).
SECC.cyanobacteria <- checkSECCdata(SECC.cyanobacteria, "SECC.cyanobacteria", CheckDuplicates = FALSE)

SampleID.length <- length(SECC.cyanobacteria$SampleID)
SampleID.unique <- length(unique(SECC.cyanobacteria$SampleID))
SampleID.first  <- match(unique(SECC.cyanobacteria$SampleID),
                        SECC.cyanobacteria$SampleID)
SampleID.counts <- table(SECC.cyanobacteria$SampleID)  # Counts of each 'factor level'
SampleID.wrong  <- SampleID.counts[SampleID.counts!=2] # IDs with counts other than 2
if( length(SampleID.wrong)>0 ) {
  Duplicates <- paste(names(SampleID.counts[SampleID.counts>2]), collapse=", ")
  Missing    <- paste(names(SampleID.counts[SampleID.counts<2]), collapse=", ")
  text.empty <- "<none>"
  Pblm.msg <- paste("\nDuplicates: ",
                    ifelse(Duplicates=="", text.empty, Duplicates),
                    "\nMissing:    ",
                    ifelse(Missing==""   , text.empty, Missing   ),
                    sep = ''
                    )
  ## This is downgraded to a warning, because the calculations can handle it:
  ##   they are all averaged together, and I'm not sure it's that big a deal.
  ## Ideally, I'd have an explicit way of including all replicates,
  ##   or keeping only the first 2, but I haven't figured out an easy / elegant way to do this yet.
  warning( paste("There are not exactly 2 subsamples of each Cyanobacteria sample.",
                 Pblm.msg, sep='' ) )
}

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
# Manually clean & prepare data for automatic checking.
names(SECC.cyanobacteria)[names(SECC.cyanobacteria)=="Count.."] <- "Count"  # rename column
SECC.cyanobacteria <- within( SECC.cyanobacteria, {
  # Count should really be a factor
  Count <- factor(Count, levels = c(1, 2) )
  ## Prepare for pooling
  ARA.dwt   <- Total.Dry.wt  # Dry weight (mg) of tube moss sample used for ARA (N-fixation)
  Cells.dwt <- CB.Dry.wt     # Dry weight (mg) of 2 moss shoots used for cyanobacteria counts.
  # Discard data from second count (keep count 1).  I need values in both, for the calculations!
  ARA.dwt[Count==2]   <- ARA.dwt[Count==1]  
  Cells.dwt[Count==2] <- Cells.dwt[Count==1]
})


##################################################
## CALCULATIONS
##################################################
sampleA <- 6	# sample Area, in cm^2: 6 for rough estimate of inner tube diameter (as used in ARA excel file), or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100/sampleA)  # scale sample area, in cm^2 to m^2
sample.ul <- 1000	          # volume of water added to sample, in ul.
sample.sq <- sample.ul / 0.1  # 1 lg hemacytometer square = 0.1 ul.

SECC.cyanobacteria <- within( SECC.cyanobacteria, {
  Calothrix   <- Calothrix.V + Calothrix.H
  Nostoc      <- Nostoc.V + Nostoc.H
  Stigonema   <- Stigonema.V + Stigonema.H
  Other.cells <- single.cells..spp.. + Other
  Total.cells <- Calothrix + Nostoc + Stigonema + Other.cells # Total cell counts
  Total.h     <- Calothrix.H + Nostoc.H + Stigonema.H         # Heterocysts
  Cells   <- (Total.cells / X..lg.squares) * sample.sq / Cells.dwt * 1000
          # Total Count / # hemacytometer squares * squares/sample
          # / Dry Weight (mg) * 1000 mg/g -> cells / g
  Cells.m <- Cells * ARA.dwt/1000 * sample.to.m2
          # scale cells/mg to cells/ARA sample -> scale ARA sample up to m^2
  H.cells <- (Total.h / X..lg.squares) * sample.sq / Cells.dwt * 1000
})

Cyanobacteria.full <- SECC.cyanobacteria  # save a copy, just in case.
SampleID.date      <- Cyanobacteria.full[SampleID.first, c('SampleID', 'Date')]  # first date
SampleID.counts    <- as.data.frame(SampleID.counts)  # counts per unique SampleID
names(SampleID.counts) <- c("SampleID", "Counts")
## Dry weight (mg) of tube moss sample used for ARA (N-fixation) (assuming same sort order)
# ARA.dwt     <- Cyanobacteria.full[SampleID.first, "Total.Dry.wt"]  # first occurence
## Dry weight (mg) of 2 moss shoots used for cyanobacteria counts. (assuming same sort order)
# Cells.dwt   <- Cyanobacteria.full[SampleID.first, "CB.Dry.wt"]     # first occurence


## Pool patch data across subsamples.
col.types <- lapply(SECC.cyanobacteria, class)
aggregate.columns <- names(SECC.cyanobacteria)[(col.types!="factor" & col.types!="character")]  # a single '&' performs element-wise 'AND' **
## aggregate data: means after calcs might be more appropriate than pooling by sum before calcs?
SECC.cyanobacteria <- with(SECC.cyanobacteria,
                           aggregate(SECC.cyanobacteria[, aggregate.columns], 
                                     by = list( SampleID = SampleID,
                                       Block   = Block,
                                       Time    = Time,
                                       Chamber = Chamber,
                                       Frag    = Frag,
                                       Pos     = Pos
                                       ),
                                     FUN = 'mean', na.rm = TRUE
                                     )
                           )  # 2 counts / sample
## NAs in data lead to NAs in the result, rather than being omitted -- without 'na.rm = TRUE' **

## aggregate **re-sorts** results
# Re-producible sort
SECC.cyanobacteria  <- sort_df( SECC.cyanobacteria,  vars=c('SampleID') )
SampleID.date       <- sort_df( SampleID.date, vars=c('SampleID') )
SampleID.counts     <- sort_df( SampleID.counts,     vars=c('SampleID') )

## Additional calculations & clean-up
SECC.cyanobacteria <- within( SECC.cyanobacteria, {
  Counts <- SampleID.counts$Counts  # Counts for each SampleID
  Date   <- SampleID.date$Date    # First Date in data frame
})


##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.cyanobacteria, "SECC columns") <- c("ARA.dwt", "Cells.dwt",
                                              "Nostoc", "Nostoc.H",
                                              "Stigonema", "Stigonema.H",
                                              "Other.cells", "Cells", "H.cells")
attr(SECC.cyanobacteria, "labels") <- list("Cells"="Cyanobacteria Cell Density")
attr(SECC.cyanobacteria, "units")  <- list("Cells"="cells/g dwt")


##################################################
## CHECK DATA
##################################################

head(SECC.cyanobacteria)  # have a peek at the first 6 rows & columns: is this what you expected?
# str(SECC.cyanobacteria)   # check structure: are the appropriate variables factors, numeric, etc.?
## Check data structure of main data for analysis: I, O; A, C
SECCstr(SECC.cyanobacteria[SECC.cyanobacteria$Chamber %in% c("A", "C") &
                           SECC.cyanobacteria$Pos     %in% c("I", "O"),
                           ])

##################################################
## SAVE DATA
##################################################
# leave in memory

##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SampleID.length', 'SampleID.unique', 'SampleID.first', 'SampleID.counts', 'SampleID.wrong', 'SampleID.date')
rm(list=rm.objects)
## Update list of Data_objects for importing
# Data_objects <- c( Data_objects[Data_objects!=rm.objects] , 'SECC.cyanobacteria' )

