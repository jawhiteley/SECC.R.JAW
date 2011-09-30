##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process Cyanobacteria data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-03-23
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
  rm(list=c('Duplicates', 'Missing', 'text.empty', 'Pblm.msg'))
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
    # t4:  Total Dry Weight of moss sample in tube
    #      = Total weight - thread weight (if present)
    #      thread weight is ignored (treated as 0) if missing (NA).
    #      This is currently calculated in Excel.  R would replace values with 'NA'
    #      If there was any NA involved in the calculation :-(
  ARA.dwt   <- Total.Dry.wt  # Dry weight (mg) of tube moss sample used for ARA (N-fixation)
  Cells.dwt <- CB.Dry.wt     # Dry weight (mg) of 2 moss shoots used for cyanobacteria counts.
  # Discard data from second count (keep count 1).  I need values in both, for the calculations!
  ARA.dwt[Count==2]   <- ARA.dwt[Count==1]  
  Cells.dwt[Count==2] <- Cells.dwt[Count==1]
})

##################################################
## CALCULATIONS
##################################################
sampleA <- 6	# sample Area, in cm^2
      #    6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      # or 6.4 for 20 shoots, based on density survey.
sample.m2 <- sampleA/(100*100)  # sample area in (cm^2 to) m^2
sample.ul <- 1000	            # volume of water added to sample: 1 ml in ul.
sample.sq <- sample.ul / 0.1    # 1 lg hemacytometer square = 0.1 ul.

SECC.cyanobacteria <- within( SECC.cyanobacteria, {
  Calothrix   <- Calothrix.V + Calothrix.H
  Nostoc      <- Nostoc.V + Nostoc.H
  Stigonema   <- Stigonema.V + Stigonema.H
  Other.cells <- single.cells..spp.. + Other
  Total.cells <- Calothrix + Nostoc + Stigonema + Other.cells # Total cell counts
  Total.h     <- Calothrix.H + Nostoc.H + Stigonema.H         # Heterocysts
  Cells   <- (Total.cells / X..lg.squares) * sample.sq / 2
          # Total Count / # hemacytometer squares * squares/sample
          # / 2 shoots/sample = Cells/shoot
  Cells.g <- Cells * 2 / (Cells.dwt / 1000)
          # Cells/shoot * 2 shoots/sample
          # / ( Dry Weight mg / 1000 mg/g ) -> cells / g
  Cells.m <- (Cells.g / 1000) * ARA.dwt / sample.m2
          # scale cells/mg to cells/ARA sample -> scale ARA sample up to m^-2
          # OR: scale to Total Patch Dry Weight (available in SECC, not here) before scaling to m^-2 based on patch Area.
          # could also take # cells /2 shoots and scale up to 20 in sample (*20)?
  Hcells   <- (Total.h / X..lg.squares) * sample.sq / 2
  Hcells.g <- Hcells * 2 / (Cells.dwt / 1000)
  Hcells.m <- (Hcells.g / 1000) * ARA.dwt / sample.m2
})

Cyanobacteria.full <- SECC.cyanobacteria  # save a copy, just in case.
SampleID.date      <- SECC.cyanobacteria[SampleID.first, c('SampleID', 'Date')]  # first date
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
## Values are averaged together (FUN = 'mean'): mean(f(x)) != f(mean(x))
## This also means that, if values are the same for counts 1 & 2, the net result will be the same ( mean(2x) == x ).
## NAs in data lead to NAs in the result, rather than being ignored -- without 'na.rm = TRUE' **

## Merge in single values for Dates & Counts, which can not be aggregated as numerics.
 # Although I find merge() a little confusing,
 # I prefer matching by SampleID rather than relying on a common sort order.
SECC.cyanobacteria <- merge( SECC.cyanobacteria, SampleID.date,   by = c("SampleID"), all.x = TRUE )
SECC.cyanobacteria <- merge( SECC.cyanobacteria, SampleID.counts, by = c("SampleID"), all.x = TRUE )

## aggregate **re-sorts** results
# Re-producible sort
SECC.cyanobacteria  <- sort_df( SECC.cyanobacteria, vars=c('SampleID') )
SampleID.date       <- sort_df( SampleID.date,      vars=c('SampleID') )  # deprecated
SampleID.counts     <- sort_df( SampleID.counts,    vars=c('SampleID') )  # deprecated


##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.cyanobacteria, "SECC columns") <- 
    c("ARA.dwt", "Cells.dwt",
      "Nostoc", "Nostoc.H", "Stigonema", "Stigonema.H", "Other.cells", 
      "Cells",  "Cells.g",  "Cells.m", "Hcells", "Hcells.g", "Hcells.m"
      )
attr(SECC.cyanobacteria, "labels") <- 
    list("ARA.dwt"      = "Dry Weight",
         "Cells.dwt"    = "Dry Weight",
         "Cells"    = "Cyanobacteria Cell Density",
         "Cells.g"  = "Cyanobacteria Cell Density",
         "Cells.m"  = "Cyanobacteria Cell Density",
         "Hcells"   = "Heterocyst Cell Density", 
         "Hcells.g" = "Heterocyst Cell Density", 
         "Hcells.m" = "Heterocyst Cell Density" 
         )
attr(SECC.cyanobacteria, "units")  <- 
    list("ARA.dwt"      = "mg",
         "Cells.dwt"    = "mg",
         "Cells"    = "cells/shoot",
         "Cells.g"  = quote("cells" %.% g^-1 * "dwt"),
         "Cells.m"  = quote("cells" %.% m^-2),
         "Hcells"   = "cells/shoot",
         "Hcells.g" = quote("cells" %.% g^-1 * "dwt"),
         "Hcells.m" = quote("cells" %.% m^-2)
         )


##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d
  head(SECC.cyanobacteria)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.cyanobacteria)   # check structure: are the appropriate variables factors, numeric, etc.?

  # Check for rows where dry weight data does not match between counts
  wt.cols <- c("SampleID", "Total.Dry.wt", "CB.Dry.wt")
  CB.c1   <- Cyanobacteria.full[Cyanobacteria.full$Count == 1, wt.cols]
  CB.c2   <- Cyanobacteria.full[Cyanobacteria.full$Count == 2, wt.cols]
  CB.wts  <- merge( CB.c1, CB.c2, by = c("SampleID"), all = TRUE, suffixes = c(".1", ".2") )
  CB.wts[CB.wts$Total.Dry.wt.1 != CB.wts$Total.Dry.wt.2, ]
  CB.wts[CB.wts$CB.Dry.wt.1    != CB.wts$CB.Dry.wt.2   , ]
  # could also just check for rows in SECC.cyanobacteria where Total.Dry.wt != ARA.dwt, etc.
  SECC.cyanobacteria[SECC.cyanobacteria$Total.Dry.wt != SECC.cyanobacteria$ARA.dwt,   ]
  SECC.cyanobacteria[SECC.cyanobacteria$CB.Dry.wt    != SECC.cyanobacteria$Cells.dwt, ]


  ## Check data structure of main data for analysis: I, O; A, C
  SECCstr(SECC.cyanobacteria[SECC.cyanobacteria$Chamber %in% c("A", "C") &
                             SECC.cyanobacteria$Pos     %in% c("I", "O"),
                             ])

  hist(SECC.cyanobacteria[SECC.cyanobacteria$Time == 1, "ARA.dwt"])
  hist(SECC.cyanobacteria[SECC.cyanobacteria$Time == 2, "ARA.dwt"])
  hist(SECC.cyanobacteria[SECC.cyanobacteria$Time == 4, "ARA.dwt"])
  hist(SECC.cyanobacteria[SECC.cyanobacteria$Time == 1, "Cells.dwt"])
  hist(SECC.cyanobacteria[SECC.cyanobacteria$Time == 2, "Cells.dwt"])
  hist(SECC.cyanobacteria[SECC.cyanobacteria$Time == 4, "Cells.dwt"])

  ## Are calculated values in the right ballpark?
  ## Early results (t1):     0 ~ 5.6e+07 cells m^-2  (56,000,000) mean
  ## DeLuca et al. (2007):   0 ~ 80,000  cells / shoot 
  ## * 20 shoots / 6.4 cm^2    = 2.5e+09 cells m^-2  (2,500,000,000)
  hist(       SECC.cyanobacteria[, "Cells"])  # Cells / shoot 
  boxplot(    SECC.cyanobacteria[, "Cells"])
  boxplot(log(SECC.cyanobacteria[, "Cells"] +1))
  boxplot( list("t1"=log(SECC.cyanobacteria[SECC.cyanobacteria$Time == 1, "Cells"] +1),
                "t2"=log(SECC.cyanobacteria[SECC.cyanobacteria$Time == 2, "Cells"] +1),
                "t4"=log(SECC.cyanobacteria[SECC.cyanobacteria$Time == 4, "Cells"] +1)
                )
          )
  mean(SECC.cyanobacteria[, "Cells"], na.rm = TRUE) 

  hist(       SECC.cyanobacteria[, "Cells.m"])
  boxplot(    SECC.cyanobacteria[, "Cells.m"])
  boxplot(log(SECC.cyanobacteria[, "Cells.m"] +1))
  boxplot( list("t1"=log(SECC.cyanobacteria[SECC.cyanobacteria$Time==1, "Cells.m"] +1),
                "t2"=log(SECC.cyanobacteria[SECC.cyanobacteria$Time==2, "Cells.m"] +1),
                "t4"=log(SECC.cyanobacteria[SECC.cyanobacteria$Time==4, "Cells.m"] +1)
                )
          )
  mean(SECC.cyanobacteria[, "Cells.m"], na.rm = TRUE) 

  ## Plots: use identify() to id rows of interest.
  t1.dwt <- SECC.cyanobacteria[SECC.cyanobacteria$Time == 1, c("Cells.dwt", "ARA.dwt")]
  plot(t1.dwt)
  abline(0, 10, lty = 3, col="#444444")
  rows.id <- identify(t1.dwt)
  SECC.cyanobacteria[as.numeric(row.names(t1.dwt)[rows.id]), ]

  t2.dwt <- SECC.cyanobacteria[SECC.cyanobacteria$Time == 2, c("Cells.dwt", "ARA.dwt")]
  plot(t2.dwt)
  abline(0, 10, lty = 3, col="#444444")
  rows.id <- identify(t2.dwt)
  SECC.cyanobacteria[as.numeric(row.names(t2.dwt)[rows.id]), ]

  t4.dwt <- SECC.cyanobacteria[SECC.cyanobacteria$Time == 4, c("Cells.dwt", "ARA.dwt")]
  plot(t4.dwt)
  abline(0, 10, lty = 3, col="#444444")
  rows.id <- identify(t4.dwt)
  SECC.cyanobacteria[as.numeric(row.names(t4.dwt)[rows.id]), ]

  t1.cells <- SECC.cyanobacteria[SECC.cyanobacteria$Time == 1, c("Cells", "H.cells")]
  plot(t1.cells)
  rows.id <- identify(t1.cells)
  SECC.cyanobacteria[as.numeric(row.names(t1.cells)[rows.id]), ]

  t2.cells <- SECC.cyanobacteria[SECC.cyanobacteria$Time == 2, c("Cells", "H.cells")]
  plot(t2.cells)
  rows.id <- identify(t2.cells)
  SECC.cyanobacteria[as.numeric(row.names(t2.cells)[rows.id]), ]

  t4.cells <- SECC.cyanobacteria[SECC.cyanobacteria$Time == 4, c("Cells", "H.cells")]
  plot(t4.cells)
  rows.id <- identify(t4.cells)
  SECC.cyanobacteria[as.numeric(row.names(t4.cells)[rows.id]), ]

### Potential Outliers:
  ## 31C-3.N  too many Heterocysts?
  ## 11A-1.E  too few  Heterocysts?
  ## 24C-2.O  no heterocysts?
  ## 84C-2.O  too few  Heterocysts?
}

##################################################
## SAVE DATA
##################################################
# leave in memory

##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SampleID.length', 'SampleID.unique', 'SampleID.first',
                'SampleID.counts', 'SampleID.wrong', 'SampleID.date',
                'sample.ul', 'sample.sq', 'col.types', 'aggregate.columns')
rm(list=rm.objects)
## Update list of Data_objects for importing
# Data_objects <- c( Data_objects[Data_objects!=rm.objects] , 'SECC.cyanobacteria' )

