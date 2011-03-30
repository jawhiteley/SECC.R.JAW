##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process Moisture data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-03-29
##################################################
## This script is run as part of `./lib/load.R`

##================================================
## CHECK DATA
##================================================
# str(SECC.H2O.t1)  # Should already be in memory.    str() produces output on source().

## Standardize ID column names & values
SECC.H2O.t1 <- checkSECCdata(SECC.H2O.t1, 'SECC.H2O.t1')
SECC.H2O.t2 <- checkSECCdata(SECC.H2O.t2, 'SECC.H2O.t2')
SECC.H2O.t4 <- checkSECCdata(SECC.H2O.t4, 'SECC.H2O.t4')

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.
## All Weights in grams (g)

###  T1 H2O calculations
SECC.H2O.t1 <- within( SECC.H2O.t1, {
  Total.WET <- tin.moss_wwt
  Total.DRY <- tin.moss_dwt
})

###  T2 H2O calculations
SECC.H2O.t2 <- within( SECC.H2O.t2, {
  Total.WET <- tin.moss_wwt
  Total.DRY <- tin.moss_dwt
})


##################################################
## MERGE SEPARATE TIME
##################################################
## merge(), or rbind() on common columns?

SECC.H2O <- SECC.H2O.t1
SECC.H2O <- merge( SECC.H2O, SECC.H2O.t2, all=TRUE, suffixes = c(".t1", ".t2") )
SECC.H2O <- merge( SECC.H2O, SECC.H2O.t4, all=TRUE, suffixes = c(".12", ".t4"), sort = TRUE )

##================================================
## Check Results
nrow(SECC.H2O)  # 1152
nrow(SECC.H2O.t1) + nrow(SECC.H2O.t2) + nrow(SECC.H2O.t4)  # 1152



##################################################
## CALCULATIONS
##################################################

### H2O calculations
SECC.H2O <- within( SECC.H2O, {
  SampleID  <- SECC_sampleID( SECC.H2O )
  wwt <- Total.WET - tin
  dwt <- Total.DRY - tin
  H2O.wt  <- wwt - dwt
  H2O     <- H2O.wt / dwt
  H2O.wwt <- H2O.wt / wwt
})



##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.H2O, "SECC columns") <- c('H2O', 'H2O.wwt')
attr(SECC.H2O, "labels") <- list("H2O"     = "Water Content",
                                 "H2O.wwt" = "% Water"
                                )
attr(SECC.H2O, "units")  <- list("H2O"     = quote("% moss dry weight"),
                                 "H2O.wwt" = quote("% moss wet weight")
                                 )

##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d
  head(SECC.H2O)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.H2O)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.H2O)

###===============================================
### CHECK Calculations
  # Histograms
  hist(SECC.H2O[SECC.H2O$Time == 1, "H2O"])
  hist(SECC.H2O[SECC.H2O$Time == 2, "H2O"])
  hist(SECC.H2O[SECC.H2O$Time == 4, "H2O"])

  boxplot( list("t1"=SECC.H2O[SECC.H2O$Time == 1, "H2O"],
                "t2"=SECC.H2O[SECC.H2O$Time == 2, "H2O"],
                "t4"=SECC.H2O[SECC.H2O$Time == 4, "H2O"]
                )
          )
  mean(SECC.H2O[, "H2O"], na.rm = TRUE)

}

##################################################
## SAVE DATA
##################################################
# leave in memory
H2O.full <- SECC.H2O
# only rows for samples, exclude controls
# SECC.H2O <- SECC.H2O[SECC.H2O$SampleControl=="Sample",]

##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.H2O.t1', 'SECC.H2O.t2', 'SECC.H2O.t4')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.H2O' )

