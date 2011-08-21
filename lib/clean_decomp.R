##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process Decomposition data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-08-21
##################################################
## This script is run as part of `./lib/load.R`

##================================================
## CHECK DATA
##================================================
# str(SECC.Decomposition)  # Should already be in memory.    str() produces output on source().

## Standardize ID column names & values
SECC.Decomposition <- checkSECCdata(SECC.Decomposition, 'SECC.Decomposition')

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.

SECC.Decomposition$SampleID <- SECC.Decomposition$bagID


##################################################
## CALCULATIONS
##################################################

SECC.decomp <- within( SECC.decomp, {
  SampleID  <- SECC_sampleID( SECC.decomp )
})



##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.decomp, "SECC columns") <- c('')
attr(SECC.decomp, "labels") <- list(""     = "",
                                 "" = ""
                                )
attr(SECC.decomp, "units")  <- list(""     = quote(""),
                                 "" = quote("")
                                 )

##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d
  head(SECC.decomp)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.decomp)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.decomp)

###===============================================
### CHECK Calculations
  # Histograms
  hist(SECC.decomp)
  boxplot(SECC.decomp)

  mean(SECC.decomp[, "H2O"], na.rm = TRUE)

}

##################################################
## SAVE DATA
##################################################
# leave in memory
Decomp.full <- SECC.decomp
# only rows for samples, exclude controls
SECC.decomp <- SECC.decomp[SECC.decomp$Type=="Experiment",]

##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.Decomposition')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.decomp' )

