##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process Decomposition data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-08-21
##################################################
## This script is run as part of `./lib/load.R`

##================================================
## CHECK DATA
##================================================
SECC.decomp.raw <- SECC.decomposition
if (FALSE){
  SECC.decomposition <- SECC.decomp.raw # in case recovery is needed
}
## strip empty rows
SECC.decomposition <- strip_empty_dims(SECC.decomposition, 1, col.class = "numeric")

SECC.decomposition <- within( SECC.decomposition, {
  ChamberID <- Chamber
  rm(Chamber)   # prevent duplication: Warming is going to get auto-renamed to this.
  SampleID <- bagID
  Type <- factor(Type)   # drop unused factor levels
})

## Standardize ID column names, but not values
SECC.decomposition <- checkSECCdata(SECC.decomposition, 'SECC.decomposition', CheckValues = FALSE)

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.



##################################################
## CALCULATIONS
##################################################

SECC.decomposition <- within( SECC.decomposition, {
  ##   SampleID <- SECC_sampleID( SECC.decomposition )
  mass.loss <- dwt.initial - dwt.final
  prop.mass.loss <- mass.loss / dwt.initial
})

##================================================
## Preliminary analysis and further calculations
##================================================
if (FALSE) {
  with(SECC.decomposition, boxplot(prop.mass.loss ~ Type) )
  with(SECC.decomposition, 
       plotMeans(prop.mass.loss, Type, error.bars = "conf.int", level=0.95) 
  )
}


##================================================
## Final Calculations
##================================================
## decomp <- prop.mass.loss - mean(prop.mass.loss[Type=="Control 2"])
## "Control" was the original controls I carried out in the field
##    Unfortunately, I think I transported them home dry (why?  I don't know), which
##    artificially inflated the values, making them higher than the experimental values!!
## "Control 2" was a redo of the controls I carried out in the summer of 2011
##    during a trip to Cape Breton: roughly the same amount of flying and time in transit,
##    but transported wet instead.
##    I'm pleased that these new controls turned out to be measureable,
##    but just below nearly all experimental values.
## Use "Control 2" instead of "Control" to correct values for mass lost due to handling.**
decomp.control <- with(SECC.decomposition, mean(prop.mass.loss[Type=="Control 2"]) )
SECC.decomposition <- within(SECC.decomposition, {
    Decomposition <- prop.mass.loss - decomp.control
    Decomposition[Type %in% c('Control', 'Control 2')] <- NA
  })

##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d

  head(SECC.decomposition)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.decomposition)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.decomposition)

###===============================================
### CHECK Calculations
  # Histograms
  hist(SECC.decomposition)
  boxplot(SECC.decomposition)

  mean(SECC.decomposition[, "Decomposition"], na.rm = TRUE)
}

##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.decomposition, "SECC columns") <- c('Decomposition')
attr(SECC.decomposition, "labels") <- list("Decomposition"     = "Decomposition"
                                )
attr(SECC.decomposition, "units")  <- list("Decomposition"     = quote("% mass loss" %.% "yr"^-1)
                                 )

##################################################
## SAVE DATA
##################################################
# leave in memory
Decomp.full <- SECC.decomposition
# only rows for samples, exclude controls
SECC.decomp <- SECC.decomposition[SECC.decomposition$Type=="Experiment",]
SECC.decomp$SampleID <- SECC_sampleID(SECC.decomp) # important for merging!


##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.decomposition', 'SECC.decomp.raw')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.decomp' )

