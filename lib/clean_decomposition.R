##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process Decomposition data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-08-21
##################################################
## This script is run as part of `./lib/load.R`

##================================================
## CHECK DATA
##================================================
SECC.decomp.raw <- SECC.Decomposition
if (FALSE){
  SECC.Decomposition <- SECC.decomp.raw # in case recovery is needed
}
## strip empty rows
SECC.Decomposition <- strip_empty_dims(SECC.Decomposition, 1, col.class = "numeric")

SECC.Decomposition <- within( SECC.Decomposition, {
  ChamberID <- Chamber
  rm(Chamber)   # prevent duplication: Warming is going to get auto-renamed to this.
  SampleID <- bagID
  Type <- factor(Type)   # drop unused factor levels
})

## Standardize ID column names, but not values
SECC.Decomposition <- checkSECCdata(SECC.Decomposition, 'SECC.Decomposition', CheckValues = FALSE)

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.



##################################################
## CALCULATIONS
##################################################

SECC.Decomposition <- within( SECC.Decomposition, {
  ##   SampleID <- SECC_sampleID( SECC.Decomposition )
  mass.loss <- dwt.initial - dwt.final
  prop.mass.loss <- mass.loss / dwt.initial
})

##================================================
## Preliminary analysis and further calculations
##================================================
if (FALSE) {
  with(SECC.Decomposition, boxplot(prop.mass.loss ~ Type) )
  with(SECC.Decomposition, 
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
decomp.control <- with(SECC.Decomposition, mean(prop.mass.loss[Type=="Control 2"]) )
SECC.Decomposition <- within(SECC.Decomposition, {
    Decomposition <- prop.mass.loss - decomp.control
    Decomposition[Type %in% c('Control', 'Control 2')] <- NA
  })

##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d
  head(SECC.Decomposition)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.Decomposition)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.Decomposition)

###===============================================
### CHECK Calculations
  # Histograms
  hist(SECC.Decomposition)
  boxplot(SECC.Decomposition)

  mean(SECC.Decomposition[, "Decomposition"], na.rm = TRUE)
}

##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.Decomposition, "SECC columns") <- c('Decomposition')
attr(SECC.Decomposition, "labels") <- list("Decomposition"     = "Decomposition"
                                )
attr(SECC.Decomposition, "units")  <- list("Decomposition"     = quote("% mass loss")
                                 )

##################################################
## SAVE DATA
##################################################
# leave in memory
Decomp.full <- SECC.Decomposition
# only rows for samples, exclude controls
SECC.decomp <- SECC.Decomposition[SECC.Decomposition$Type=="Experiment",]


##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.Decomposition', 'SECC.decomp.raw')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.decomp' )

