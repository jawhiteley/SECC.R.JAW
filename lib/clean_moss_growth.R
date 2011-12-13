##################################################
## Schefferville Experiment on Climate Change (SEC-C)
## Load, clean & process Moss Growth data loaded from "./data/"
## Jonathan Whiteley		R v2.12		2011-12-12
##################################################
## This script is run as part of `./lib/load.R`

##================================================
## CHECK DATA
##================================================
SECC.moss.raw <- SECC.moss.growth
if (FALSE){
  SECC.moss.growth <- SECC.moss.raw # in case recovery is needed
}

## Standardize ID column names, but not values
SECC.moss.growth <- checkSECCdata(SECC.moss.growth, 'SECC.moss.growth', CheckValues = FALSE)

##################################################
## CALCULATIONS
##################################################
## raw data is in mm of vertical extension from previous time point.
## incorporate biomass data to convert mm to productivity (g dwt)?
## * as sample-specific as possible, with available data.
## + Time 2,4 ; Block 3, 5, 6, 7, 8 ; Chamber A, C ; Frag 1, 4 ; Position I, O
## + Some preliminary data from: Time 4, Chamber B ; Block 1, 2, 3, 6, 7, 8 ; Frag 1, 2, 3, 4 ; Position I, O, S
SECC.moss.growth <- within( SECC.moss.growth, 
                 {
                   grow01[grow01=="."] <- NA
                   grow01 <- as.numeric(as.character(grow01))
                   grow12[grow12=="."] <- NA
                   grow12 <- as.numeric(as.character(grow12))
                   grow23[grow23=="."] <- NA
                   grow23 <- as.numeric(as.character(grow23))
                 })


##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d

  head(SECC.moss.growth)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.moss.growth)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.moss.growth)

###===============================================
### CHECK Calculations
  ## Histograms
  hist(SECC.moss.growth$grow01)
  hist(SECC.moss.growth$grow12)
  hist(SECC.moss.growth$grow23)
  boxplot(SECC.moss.growth[, c("grow01", "grow12", "grow23")])

  summary(SECC.moss.growth, na.rm = TRUE)
}

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.

SECC.moss.growth <- within( SECC.moss.growth, {
})


##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.moss.growth, "SECC columns") <- c('grow01', 'grow12', 'grow23')
attr(SECC.moss.growth, "labels") <- list("grow01"     = "Moss growth  0-12 months",
                                         "grow12"     = "Moss growth 12-22 months",
                                         "grow23"     = "Moss growth 22-24 months"
                                         )
attr(SECC.moss.growth, "units")  <- list("grow01"     = "mm",
                                         "grow12"     = "mm",
                                         "grow23"     = "mm"
                                         )

##################################################
## SAVE DATA
##################################################
# leave in memory
SECC.moss.full <- SECC.moss.growth
# Keep rows for experiment samples, drop other samples (for merging)
CorridorRows <- grep("c.*", SECC.moss.growth$Pos)
ExternalRows <- grep("[^1-4]", SECC.moss.growth$Frag)
SECC.moss.growth <- SECC.moss.growth[-c(CorridorRows, ExternalRows), ]
SECC.moss.growth$SampleID <- SECC_sampleID(SECC.moss.growth) # important for merging!


##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.moss.raw')
rm(list=rm.objects)
## Update list of Data_objects for importing
## Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.moss.growth' )
## Data_objects <- Data_objects[!(Data_objects %in% rm.objects)]

