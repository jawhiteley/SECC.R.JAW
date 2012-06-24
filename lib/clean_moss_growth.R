##################################################
## Schefferville Experiment on Climate Change (SEC-C)
## Load, clean & process Moss Growth data loaded from "./data/"
## Jonathan Whiteley		R v2.12		2011-12-14
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

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.

SECC.moss.growth <- within( SECC.moss.growth, {
						 ## manually fix ID codes
						 SampleID <- gsub("\\.0", ".O", SampleID)
						 SampleID <- gsub("\\.1", ".I", SampleID)
						 ## fix growth data
						 grow01[grow01=="."] <- NA
						 grow01 <- as.numeric(as.character(grow01))
						 grow12[grow12=="."] <- NA
						 grow12 <- as.numeric(as.character(grow12))
						 grow23[grow23=="."] <- NA
						 grow23 <- as.numeric(as.character(grow23))
})


##################################################
## CALCULATIONS
##################################################
## raw data is in mm of vertical extension from previous time point.
## incorporate biomass data to convert mm to productivity (g dwt)?
## * as sample-specific as possible, with available data.
## + Time 2,4 ; Block 3, 5, 6, 7, 8 ; Chamber A, C ; Frag 1, 4 ; Position I, O
## - Some preliminary data from: Time 4, Chamber B ; Block 1, 2, 3, 6, 7, 8 ; Frag 1, 2, 3, 4 ; Position I, O, S
load(sprintf("%sMoss_mg.R", SaveDir.obj()))
SECC.moss <- merge(SECC.moss.growth, subset(Moss.mg, select=c("SampleID", "mg.cm")), 
				   by="SampleID", all.x=TRUE, all.y=FALSE)
## Locate shoots with misleading growth values: skinny / thin / wiry, or growth > 20 mm
## It's not very many, but helps reduce outliers :D
Skinny01 <- union( grep("skinny|thin|wiry", SECC.moss$notes01), which(SECC.moss$grow01 > 20) ) 
Skinny12 <- union( grep("skinny|thin|wiry", SECC.moss$notes12), which(SECC.moss$grow12 > 20) )
Skinny23 <- union( grep("skinny|thin|wiry", SECC.moss$notes23), which(SECC.moss$grow23 > 20) )

SECC.moss <- within( SECC.moss, 
					{
					  if (F) {  ## convert growth in mm to growth in cm
						grow01 <- grow01 / 10           # growth in mm /10 mm/cm -> cm
						grow12 <- grow12 / 10           # growth in mm /10 mm/cm -> cm
						grow23 <- grow23 / 10           # growth in mm /10 mm/cm -> cm
					  }
					  mg.mm  <- mg.cm  / 10           # mg/cm / 10 mm/cm -> mg/mm
					  ## Calculate Biomass Production in mg / shoot
					  Prod23 <- grow23 * mg.mm
					  Prod12 <- grow12 * mg.mm
					  Prod01 <- grow01 * mg.mm        # growth in mm * mg/mm -> mg
					  ## correct for 'skinny' shoots by using 1/2 estimated mg/cm
					  Prod01[Skinny01] <- (grow01 * mg.mm/2)[Skinny01]
					  Prod12[Skinny12] <- (grow12 * mg.mm/2)[Skinny12]
					  Prod23[Skinny23] <- (grow23 * mg.mm/2)[Skinny23]
                      ## convert Productivity to mg / patch (or square meter)
					})
SECC.moss <- SECC.moss[, c("SampleID", "Block", "Time", "Chamber", "Frag", "Pos",
						   'grow01', 'Prod01', 'notes01', 
						   'grow12', 'Prod12', 'notes12', 
						   'grow23', 'Prod23', 'notes23',
						   'mg.cm', 'mg.mm')]


##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d

  head(SECC.moss)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.moss)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.moss)

###===============================================
### CHECK Calculations
  ## Histograms
  op <- par(mfrow=c(3,3))
  hist(SECC.moss$grow01)
  hist(SECC.moss$grow12)
  hist(SECC.moss$grow23)
  hist(SECC.moss$Prod01)
  hist(SECC.moss$Prod12)
  hist(SECC.moss$Prod23)
  boxplot(SECC.moss[, c("grow01", "grow12", "grow23")])
  boxplot(SECC.moss[, c("Prod01", "Prod12", "Prod23")])
  par(op)

  summary(SECC.moss, na.rm = TRUE)
}

##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.moss, "SECC columns") <- c('grow01', 'grow12', 'grow23',  # 'mg.cm', 'mg.mm'
									 'Prod01', 'Prod12', 'Prod23')
attr(SECC.moss, "labels") <- list("Grow"       = "Moss growth",
								  "grow01"     = "Moss growth  0-12 months",
								  "grow12"     = "Moss growth 12-22 months",
								  "grow23"     = "Moss growth 22-24 months",
								  "Prod"       = "Moss Biomass Production",
								  "Prod01"     = "Moss Productivity  0-12 months",
								  "Prod12"     = "Moss Productivity 12-22 months",
								  "Prod23"     = "Moss Productivity 22-24 months"
								  )
attr(SECC.moss, "units")  <- list("Grow"       = "mm",
								  "grow01"     = "mm",
								  "grow12"     = "mm",
								  "grow23"     = "mm",
								  "Prod"       = quote("mg" %.% "shoot"^-1),
								  "Prod01"     = quote("mg" %.% "shoot"^-1),
								  "Prod12"     = quote("mg" %.% "shoot"^-1),
								  "Prod23"     = quote("mg" %.% "shoot"^-1)
								  )

##################################################
## SAVE DATA
##################################################
# leave in memory
SECC.moss.full <- SECC.moss
# Keep rows for experiment samples, drop other samples: for merging
CorridorRows <- grep("c.*", SECC.moss$Pos)
ExternalRows <- grep("[^1-4]", SECC.moss$Frag)
SECC.moss <- SECC.moss[-c(CorridorRows, ExternalRows), ]
SECC.moss$SampleID <- SECC_sampleID(SECC.moss) # important for merging?


##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.moss.growth', 'SECC.moss.raw')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.moss' )

