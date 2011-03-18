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

# SECC.cyanobacteria <- SECCcolumnNames(SECC.cyanobacteria)  # Standardize ID column names (but not values)
SECC.cyanobacteria <- checkSECCdata(SECC.cyanobacteria)  # Standardize column names & types (needed for aggregate & other processing based on column types).

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
# Manually clean & prepare data for automatic checking.
SECC.cyanobacteria <- SECC.cyanobacteria[!is.na(SECC.cyanobacteria$Time),]  # strip empty rows.
names(SECC.cyanobacteria)[names(SECC.cyanobacteria)=="Count.."] <- "Count"  # rename column
SECC.cyanobacteria <- within( SECC.cyanobacteria, {
  # Count should really be a factor
  Count <- factor(Count, levels = c(1, 2) )
})


##################################################
## CALCULATIONS
##################################################
## Prepare for pooling
SECC.cyanobacteria <- within( SECC.cyanobacteria, {
  ARA.dwt   <- Total.Dry.wt  # Dry weight (mg) of tube moss sample used for ARA (N-fixation)
  Cells.dwt <- Dry.wt        # Dry weight (mg) of 2 moss shoots used for cyanobacteria counts.
  ARA.dwt[Count==2]   <- NA  # discard data from second count (keep count 1)
  Cells.dwt[Count==2] <- NA  # discard data from second count (keep count 1)
})
Cyanobacteria.full <- SECC.cyanobacteria  # save a copy, just in case.
## Pool patch data across subsamples.
col.types <- lapply(SECC.cyanobacteria, class)
aggregate.columns <- names(SECC.cyanobacteria)[(col.types!="factor" & col.types!="character")]  # a single '&' performs element-wise 'AND' **
SECC.cyanobacteria <- with(SECC.cyanobacteria,
                           aggregate(SECC.cyanobacteria[, aggregate.columns], 
                                     by = list( SampleID = SampleID,
                                       Block   = Block,
                                       Time    = Time,
                                       Chamber = Chamber,
                                       Frag    = Frag,
                                       Pos     = Pos
                                       ),
                                     FUN = 'sum', na.rm = TRUE
                                     )
                           )  # 2 counts / sample
## NAs in data lead to NAs in the result, rather than being omitted -- without 'na.rm = TRUE' **
SECC.cyanobacteria[["Date"]] <- Cyanobacteria.full[Cyanobacteria.full$Count==1, "Date"]

## Additional calculations
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

##################################################
## SAVE DATA
##################################################
# leave in memory
## Remove old objects from memory
# rm.objects <- c()
# rm(list=rm.objects)
## Update list of Data_objects for importing
# Data_objects <- c( Data_objects[Data_objects!=rm.objects] , 'SECC.cyanobacteria' )

