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
# str(SECC.ARA.t1)  # Should already be in memory.    str() produces output on source().

SECC.ARA.t1 <- SECCcolumnNames(SECC.ARA.t1)  # Standardize ID column names (but not values)

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
# Manually clean & prepare data for automatic checking.

SECC.ARA.t1 <- within( SECC.ARA.t1, {
  # Some Fragmentation entries are '-' and need to be recoded as 'NA'
# Frag <- recode(Frag,
#                "'1'='1' ; '2'='2' ; '3'='3' ; '4'='4' ; else=''",
#                as.factor.result = TRUE,
#                levels = c('1', '2', '3', '4')
#               )
})


##################################################
## CALCULATIONS
##################################################
##
GAS.CONSTANT <- 8.314472 / 100000  # gas constant J / mol K -> umol
Vol.injection <- 1/1000000  # 1 ml injection in m^3

sampleA      <- 6	# sample Area, in cm^2: 6 for rough estimate of inner tube diameter (as used in ARA excel file), or 6.4 for 20 shoots, based on density survey.
sample_to_m2 <- (100*100/sampleA)	# scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample

SECC.ARA.t1 <- within( SECC.ARA.t1, {
  umol.ml <- ( (Pressure..kPa..lab * 1000)*Vol.injection ) /
             ( GAS.CONSTANT * (Temperature...C..lab + 273) )
             # Pressure in pa (*1000)
             # Temperature in C -> K (+273)
  Eth.umol <- X_mol / umol.ml  # % umol C2H4 / umol in 1ml injection.
  Control <- 0  # Control % umol C2H4 corresponding to sample. 
  Blank   <- 0  # Blank   % umol C2H4 corresponding to sample.
  for (i in 1:dim(SECC.ARA.t1)[1])
  {
    Sample.id <- SampleID[i]
    ## GET CONTROL
    Control.id <- paste( substr(Sample.id, 1, 3), "-C", sep="" )   # try for an exact match for ID
    if (Control.id %in% SampleID == FALSE)
      Control.id <- paste( substr(Sample.id, 1, 2), ".-C", sep="" ) # try a control from any chamber in the corresponding block & time.
    # get metching values of SampleID
    ID.control   <- grep( Control.id, SampleID , perl = TRUE, value = TRUE)
    # get values of C2H4 / umol
    Control.umol <- Eth.umol[SampleID==ID.control]
    # take the mean if there is more than one.
    if (length(Control.umol) > 1)
      Control.umol <- mean(Control.umol)
    else if (length(Control.umol) < 1)
      Control.umol <- NA
    ## Store Value
    Control[i] <- Control.umol
    ## GET BLANK
    Blank.id <- paste( substr(Sample.id, 1, 7), "c", sep="" )   # try for an exact match for ID
    if (Blank.id %in% SampleID == FALSE)
      Blank.id <- paste( substr(Sample.id, 1, 6), ".c", sep="" ) # try a from any patch in the corresponding fragmentation treatment.
    # get metching values of SampleID
    ID.blank   <- grep( Blank.id, SampleID , perl = TRUE, value = TRUE)
    # get values of C2H4 / umol
    Blank.umol <- Eth.umol[SampleID==ID.blank]
    # take the mean if there is more than one.
    if (length(Blank.umol) > 1)
      Blank.umol <- mean(Blank.umol)
    else if (length(Blank.umol) < 1)
      Blank.umol <- NA
    ## Store Value
    Blank[i] <- Blank.umol
  }
  ARA  <- X_mol - ( (Blank + Control) * umol.ml )  # umol C2H4 in sample
          #     - (% from Blank and Control) * Total umol
  ARA[SampleControl != "Sample"] <- NA  # no values for non-samples (they would be meaningless).
  ARA.m <- ARA * sample_ml * sample_to_m2  # ARA in umol / m^2 / day

  rm(list=c('i', 'Sample.id', 'Control.id', 'ID.control', 'Control.umol',
            'Blank.id', 'ID.blank', 'Blank.umol'))  # housecleaning
})

##################################################
## MERGE TIME POINTS
##################################################
## 

SECC.ARA <- SECC.ARA.t1

## Remove old objects from memory
rm.objects <- c('SECC.ARA.t1')
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[Data_objects!=rm.objects] , 'SECC.ARA' )

##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.ARA, "SECC columns") <- c('ARA', 'ARA.m')
attr(SECC.ARA, "labels") <- list(
                                 "ARA" ="Acetylene Reduction Assay",
                                 "ARA.m"="Acetylene Reduction Assay",
                                 "ARA.g"="Acetylene Reduction Assay"
                                 )
attr(SECC.ARA, "units")  <- list(
                                 "ARA" =expression(mu * "mol " * ml^-1* d^-1),
                                 "ARA.m"=expression(mu * "mol " * g^-1 * d^-1),
                                 "ARA.g"=expression(mu * "mol " * m^-2 * d^-1)
                                 )

##################################################
## CHECK DATA
##################################################

head(SECC.ARA)  # have a peek at the first 6 rows & columns: is this what you expected?
# str(SECC.ARA.t1)   # check structure: are the appropriate variables factors, numeric, etc.?

##################################################
## SAVE DATA
##################################################
# leave in memory
# only rows for samples, exclude controls
# SECC.ARA.t1 <- SECC.ARA.t1[SECC.ARA.t1$SampleControl=="Sample",]
