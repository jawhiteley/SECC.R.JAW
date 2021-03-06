##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process ARA data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-03-23
##################################################
## This script is run as part of `./lib/load.R`
# require(car)		# load external package 'car', for recode()

##================================================
## CHECK DATA
##================================================
# str(SECC.ARA.t1)  # Should already be in memory.    str() produces output on source().

## Standardize ID column names (but not values) 
SECC.ARA.t1 <- SECCcolumnNames(SECC.ARA.t1)
SECC.ARA.t2 <- SECCcolumnNames(SECC.ARA.t2)
SECC.ARA.t4 <- SECCcolumnNames(SECC.ARA.t4)
## Store Raw (non-standard, but accurate) IDs for reference later
SampleID.t1 <- SECC.ARA.t1[,c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos')]
SampleID.t2 <- SECC.ARA.t2[,c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos')]
SampleID.t4 <- SECC.ARA.t4[,c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos')]
SECC.ARA.t4 <- within( SECC.ARA.t4, {
  Frag[SampleControl=="blank"]  <- NA  # Remove fake / nonsensical Fragmentation IDs. 
  Pos[ SampleControl=="blank"]  <- NA  # Remove fake / nonsensical Position IDs.
  ## Frag column in SamplID.t4 (original IDs) indicates which
   # fragmentation treatments are included in the blank from that particular Chamber.
})
## Standardize ID column values
SECC.ARA.t1 <- checkSECCdata(SECC.ARA.t1, 'SECC.ARA.t1', CheckDuplicates = FALSE)
SECC.ARA.t2 <- checkSECCdata(SECC.ARA.t2, 'SECC.ARA.t2', CheckDuplicates = FALSE)
SECC.ARA.t4 <- checkSECCdata(SECC.ARA.t4, 'SECC.ARA.t4', CheckDuplicates = FALSE)

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.
Blank.code   <- "G"  # code for Gas Blanks (gass with no moss)
Control.code <- "c"  # code for moss controls (moss with no gas)

SECC.ARA.t1 <- within( SECC.ARA.t1, {
  ## Re-name "Controls" to be less confusing
  levels(SampleControl)[levels(SampleControl)=="Control2"] <- "blank"
  levels(SampleControl)[levels(SampleControl)=="control1"] <- "control"
  ## re-order factor levels.
  SampleControl <- factor(SampleControl, levels=c("blank", "control", "Sample")) 
  ## Re-calculate SampleID to indicate controls
  SampleID <- paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep="")
  SampleID[SampleControl=="blank"] <- paste(Block, Time, Chamber, "-", Blank.code, sep="")[SampleControl=="blank"]
                                         #'G' for gas Blanks?
  SampleID[SampleControl=="control"] <- paste(SampleID, Control.code, sep="")[SampleControl=="control"]
})

SECC.ARA.t2 <- within( SECC.ARA.t2, {
  ## re-order factor levels.
  SampleControl <- factor(SampleControl, levels=c("blank", "control", "Sample")) 
  ## Re-calculate SampleID to indicate controls
  SampleID <- paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep="")
  SampleID[SampleControl=="blank"]   <- paste(Block, Time, Chamber, "-", Blank.code, sep="")[SampleControl=="blank"]
  SampleID[SampleControl=="control"] <- paste(SampleID, Control.code, sep="")[SampleControl=="control"]
})


SECC.ARA.t4 <- within( SECC.ARA.t4, {
  ## re-order factor levels.
  SampleControl <- factor(SampleControl, levels=c("blank", "control", "Sample")) 
  ## Re-calculate SampleID to indicate controls
  SampleID <- paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep="")
  SampleID[SampleControl=="control"] <- paste(SampleID, Control.code, sep="")[SampleControl=="control"]
  SampleID[SampleControl=="blank"]   <- paste(Block, Time, Chamber, "-", Blank.code, sep="")[SampleControl=="blank"]  # not unique: need to append Frag treatments to make them unique.
  Raw.blanks <- grep( "G\\d\\b" , SampleID.t4$SampleID, perl = TRUE)  # indices of RAW blank SampleIDs
  Blank.ids <- SampleID[SampleControl=="blank"]  # temporary container
  Blank.Frags <- SampleID.t4$Frag[Raw.blanks]    # Frag levels: comma-delimited list
  Blank.Frags <- gsub( ",\\s?", "", Blank.Frags )# remove commas
  Blank.ids <- paste( Blank.ids, Blank.Frags, sep = '' )  # append Frag labels
  SampleID[SampleControl=="blank"] <- Blank.ids  # update ID values
  rm(list=c( 'Raw.blanks', 'Blank.Frags', 'Blank.ids' ))  # housekeeping
})



##################################################
## CALCULATIONS
##################################################
## Leave the original calculated columns in, to act as an internal reference
## for missing values in calculations performed here.

GAS.CONSTANT <- 8.314472 / 1e+06  # gas constant J / mol K -> umol
Vol.injection <- 1e-06  # 1 ml injection in m^3 : 1/(1000*1000)

sampleA      <- 6	# sample Area, in cm^2
           #    6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
           # or 6.4 for 20 shoots, based on density survey: 31440 shoots /m^2.
sample.to.m2 <- (100*100)/sampleA	# scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample

###  T1 ARA calculations
SECC.ARA.t1 <- within( SECC.ARA.t1, {
  C2H4._mol <- X_mol  # for compatibility with data from other time points.
  umol.ml <- ( (Pressure..kPa..lab * 1000)*Vol.injection ) /
             ( GAS.CONSTANT * (Temperature...C..lab + 273) )
             # Pressure in pa (*1000)
             # Temperature in C -> K (+273)
  Eth.umol <- X_mol / umol.ml  # % umol C2H4 / umol in 1ml injection.
  Control <- NA  # Control % umol C2H4 corresponding to sample. 
  Blank   <- NA  # Blank   % umol C2H4 corresponding to sample.
  for (i in 1:NROW(SECC.ARA.t1)) {
    Sample.id <- SampleID[i]
    ## GET BLANK
    Blank.id <- paste( substr(Sample.id, 1, 3), "-", Blank.code, sep="" )   # try for an exact match for ID
    if (Blank.id %in% SampleID == FALSE)
      Blank.id <- paste( substr(Sample.id, 1, 2), ".-", Blank.code, sep="" ) # try a control from any chamber in the corresponding block & time.
    # get metching values of SampleID
    ID.blank   <- grep( Blank.id, SampleID , perl = TRUE )  # get indices
    # get values of C2H4 / umol
    Blank.umol <- Eth.umol[ID.blank]
    # take the mean if there is more than one.
    if (length(Blank.umol) > 1) {
      Blank.umol <- mean(Blank.umol, na.rm = TRUE)
    } else if (length(Blank.umol) < 1) {
      Blank.umol <- NA
    } else {
      Blank.umol <- Blank.umol
    }
    ## Store Value
    Blank[i] <- Blank.umol
    ## GET CONTROL
    Control.id <- paste( substr(Sample.id, 1, 7), Control.code, sep="" )   # try for an exact match for ID
    if (Control.id %in% SampleID == FALSE)
      Control.id <- paste( substr(Sample.id, 1, 6), ".", Control.code, sep="" ) # try a from any patch in the corresponding fragmentation treatment.
    # get metching values of SampleID
    ID.control   <- grep( Control.id, SampleID , perl = TRUE )
    # get values of C2H4 / umol
    Control.umol <- Eth.umol[ID.control]
    # take the mean if there is more than one.
    if (length(Control.umol) > 1) {
      Control.umol <- mean(Control.umol, na.rm = TRUE)
    } else if (length(Control.umol) < 1) {
      Control.umol <- NA
    } else {
      Control.umol <- Control.umol
    }
    ## Store Value
    Control[i] <- Control.umol
  }
  rm(list=c('i', 'Sample.id', 'Blank.id', 'ID.blank', 'Blank.umol',
                 'Control.id', 'ID.control', 'Control.umol' ))  # housecleaning

  ARA.ml <- X_mol - ( (Blank + Control) * umol.ml )  # umol C2H4 in 1 ml sample
          #       - (% from Blank and Control) * Total umol
  ARA.ml[SampleControl != "Sample"] <- NA  # no values for non-samples (they would be meaningless).
  ARA.m <- ARA.ml * sample_ml * sample.to.m2  # ARA in umol / m^2 / day
})


###  T2 ARA calculations
SECC.ARA.t2 <- within( SECC.ARA.t2, {
  umol.ml <- ( (Pressure..kPa..lab * 1000)*Vol.injection ) /
             ( GAS.CONSTANT * (Temperature...C..lab + 273) )
             # Pressure in pa (*1000)
             # Temperature in C -> K (+273)
  Eth.umol <- C2H4._mol / umol.ml  # % umol C2H4 / umol in 1ml injection.
  Ac.umol  <- C2H2._mol / umol.ml  # % umol C2H2 / umol in 1ml injection.
  ARA.umol <- C2H2._mol + C2H4._mol  # Total umol of ARA gas in injection.
  C2H4.ARA <- C2H4._mol / ARA.umol   # % C2H4 of ARA gases. *****
  C2H2.ARA <- C2H2._mol / ARA.umol   # % C2H2 of ARA gases (not used).
  Control  <- NA  # Control % umol C2H4 corresponding to sample. 
  Blank    <- NA  # Blank   % umol C2H4 corresponding to sample.
  for (i in 1:nrow(SECC.ARA.t2)) {
    Sample.id <- SampleID[i]
    ## GET BLANK
    Blank.id <- paste( substr(Sample.id, 1, 3), "-", Blank.code, sep="" )   # try for an exact match for ID
    if (Blank.id %in% SampleID == FALSE)
      Blank.id <- paste( substr(Sample.id, 1, 2), ".-", Blank.code, sep="" ) # try a control from any chamber in the corresponding block & time.
    # get metching values of SampleID
    ID.blank   <- grep( Blank.id, SampleID , perl = TRUE )  # get indices
    # get values of [ C2H4 / (C2H4+C2H2) ] in Blank: %C2H4 of (C2H4+C2H2)
    Blank.ARA <- C2H4.ARA[ID.blank]
    # take the mean if there is more than one.
    if (length(Blank.ARA) > 1) {
      Blank.ARA <- mean(Blank.ARA, na.rm = TRUE)
    } else if (length(Blank.ARA) < 1) {
      Blank.ARA <- NA
    } else {
      Blank.ARA <- Blank.ARA
    }
    ## Store Value
    Blank[i] <- Blank.ARA
    ## CALCULATE CONTROL
    Control.id <- paste( substr(Sample.id, 1, 7), Control.code, sep="" )   # try for an exact match for ID
    if (Control.id %in% SampleID == FALSE)
      Control.id <- paste( substr(Sample.id, 1, 6), ".", Control.code, sep="" ) # try a from any patch in the corresponding fragmentation treatment.
    # get metching values of SampleID
    ID.control   <- grep( Control.id, SampleID , perl = TRUE )
    # get values of C2H4 / umol: %C2H4 of (Total umol of gas in 1ml injection)
    Control.ARA <- Eth.umol[ID.control]  # cheap approximation by volume
                   # using this for t2 instead of formulas below
                   # increases # net values < 0 from 95 to 99.
    # umol C2H4 in Control from contamination (Blank):
    Control.Blank <- (Blank.ARA * C2H2._mol[ID.control]) / (1 - Blank.ARA)
    # umol C2H4 in Control from non-ARA sources:
    Control.ARA   <- C2H4._mol[ID.control] - Control.Blank
    # % C2H4 in Total umol of gas in Control
    Control.ARA   <- Control.ARA / umol.ml[ID.control]
    # take the mean if there is more than one.
    if (length(Control.ARA) > 1) {
      Control.ARA <- mean(Control.ARA, na.rm = TRUE)
    } else if (length(Control.ARA) < 1) {
      Control.ARA <- NA
    } else {
      Control.ARA <- Control.ARA
    }
    ## Store Value
    Control[i] <- Control.ARA
  }
  rm(list=c('i', 'Sample.id', 'Blank.id', 'ID.blank', 'Blank.ARA',
                 'Control.id', 'ID.control', 'Control.ARA', 'Control.Blank' ))  # housecleaning

  ## CALCULATE ARA
  ARA.control <- C2H4._mol - (Control * umol.ml)  # umol C2H4 after removing Control
  ARA.ml <- ARA.control - (Blank * (C2H2._mol + ARA.control))  # umol C2H4 / 1 ml sample
          # umol C2H4 without Control (by volume)
          # - (% Blank * (Total ARA umol without Control))
  ARA.ml[SampleControl != "Sample"] <- NA  # no values for non-samples (they would be meaningless).
  ARA.m <- ARA.ml * sample_ml * sample.to.m2  # ARA in umol / m^2 / day
})


###  T4 ARA calculations
SECC.ARA.t4 <- within( SECC.ARA.t4, {
  umol.ml <- ( (Pressure..kPa..lab * 1000)*Vol.injection ) /
             ( GAS.CONSTANT * (Temperature...C..lab + 273) )
             # Pressure in pa (*1000)
             # Temperature in C -> K (+273)
  Eth.umol <- C2H4._mol / umol.ml  # % umol C2H4 / umol in 1ml injection.
  Ac.umol  <- C2H2._mol / umol.ml  # % umol C2H2 / umol in 1ml injection.
  ARA.umol <- C2H2._mol + C2H4._mol  # Total umol of ARA gas in injection.
  C2H4.ARA <- C2H4._mol / ARA.umol   # % C2H4 of ARA gases. *****
  C2H2.ARA <- C2H2._mol / ARA.umol   # % C2H2 of ARA gases (not used).
  Control  <- NA  # Control % umol C2H4 corresponding to sample. 
  Blank    <- NA  # Blank   % umol C2H4 corresponding to sample.
  for (i in 1:nrow(SECC.ARA.t4)) {
    Sample.id <- SampleID[i]
    ## GET BLANK
    Blank.id <- paste( substr(Sample.id, 1, 3), "-", Blank.code, sep="" )   # try for an exact match for ID
#   if (Blank.id %in% SampleID == FALSE)
#     Blank.id <- paste( substr(Sample.id, 1, 2), ".-", Blank.code, sep="" ) # try a control from any chamber in the corresponding block & time.
    ## get metching values of SampleID
    ID.blank   <- grep( Blank.id, SampleID , perl = TRUE )  # get indices
    ## get values of [ C2H4 / (C2H4+C2H2) ] in Blank: %C2H4 of (C2H4+C2H2)
    Blank.ARA <- C2H4.ARA[ID.blank]
    ## further matching or take the mean if there is more than one.
    if (length(Blank.ARA) > 1) {
      ## browser()  # debug
      Sample.Frag <- Frag[i]
      Blank.id    <- SampleID[ID.blank]
      Blank.Frag  <- strsplit(as.character(Blank.id), "")  # split each character
      Blank.Frag  <- lapply( Blank.Frag, function (x) x[-(1:5)] )  # drop first 5 code digits, leaving only Frag codes.
      for ( blank in 1:length(Blank.Frag) ) {
        if  (Sample.Frag %in% unlist(Blank.Frag[blank]) ) {
          ## This is the correct blank for this sample.
          Blank.ARA <- C2H4.ARA[SampleID == Blank.id[blank]]
                    # assumes blanks are in same relative order as returned by grep()
        }
      }
      ## Take the mean if there is still more than 1
      if (length(Blank.ARA) > 1) {
        Blank.ARA <- mean(Blank.ARA, na.rm = TRUE)
      } else {
        Blank.ARA <- Blank.ARA
      }
    } else if (length(Blank.ARA) < 1) {
      Blank.ARA <- NA
    } else {
      Blank.ARA <- Blank.ARA
    }
    ## Store Value
    Blank[i] <- Blank.ARA
    ## CALCULATE CONTROL
    Control.id <- paste( substr(Sample.id, 1, 7), Control.code, sep="" )   # try for an exact match for ID
    if (Control.id %in% SampleID == FALSE)
      Control.id <- paste( substr(Sample.id, 1, 6), ".", Control.code, sep="" ) # try a from any patch in the corresponding fragmentation treatment.
    # get metching values of SampleID
    ID.control   <- grep( Control.id, SampleID , perl = TRUE )
    # get values of C2H4 / umol: %C2H4 of (Total umol of gas in 1ml injection)
    Control.ARA   <- Eth.umol[ID.control]  # cheap approximation by volume
    # umol C2H4 in Control from contamination (Blank):
    Control.Blank <- (Blank.ARA * C2H2._mol[ID.control]) / (1 - Blank.ARA)
    # umol C2H4 in Control from non-ARA sources:
    Control.ARA   <- C2H4._mol[ID.control] - Control.Blank
    # % C2H4 in Total umol of gas in Control
    Control.ARA   <- Control.ARA / umol.ml[ID.control]
    # take the mean if there is more than one.
    if (length(Control.ARA) > 1) {
      Control.ARA <- mean(Control.ARA, na.rm = TRUE)
    } else if (length(Control.ARA) < 1) {
      Control.ARA <- NA
    } else {
      Control.ARA <- Control.ARA
    }
    ## Store Value
    Control[i] <- Control.ARA
  }

  ## Housekeeping
  rm(list=c('i', 'Sample.id', 'Blank.id', 'ID.blank', 'Blank.ARA',
                 'Blank.Frag', 'blank', 'Sample.Frag',
                 'Control.id', 'ID.control', 'Control.ARA', 'Control.Blank' ))
  Blank[SampleControl == "blank"] <- NA  # no values for blanks (they would be meaningless).
  Control[SampleControl != "Sample"] <- NA  # no values for non-samples (they would be meaningless).

  ## CALCULATE ARA
  ARA.control <- C2H4._mol - (Control * umol.ml)  # umol C2H4 after removing Control
  ARA.ml <- ARA.control - (Blank * (C2H2._mol + ARA.control))  # umol C2H4 / 1 ml sample
          # umol C2H4 without Control (by volume)
          # - (% Blank * (Total ARA umol without Control))
  ARA.ml[SampleControl != "Sample"] <- NA  # no values for non-samples (they would be meaningless).
  ARA.m <- ARA.ml * sample_ml * sample.to.m2  # ARA in umol / m^2 / day
})



##################################################
## MERGE TIME POINTS
##################################################
## merge(), or rbind() on common columns?

SECC.ARA <- SECC.ARA.t1
SECC.ARA <- merge( SECC.ARA, SECC.ARA.t2, all=TRUE, suffixes = c(".t1", ".t2") )
SECC.ARA <- merge( SECC.ARA, SECC.ARA.t4, all=TRUE, suffixes = c(".12", ".t4"), sort = TRUE )

##================================================
## Check Results
nrow(SECC.ARA)  # 1780
nrow(SECC.ARA.t1) + nrow(SECC.ARA.t2) + nrow(SECC.ARA.t4)  # 1780

##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.ARA, "SECC columns") <- c('ARA.ml', 'ARA.m')
attr(SECC.ARA, "labels") <- list(
                                 "ARA.ml" = "Acetylene Reduction",
                                 "ARA.m"  = "Acetylene Reduction",
                                 "ARA.g"  = "Acetylene Reduction"
                                 )
attr(SECC.ARA, "units")  <- list(
                                 "ARA.ml" = quote(mu*"mol" %.% ml^-1 %.% d^-1),
                                 "ARA.g"  = quote(mu*"mol" %.% g^-1  %.% d^-1),
                                 "ARA.m"  = quote(mu*"mol" %.% m^-2  %.% d^-1)
                                 )

##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d
  head(SECC.ARA)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.ARA)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.ARA[SECC.ARA$SampleControl == "Sample", ARAcalc.cols])
  SECCstr(SECC.ARA[SECC.ARA$SampleControl == "Sample", ])

###===============================================
### CHECK Calculations
  # Histograms
  hist(SECC.ARA[SECC.ARA$Time == 1, "ARA.ml"])
  hist(SECC.ARA[SECC.ARA$Time == 2, "ARA.ml"])
  # compare to calculations using standardization by total umol (in Excel file).
  hist(SECC.ARA[SECC.ARA$Time == 2, "ARA.._mol.ml."])
  hist(SECC.ARA[SECC.ARA$Time == 4, "ARA.ml"])
  hist(SECC.ARA[SECC.ARA$Time == 4, "ARA.._mol.ml."])

  ## Are calculated values in the right ballpark?
  ## Early results (t1):       0 ~  300 umol m^-2 day^-1 ARA
  ## DeLuca et al. (2002):     0 ~  700 umol m^-2 day^-1 ARA
  ## Zackrisson et al. (2004): 0 ~  400 umol m^-2 day^-1 ARA
  ## DeLuca et al. (2007):     0 ~ 1000 umol m^-2 day^-1 ARA
  ## Gentili et al. (2005):    0 ~ 3500 umol g^-1 day^-1 ARA (dwt Moss)
  hist(SECC.ARA[, "ARA.m"])
  boxplot( list("t1"=SECC.ARA[SECC.ARA$Time == 1, "ARA.m"],
                "t2"=SECC.ARA[SECC.ARA$Time == 2, "ARA.m"],
                "t4"=SECC.ARA[SECC.ARA$Time == 4, "ARA.m"]
                )
          )
  mean(SECC.ARA[, "ARA.m"], na.rm = TRUE)

  # Number of values > 0?
  length(na.omit(SECC.ARA[SECC.ARA$Time == 1 & SECC.ARA$ARA.ml > 0, "ARA.ml"]))
  length(na.omit(SECC.ARA[SECC.ARA$Time == 2 & SECC.ARA$ARA.ml > 0, "ARA.ml"]))  # low
  length(na.omit(SECC.ARA[SECC.ARA$Time == 4 & SECC.ARA$ARA.ml > 0, "ARA.ml"]))
  # compare to calculations using standardization by total umol (in Excel file).
  length(na.omit(SECC.ARA[SECC.ARA$Time == 2 & SECC.ARA$ARA.._mol.ml. > 0, "ARA.._mol.ml."]))
  length(na.omit(SECC.ARA[SECC.ARA$Time == 4 & SECC.ARA$ARA.._mol.ml. > 0, "ARA.._mol.ml."]))

  # Compare raw values, blanks & controls
  hist(SECC.ARA[SECC.ARA$SampleControl == "Sample",  "C2H4._mol"])
  hist(SECC.ARA[SECC.ARA$SampleControl == "control", "C2H4._mol"])
  hist(SECC.ARA[SECC.ARA$SampleControl == "blank",   "C2H4._mol"])


  ARA.calcs <- c("SampleID", # "Block", "Time", "Chamber", "Frag", "Pos",
                 "SampleControl", "C2H4._mol", "umol.ml", "Eth.umol", "C2H4.ARA", 
                 "Control", "Blank", "ARA.ml", "ARA.._mol.ml.",
                 "ARA.m", "ARA.._mol..m.2..d.1.")
  invisible(edit(SECC.ARA[, ARA.calcs]))

  ## largest ARA values (are they really outliers?)
  ARA.big <- na.omit( SECC.ARA$SampleID[SECC.ARA$ARA.ml > 0.01] )
  SECC.ARA[SECC.ARA$SampleID %in% ARA.big, ARA.calcs]

  ## major differences between calcualted values here, and Excel Spreadsheet?
  Big.diffs <- with( SECC.ARA, SampleID[abs(ARA.m - ARA.._mol..m.2..d.1.) > 1])
  length(na.omit(Big.diffs))  # check here before proceeding with next line.
  SECC.ARA[SECC.ARA$SampleID %in% Big.diffs, ARA.calcs]
  invisible(edit(SECC.ARA[SECC.ARA$SampleID %in% Big.diffs, ARA.calcs]))
  

### t4 Control & Blank in Excel spreadsheet are incorrect approximations
### (used to be empty)
}

##################################################
## SAVE DATA
##################################################
# leave in memory
ARA.full <- SECC.ARA
# only rows for samples, exclude controls
# SECC.ARA <- SECC.ARA[SECC.ARA$SampleControl=="Sample",]

##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.ARA.t1', 'SECC.ARA.t2', 'SECC.ARA.t4',
                'SampleID.t1', 'SampleID.t2', 'SampleID.t4',
                'Blank.code', 'Control.code',
                'GAS.CONSTANT', 'Vol.injection', 'sample_ml'
                )
rm(list=rm.objects)
## Update list of Data_objects for importing
Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.ARA' )

