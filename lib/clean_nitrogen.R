##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process Decomposition data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-12-11
##################################################
## This script is run as part of `./lib/load.R`

##================================================
## CHECK DATA
##================================================
SECC.N.raw <- SECC.N
if (FALSE){
  SECC.N <- SECC.N.raw # in case recovery is needed
}

## Standardize ID column names, but not values
SECC.N <- checkSECCdata(SECC.N, 'SECC.N', CheckValues = FALSE)

##################################################
## CALCULATIONS
##################################################
## convert to useful units
## total amount / ... g of resin?  capsule? surface area?
## raw data is in mg / L
mg.g <- 1/1000                         # 1000 mg / g
Capsule.L     <- 3*16 / 1000           # 3 aliquots of 16 ml / 1000 ml/L (per capsule); sometimes a little less (-5 ml)
CapsuleDiam   <- 1.9                   # 0.75 inches (1.9 -- 2 cm) diameter
CapsuleRadius <- CapsuleDiam/2         # Capsule Radius (for SA calculation)
CapsuleVol    <- (4/3)*pi*(CapsuleRadius^3) # VOLUME of a sphere = 4/3*pi*r^3
CapsuleSA     <- 4*pi*(CapsuleRadius^2)     # Surface area of a sphere = 4*pi*r^2
Capsule.m2    <- CapsuleSA / (100*100) # cm^2 -> m^2
Blank1 <- which(SECC.N$Label=="B")
Blank2 <- which(SECC.N$Label=="B2")
Blank3 <- which(SECC.N$Label=="B3")
BlankRows  <- c(Blank1, Blank2, Blank3)
Blank1rows <- which(SECC.N$Label %in% c(1:96, 289:300))
Blank2rows <- which(SECC.N$Label %in% 97:192)
Blank3rows <- which(SECC.N$Label %in% 193:288) 
SECC.N <- within( SECC.N, 
                 {
                   Years <- Days / 365   # Length of time capsules were in the field
                   Label <- as.character(Label) # I want raw values, not factor levels
                   ## Lookup Blanks
                   NH4Blanks <- rep(NA, length(NH4))
                   NH4Blanks[Blank1rows] <- NH4[Blank1]
                   NH4Blanks[Blank2rows] <- NH4[Blank2]
                   NH4Blanks[Blank3rows] <- NH4[Blank3]
                   NO3Blanks <- rep(NA, length(NO3))
                   NO3Blanks[Blank1rows] <- NO3[Blank1]
                   NO3Blanks[Blank2rows] <- NO3[Blank2]
                   NO3Blanks[Blank3rows] <- NO3[Blank3]
                   ## Subtract Blanks
                   NH4 <- NH4 - NH4Blanks
                   NO3 <- NO3 - NO3Blanks
                   ##                    browser()
                   ## Convert units
                   ## raw data is in mg / L
                   ## after DeLuca et al. 2008 Science (Supplementary Methods) 
                   ## - N *deposition* (not available in soil?)
                   ## - scaled up from SA of capsule (4.9 cm^2 ?) to m^2
                   NH4 <- (NH4 * mg.g * Capsule.L) / Capsule.m2 # mg/L -> g/capsule -> g/m^2
                   NO3 <- (NO3 * mg.g * Capsule.L) / Capsule.m2 # mg/L -> g/capsule -> g/m^2
                   NH4 <- NH4 / Years    # g / m^2 / Year
                   NO3 <- NO3 / Years    # g / m^2 / Year
                   ## Some samples were processed with dummy samples, to simplify labels & tracking
                   ## Replace dummy values with real NAs
                   NArows <- which(is.na(Total.N))
                   NH4[NArows] <- NA
                   NO3[NArows] <- NA
				   ## Calculate Total N
                   TAN <- NH4 + NO3      # Total Available N: In/Organic?
                   ## Clean-up
                   rm(NH4Blanks, NO3Blanks, NArows)
                 })


##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run when source()'d

  head(SECC.N)  # have a peek at the first 6 rows & columns: is this what you expected?
  str(SECC.N)   # check structure: are the appropriate variables factors, numeric, etc.?
  ## Check structure
  SECCstr(SECC.N)

###===============================================
### CHECK Calculations
  # Histograms
  hist(SECC.N[, c("NO3", "NH4", "TAN")])
  boxplot(SECC.N[, c("NO3", "NH4", "TAN")])
  boxplot(log10(SECC.N[, c("NO3", "NH4", "TAN")] +1))

  summary(SECC.N, na.rm = TRUE)
}


##================================================
## Assign Attributes
##================================================
# "SECC columns" determines which response variable columns will be merged into final data frame.

attr(SECC.N, "SECC columns") <- c('NH4', 'NO3', 'TAN')
attr(SECC.N, "labels") <- list("NH4"     = quote("Available"~NH[4]^"+"~"-N"),
                               "NO3"     = quote("Available"~NO[3]^"-"~"-N"),
                               "TAN"     = "Total Available Nitrogen" # In/Organic?
                               )
attr(SECC.N, "units")  <- list("NH4"     = quote(g %.% m^-2 %.% yr^-1),
                               "NO3"     = quote(g %.% m^-2 %.% yr^-1),
                               "TAN"     = quote(g %.% m^-2 %.% yr^-1)
                               )

##################################################
## SAVE DATA
##################################################
# leave in memory
N.full <- SECC.N
##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
## Manually clean & prepare data for automatic checking.
# Keep rows for experiment samples, drop blanks, other samples
SECC.N <- SECC.N[-BlankRows, ]
SECC.N <- SECC.N[-which(SECC.N$Chamber %in% c("", "X")), ]
SECC.N$SampleID <- SECC_sampleID(SECC.N) # important for merging!


##================================================
## Housekeeping
##================================================
## Remove old objects from memory
rm.objects <- c('SECC.N.raw', 'N.full')
rm(list=rm.objects)
## Update list of Data_objects for importing
## Data_objects <- c( Data_objects[!(Data_objects %in% rm.objects)] , 'SECC.N' )

