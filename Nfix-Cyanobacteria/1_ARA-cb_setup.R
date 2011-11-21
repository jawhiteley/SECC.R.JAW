################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-11-02
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (MBP).
  getwd()  # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')

library(car)                           # diagnostic plots & tools

################################################################
## CONFIGURE BASIC ANALYSIS
################################################################

### Response Variable *****
X.col <- 'Cells.m'   # Column to analyze as explanatory variable        *****
Y.col <- 'ARA.m'     # Column to analyze as response variable           *****
## I would prefer to use measurements on a per-gram dwt of moss basis ('.g'),
## but there aren't enough of such measurements for all ARA values (univariate).
## There are certainly enough for this analysis, but I'm worried about inconsistency
## for readers if N-fixation is analyzed on its own on a per-m^2 basis, 
## while everything else is on a (preferred) per-g dwt basis.

# explanatory vars for data exploration?
vars.ls   <- c("ARA.m", "Cells.m", "Hcells.m", "H2O")  

##==============================================================
## SETTINGS 
##==============================================================
## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]      # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)         # Frag treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)]     # Patch Positions to include


## Output Results?
Save.results  <- FALSE  


##==============================================================
## CALCULATIONS 
##==============================================================
SECC.prime <- SECC    # save a copy of the original for reference.

# str(SECC)
sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample
ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * (12.5^2)      # patch area
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2
Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.

SECC <- within( SECC, { 
  ## change negative ARA values to 0 - should I wait until after aggregation?
  ARA.ml[ARA.ml < 0] <- 0
  ARA.m[ ARA.m  < 0] <- 0
  ARA.g[ ARA.g  < 0] <- 0
  Nfix    <- ARA.m * Nfix.ARA.ratio
  H2O     <- H2O * 100
  H2O.wwt <- H2O.wwt * 100
})


##==============================================================
## LABELS
##==============================================================

X.label <- attr(SECC, "labels")[[X.col]]  # explanatory variable label
X.units <- attr(SECC, "units" )[[X.col]]  # explanatory variable units

Y.label <- attr(SECC, "labels")[[Y.col]]  # response variable label
Y.units <- attr(SECC, "units" )[[Y.col]]  # response variable units

X.plotlab <- bquote( .(X.label) * "  " * .(X.units) *  "" )
Y.plotlab <- bquote( .(Y.label) * "  " * .(Y.units) *  "" )

labels.ls <- c()
for(i in 1:length(vars.ls) ){
  var <- vars.ls[i]
  labels.ls[i] <- attr(SECC, "labels")[[var]]
}


## Save Output to Files - set to NULL to prevent output.
Save.filename <- paste("Results - ", Y.col, "~" , X.col, " - ",
                       paste(which(levels(SECC$Time) == Time.use), collapse=""),
                       sep = ""
                   )
Save.text  <- paste("./output/", Save.filename, ".txt", sep = "")
Save.plots <- paste("./graphs/", Save.filename, ".pdf", sep = "")
Save.final <- Save.plots              # Destination for final plots.

## Output text
Save.divider <-        "================================================================\n"
Save.header  <- paste( "GLM Results for:", Y.label, "(", Y.col, ")",
                     "\n               ~", X.label, "(", X.col, ")",
                     "\nExpt. Time:   ", paste(Time.use,     collapse = ", "),
                     "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                     "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                     "\nPatches:      ", paste(Position.use, collapse = ", "),
                     paste("\n\n", date(), "\n\n", Save.divider, sep = "")
                     )
Save.end      <- paste("\n", 
					   "<============================= END ============================>",
						sep = "\n"
					  )

###===============================================
### Include Time as a factor 
###===============================================
## Note that Samples at different times are actually independent
## in this design, due to destructive sampling.

## Time.use     <- levels(SECC$Time)      # Include *ALL* Times (as a Treatment)



################################################################
## PROCESS DATA: planned
################################################################
## filter & clean SECC data frame.
SECC.use <- SECCclean(SECC, Time.use, Chamber.use, Frag.use, Position.use)

## Merge in spatial coordinates, and retrieve attributes
SECC.use <- merge(SECC.use, SECC.xy[, c("SampleID", "xE", "yN")], 
                  by="SampleID", all.x=TRUE, all.y=FALSE, sort=TRUE)
attr(SECC.use, "labels") <- attr(SECC, "labels")
attr(SECC.use, "units" ) <- attr(SECC, "units" )

## Summarize data by means across different (or all) positions to prevent unbalanced effects?
## aggregate 'other' patch Positions?
##   unnecessary: only 'inner' & 'outer' included.
SECCp  <- if (FALSE)  SECC_aggregate( SECC.use, trt = 'Position' )  else SECC.use

## Meta-community (regional) -level analyses (ignoring position):
SECCmc <- SECC_aggregate( SECC.use, trt = 'Frag' )

SECC.scale  <- "patch"  # c("patch", "mc")
SECCa <- if(SECC.scale == "patch") SECCp else if (SECC.scale == "mc") SECCmc else SECC


SECCa <- within( SECCa, {
                X <- as.numeric( get(X.col) )
                Y <- as.numeric( get(Y.col) )
                X.log <- log10(X)
                X.log[X <= 0] <- 0
                Y.log <- log10(Y)
                Y.log[Y <= 0] <- 0
                Y.trans <- Y.log  # convenience
                X.trans <- X.log  # convenience
                Climate <- factor( paste(Position, Chamber) ) # psuedo-factor to simplify modelling: fewer interactions to deal with.
})




################################################################
## CHECK DATA
################################################################
if (FALSE) {  # do not run if source()d
  head(SECCa)     # have a peek at the first 6 rows & columns: is this what you expected?
  str( SECCa)     # check structure: are the appropriate variables factors, numeric, etc.?
  summary(SECCa)  # summary statistics
}

