################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Jonathan Whiteley     R v2.12     2012-07-31
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("..")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}


##==============================================================
## CALCULATIONS 
##==============================================================
SECC.prime <- SECC    # save a copy of the original for reference.

# merge in Fauna statistics
source("./CH4-model-fitting/1_Fauna_proc.R")                 # Prepare Fauna data & stats for merging.

## Merge with SECC data for direct comparisons
SECC <- merge(SECC, recodeSECC(SECC.sp.sum), all = TRUE)
## merge attributes: this leaves duplicates, but at least it's all there.
attr(SECC, "labels") <- c( attr(SECC.prime, "labels"), attr(SECC.sp.sum, "labels") )
attr(SECC, "units")  <- c( attr(SECC.prime, "units"),  attr(SECC.sp.sum, "units") )
SECC.labels <- unique( names( attr(SECC, "labels") ) )
SECC.units  <- unique( names( attr(SECC, "units" ) ) )
attr(SECC, "labels") <- attr(SECC, "labels")[SECC.labels]
attr(SECC, "units")  <- attr(SECC, "units" )[SECC.units]
cat("Fauna data merged into main data frame.\n")


## Conversion factors defined in /lib/SECC.functions.R
SECC <- within( SECC, { 
  ## change negative ARA values to 0 - should I wait until after aggregation?
  ARA.ml[ARA.ml < 0] <- 0
  ARA.m[ ARA.m  < 0] <- 0
  ARA.g[ ARA.g  < 0] <- 0
  Nfix    <- ARA.m * Nfix.ARA.ratio
  logNfix <- log10(Nfix)               # There's a value of Nfix <1 :(
  logNfix[Nfix < 1] <- 0               # log10(<1) will be negative ... log10(<0) = NaN
  H2O     <- H2O * 100
  H2O.wwt <- H2O.wwt * 100
  Growth  <- grow12 + grow23           # moss growth during second year **
  ## Should I be log-transforming moss Growth, or just use non-linear GLMM (log-link?)?  Note negative values!
  logCells <- log10(Cells.m +1)        # log-transform of Cyanobacteria;  log10(Cells +1) 
  ##   logCells[Cells.m <= 0] <- 0
  Decomp.asq <- asin(sqrt(Decomposition)) # proportions 0-1; probably better off using glm() anyway.
  logTAN <- log10(TAN)                 # log-transform of Total (Available) Nitrogen
  ## recoded factors / new explanatory variables
  ## Don't use these: predicted values are weird and unconvincing.
  ## Chamber treatments as degrees of warming (I'm interpolating for Partial chambers for now, but I probably have the real values somewhere - will probably never matter, as these are unlikely to be included in the analyses)
  Warming <- factor( Chamber, levels = c("Ambient", "Partial Chamber", "Full Chamber"), 
                    labels = c(0, 0.26, 0.52) )
  Warming <- as.numeric(as.character(Warming))
  TempC   <- factor( Chamber, levels = c("Ambient", "Partial Chamber", "Full Chamber"), 
                    labels = c(6.65, 6.98, 7.32) )
  TempC   <- as.numeric(as.character(TempC))
  Climate <- paste(Position, Chamber) # pseudo-factor to simplify modelling: fewer interactions to deal with.
  Climate <- gsub(".*Ambient", "Ambient", Climate)
  Climate <- gsub("(.*)Full Chamber", "\\1Chamber", Climate)
  Climate <- factor(Climate)
})

attr(SECC, "labels")[["logCells"]] <- attr(SECC, "labels")[["Cells.m"]] # Cells ?
attr(SECC, "units" )[["logCells"]] <- quote( log("cells" %.% m^-2) )    # cells / shoot ?
attr(SECC, "labels")[["Growth"]] <- "Moss growth"
attr(SECC, "units" )[["Growth"]] <- quote("mm" %.% "yr"^-1)
attr(SECC, "labels")[["Decomp.asq"]] <- attr(SECC, "labels")[["Decomposition"]]
attr(SECC, "units" )[["Decomp.asq"]] <- quote(asin(sqrt("% mass loss" %.% "yr"^-1)))
attr(SECC, "labels")[["Warming"]] <- "Warming"
attr(SECC, "units" )[["Warming"]] <- quote(delta ~ "°C vs. ambient")
attr(SECC, "labels")[["TempC"]] <- "Mean annual T"
attr(SECC, "units" )[["TempC"]] <- quote("°C")



##==============================================================
## LABELS
##==============================================================

Y.label <- attr(SECC, "labels")[[Y.col]]  # response variable label
Y.units <- attr(SECC, "units" )[[Y.col]]  # response variable units

## Y.plotlab <- bquote( .(Y.label) * "  " * log[10](.(Y.units)) *  "" )
Y.plotlab <- SECC.axislab(SECC, Y.col)


## Save Output to Files - set to NULL to prevent output.
Subdir <- "synth/"
Save.plot.dir <- sprintf("%s%s", SaveDir.plots(), Subdir)
Save.filename <- paste("Results-", Y.col, "~",
                       paste(which(levels(SECC$Time) == Time.use), collapse=""),
                       sep = ""
                   )
Save.filename <- sprintf("Results-%s~%s", Y.col,
                         paste(which(levels(SECC$Time) == Time.use), collapse="") )
Save.text  <- paste(SaveDir.text(),  Save.filename, ".txt", sep = "")
Save.plots <- paste(SaveDir.plots(), Subdir, Save.filename, ".pdf", sep = "")
Save.final <- Save.plots              # Destination for final plots.

## Output text
Save.divider <- Save.div()
Save.head    <- paste(  "Results for:",   Y.label, "(", Y.col, ") ~", "(others)",
                      "\nExpt. Time:   ", paste(Time.use,     collapse = ", "),
                      "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                      "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                      "\nPatches:      ", paste(Position.use, collapse = ", ")
                      )
Save.head.txt <- Save.header(Save.head)
Save.end.txt  <- Save.end()


Fig.filename <- sprintf("%s%sFigure-%s~", SaveDir.plots(), Subdir, Y.col
###                         , paste(which(levels(SECC$Time) == Time.use), collapse="") 
)
Suppl.filename <- sprintf("%s%sSupplemental-%s~", SaveDir.plots(), Subdir, Y.col)



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


## Filter Outliers (see Exploration Graphs)
if (ExcludeOutliers)
{
  SECCa$Cells.m[SECCa$Cells.m > 5e+09] <- NA # moot after transformation
  SECCa$Growth[SECCa$Growth > 30] <- NA
  SECCa$H2O[SECCa$H2O > 800] <- NA     #  I believe it's real, but it may be highly influential
  ##   SECCa$Decomposition[SECCa$Decomposition > 0.30] <- NA
}


SECCa <- within( SECCa, 
                {
                  Y <- as.numeric( get(Y.col) )
                  Y.log <- log10(Y)    # Generic: some variables (Nfix) will use a special version
                  Y.trans <- Y         # default
                })




################################################################
## CHECK DATA
################################################################
if (FALSE) {  # do not run if source()d
  head(SECCa)     # have a peek at the first 6 rows & columns: is this what you expected?
  str( SECCa)     # check structure: are the appropriate variables factors, numeric, etc.?
  summary(SECCa)  # summary statistics
}

cat("== SECC data loaded & processed for multivariate / regression analyses ==\n")
