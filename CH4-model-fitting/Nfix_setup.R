################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Acetylene Reduction Assay (ARA: N-fixation)
### Jonathan Whiteley     R v2.12     2012-07-15
################################################################
## INITIALISE (wth default settings)
################################################################
source('./CH4-model-fitting/1_reg_settings.R')

Y.col <- 'Nfix'      # Column to analyze as response variable           *****
# explanatory vars for data exploration (and labels)
X.cols <- c("Cells.m", "Hcells.m", "H2O", "TAN")  

##==============================================================
## CUSTOM SETTINGS
##==============================================================
Save.results  <- TRUE                 # Output Results?
Save.glmulti  <- paste(SaveDir.obj(), "Nfix.glmulti.R", sep="")

## Default labels & calcs
source('./CH4-model-fitting/1_reg_setup.R')



##==============================================================
## PROCESS DATA: planned
##==============================================================
SECCa <- within( SECCa, 
				{
				  Y.log <- log10(Y)                # There's a value of Nfix <1 :(
				  Y.log[Y < 1] <- 0                # log10(<1) will be negative ... log10(<0) = NaN
				  Y.trans <- Y.log                 # 
				})


##==============================================================
## LABELS
##==============================================================

## Y.plotlab <- bquote( .(Y.label) * "  " * log[10](.(Y.units)) *  "" )
Y.plotlab <- SECC.axislab(SECCa, Y.col)

Fig.filename <- sprintf("%sFigure-%s~", SaveDir.plots(), Y.col
###                         , paste(which(levels(SECC$Time) == Time.use), collapse="") 
)
Suppl.filename <- sprintf("%sSupplemental-%s~", SaveDir.plots(), Y.col)

cat("== Setup complete ==\n")
cat("--< Ready for N-fixation analysis >--\n")
