################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Cyanobacteria density (Cells.m)
### Jonathan Whiteley     R v2.12     2012-08-04
################################################################
## INITIALISE (wth default settings)
################################################################
source('./CH4-model-fitting/1_reg_settings.R')

Y.col <- 'Cells'      # Column to analyze as response variable           *****
# explanatory vars for data exploration (and labels)
X.cols <- c("H2O", "TAN")  

##==============================================================
## CUSTOM SETTINGS
##==============================================================
Save.results  <- TRUE                 # Output Results?
Save.glmulti  <- paste(SaveDir.obj(), "Cells.glmulti.R", sep="")

## Default labels & calcs
source('./CH4-model-fitting/1_reg_setup.R')



##==============================================================
## PROCESS DATA: planned
##==============================================================
SECCa <- within( SECCa, 
				{
				  Y.log <- logCells                # pre-transform (custom)
				  Y.trans <- Y.log                 # 
				})


##==============================================================
## CUSTOM LABELS
##==============================================================


cat("== Setup complete ==\n")
cat("--< Ready for Cyanobacteria density analysis >--\n")
