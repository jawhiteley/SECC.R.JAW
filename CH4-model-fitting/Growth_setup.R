################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Moss Growth (mm)
### Jonathan Whiteley     R v2.12     2012-07-12
################################################################
## INITIALISE (wth default settings)
################################################################
source('./CH4-model-fitting/1_reg_settings.R')

Y.col <- 'Growth'      # Column to analyze as response variable           *****
# explanatory vars for data exploration (and labels)
X.cols <- c("Nfix", "H2O", "TAN")  

##==============================================================
## CUSTOM SETTINGS
##==============================================================
Save.results  <- TRUE                 # Output Results?
Save.glmulti  <- paste(SaveDir.obj(), "Growth.glmulti.R", sep="")

## Default labels & calcs
source('./CH4-model-fitting/1_reg_setup.R')



##==============================================================
## PROCESS DATA: planned
##==============================================================
## Some values of Growth <0 (legitimately), so log-transform may not be viable, or would have to ne handled specially
SECCa <- within( SECCa, 
                {
                  Y.trans <- Y  # convenience
                })


##==============================================================
## LABELS
##==============================================================


cat("== Setup complete ==\n")
cat("--< Ready for Moss Growth analysis >--\n")
