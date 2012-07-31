################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Moss Growth (mm)
### Jonathan Whiteley     R v2.12     2012-07-31
################################################################
## INITIALISE (wth default settings)
################################################################
source('./CH4-model-fitting/1_reg_settings.R')

Y.col <- 'Decomp.asq'      # Column to analyze as response variable           *****
# explanatory vars for data exploration (and labels)
X.cols <- c("Richness", "Grazers", "H2O", "Nfix", "TAN")  
## including Richness severely reduces the available samples, and treatment coverage (no pseudo-corridors)

##==============================================================
## CUSTOM SETTINGS
##==============================================================
Save.results  <- TRUE                 # Output Results?
Save.glmulti  <- paste(SaveDir.obj(), "Decomp.glmulti.R", sep="")

## Default labels & calcs
source('./CH4-model-fitting/1_reg_setup.R')



##==============================================================
## PROCESS DATA: planned
##==============================================================
## arcsin-square-root transformed already
SECCa <- within( SECCa, 
                {
                  Y.trans <- Y  # convenience
                })


##==============================================================
## LABELS
##==============================================================


cat("== Setup complete ==\n")
cat("--< Ready for Decomposition analysis >--\n")
