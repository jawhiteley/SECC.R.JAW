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
X.cols <- c("Nfix", "H2O")  

##==============================================================
## CUSTOM SETTINGS
##==============================================================
Save.results  <- FALSE                 # Output Results?
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

## Y.plotlab <- bquote( .(Y.label) * "  " * log[10](.(Y.units)) *  "" )
Y.plotlab <- SECC.axislab(SECCa, Y.col)

Fig.filename <- sprintf("%sFigure-%s~", SaveDir.plots(), Y.col
###                         , paste(which(levels(SECC$Time) == Time.use), collapse="") 
)
Suppl.filename <- sprintf("%sSupplemental-%s~", SaveDir.plots(), Y.col)

cat("== Setup complete ==\n")
cat("--< Ready for N-fixation analysis >--\n")
