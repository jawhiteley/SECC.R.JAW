################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-12-06
################################################################
## INITIALISE (wth default settings)
################################################################
source('./Nfix-Cyanobacteria/1_ARA-cb_setup.R')



################################################################
## PROCESS DATA: planned
################################################################
## Repeated here for assurance and easy reference
SECCa <- within( SECCa, {
                Y.trans <- Y.log  # convenience: log10
                X.trans <- X.log  # convenience: log10
})



################################################################
## CUSTOM SETTINGS
################################################################
Save.results  <- TRUE                  # Output Results?
Save.glmulti  <- paste(SaveDir.obj(), "ARA-cb.glmulti.R", sep="")
UseClimateFac <- FALSE                 # Combine Chamber & Position into pseudo-factor?
UseMM         <- FALSE                 # Use Mixed Models?



##==============================================================
## LABELS
##==============================================================

X.plotlab <- bquote( .(X.label) * "  " * log[10](.(X.units)) *  "" )
Y.plotlab <- bquote( .(Y.label) * "  " * log[10](.(Y.units)) *  "" )

Fig.filename <- sprintf("%sFigure-%s~%s", SaveDir.plots(), Y.col, X.col
###                         , paste(which(levels(SECC$Time) == Time.use), collapse="") 
)
Suppl.filename <- sprintf("%sSupplemental-%s~%s", SaveDir.plots(), Y.col, X.col)

## Initialize Output file (following scripts *append* to this file)
## if (Save.results == TRUE && is.null(Save.text) == FALSE) capture.output(cat(""), file=Save.text)

