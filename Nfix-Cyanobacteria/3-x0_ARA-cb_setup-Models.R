################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Initialize, load & process data, configure analysis options
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-12-06
################################################################
## INITIALISE (wth default settings)
################################################################
source('./Nfix-Cyanobacteria/3_ARA-cb_setup-Models.R')



################################################################
## PROCESS DATA: planned
################################################################
## drop values of X == 0 
## - detection errors where I didn't count any cells 
##   (doesn't mean there were none in the sample)
## Unfortunately, this happens too often: errors in model fitting.  
## May have to drop some variables?
SECCa <- SECCa[SECCa$X.trans != 0, ]
## SECCa <- SECCa[SECCa$Y.trans != 0, ] # effect() won't work (too unbalanced?)



################################################################
## CUSTOM SETTINGS
################################################################
Save.results  <- TRUE                  # Output Results?
Save.glmulti  <- paste(SaveDir.obj(), "ARA-cb.glmulti_x0.R", sep="")
UseClimateFac <- FALSE                 # Combine Chamber & Position into pseudo-factor?
UseMM         <- FALSE                 # Use Mixed Models?



##==============================================================
## LABELS
##==============================================================
## Save Output to Files - set to NULL to prevent output.
Save.filename <- paste("Results - ", Y.col, "~", X.col, "-x0 - ",
                       paste(which(levels(SECC$Time) == Time.use), collapse=""),
                       sep = ""
                   )
Save.text  <- paste(SaveDir.text(),  Save.filename, ".txt", sep = "")
Save.plots <- paste(SaveDir.plots(), Save.filename, ".pdf", sep = "")
Save.final <- Save.plots              # Destination for final plots.

## Output text
Save.head     <- paste(Save.head, "** Rows with Cell densities of 0 removed **", sep="\n")
Save.head.txt <- Save.header(Save.head)

## Initialize Output file (following scripts *append* to this file)
## if (Save.results == TRUE && is.null(Save.text) == FALSE) capture.output(cat(""), file=Save.text)

