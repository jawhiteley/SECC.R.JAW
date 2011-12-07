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
## CUSTOM SETTINGS
################################################################
Save.results  <- TRUE                  # Output Results?
UseClimateFac <- FALSE                 # Combine Chamber & Position into pseudo-factor?



##==============================================================
## LABELS
##==============================================================

X.plotlab <- bquote( .(X.label) * "  " * .(X.units) *  "" )
Y.plotlab <- bquote( .(Y.label) * "  " * .(Y.units) *  "" )


## Save Output to Files - set to NULL to prevent output.
Save.filename <- paste("Results-", Y.col, "~", X.col, "-RegTrees-",
                       paste(which(levels(SECC$Time) == Time.use), collapse=""),
                       sep = ""
                   )
Save.text  <- paste(SaveDir.text(),  Save.filename, ".txt", sep = "")
Save.plots <- paste(SaveDir.plots(), Save.filename, ".pdf", sep = "")
Save.final <- Save.plots              # Destination for final plots.


## Output text
Save.head  <- paste("Analysis:    Regression Trees", Save.head, sep="\n")
Save.head.txt <- Save.header(Save.head)


## Initialize Output file (following scripts *append* to this file)
## if (Save.results == TRUE && is.null(Save.text) == FALSE) capture.output(cat(""), file=Save.text)

