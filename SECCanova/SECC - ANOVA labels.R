##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Default labels for standardized SECC Nested ANOVA
### Jonathan Whiteley     R v2.12     2011-03-26
##################################################

Y.plotlab <- bquote( .(Y.label) * "  " * .(Y.units) *  "" )

Dataset.labels <- c( "Patch scale data", "Meta-Community scale data" )
Dataset.list   <- c("SECCp", "SECCmc")
ID.cols    <- c('SampleID', 'Time', 'Block', 'Chamber', 'Frag', 'Pos', 'Position')
Trt.nested <- c('Time', 'Block', 'Chamber', 'Frag', 'Position')

## Save Output to Files - set to NULL to prevent output.
Save.filename <- paste("Results - ", Y.col, " - ",
                       paste(which(levels(SECC$Time) == Time.use), collapse=""),
                       sep = ""
                   )
Save.text  <- paste("./output/", Save.filename, ".txt", sep = "")
Save.plots <- paste("./graphs/", Save.filename, ".pdf", sep = "")
Save.final <- paste("./graphs/", "Figure - ", Y.col, sep = "")    # Destination for final plots.


## Output text
Save.divider <-        "================================================================\n"
Save.header  <- paste( "Nested ANOVA Results for:", Y.label, "(", Y.col, ")",
                     "\nTransformation used:     ", Y.use,
                     "\nExpt. Time:   ", paste(Time.use,     collapse = ", "),
                     "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                     "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                     "\nPatches:      ", paste(Position.use, collapse = ", "),
                     paste("\n\n", date(), "\n\n", Save.divider, sep = "")
                     )
Save.patch.header   <- "================  Patch scale Results  =========================\n\n"
Save.mc.header <- paste("\n",
					   "================================================================", 
					   "================  Meta-Community scale Results  ================",
					   "================================================================\n",
					   sep = "\n"
					   )

Save.end      <- paste("\n", 
					   "<============================= END ============================>",
						sep = "\n"
					  )

