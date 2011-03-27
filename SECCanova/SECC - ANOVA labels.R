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

## Output Files - set to NULL to prevent output.
Out.filename  <- paste("Results - ", Y.col, " - ",
                   which(levels(SECC$Time) == Time.use), sep = ""
                   )
Out.text   <- paste("./output/", Out.filename, ".txt", sep = "")
Out.plots  <- paste("./graphs/", Out.filename, ".pdf", sep = "")
Out.final  <- Out.plots              # Destination for final plots.


## Output text
Out.divider <-        "================================================================\n"
Out.header  <- paste( "Nested ANOVA Results for:", Y.label, "(", Y.col, ")",
                    "\nTransformation used:     ", Y.use,
                    "\nSample Time:  ", paste(Time.use,     collapse = ", "),
                    "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                    "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                    "\nPatches:      ", paste(Position.use, collapse = ", "),
                    "\n\n", date(), paste("\n\n", Out.divider, sep = "")
                    )
Out.patch.header    <- "================  Patch scale Results  =========================\n\n"
Out.mc.header <- paste("\n",
					   "================================================================", 
					   "================  Meta-Community scale Results  ================",
					   "================================================================\n",
					   sep = "\n"
					   )

Out.end       <- paste( "\n", 
						"<============================= END ============================>",
						sep = "\n"
					  )

