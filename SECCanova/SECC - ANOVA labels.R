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

## Load eps graphics for Plots
library(grImport)
setwd("./save/")
if ( !file.exists("hexagon.eps.xml") & !file.exists("Frag-Black-icons.eps.xml") )
{
  FragIcons <- PostScriptTrace("Frag-Black-icons.eps")
  Hex       <- PostScriptTrace("hexagon.eps")
}
FragIcons   <- readPicture("Frag-Black-icons.eps.xml")
Hex         <- readPicture("hexagon.eps.xml")
setwd("..")
FragIcon1 <- FragIcons[49:50]          # 1. Continuous 
FragIcon2 <- FragIcons[29:48]          # 2. Corridors
FragIcon3 <- FragIcons[ 9:28]          # 3. Pseudo-corridors
FragIcon4 <- FragIcons[ 1:8 ]          # 4. Isolated

