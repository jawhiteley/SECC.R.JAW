################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Regression trees: Moss Growth (mm)
### Jonathan Whiteley     R v2.12     2012-07-31
################################################################
## INITIALISE
################################################################
if (FALSE) 
{  ## Working Directory: see lib/init.R below [\rd in Vim]
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("../")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}
## Clear memory, load data, functions, etc.  Process data & setup config. settings
source('./CH4-model-fitting/Decomp_setup.R')


## library(lattice)    # ggplot2 with faceting is easier!
library(ggplot2)
library(ggdendro)                      # ggplot functions for dendrograms :D
theme_set(theme_bw())                  # change global ggplot2 theme
library(rpart)                         # Recursive Partitioning & Regression Trees 
                                       # (used by Zuur et al. 2007) 
## library(tree)                          # Classification & Regression Trees

################################################################
## PROCESS DATA: planned
################################################################
## Regression trees are not affected by transformations in the EXPLANATORY variables,
##   but they are sensitive to transformations in the RESPONSE variable...
## log-transformation DOES change the resulting tree in subtle ways, 
##   but qualitative results are similar.
##==============================================================
cat("- Processing Data (excluding NAs)\n")
## Remove NA rows from data, but only for variables in the model ...
## Easier to do it once here, than in every single fitting function in the script ;)
SECCf <- SECCa
Mod.cols <- unlist( strsplit("SampleID, Block, Time, Chamber, Frag, Position, TempC, H2O, ARA.m, Nfix, logNfix, TAN, logTAN, Growth, Richness, Grazers, Y, Y.trans", 
                             ", ", fixed = TRUE) )
SECCa <- na.exclude( SECCa[, Mod.cols] ) # only exclude NAs from relevant (used) columns
attr(SECCa, "labels") <- attr(SECCf, "labels")
attr(SECCa, "units")  <- attr(SECCf, "units")
levels(SECCa$Chamber)[levels(SECCa$Chamber) == "Full Chamber"] <- "Chamber"




################################################################
## ANALYSIS
################################################################
## Goal: Identify which variables have most relative importance,
##       and identify any interactions with experimental treatments.
##==============================================================
## Model Formula
##==============================================================
## Main effects only
##  include Block just to see relative importance
##  can also include Block as an 'offset' variable: offset(Block)
Y.main   <- Y ~ Block + Chamber + Frag + H2O + Richness + Grazers + logTAN + logNfix


##==============================================================
## Regression Trees
##==============================================================
minsplit <- 10
rpc <- rpart.control(minsplit=minsplit, cp=0.001)      
## default minsplit is 20; smaller number = more (deeper) splits!
##  more splits => different tree structure, and larger optimal tree?
##  deeper splits: may end up pruning less 'significant' branches elsewhere
##  - deep branches end up with large groups, other leaves split to smaller groups.
## Do I care more about big differences within small groups (smaller minsplit)
## or small differences within large groups (larger minsplit)?
## default cp = 0.01; smaller number = longer cptable (more candidate splits)
Y.tree  <- rpart(Y.main, data=SECCa, control=rpc)  # no interaction terms
Y.treeD <- rpart(Y.main, data=SECCa) # using default settings; no pruning

## Pruning
if (FALSE) {                           # repeated tree-fitting for 'optimal' size
  ## cptable values are unstable, due to cross-validation process; better to choose a consistent value, based on multiple runs...

  cpMinError <- c()                      # empty vector
  cpMin      <- c()                      # empty vector
  for (i in 1:1000) {
    ThisTree <- rpart(Y.main, data=SECCa, control=rpc)  # no interaction terms
    Y.cptable <- ThisTree$cptable             # matrix, not a data.frame!
    Y.cpMinError <- Y.cptable[which.min(Y.cptable[,"xerror"]),"CP"]
    MaxXerror <- mean(Y.cptable[, "xerror"]) + Y.cptable[nrow(Y.cptable),"xstd"]
    SE1 <- min(Y.cptable[, "xerror"]) + Y.cptable[which.min( Y.cptable[, "xerror"] ), "xstd"]
    Y.cpMin <- Y.cptable[ which(Y.cptable[, "xerror"] < SE1)[1], "CP"]
    cpMinError <- c(cpMinError, Y.cpMinError)
    cpMin      <- c(cpMin     , Y.cpMin     ) # less variable over multiple repetitions
  }

  Y.cp  <- mean(cpMin)
  Y.cpE <- mean(cpMinError)

  ## Approx. long-term averages:
  ## Block as a Factor
  ## cpMin:           minsplit 20: 0. ; minsplit 10: 0.24
  ## cpMinError:      minsplit 20: 0. ; minsplit 10: 0.074
  ## Block as an Offset
  ## cpMin:           minsplit 20: 0. ; minsplit 10: 0.
  ## cpMinError:      minsplit 20: 0. ; minsplit 10: 0.
}

if (minsplit==20) {
  Y.cp <- 0.01
} else {
  Y.cp <- 0.01                        # no pruning - it's pretty thin, anyway
}
Y.treeP  <- prune(Y.tree , cp=Y.cp)



##==============================================================
## Diagnostics
##==============================================================
op <- par(mfrow=c(2,2))
for (i in c(0, "D")) {
  thisTree <- paste("Y.tree", ifelse(i==0, "", i), sep="")
  thisTree <- get(thisTree)
  if (FALSE) {
    summary(thisTree)                  # verbose!
    printcp(thisTree)                  # included in rsq.rpart()
  }
  ## plots
  meanvar(thisTree)
  rsq.rpart(thisTree)                   # 2 plots + printcp()
  plotcp(thisTree)
  mtext(thisTree$call, side=3, padj=2, outer=TRUE)
}
par(op)



##################################################
## ANALYSIS: GET RESULTS
##################################################
print(Y.tree)
print(Y.treeP)

par.label <- paste("minsplit =", minsplit)
plot(Y.tree, main=paste("Full", par.label), 
                 minbranch = 0, margin=0.1, compress=FALSE)
text(Y.tree, use.n = TRUE, cex=0.8)

Y.plot <- plot(Y.treeP, main=paste("Pruned", par.label), 
                 minbranch = 0, margin=0.1, compress=FALSE)
text(Y.treeP, use.n = TRUE, cex=0.8)



################################################################
## SAVE OUTPUT
################################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  Save.text <- gsub("\\.txt", "-trees.txt", Save.text)
  capture.output(cat(Save.head.txt), 
				 print(Y.treeP),                 # model
				 summary(Y.treeP),               # model summary
				 cat("\n\n"),                      # for output
				 cat(Save.end.txt),                # END OUTPUT #
				 file = Save.text
				)
}




################################################################
## FINAL GRAPHICS
################################################################
## I want a HORIZONTAL dendrogram, with criteria ABOVE each branch,
##   and mean values & n's BELOW each branch
##   Mean values & n's at each terminal leaf
##   + histograms of y-values at each leaf?

Y.tree.plot  <- RegTreePlot.SECC(Y.treeP , minsplit)

if (Save.results == TRUE) {
  FileName = paste(Save.plot.dir, gsub("Results", "Figure", Save.filename), "-tree.eps", sep="")
  ggsave(filename = FileName, plot = Y.tree.plot, width=6, height=6, scale=1.5)
} else {
  print(Y.tree.plot)
}
