################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Regression trees: N-fixation (ARA)
### Jonathan Whiteley     R v2.12     2012-07-12
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
## Load data, functions, etc.  Process data & setup config. values.  
source('./CH4-model-fitting/Nfix_setup.R')


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
## My preference is to use untransformed data, with no strong reason to use transformed data;
## log-transformed ARA data may be used as a response in regression,
## in order to linearize the relationship.
## log-transformation should not be necessary here, although might be more consistent?




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
##  - not 100% sure what this does, but should model a Block 'effect', without including it in the tree (or calculating parameters).
##  - splits do not change, but group means do.
##  - The position of the offset variable **does** matter
##    - I think it should be first, to account for Blocks *before* any other factors?
Y.main   <- Y ~ Block + Warming + Frag + H2O + log(Cells.m) + TAN
Y.mainB  <- Y ~ offset(Block) + Warming + Frag + H2O + log(Cells.m) + TAN
 
BlockOffset <- FALSE
if (BlockOffset) Y.main <- Y.mainB


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
  ## cpMin:           minsplit 20: 0.024 ; minsplit 10: 0.040
  ## cpMinError:      minsplit 20: 0.006 ; minsplit 10: 0.015
  ## Block as an Offset
  ## cpMin:           minsplit 20: 0. ; minsplit 10: 0.
  ## cpMinError:      minsplit 20: 0. ; minsplit 10: 0.
  ## Block as a Factor (with Growth)
  ## cpMin:           minsplit 20: 0.060 ; minsplit 10: 0.062
  ## cpMinError:      minsplit 20: 0.015 ; minsplit 10: 0.015
}

if (minsplit==20) {
  Y.cp <- 0.006                    # actual means are a little lower, but the resulting cutoff is effectively the same
} else {
  Y.cp <- 0.015                    # reduce smaller branches (over-fitting)
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
Y.plot <- plot(Y.treeP, main=paste("Pruned", par.label), 
                 minbranch = 0, margin=0.1, compress=FALSE)
text(Y.treeP, use.n = TRUE, cex=0.8)

if (FALSE) {                           # Full trees: normally not needed, once pruning criteria established
  plot(Y.treeD, main="rpart() Default settings")
  text(Y.treeD, use.n = TRUE)
  plot(Y.tree, main=par.label, compress=TRUE)
  text(Y.tree, cex=0.6)
}

if (FALSE) {                           # tree package
  print(ARA.ttree)
  plot(ARA.ttree)
  text(ARA.ttree)
}





################################################################
## SAVE OUTPUT
################################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  capture.output(cat(Save.head.txt), 
				 print(ARA.treeP),                 # model
				 summary(ARA.treeP),               # model summary
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

RegTreePlot.SECC <- function(Reg.tree, minsplit) {
  require(ggplot2)
  require(ggdendro)
  ARA.cp <- min(Reg.tree$cptable[, "CP"])
  Suppl.text  <- paste("cp = ", ARA.cp, ", split min. n = ", minsplit, sep="")

  Yvals       <- as.numeric(Reg.tree$frame[, "yval"])
  Ylabels     <- paste(" ", formatC(Yvals, digits=1, flag="-", format="f"), sep="")
  Nvals       <- as.numeric(Reg.tree$frame[, "n"])
  Nlabels     <- paste("", formatC(Nvals), sep="") # "n="
  Droot       <- as.numeric(Reg.tree$frame[1, "dev"]) 
  Dvals       <- as.numeric(Reg.tree$frame[, "dev"]) / Droot # proportion of total deviance
  Dlabels     <- paste("", formatC(Dvals, digits=2, big.mark=",", flag="-", format="f"), sep="")
  Leaf.labels <- paste(Ylabels, " (", Dlabels, ")", sep="") #
  Info.labels <- paste(" n=", Nlabels, sep="") # ARA.tree.labels
  Var.labels  <- labels(Reg.tree, pretty=F) # includes full factor levels, decision criteria, etc.
  ## clean up labels
  Var.labels <- gsub(" months,", ",", Var.labels) #  duplicate month labels
  Var.labels <- gsub("H2O(..[0-9.]+)", "Moisture\\1%", Var.labels)
  Var.labels <- gsub("X.trans", "Cells", Var.labels)
  Var.labels <- gsub("Full Corridors", "Corridors", Var.labels)
  Var.labels <- gsub("Pseudo-Corridors", "Pseudo-Cs", Var.labels)
  if (FALSE) { # ggplot2 doesn't support math expressions (yet)?
    Var.labels <- gsub("H2O(..[0-9.]+)", "H[2]*O\\1%", Var.labels)
    Var.labels <- gsub("\\de\\+(\\d+)", " %*% 10^{\\1}", Var.labels)
    Var.labels <- expression(Var.labels)      # expressions
  } else {
  }
  Var.labels  <- paste(" ", Var.labels, sep="")
  Node.labels <- paste(Var.labels, "\n", Dlabels, sep="") # space in front of each line for padding (hack) ; blank line for leaf labels in between
  Full.labels <- paste(" ", Var.labels, "\n ", Ylabels, " (n=", Nlabels, ") ", Dlabels, sep="")
  ARA.tree.labels <- data.frame(label=Full.labels, 
                                node=Node.labels, leaf=Leaf.labels, info=Info.labels,
                                var=Var.labels, dev=Dlabels, n=Nlabels, yval=Ylabels)
  ## x=ARA.plot$x, y=ARA.plot$y, 

  class(Reg.tree) <- c("rpart", "tree") # largely the same, but no methods for "rpart"
  dend_data <- dendro_data(Reg.tree)
  ## dend.labels <- ARA.tree.labels$label[as.numeric(row.names(dend_data$labels)) +1]
  numLabels <- nrow(ARA.tree.labels) -1  # drop "root"
  branch_labels <- segment(dend_data)[1:numLabels, c("x", "y", "xend", "yend")]     # coordinates for segments / branches
  ## names(branch_labels) <- c("x", "y")    # rename columns ;)
  branch_labels <- cbind(branch_labels, ARA.tree.labels[-1, 3:ncol(ARA.tree.labels)])
  dend_data$branch_label <- branch_labels # package it all up
  dend.leaf   <- ARA.tree.labels$leaf[ as.numeric(row.names(dend_data$leaf_label)) ]
  leaf_labels <- dend_data$leaf_label
  leaf_labels$label <- dend.leaf
  dend_data$leaf_label <- leaf_labels
  ## move leaves from edge to leave space for labels?
  ## ylim doesn't work too well, but I did discover that the plotting area will expand
  ## to include all coordinates (on all layers), so I can just add some (blank)
  ## text at the limits to ensure the axes go at least that far.
  Y.lim <- c( min(segment(dend_data)$y), max(segment(dend_data)$y) )
  Y.lim[1] <- Y.lim[1] *0
  Suppl.label <- data.frame(label=Suppl.text, x=(max(segment(dend_data)$x) +0.5), 
                            y=max(segment(dend_data)$y), stringsAsFactors=FALSE)
  Suppl.label <- rbind(Suppl.label, 
                       data.frame(label=" ", x=1, y=Y.lim[1], stringsAsFactors=FALSE)
  ) # add some space 'below' bottom leaves for labels.

  RegTree.plot <- ggplot(segment(dend_data)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=1, color="grey50") +
  scale_size_identity() + scale_color_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=branch_labels, 
                                             aes(label=var, x=xend, y=yend),
                                             hjust=0, vjust=-0.7, size=3.1, colour="grey20") + 
                                   scale_size_identity() + scale_colour_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=branch_labels, 
                                             aes(label=info, x=xend, y=yend),
                                             hjust=0, vjust=1.75, size=3.1, colour="grey20") + 
                                   scale_size_identity() + scale_colour_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=branch_labels, 
                                             aes(label=leaf, x=x, y=y),
                                             hjust=0, vjust=0.5, size=3) + 
                                   scale_size_identity()
  ## RegTree.plot <- RegTree.plot + geom_text(data=leaf_labels, aes(label=label, x=x, y=y), 
  ##                                            hjust=-0.1, size=4) + scale_size_identity()
  RegTree.plot <- RegTree.plot + geom_text(data=Suppl.label, aes(label=label, x=x, y=y), 
                                             hjust=0, size=4) + scale_size_identity()
  RegTree.plot <- RegTree.plot + coord_flip() + scale_y_reverse() + 
  theme_dendro() + opts(panel.border = theme_blank())
  ## expand=c(0.2, 0) within the sacle_y_reverse() adds space to *both* sides, to make room for long labels.  Not ideal, but it works
  ## expand=c(0.2, 0)

  ## these ggplot objects are not entirely stand-alone if they depend on the dendro objects!!
  ## dendro_data()
  RegTree.plot                          # return
}

Y.tree.plot  <- RegTreePlot.SECC(Y.treeP , minsplit)

if (Save.results == TRUE) {
  FileName = paste(Save.plot.dir, gsub("Results", "Figure", Save.filename), ".eps", sep="")
  ggsave(filename = FileName, plot = Y.tree.plot, width=6, height=6, scale=1.5)
} else {
  print(Y.tree.plot)
}
