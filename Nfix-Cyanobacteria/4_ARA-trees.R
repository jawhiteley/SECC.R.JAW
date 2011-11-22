################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Regression trees:
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density (& other explanatory variables)
### Jonathan Whiteley     R v2.12     2011-11-05
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (MBP).
  getwd()  # Check that we're in the right place

  ## Load data, functions, etc.  Process data & setup config. values.  
  ## Includes rm(list=ls()) to clear memory
  source('./Nfix-Cyanobacteria/1_ARA-cb_setup.R')
}


## library(lattice)    # ggplot2 with faceting is easier!
library(ggplot2)
library(ggdendro)                      # ggplot functions for dendrograms :D
theme_set(theme_bw())                  # change global ggplot2 theme
library(rpart)                         # Recursive Partitioning & Regression Trees 
                                       # (used by Zuur et al. 2007) 
library(tree)                          # Classification & Regression Trees
## rpart or tree?
## tree appears to be "newer", being an alternative to rpart written by B.D. Ripley
## (rpart was ported to R, from S, by B.D. Ripley)
## but I find the output of rpart() easier to read, clearer, and easier to control.
## rpart seems to offer more options, information, diagnostic plotting tools, etc.

################################################################
## PROCESS DATA: planned
################################################################
## Repeated here for assurance and easy reference
SECCa <- within( SECCa, {
                Y.trans <- Y           # convenience
                X.trans <- X           # convenience
})
## Regression trees are not affected by transformations in the EXPLANATORY variables,
##   but they are sensitive to transformations in the RESPONSE variable...
## log-transformation DOES change the resulting tree in subtle ways, 
##   but qualitative results are similar.
## My preference is to use untransformed data, with no strong reason to use transformed data;
## log-transformed ARA data may be used as a response in regression,
## in order to linearize the relationship.
## log-transformation should not be necessary here, although might be more consistent?

## drop values of X == 0 ?
## - detection errors where I didn't count any cells 
##   (doesn't mean there were none in the sample)
## SECC.X0 <- SECCa[SECCa$X.trans != 0, ]
UseClimateFac <- FALSE                 # appears to make no difference! (equally unimportant)

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
Y.main   <- Y.trans ~ Block + X.trans + H2O + Time + Chamber + Frag + Position
Y.mainCl <- Y.trans ~ Block + X.trans + H2O + Time + Climate + Frag

if (UseClimateFac) {
  Y.main  <- Y.mainCl
}

## not strictly necessary: even if only 1 level of Time is included, it just becomes a variable that is useless for splitting (all observations are already in the same group)
## included for comparison, to verify the above.
Y.mainT <- update(Y.main, .~. -Time)   # this changes the order!! puts offset() @ end :(
Y.mainT <- Y.trans ~ offset(Block) + X.trans + H2O + Chamber + Frag + Position
 
BlockOffset <- if (length(grep("offset\\(Block\\)", paste(Y.main, collapse=" "))) > 0) TRUE else FALSE


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
ARA.tree  <- rpart(Y.main, data=SECCa, control=rpc)  # no interaction terms
ARA.tree1 <- rpart(Y.main, data=SECCa, subset=Time==levels(Time)[1], control=rpc)
ARA.tree2 <- rpart(Y.main, data=SECCa, subset=Time==levels(Time)[2], control=rpc)
ARA.tree3 <- rpart(Y.main, data=SECCa, subset=Time==levels(Time)[3], control=rpc)
ARA.treeD <- rpart(Y.main, data=SECCa) # using default settings; no pruning

## Pruning
if (FALSE) {                           # repeated tree-fitting for 'optimal' size
  ## cptable values are unstable, due to cross-validation process; better to choose a consistent value, based on multiple runs...

  cpMinError <- c()                      # empty vector
  cpMin      <- c()                      # empty vector
  for (i in 1:1000) {
    ThisTree <- rpart(Y.main, data=SECCa, control=rpc)  # no interaction terms
    ARA.cptable <- ThisTree$cptable             # matrix, not a data.frame!
    ARA.cpMinError <- ARA.cptable[which.min(ARA.cptable[,"xerror"]),"CP"]
    MaxXerror <- mean(ARA.cptable[, "xerror"]) + ARA.cptable[nrow(ARA.cptable),"xstd"]
    SE1 <- min(ARA.cptable[, "xerror"]) + ARA.cptable[which.min( ARA.cptable[, "xerror"] ), "xstd"]
    ARA.cpMin <- ARA.cptable[ which(ARA.cptable[, "xerror"] < SE1)[1], "CP"]
    cpMinError <- c(cpMinError, ARA.cpMinError)
    cpMin      <- c(cpMin     , ARA.cpMin     ) # less variable over multiple repetitions
  }

  ARA.cp <- mean(cpMin)
  ARA.cp <- mean(cpMinError)
  ## Approx. long-term averages:
  ## Block as a Factor
  ## cpMin:           minsplit 20: 0. ; minsplit 10: 0.038
  ## cpMinError:      minsplit 20: 0. ; minsplit 10: 0.0105
  ## minsplit 10, t1: cpMin 0.147 ; cpMinError 0.0485
  ## minsplit 10, t2: cpMin 0.183 ; cpMinError 0.176
  ## minsplit 10, t3: cpMin 0.113 ; cpMinError 0.038
  ## Block as an Offset
  ## cpMin:           minsplit 20: 0.039  ; minsplit 10: 0.056
  ## cpMinError:      minsplit 20: 0.0088 ; minsplit 10: 0.023
  ## minsplit 10, t1: cpMin 0.17  ; cpMinError 0.085
  ## minsplit 10, t2: cpMin 0.136 ; cpMinError 0.126
  ## minsplit 10, t3: cpMin 0.100 ; cpMinError 0.047
}

if (BlockOffset) {
  if (minsplit==20) {
    ARA.cp <- 0.01                       # actual means are a little lower, but the resulting cutoff is effectively the same
  } else {
    ARA.cp <- 0.02                       # 0.02
  }
  ARA.treeP  <- prune(ARA.tree , cp=ARA.cp)
  ARA.treeP1 <- prune(ARA.tree1, cp=0.05) # 0.08 ; but this allows an extra split
  ARA.treeP2 <- prune(ARA.tree2, cp=0.02) # 0.12 ; strange tree: # splits accelerates rapidly
  ARA.treeP3 <- prune(ARA.tree3, cp=0.03) # 0.045 ; a couple of extra splits?
} else {
  if (minsplit==20) {
    ARA.cp <- 0.01                     # actual means are a little lower, but the resulting cutoff is effectively the same
  } else {
    ARA.cp <- 0.015                    # reduce smaller branches (over-fitting)
  }
  ARA.treeP  <- prune(ARA.tree , cp=ARA.cp)
  ARA.treeP1 <- prune(ARA.tree1, cp=0.04) # 0.04 ; but this allows an extra split
  ARA.treeP2 <- prune(ARA.tree2, cp=0.17) # 0.17 ; strange tree: # splits accelerates rapidly
  ARA.treeP3 <- prune(ARA.tree3, cp=0.04) # 0.04 ; a couple of extra splits?
}

if (FALSE) {                           # tree package
  ## using the tree package - pretty much the same, but with minor differences
  ## doesn't understand offset variables?
  ARA.ttree <-  tree(Y.main, data=SECCa)
}



##==============================================================
## Diagnostics
##==============================================================
op <- par(mfrow=c(2,2))
for (i in 0:3) {
  thisTree <- paste("ARA.tree", ifelse(i==0, "", i), sep="")
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
print(ARA.tree)
print(ARA.treeP)
print(ARA.treeP1)
print(ARA.treeP2)
print(ARA.treeP3)

par.label <- paste("minsplit =", minsplit)
ARA.plot <- plot(ARA.treeP, main=paste("Pruned", par.label), 
                 minbranch = 0, margin=0.1, compress=FALSE)
text(ARA.treeP, use.n = TRUE, cex=0.8)
## trees for separate time points: t2 is essentially the same as what would appear under the t2 branch of the full tree: the interesting comparison is t1 & t3 (the August samples)
plot(ARA.treeP2, main=ARA.treeP2$call, 
     minbranch = 0, margin=0.1, compress=FALSE)
text(ARA.treeP2, use.n = TRUE, cex=0.8)
plot(ARA.treeP1, main=ARA.treeP1$call, 
     minbranch = 0, margin=0.1, compress=FALSE)
text(ARA.treeP1, use.n = TRUE, cex=0.8)
plot(ARA.treeP3, main=ARA.treeP3$call, 
     minbranch = 0, margin=0.1, compress=FALSE)
text(ARA.treeP3, use.n = TRUE, cex=0.8)

if (FALSE) {                           # Full trees: normally not needed, once pruning criteria established
  plot(ARA.treeD, main="rpart() Default settings")
  text(ARA.treeD, use.n = TRUE)
  plot(ARA.tree, main=par.label, compress=TRUE)
  text(ARA.tree, cex=0.6)
  plot(ARA.tree2, main=ARA.tree2$call, 
       minbranch = 0)
  text(ARA.tree2, use.n = TRUE, cex=0.8)
  plot(ARA.tree1, main=ARA.tree1$call, 
       minbranch = 0)
  text(ARA.tree1, use.n = TRUE, cex=0.8)
  plot(ARA.tree3, main=ARA.tree3$call, 
       minbranch = 0)
  text(ARA.tree3, use.n = TRUE, cex=0.8)
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
  capture.output(cat(Save.header), 
				 print(ARA.treeP),                 # model
				 summary(ARA.treeP),               # model summary
				 cat("\n\n"),                      # for output
				 cat(Save.end),                    # END OUTPUT #
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

## post(ARA.treeP, filename="", horizontal=FALSE)   # ugly
## plot.rpart returns a list of x & y coordinates for *nodes & leaves* (ARA.plot)
## labels as mathematical expressions?
## grid.text(parse(text=Node.labels), x=ARA.tree.labels$x, y=ARA.tree.labels$y)

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
  Y.lim[1] <- Y.lim[1] *0.8
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

ARA.tree.plot  <- RegTreePlot.SECC(ARA.treeP , minsplit)
ARA.tree1.plot <- RegTreePlot.SECC(ARA.treeP1, minsplit) + opts(title = levels(SECCa$Time)[1])
ARA.tree2.plot <- RegTreePlot.SECC(ARA.treeP2, minsplit) + opts(title = levels(SECCa$Time)[2])
ARA.tree3.plot <- RegTreePlot.SECC(ARA.treeP3, minsplit) + opts(title = levels(SECCa$Time)[3])


if (Save.results == TRUE && is.null(Save.final) == FALSE && Save.plots != Save.final) {
  FileName = sub(".pdf$", ".eps", Save.final)
  FileName = paste("./graphs/Figure -", Y.col, "- RegTree - 123.eps")
  ggsave(filename = FileName, plot = ARA.tree.plot )
} else {
  print(ARA.tree.plot)
  print(ARA.tree1.plot)
  print(ARA.tree2.plot)
  print(ARA.tree3.plot)
}
