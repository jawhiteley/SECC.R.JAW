################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Regression trees:
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
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

## drop values of X == 0 
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
### Fixed effects
if ( length(Time.use) > 1 ) {
} else {
  Y.fixed <- Y.trans ~ X.trans * H2O * Chamber * Frag * Position
  Y.fixCl <- Y.trans ~ X.trans * H2O * Climate * Frag
  ## Main effects only
  Y.main  <- Y.trans ~ X.trans + H2O + Chamber + Frag + Position
  Y.mainCl<- Y.trans ~ X.trans + H2O + Climate + Frag
}
Y.fixed <- Y.trans ~ X.trans * H2O * Time * Chamber * Frag * Position
Y.fixCl <- Y.trans ~ X.trans * H2O * Time * Climate * Frag
## Main effects only
##  include Block just to see relative importance
##  can also include Block as an 'offset' variable
##  - not 100% sure what this does, but should model a Block 'effect', without including it in the tree (or calculating parameters).
##  - splits do not change, but group means do.
##  - The position of the offset variable **does** matter
##    - I think it should be first, to account for Blocks *before* any other factors?
Y.main   <- Y.trans ~ offset(Block) + X.trans + H2O + Time + Chamber + Frag + Position
Y.mainCl <- Y.trans ~ offset(Block) + X.trans + H2O + Time + Climate + Frag

if (UseClimateFac) {
  Y.fixed <- Y.fixCl
  Y.main  <- Y.mainCl
}


##==============================================================
## Regression Trees
##==============================================================
minsplit <- 20
rpc <- rpart.control(minsplit=minsplit, cp=0.001)      
## default minsplit is 20; smaller number = more (deeper) splits!
##  more splits => different tree structure, and larger optimal tree?
##  deeper splits: may end up pruning less 'significant' branches elsewhere
##  - deep branches end up with large groups, other leaves split to smaller groups.
## default cp = 0.01; smaller number = longer cptable (more candidate splits)
ARA.tree <- rpart(Y.main, data=SECCa, control=rpc)  # no interaction terms
ARA.treeD <- rpart(Y.main, data=SECCa) # using default settings; no pruning

## Pruning
ARA.cptable <- ARA.tree$cptable             # matrix, not a data.frame!
## Optimal cp, according to Quick-R http://www.statmethods.net/advstats/cart.html
## cross-validation is not always reproducible, hence this value may be somewhat unstable.
ARA.cpMinError <- ARA.tree$cptable[which.min(ARA.tree$cptable[,"xerror"]),"CP"]
## optimal cp, according to Zuur et al. 2007
MaxXerror <- mean(ARA.cptable[, "xerror"]) + ARA.cptable[nrow(ARA.cptable),"xstd"]
## actually, rpart calculates the line differently than Zuur et al. describe:
SE1 <- min(ARA.cptable[, "xerror"]) + ARA.cptable[which.min( ARA.cptable[, "xerror"] ), "xstd"]
ARA.cpMin <- ARA.cptable[ which(ARA.cptable[, "xerror"] < SE1)[1], "CP"]
ARA.cp <- ARA.cpMinError

if (FALSE) {                           # repeated tree-fitting for 'optimal' size
  ## cptable values are unstable, due to cross-validation process; better to choose a consistent value, based on multiple runs...
  cpMinError <- c()                      # empty vector
  cpMin      <- c()                      # empty vector
  for (i in 1:1000) {
    ARA.tree <- rpart(Y.main, data=SECCa, control=rpc)  # no interaction terms
    ARA.cptable <- ARA.tree$cptable             # matrix, not a data.frame!
    ARA.cpMinError <- ARA.cptable[which.min(ARA.cptable[,"xerror"]),"CP"]
    MaxXerror <- mean(ARA.cptable[, "xerror"]) + ARA.cptable[nrow(ARA.cptable),"xstd"]
    SE1 <- min(ARA.cptable[, "xerror"]) + ARA.cptable[which.min( ARA.cptable[, "xerror"] ), "xstd"]
    ARA.cpMin <- ARA.cptable[ which(ARA.cptable[, "xerror"] < SE1)[1], "CP"]
    cpMinError <- c(cpMinError, ARA.cpMinError)
    cpMin      <- c(cpMin     , ARA.cpMin     ) # less variable over multiple repetitions
  }
                                       # approx. long-term averages:
  ARA.cp <- mean(cpMin)                # minsize 20: 0.039  ; minsize 10: 0.056
  ARA.cp <- mean(cpMinError)           # minsize 20: 0.0088 ; minsize 10: 0.023
}

if (minsplit==20) {
  ARA.cp <- 0.01                         # actual means are a little lower, but the resulting cutoff is effectively the same
} else {
  ARA.cp <- 0.02
}
ARA.treeP <- prune(ARA.tree, cp=ARA.cp)


if (FALSE) {                           # tree package
  ## using the tree package - pretty much the same, but with minor differences
  ## doesn't understand offset variables?
  ARA.ttree <-  tree(Y.main, data=SECCa)
}



##==============================================================
## Diagnostics
##==============================================================
## summary(ARA.tree)
printcp(ARA.tree)

op <- par(mfrow=c(2,2))
meanvar(ARA.tree)
rsq.rpart(ARA.tree)                   # 2 plots
plotcp(ARA.tree)
par(op)




##################################################
## ANALYSIS: GET RESULTS
##################################################
print(ARA.tree)
print(ARA.treeP)

par.label <- paste("minsplit =", minsplit)
ARA.plot <- plot(ARA.treeP, main=paste("Pruned", par.label), 
                 minbranch = 0, margin=0.1, compress=FALSE)
text(ARA.treeP, use.n = TRUE, cex=0.8)
plot(ARA.tree, main=par.label, compress=TRUE)
text(ARA.tree, cex=0.6)
plot(ARA.treeD, main="rpart() Default settings")
text(ARA.treeD, use.n = TRUE)

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
				 print(Y.formula),              # model
				 anova(Y.model),                # model summary
				 summary(Y.model),              # model summary
				 cat("\n\n"),                   # for output
				 cat(Save.end),                 # END OUTPUT #
				 file = Save.text
				)
}

if (Save.results == TRUE && is.null(Save.plots) == FALSE && Save.plots != Save.final) dev.off()



################################################################
## FINAL GRAPHICS
################################################################
if (Save.results == TRUE && is.null(Save.final) == FALSE && Save.plots != Save.final) pdf( file = Save.final )

## I want a HORIZONTAL dendrogram, with criteria ABOVE each branch,
##   and mean values & n's BELOW each branch
##   Mean values & n's at each terminal leaf

## post(ARA.treeP, filename="", horizontal=FALSE)   # ugly
## ggplot(ARA.tree)
## plot.rpart returns a list of x & y coordinates for *nodes & leaves* (ARA.plot)
## ARA.dend <- as.dendro(ARA.plot)

Suppl.text  <- paste("cp=", ARA.cp, ", split min. n = ", minsplit, sep="")
Yvals       <- as.numeric(ARA.treeP$frame[, "yval"])
Ylabels     <- paste("", formatC(Yvals, digits=1, flag="-", format="f"), sep="")
Nvals       <- as.numeric(ARA.treeP$frame[, "n"])
Nlabels     <- paste("n=", formatC(Nvals), sep="")
Info.labels <- paste(Ylabels, ", ", Nlabels, sep="")
Node.labels <- labels(ARA.treeP, pretty=F) # includes full factor levels, decision criteria, etc.
Full.labels <- paste(" ", Node.labels, "\n", " ", Info.labels, sep="") #space in front of each line for padding (hack)
ARA.dend.data <- data.frame(x=ARA.plot$x, y=ARA.plot$y, label=Full.labels, node=Node.labels, leaf=Info.labels)

class(ARA.treeP) <- c("rpart", "tree") # largely the same, but no methods for "rpart"
ARA.dend <- dendro_data(ARA.treeP)
## ARA.dend.labels <- ARA.dend.data$label[as.numeric(row.names(ARA.dend$labels)) +1]
ARA.dend.labels <- ARA.dend.data$label[-1] # drop "root"
numLabels <- length(ARA.dend.labels)
branch_labels <- segment(ARA.dend)[1:numLabels, c("xend", "yend")]     # coordinates for segments / branches
names(branch_labels) <- c("x", "y")    # rename columns ;)
branch_labels$label <- ARA.dend.labels
ARA.dend.leaf   <- ARA.dend.data$leaf[ as.numeric(row.names(ARA.dend$leaf_label)) ]
leaf_labels <- ARA.dend$leaf_label
leaf_labels$label <- ARA.dend.leaf

Y.lim <- c( min(segment(ARA.dend)$y), max(segment(ARA.dend)$y) )
Y.lim[1] <- 0

ARA.tree.plot <- ggplot(segment(ARA.dend)) +
                 geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
ARA.tree.plot <- ARA.tree.plot + geom_text(data=branch_labels, aes(label=label, x=x, y=y),
                                           hjust = 0, size=3) + scale_size_identity()
## ARA.tree.plot <- ARA.tree.plot + geom_text(data=leaf_labels, aes(label=label, x=x, y=y), 
##                                            hjust=-0.1, size=4) + scale_size_identity()
ARA.tree.plot <- ARA.tree.plot + geom_text(aes(label=Suppl.text, x=max(segment(ARA.dend)$x), y=max(Y.lim[2])), hjust = 0, size=4) + scale_size_identity()
## ARA.tree.plot <- ARA.tree.plot + coord_cartesian(ylim=Y.lim)
ARA.tree.plot <- ARA.tree.plot + coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
                 theme_dendro() + opts(panel.border = theme_blank())
## expand() within the sacle_y_reverse() adds space to *both* sides, to make room for long labels.  Not ideal, but it works
print(ARA.tree.plot)

## ggdendrogram(dendro_data(ARA.ttree))   # still get "Error in eval(expr, envir, enclos) : object 'x0' not found"   wtf?


if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()

