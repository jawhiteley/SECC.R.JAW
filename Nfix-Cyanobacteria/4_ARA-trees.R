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
theme_set(theme_bw())                  # change global ggplot2 theme
library(rpart)                         # Recursive Partitioning & Regression Trees 
                                       # (used by Zuur et al. 2007) 
library(tree)                          # Classification & Regression Trees
## library(MASS)                          # Modern Applied Statistics in S (Venables & Ripley 2002)

################################################################
## PROCESS DATA: planned
################################################################
## Repeated here for assurance and easy reference
SECCa <- within( SECCa, {
                Y.trans <- Y  # convenience
                X.trans <- X  # convenience
})
## drop values of X == 0 
## - detection errors where I didn't count any cells 
##   (doesn't mean there were none in the sample)
## Unfortunately, this happens too often: errors in model fitting.  
## May have to drop some variables?
## SECC.X0 <- SECCa[SECCa$X.trans != 0, ]
UseClimateFac <- FALSE

################################################################
## ANALYSIS
################################################################
## Goal: Identify which variables have most relative importance,
##       and identify any interactions with experimental treatments.
##==============================================================
## Model Formula
##==============================================================
### Fixed effects
## Include H2O^2 to test for unimodal relationship?
## Including Time as a factor?
## - It might be better to analyse each time period separately, and avoid the higher-order interactions that I know exist.
if ( length(Time.use) > 1 ) {
  Y.fixed <- Y.trans ~ X.trans * H2O * Time * Chamber * Frag * Position
  Y.fixCl <- Y.trans ~ X.trans * H2O * Time * Climate * Frag
  ## Main effects only
  Y.main  <- Y.trans ~ X.trans + H2O + Time + Chamber + Frag + Position
  Y.mainCl<- Y.trans ~ X.trans + H2O + Time + Climate + Frag
} else {
  Y.fixed <- Y.trans ~ X.trans * H2O * Chamber * Frag * Position
  Y.fixCl <- Y.trans ~ X.trans * H2O * Climate * Frag
  ## Main effects only
  Y.main  <- Y.trans ~ X.trans + H2O + Chamber + Frag + Position
  Y.mainCl<- Y.trans ~ X.trans + H2O + Climate + Frag
}
if (UseClimateFac) {
  Y.fixed <- Y.fixCl
  Y.main  <- Y.mainCl
}


##==============================================================
## Regression Trees
##==============================================================
ARA.tree1 <- rpart(Y.main, data=SECCa)  # no interaction terms
ARA.tree2 <-  tree(Y.main, data=SECCa)  # no interaction terms



##################################################
## ANALYSIS: GET RESULTS
##################################################
print(ARA.tree1)
print(ARA.tree2)

plot(ARA.tree1)
text(ARA.tree1)
plot(ARA.tree2)
text(ARA.tree2)






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



if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()

