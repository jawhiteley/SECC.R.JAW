################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main Control script for Final Analyses
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-11-02
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("..")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}


################################################################
## CONFIGURE BASIC ANALYSIS
################################################################
## Initialize (clear memory), load & process data, load functions, etc.
## Configure analysis & run-time options
cat("N-fixation ~ Cyanobacteria: Set-up\n")
source('./Nfix-Cyanobacteria/1_ARA-cb_setup.R')




################################################################
## EXPLORE: PLOTS
################################################################
## make some meaningful plots of data to check for predicted (expected) patterns.
cat("N-fixation ~ Cyanobacteria: Data Exploration\n")
source('./Nfix-Cyanobacteria/2_ARA-cb_Explore.R')




################################################################
## ANALYSIS
################################################################
## Goal: Identify which variables have most relative importance,
##       and identify any interactions with experimental treatments.
## Challenges:
## - Complex Model: 7 explanatory variables!
## - Fixed X violation: all explanatory variables were measured, with associated error
## - Violations of Independence: 
##   - plots are arranged in a nested (multilevel) design, in a field experiment
##     with potential spatial patterns in residuals, etc.
##   - Patterns in Residuals vs. Predicted values.
## - Low replication: 8 replications of each combination of experimental treatments, each from a separate Block
## - Violations of Heterogeneity: error variance differs across Blocks, Chambers, H2O, and possibly other explanatory variables
## - Violations of Normality: the more interaction terms are included, the less normal the residuals become.
## - Correlations: among explantory variables (Cells & H2O)?
##   - Experimental treatments are orthogonal
##   - ANOVA revealed virtually no patterns in Cells vs. experimental treatments or H2O.
##==============================================================
## Methods & Approaches to use
##==============================================================
## *GLMM: fit model using Cells, H2O, and experimental treatments
##  - account for nesting of fixed factors?  How??!
##  - Include spatial autocorrelation instead?
if (F) source('./Nfix-Cyanobacteria/3_0_ARA~cb_MixedModel.R')

## Simple Linear Models with multi-model selection
cat("N-fixation ~ Cyanobacteria: Modelling\n")
if (T) {                               # full dataset **
  source('./Nfix-Cyanobacteria/3_ARA-cb_setup-Models.R')
} else {                               # remove 0 Cell densities (non-detection)
  source('./Nfix-Cyanobacteria/3-x0_ARA-cb_setup-Models.R')
}
source('./Nfix-Cyanobacteria/3_ARA~cb_Models.R')

## Regression trees: higher-order interactions & relative importance
cat("N-fixation ~ Cyanobacteria: Regression trees\n")
source('./Nfix-Cyanobacteria/4_ARA-cb_setup-Trees.R')
source('./Nfix-Cyanobacteria/4_ARA-trees.R')

##  - Mantel Test: Is N-fixation similarity related to spatial distance?
##    ~ Partial Mantel Test: Effects of treatment groups on response after removing effect of distance

## Partial Regression (Legendre & Legendre): Variance Decomposition
## - What is *pure* effect of cyanobacteria cell density on N-fixation, after removing effects of other variables?
## - What is the *pure* effect of moisture?
##   - moisture is highly related to other experimental factors, however
##     If I partial these out, the effect of moisture will be minimal.
##   - What would be more relevant is to partial out effect of moisture *before* fitting the experimental factors, to see what their 'pure' effects really are.
## - Partial out effects of moisture (using GAM) before fitting models with cell density and experimental factors...
##   - What about possible Cells:H2O interaction?

## Quantile regression: Upper or lower *limits* imposed by explanatory variables, rather than precise values. (regresstion through quantile of 0.5 == linear regression).

## Subsume Chamber & Position into "Climate" pseudo-treatment?

## Time as a fixed factor, or separtae analysis on each Time?
## - The analyses certainly support having separate analyses in each time point,
##   but I just didn't have time to coordinate consistent or different results across each one.
if (F) source('./Nfix-Cyanobacteria/3_1_ARA~cb_Model_t1.R')
if (F) source('./Nfix-Cyanobacteria/3_3_ARA~cb_Model_t3.R')



################################################################
## END
################################################################
cat("N-fixation ~ Cyanobacteria: FINISHED\n\n")
