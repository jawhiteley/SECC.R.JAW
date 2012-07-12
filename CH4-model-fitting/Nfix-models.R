################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Regression models for N-fixation (ARA)
### Jonathan Whiteley     R v2.12     2012-07-12
################################################################
if (FALSE) {  ## Working Directory: see lib/init.R below [\rd in Vim]
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("../")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}
## Clear memory, load data, functions, etc.  Process data & setup config. settings
source('./CH4-model-fitting/Nfix_setup.R')


################################################################
## MODEL FORMULA
################################################################
Y.main   <- Y.trans ~ Frag + TempC + H2O + I(H2O^2) + logCells + TAN + I(TAN^2)
Y.fixed  <- Y.trans ~ Frag * TempC * H2O * I(H2O^2) * logCells * TAN * I(TAN^2)
Y.random <-  ~ 1 | Block/TempC/Frag 


################################################################
## MODEL FITTING
################################################################
## Transformations & Linearity
##==============================================================
## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)
library(mgcv)
Y.gam <- gam(Y.trans ~ TempC + s(H2O) + s(Cells.m) + s(TAN) + Frag, data = SECCa)
Y.log.gam <- gam(Y.trans ~ TempC + s(H2O) + s(logCells) + s(log10(TAN)) + Frag, data = SECCa)
anova(Y.gam)
anova(Y.log.gam)
AIC(Y.gam, Y.log.gam)
## Cells.m is linear (no need to transform?), TAN might be, H2O is NOT
## I log-transformed Cells in original analysis: I should probably continue to do so, unless I have a really good reason not to.
## The non-linearity in H2O is driven largely by a big gap between very dry and other patches.
## log(TAN) looks a little smoother, but I'm not sure I can justify it biologically?

op <- par(mfrow=c(2,2))
plot(Y.gam, ask = FALSE)
plot(Y.log.gam, ask = FALSE)
par(op)


##==============================================================
## MODEL STRUCTURE (including mixed effects)
##==============================================================
## I will likely need to account for heterogeneity among Blocks, maybe Chambers as well (TempC/Warming)



##==============================================================
## MODEL SELECTION: Which TERMS do I need?
##==============================================================




##==============================================================
## MODEL FITTING
##==============================================================




################################################################
## CHECK ASSUMPTIONS: MODEL VALIDATION
################################################################



################################################################
## GET RESULTS
################################################################



################################################################
## OUTPUT
################################################################
## Save Text Output


##==============================================================
## PUBLICATION GRAPHS
##==============================================================



