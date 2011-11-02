################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### GLMM (regression)
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-11-02
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (MBP).
  getwd()  # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
  source('1_ARA-cb_setup.R')

## library(lattice)    # ggplot2 with faceting is easier!
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme
library(car)                           # diagnostic plots & tools




## Model II Regression: nice idea in theory
################################################################
## Model II Regression?                                       **
################################################################
## Given that both continuous explanatory variables (cyanobacteria, H2O)
## were measured and not controlled, this is most likely a "Model II regression"
##  Model 1 regression *underestimates the slope*, because both variables contain error
##  Model II regression allows X to vary (no "Fixed X" assumption),
##  but different methods for *prediction* (OLS) or
##  describing the *functional relationship* (structural relation, major axis, etc.)
##  - See Sokal & Rohlf p. 541, 545
##  - See Legendre & Legendre p.504
##  Honestly, I think I'm more interested in the latter (functional relationship) than the former (prediction), although I'm not aware of any methods for model 2 GLMMs with multiple variables and nested data, and I'm not sure I have the time to learn/develop them.

library(lmodel2)

Y.m2 <- lmodel2(Y.log ~ X.log, data=SECCa, nperm=999) # Model II regression (Legendre & Legendre)
print(Y.m2)                            # What does this tell me, exactly?
plot(Y.m2)

Y.m2 <- lmodel2(Y.log ~ X.log * H2O, data=SECCa, nperm=999) # Model II regression (Legendre & Legendre)
print(Y.m2)                            # What does this tell me, exactly?
plot(Y.m2)

## lmodel2 seems to drop all but the first explanatory variable!
## I would have to fit models for each pair of variables across all combinations
## of treatment factors?  Ouch.  How would I compare them?

## MA (Major Axis) is also called "geometric mean regression" (Sokal & Rohlf p.544)
## and can easily be computed as the geometric mean of the slopes of y~x and x~y
## <http://tolstoy.newcastle.edu.au/R/help/05/06/5992.html>
## I'm guessing that the MA also passes through the intersection of the y~x & x~y regression lines (the "centroid" of the data points?)

## Legendre & Legendre (p.515) suggest: 
## "MA may also be used with dimensionally heterogeneous variables when the purpose of the analysis is 
## **(1) to compare the slopes of the relationships between the same two variables measured under different conditions (e.g. at two or more sampling sites), or 
##   (2) to test the hypothesis that the major axis does not significantly differ from a value given by hypothesis
## "
## AND
## "If a straight line is not an appropriate model, polynomial or nonlinear regression should be considered."
## 
## see also the User Guide to the 'lmodel2' R package:
## vignette("mod2user", package="lmodel2")

## Therefore:
## In theory, I could use MA for a simple relationship between ARA & Cells,
## for each combination of experimental treatments (Chamber x2, Frag x4, Position x2).
## and compare them using some sort of non-parametric, 
## bootstapped, or permutation-based test statistic,
## while controlling for experiment-wise multiple comparisons.
## OR
## I could just proceed with non-linear regression (GLMM/GAMM), 
## and accept that the regressions parameters (slopes) are likely underestimated
## (I call this "conservative").
## Furthermore, I'm not convinced the relationship *is* linear,
## or that is dimensionally homogeneous or bivariate normal,
## and I want to include multiple fixed and random variables, with a nested structure,
## which would be rather difficult for me to do in a reasonable amount of time
## within a model 2 approach.

