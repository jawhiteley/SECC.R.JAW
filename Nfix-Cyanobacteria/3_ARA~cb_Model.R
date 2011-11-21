################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Modelling: GLMM (regression)
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

  ## Load data, functions, etc.  Process data & setup config. values.  
  ## Includes rm(list=ls()) to clear memory
  source('./Nfix-Cyanobacteria/1_ARA-cb_setup.R')
}


## library(lattice)    # ggplot2 with faceting is easier!
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme
library(rgl)                           # 3D plots
library(car)                           # diagnostic plots & tools
library(gvlma)                         # Global Validation of Linear Model Assumptions
library(nlme)                          # GLMMs (older, but still works)
## library(lme4)    # I'd rather use lme4 for GLMMs, but all my references use nlme
## library\(mgcv)   # Additive Modelling, Cross-Validation

################################################################
## PROCESS DATA: planned
################################################################
## Repeated here for assurance and easy reference
SECCa <- within( SECCa, {
                Y.trans <- Y.log  # convenience
                X.trans <- X.log  # convenience
})
## drop values of X == 0 
## - detection errors where I didn't count any cells 
##   (doesn't mean there were none in the sample)
## Unfortunately, this happens too often: errors in model fitting.  
## May have to drop some variables?
SECC.X0 <- SECCa[SECCa$X.trans != 0, ]
UseClimateFac <- FALSE

################################################################
## ANALYSIS
################################################################
## Goal: Identify which variables have most relative importance,
##       and identify any interactions with experimental treatments.
##==============================================================
## Alternative outline (based on Zuur et al. examples):
## progress from simple to "ultimate" (final) model / approach
## Cycle through:
##  1. Model fitting
##  2. Test Assumptions (graphical analysis of residuals)
##  2a) Identify violated assumptions, and possible solutions
##  3. Model Selection on *valid* models
##  3a) Check Model output, ANOVA tables, parameter values, AIC, R^2, etc.
##  3b) Identify possible improvements
##  END:    choose optimal model; produce model output: tables and graphs.
## Start with GLMM, to incorporate known nested treatments.
## Build in extensions as necessary, time permitting.
##==============================================================
## Methods & Approaches to use
##==============================================================
## *GLMM: fit model using Cells, H2O, and experimental treatments
##  - account for nesting of fixed factors?  How??!
##  - Include spatial autocorrelation instead?
##  - Zero-inflated model to account for exessive 0s in Cells AND ARA?

## H2O
## - include quadratic term to allow for unimodal effect of moisture?
##   - This term was NS in the glm, 
##      but should probably be included anyway for theoretical reasons, 
##      and in case it helps remove patterns in the residuals.
## - Partial regression to remove effect of moisture FIRST, then model the resulting residuals
## - does it affect N-fixation directly, or mediate the effect of cyanobacteria?

## * Subsume Chamber & Position into "Climate" pseudo-treatment?
## * Time as a fixed factor, or separate analysis on each Time?
## Include Temperature values as explanatory variables?
## - Average Temp. during sample collection
## - Avg. Temp. the week before; the month before; since last sample.
## - I don't have Temperature data at the plot resolution ...
## Try other transformations (for residual patterns)
## Try GAM(M)s to see if relationships are really linear?
## Remove data where Cells == 0? (detection failure)
## Bootstrapping to overcome "Fixed X" violations?

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

## Random effects for GLMM
##  Although it's hard to find examples of this type of nesting, 
##  the slash '/' does appear to be the correct operator 
##  for this type of serially nested treatment structure 
##  (based on examples found online).
##  - http://r.789695.n4.nabble.com/lmer-random-factor-nested-in-a-fixed-factor-td865029.html
## The only truly random factor I have is 'Block', 
## which is also the only effective level of replication.  
## Including Block as the highest-level random factor causes 
## nlme to take a *LONG* time to converge.
##   i.e. R hangs for several minutes / hours
## It might be more practical to exclude Block entirely, and just chalk it up to replication.
##  If I do, it might also make the whole 'mixed effects' side of things unnecessary?
## The nesting structure I have is really a violation of independence: 
## patches in the same block are closer together 
## and likely to experience more similar environmental conditions 
## (apart from the chamber effect).
## With Block as largest nested unit: > 12 hours to fit ris model (random intercet + slope), with NaNs in some p-values (Time, Climate, and their interaction)
if (!UseClimateFac) {
  Y.RiN  <- ~ 1|Block/Time/Chamber/Frag                 # Random intercept across Blocks & nested factors
  Y.RisN <- ~ 1 + X.trans * H2O | Block/Time/Chamber/Frag # Random intercept + slopes
} else {
  Y.RiN  <- ~ 1|Block/Time/Climate
  Y.RisN <- ~ 1 + X.trans * H2O | Block/Time/Climate # Random intercept + slopes
}

## This might be the wrong approach: Block is the only real "random" factor.
## Fixed factors happen to be nested within it.
## Maybe what I should really be doing is accounting for heterogeneity across nested fixed factors?
Y.Ri  <- ~ 1 | Block
Y.Ris <- ~ 1 + X.trans * H2O | Block     # Random intercept + slopes

## Block / Chamber / Frag / Position
## - using a correlation structure?
## correlation = corAR1(form = ~ 1 | Block/Time/Chamber/Frag)
## - using variance-covariance structure?
## weights = varIdent(form = ~ 1 | Block/Time/Chamber/Frag)
## - equivalent: correlation = corCompSymm(form = ~ 1 | Block/Time/Chamber/Frag)




if (FALSE) {    # GLM: deprecated. wrapped to allow folding
##################################################
## BASIC GLM - soon to be deprecated by GLMM (below)
##################################################
##================================================
## Model Fitting
##================================================
### "basic" GLM: initial foray into analysis.  Soon to be deprecated by GLMM.
Y.model <- glm( Y.fixed, data=SECCa, family="gaussian" )
Y.model.full <- Y.model
Y.model.main <- glm(Y.trans ~ X.trans + Time + Chamber + Frag + Position + H2O + I(H2O^2), data = SECCa)  # Main factors only; no interaction terms.

# Should I not be using a mixed-effects model to account for treatments of different sizes? Chamber / Frag / Position
# YES (See below)

##================================================
## MODEL SELECTION
##================================================
drop1(Y.model)
Y.model.selected <- step(Y.model, direction = "backward")
Y.model.main.selected <- step(Y.model.main, direction = "backward")

Y.model  <- Y.model.selected


##################################################
## CHECK ASSUMPTIONS (MODEL VALIDATION)
##################################################
## analyse residuals, standard diagnostic plots
op <- par(mfrow=c(2,2))	 # panel of figures: 3 rows & 2 columns
plot(Y.model)
par(op)

## Residuals
Model.resid <- resid(Y.model)
## par(mfrow=c(1,1))	 # panel of figures: 1 rows & 1 columns
hist(Model.resid)    # plot residuals
plot(SECCa$X, Model.resid)
plot(SECCa$H2O, Model.resid)
qplot(Frag, Model.resid, data = SECCa, facets = Chamber * Position ~ Time ) + theme_bw()

## Plot Residuals: see Zuur et al. 2007, pg. 61-63
coplot( Model.resid ~ X | Frag * Position, data=SECCa, 
       pch=SECCa$pt, col=SECCa$colr,
       ylab = "Residuals"
)
## much the same breakdown, using TRELLIS xyplot over 3 factors.
print( xyplot( Model.resid ~ X | Chamber * Frag * Position, data=SECCa, 
       pch=point, col=SECCa$colr, bg = SECCa$fill,
       ylab = "Residuals"
) )

## Plot Residuals: see Zuur et al. 2007, pg. 131-133
plot(Y.model$fitted[,2],resid(Y.model,type="p"),xlab="Fitted values",ylab="Residuals")
qqnorm(resid(Y.model,type="p"),cex=1,main="")



##################################################
## ANALYSIS: GET RESULTS
##################################################
anova(Y.model)
summary(Y.model)

}




################################################################
## LMM - Hierarchical / Multilevel Mixed Model                **
################################################################
## Control options for lme
## - the default optimizer since R 2.2.x is 'nlminb', which sometimes causes false convergence errors.  'optim' sometimes avoids the issue, but may change results?
##   <http://r.789695.n4.nabble.com/unexpected-quot-false-convergence-quot-td790625.html>
## - adding more iterations allows for convergence for complex models, but takes longer.
## - returnObject=TRUE may also avoid false convergence errors
##   <http://r.789695.n4.nabble.com/quot-False-convergence-quot-in-LME-td860675.html>
lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 500, msMaxIter = 100, opt="optim")

##==============================================================
## Model Fitting
##==============================================================
## fixed effects only for assessment (How much variation is explained?
Y.lm <- lm( Y.fixed , data=SECCa)
summary(Y.lm)
Y.lm <- lm( Y.trans ~ X.trans + H2O + Block + Time * Chamber * Frag * Position , data=SECCa)
summary(Y.lm)

## Using a mixed-effects model to account for treatments of different sizes? 
##  Block / Chamber / Frag / Position
Y.fm <- gls( Y.fixed, data=SECCa, method="REML") # fixed effects only for comparison
Y.rim  <- lme(Y.fixed, random=Y.Ri,  data=SECCa, control=lmd, method="REML")
## SLOW! :(
## more complex models take a few minutes (or hours) to fit.  Use with discretion.
Y.rism <- lme(Y.fixed, random=Y.Ris, data=SECCa, control=lmc, method="REML")
## Allow different error variances within treatment groups (see Zuur et al. 2009, pg. 188)
##   too many groups causes lme errors: "false convergence (8)"; see notes on lmeControl, above
## Nesting causes Error in getGroups.data.frame(...): Invalid formula for groups
Y.rie  <- lme(Y.fixed, random=Y.Ri,  
              weights=varIdent(form=~ 1 | Block * Time), # Block * Time?
              data=SECCa, control=lmc, method="REML")
Y.rise <- lme(Y.fixed, random=Y.Ris,  
              weights=varIdent(form=~ 1 | Time),
              data=SECCa, control=lmc, method="REML")
## Nesting with correlation structure?  Not sure an auto-regression model is the most appropriate
if (UseClimateFac) {
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block),
                 correlation=corAR1(form = ~ 1 | Block/Time/Frag/Climate),
                 data=SECCa, control=lmc, method="REML")
} else {
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block),
                 correlation=corAR1(form = ~ 1 | Block/Time/Chamber/Frag),
                 data=SECCa, control=lmc, method="REML")
}

if (FALSE) {                           # these attempts are not working :(
  ## Error in getGroups.data.frame(data, form, level = length(splitFormula(grpForm,  : Invalid formula for groups
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block/Time/Chamber/Frag),
                 data=SECCa, control=lmc, method="REML")
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block/Time/Frag/Climate),
                 data=SECCa, control=lmc, method="REML")

  ## I want to use varConstPower (variance covariate with values of 0), but that produces an error (false convergence) :(
  Y.rieXP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varConstPower(form = ~ X.trans),
                 data=SECCa, control=lmd, method="REML")

  ## varFixed gives no errors, but takes >12 hours to fit (if at all)
  Y.rieXP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varFixed(~ X.trans),
                 data=SECCa, control=lmd, method="REML")
  ## varPower works, but not for X.trans, which has 0's!   :/
  Y.rieXP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varPower(form = ~ X.trans),
                 data=SECCa, control=lmd, method="REML")
  Y.rieHP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varPower(form = ~ H2O),
                 data=SECCa, control=lmd, method="REML")

  Y.riceX <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ X.trans)),
                 data=SECCa, control=lmc, method="REML")
  Y.riceH <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ H2O)),
                 data=SECCa, control=lmc, method="REML")
  Y.rice  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ X.trans),
                                 varConstPower(form=~ H2O)),
                 data=SECCa, control=lmc, method="REML")
}

## Main effects only: ignore interactions?
Y.mainMM <- lme(Y.main, data = SECCa, random = Y.Ri, method="REML")
Y.mainML <- update(Y.mainMM, method="ML")


## Spatial autocorrelation?
Y.ce <- update(Y.fm, correlation = corExp(form=~ xE + yN, nugget = TRUE) )
Y.cg <- update(Y.fm, correlation = corGaus(form=~ xE + yN, nugget = TRUE) )
Y.cs <- update(Y.fm, correlation = corSpher(form=~ xE + yN, nugget = TRUE) )
## Y.cl <- update(Y.fm, correlation = corLin(form=~ xE + yN, nugget = TRUE) )




##==============================================================
## MODEL SELECTION
##==============================================================
## RANDOM structure
anova(Y.fm, Y.ce, Y.cg, Y.cs)          # does spatial autocorrelation improve the model?
anova(Y.fm, Y.rim)                     # do random effects improve the model?
anova(Y.fm, Y.rism)                    # do random effects improve the model?
anova(Y.fm, Y.rie)                     # do random effects improve the model?
anova(Y.fm, Y.rise)                    # do random effects improve the model?
anova(Y.rieN, Y.rie, Y.rim, Y.rism, Y.rise) # do we need random slopes or error terms?
anova(Y.rieN, Y.rie, Y.rim)                 # do we need nested error terms?

Y.mm <- Y.rim                          # Optimal random structure
Y.mm <- update(Y.rim, correlation = corExp(form=~ xE + yN, nugget = TRUE) ) # Optimal random structure
## The biggest improvement seems to come from random intercepts across Blocks.
## However, a model with only this random effect shows major heterogeneity in the residuals.
## I may need other random factors to maintain a valid model fit.
## Some models with different error structures can't be refit with ML :(
##   despite having lower AICs (e.g Y.rie)

## optimize FIXED factors
if (FALSE) {
  Y.ml  <- lme(Y.fixed, data=SECCa, random=Y.Ri, method="ML") # re-fit with ML
  Y.rieML <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block),
                 data=SECCa, control=lmc, method="ML")
}
Y.ml  <- update(Y.mm, method="ML")     # re-fit with ML; some models produce errors :(

drop1(Y.ml)                            # not encouraging
## Y.step <- step(Y.ml)                # stepwise back & forward model selection?  Not for lme
if (!UseClimateFac) {                  # All factors, or Climate pseudo-factor
  ## model selection for Y.rim **
  Y.ml1 <- update(Y.ml, .~. - X.trans:H2O:Time:Chamber:Frag:Position) # * sig. WORSE!
  anova(Y.ml, Y.ml1)

  ## Main Effects only: No interactions?
  Y.main1 <- update(Y.mainML, .~. - X.trans)
  Y.main2 <- update(Y.mainML, .~. - H2O)
  Y.main3 <- update(Y.mainML, .~. - Time)
  Y.main4 <- update(Y.mainML, .~. - Chamber)
  Y.main5 <- update(Y.mainML, .~. - Frag)
  Y.main6 <- update(Y.mainML, .~. - Position)
  anova(Y.mainML, Y.main1)
  anova(Y.mainML, Y.main2)
  anova(Y.mainML, Y.main3)
  anova(Y.mainML, Y.main4)
  anova(Y.mainML, Y.main5)
  anova(Y.mainML, Y.main6)
  ## Full model is still the best :(
  ## Forward model selection to add interactions?
  Y.m1     <- update(Y.mainML, .~. + X.trans:H2O)
  Y.m2     <- update(Y.mainML, .~. + Chamber:Position)         # *
  Y.m3     <- update(Y.mainML, .~. + Chamber:Frag)
  Y.m4     <- update(Y.mainML, .~. + Frag:Position)
  Y.m5     <- update(Y.mainML, .~. + Time:X.trans)
  Y.m6     <- update(Y.mainML, .~. + Time:H2O)                 # **
  Y.m7     <- update(Y.mainML, .~. + X.trans:Chamber)
  Y.m8     <- update(Y.mainML, .~. + H2O:Chamber)              # *
  Y.m1.1   <- update(Y.mainML, .~. + X.trans:H2O
                     + X.trans:Time + H2O:Time 
                     + X.trans:H2O:Time)         # *
  Y.m1.1.1 <- update(Y.mainML, .~. + X.trans:H2O + X.trans:Time + H2O:Time 
                     + X.trans:Chamber + H2O:Chamber + Time:Chamber
                     + X.trans:H2O:Time + X.trans:H2O:Chamber 
                     + X.trans:Time:Chamber + H2O:Time:Chamber
                     + X.trans:H2O:Time:Chamber) # *
  anova(Y.mainML, Y.m1)
  anova(Y.mainML, Y.m2)
  anova(Y.mainML, Y.m3)
  anova(Y.mainML, Y.m4)
  anova(Y.mainML, Y.m5)
  anova(Y.mainML, Y.m6)
  anova(Y.mainML, Y.m7)
  anova(Y.mainML, Y.m8)
  anova(Y.mainML, Y.m1.1)
  anova(Y.mainML, Y.m1.1.1)

  Y.mainMM <- update(Y.mainML, method="REML")

  AIC(Y.ml, Y.mainML)
  AIC(Y.mm, Y.mainMM)

} else {                               # Chamber & position lumped into Climate pseudo-factor
  Y.ml1 <- update(Y.ml, .~. - X.trans:H2O:Time:Climate:Frag)
  anova(Y.ml, Y.ml1)
  ## drop 4-way interactions?
  Y.ml2 <- update(Y.ml1, .~. - X.trans:H2O:Time:Climate)  # *
  Y.ml3 <- update(Y.ml1, .~. - X.trans:H2O:Time:Frag)
  Y.ml4 <- update(Y.ml1, .~. - X.trans:H2O:Climate:Frag)  # <-
  Y.ml5 <- update(Y.ml1, .~. - X.trans:Time:Climate:Frag) # **
  Y.ml6 <- update(Y.ml1, .~. - H2O:Time:Climate:Frag)
  anova(Y.ml1, Y.ml2)
  anova(Y.ml1, Y.ml3)
  anova(Y.ml1, Y.ml4)
  anova(Y.ml1, Y.ml5)
  anova(Y.ml1, Y.ml6)
  ## dropped least significant 4-way interaction: next
  Y.ml2 <- update(Y.ml4, .~. - X.trans:H2O:Time:Climate)  # *
  Y.ml3 <- update(Y.ml4, .~. - X.trans:H2O:Time:Frag)     # <-
  Y.ml5 <- update(Y.ml4, .~. - X.trans:Time:Climate:Frag) # **
  Y.ml6 <- update(Y.ml4, .~. - H2O:Time:Climate:Frag)
  anova(Y.ml4, Y.ml2)
  anova(Y.ml4, Y.ml3)
  anova(Y.ml4, Y.ml5)
  anova(Y.ml4, Y.ml6)
  ## dropped least significant 4-way interaction: next
  Y.ml2 <- update(Y.ml3, .~. - X.trans:H2O:Time:Climate)  # *
  Y.ml5 <- update(Y.ml3, .~. - X.trans:Time:Climate:Frag) # *
  Y.ml6 <- update(Y.ml3, .~. - H2O:Time:Climate:Frag)     # <-
  anova(Y.ml3, Y.ml2)
  anova(Y.ml3, Y.ml5)
  anova(Y.ml3, Y.ml6)
  ## dropped least significant 4-way interaction: next
  Y.ml2 <- update(Y.ml6, .~. - X.trans:H2O:Time:Climate)  # *
  Y.ml5 <- update(Y.ml6, .~. - X.trans:Time:Climate:Frag) # *
  anova(Y.ml6, Y.ml2)
  anova(Y.ml6, Y.ml5)
  Y.ml4 <- Y.ml6                       # optimal model with 4-way interactions
  anova(Y.ml, Y.ml4)
  Y.mm <- update(Y.ml4, method="REML") # re-fit with REML

  ## 3-way interactions?

  ## Ignore interactions?
  Y.main1 <- update(Y.mainML, .~. - X.trans)
  Y.main2 <- update(Y.mainML, .~. - H2O)
  Y.main3 <- update(Y.mainML, .~. - Time)
  Y.main4 <- update(Y.mainML, .~. - Climate)
  Y.main5 <- update(Y.mainML, .~. - Frag) # <-
  anova(Y.mainML, Y.main1)
  anova(Y.mainML, Y.main2)
  anova(Y.mainML, Y.main3)
  anova(Y.mainML, Y.main4)
  anova(Y.mainML, Y.main5)
  Y.mainML <- Y.main5
  Y.mainMM <- update(Y.mainML, method="REML")

  AIC(Y.ml, Y.ml4, Y.mainML)
  AIC(Y.rim, Y.mm, Y.mainMM)
}
##   Y.mm <- Y.mainMM
## The model with only main effects actually shows more Normal residuals, and less heterogeneity!


################################################################
## CHECK ASSUMPTIONS: MODEL VALIDATION
################################################################
## see Zuur et al. 2007, 2009 books, 
##  and [Quick-R](http://www.statmethods.net/stats/rdiagnostics.html)
## Regression (ANOVA): Normality; Homogeneity; Independence; Fixed X*
## Fixed X (Model I): *think about how X data was generated*
##                    - fixed values [I] or large mesurement error [II] ?
## Check Normality: histogram, qqnorm of Residuals
## Check Homogeneity: (standardized) Residuals vs. Fitted / vs. X
## Check Independence: (standardized) Residuals vs. X
## Assess Model fit, specification: Look for patterns in graphs of Residuals
## - Should be no areas of residuals with similar values, gaps in the cloud, etc.
## Residuals should ideally be spread out equally across all graphs (vs. X / Fitted).

diagnostics(Y.lm, more=TRUE)           # Simple linear model, with a few interactions
diagnostics(Y.fm)                      # No random effects

diagnostics(Y.rim)                     # Random Intercept
diagnostics(Y.rism)                    # Random Intercept + slope
diagnostics(Y.rise)                    # Random Intercept + slope
diagnostics(Y.rie, resType="n")        # Do more random effects help?
diagnostics(Y.rieN, resType="n")       # Do more random effects help?
diagnostics(Y.mm)                      # Optimal mixed model (after model selection)
diagnostics(Y.mainMM)                  # Main effects only: tend to violate fewer assumptions :/
diagnostics(Y.m1.1.1)                  # Main effects + interactions: better or worse?

Y.model <- Y.mm
RE <- diagnostics(Y.model, resType="normalized", more=TRUE) # full diagnostics, just to be sure

## A Model with Main Effects only has the most valid fit (despite a low AIC): 
## Normally-distributed residuals, low heterogeneity (not perfect, though)
## BUT, strong patterns in the residuals vs. fitted values (negative trend)
## Models with random terms tend to have better AICs, but also violate more assumptions, especially Normality.
## Allowing different variances per Block or Time reduces patterns in the residuals, but is also highly non-Normal, and there is still heterogeneity in the residuals across H2O

## Check for spatial patterns in residuals?
library (gstat)
RE.df <- data.frame(RE=resid(Y.model, type="normalized"), x=SECCa$xE, y=SECCa$yN) # I hope the order is the same!
coordinates(RE.df) <- c('x', 'y')
Y.bubble <- bubble(RE.df, "RE", col=c("black", "grey"), # cex=0.1,
                   main = attr(RE, "label"),
                   xlab = attr(SECC.xy, "labels")$xE,
                   ylab = attr(SECC.xy, "labels")$yN)
print(Y.bubble)                        # Evidence of spatial patterns in residuals?
## A little hard to tell how prominent the spatial patterns are.
## Patches are so close together, it's difficult to separate the large layered bubbles.
## This proximity may be more justification for including spatial autocorrelation, however.
## Variogram (Zuur et al. 2009 Ch. 7.2, pg. 167-169)
Y.vario <- variogram(RE ~ 1, RE.df)
plot(Y.vario)
Y.variogls <- Variogram(Y.fm, form=~ xE + yN,
                        robust=TRUE, maxDist=200,
                        resType="pearson")
plot(Y.variogls, smooth=TRUE)
Y.varioMM <- Variogram(Y.model, form=~ xE + yN,
                        robust=TRUE, maxDist=200,
                        resType="pearson")
plot(Y.varioMM, smooth=TRUE)

## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)



################################################################
## ANALYSIS: GET RESULTS
################################################################
anova(Y.model)
summary(Y.model)

## intervals(Fitted.Model.Object)   # Approx. Confidence intervals.  see ?intervals
## see multcomp package for multiple comparisons (Tukey's HSD on mixed effects model?)


##==============================================================
## Variance Components Analysis (variance decomposition)?
## No R^2 (ML, not LS methods), but maybe some analogs?
## - https://stat.ethz.ch/pipermail/r-help/2003-March/031162.html
## - http://www.biostat.wustl.edu/archives/html/s-news/2002-04/msg00075.html
## - http://www.biostat.wustl.edu/archives/html/s-news/2002-04/msg00076.html
## VarCorr() and then calculate percentage variance for each component from output of VarCorr 
## - http://r.789695.n4.nabble.com/Components-of-variance-with-lme-td3329753.html
## - http://r.789695.n4.nabble.com/extract-variance-components-td866971.html
## ?VarCorr
## ?getVarCov
Y.varcor <- VarCorr(Y.model)           # random components only :(






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

## generate grid to add predicted values to (X-values in all combinations of factors).
## - watch length.out: if it's too long, R will choke.
Y.pred <- expand.grid(Block    = levels(SECCa$Block) , 
                      Time     = levels(SECCa$Time) , 
                      Chamber  = levels(SECCa$Chamber) , 
                      Frag     = levels(SECCa$Frag), 
                      Position = levels(SECCa$Position), 
                      X.trans  =seq(0, max(SECCa$X.trans), length.out=20 ),
                      H2O      =seq(0, max(SECCa$H2O), length.out=20 ) 
                      )
Y.pred$predicted <- predict(Y.model, newdata=Y.pred, type="response" )  # newdata must have same explanatory variable names for predict to work.

if (FALSE) {
  pred.Chamber <- expand.grid(Chamber = levels(SECCa$Chamber) , 
                              X=seq(0, max(SECCa$X), length.out=100 ) 
  )
  pred.Chamber$predicted <- predict(Y.model, newdata=pred.Chamber, type="response" )
  pred.Frag <- expand.grid(Frag = levels(SECCa$Frag) , 
                           X=seq(0, max(SECCa$X), length.out=100 )
  )
  pred.Frag$predicted <- predict(Y.model, newdata=pred.Frag, type="response" )
  pred.Position <- expand.grid(Position = levels(SECCa$Position) , 
                               X=seq(0, max(SECCa$X), length.out=100 ) 
  )
  pred.Position$predicted <- predict(Y.model, newdata=pred.Position, type="response" )
  pred.FxP <- expand.grid(Frag     = levels(SECCa$Frag), 
                          Position = levels(SECCa$Position), 
                          X=seq(0, max(SECCa$X), length.out=100 ) 
                          )
  pred.FxP$predicted <- predict(Y.model, newdata=pred.Position, type="response" )
}


Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
point <- 21	# 21 for circles with col & bg ; 16 for solid circles
Chamber.map$pch <- c(21, 16)  # use circles for both treatments

SECCa <- within( SECCa,{
	colr = as.character(ifelse(Chamber == Chamber.map$label[1], 
                               Chamber.map$col[1], 
                               Chamber.map$col[2] 
                               )
    )
	fill = as.character(ifelse(Chamber == Chamber.map$label[1], 
                               Chamber.map$bg[1],
                               Chamber.map$bg[2]
                               )
    )
	pt = ifelse(Chamber == Chamber.map$label[1], 
			Chamber.map$pch[1], 
			Chamber.map$pch[2]
		)
})


par(mfrow=c(1,1))
pred.Y <- with( Y.pred, 
               aggregate(cbind(predicted), list(Chamber = Chamber, X = X.trans), mean)
)  # I should be getting direct predictions, not means of predictions. *****
with( SECCa,{
	# pred | augpred | ?
	plot(X, Y, type="p",
		ylab=Y.plotlab, xlab=X.plotlab,
		pch=pt, col=colr, bg=fill
        )
	lines(predicted ~ X, data=subset(pred.Y, Chamber == "Ambient"), 
          col = Chamber.map$col[1], 
          lty = Chamber.map$lty[1]
          )
	lines(predicted ~ X, data=subset(pred.Y, Chamber == "Full Chamber"), 
          col = as.character(Chamber.map$col[2]), 
          lty = Chamber.map$lty[2]
          )
	legend( "topright", legend=Chamber.map$label, pch=point, col=as.character(Chamber.map$col), pt.bg=as.character(Chamber.map$bg) )
})



##==============================================================
## Plot fitted on observed, by factor?
##==============================================================
library(lattice)
## lattice panels
print( xyplot( Y ~ X | Frag * Position , data=SECCa, 
              pch = SECCa$pt, col = SECCa$colr, 
              xlab = quote(X.plotlab), ylab = quote(Y.plotlab), 
              panel = function(..., data, subscripts) {
                panel.xyplot(...)  # regular plot of data points
                Frag.lvl <- unique(SECCa$Frag[subscripts]) # get current factor levels
                Pos.lvl  <- unique(SECCa$Position[subscripts])
                preds    <- Y.pred[which(Y.pred$Frag %in% Frag.lvl 
                                       & Y.pred$Position %in% Pos.lvl), ]
##                  browser()
                for( lvl in levels(preds$Chamber) ) {
                  preds.lvl <- subset(preds, preds$Chamber == lvl)
                  panel.xyplot(preds.lvl$X, preds.lvl$predicted, 
                               type = 'l', 
                               col = Chamber.map$col[Chamber.map$label == lvl]
                               )
                }
              },
              subscripts = T
              )
)




if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
