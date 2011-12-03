################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Modelling: GLMM (regression)
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-11-23
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (MBP).
  setwd("..")  # relative to this file (\rd in Vim-R)
  getwd()  # Check that we're in the right place

  ## Load data, functions, etc.  Process data & setup config. values.  
  ## Includes rm(list=ls()) to clear memory
  source('./Nfix-Cyanobacteria/1_ARA-cb_setup.R')
}


library(car)
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme

################################################################
## PROCESS DATA: planned
################################################################
## Repeated here for assurance and easy reference
SECCa <- within( SECCa, {
                Y.trans <- Y.log  # convenience: log10
                X.trans <- X.log  # convenience: log10
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
## Start simple and add complexity, as necessary or appropriate

## Main effects only
Y.main  <- Y.trans ~ X.trans + H2O + I(H2O^2) + Block + Time + Chamber + Frag + Position
Y.mainCl<- Y.trans ~ X.trans + H2O + I(H2O^2) + Block + Time + Climate + Frag
### Fixed effects: All?  No Replication for all combinations!
Y.fixed <- Y.trans ~ X.trans * H2O * I(H2O^2) * Block * Time * Chamber * Frag * Position
Y.fixCl <- Y.trans ~ X.trans * H2O * I(H2O^2) * Block * Time * Climate * Frag
### Main Effects + Selected Interactions (the most important / significant ones)
Y.form  <- Y.trans ~ X.trans * H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag
Y.formCl<- Y.trans ~ X.trans * H2O * I(H2O^2) + Block * Time * Climate + Frag




##==============================================================
## GAM: Is the relationship linear or not?
##==============================================================
## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)
library(mgcv)
ARA.gam <- gam(Y ~ s(X) + s(H2O), data=SECCa)
ARA.gam1 <- gam(Y.log ~ s(X) + s(H2O) + Block * Time*Chamber*Position + Frag, data=SECCa)
ARA.loggam <- gam(Y.log ~ s(X.log) + s(H2O), data=SECCa)
ARA.loggamF <- gam(Y.log ~ s(X.log) + s(H2O) + Block + Time*Chamber*Position * Frag, data=SECCa)
ARA.loggam1 <- gam(Y.log ~ s(X.log) + s(H2O) + Block * Time*Chamber*Position + Frag, data=SECCa)
## Frag-interactions are less significant than Block-interactions
AIC(ARA.gam, ARA.gam1, ARA.loggam, ARA.loggamF, ARA.loggam1)
## It *might* be linear for log(Cells) (X.log; edf = 2.18), but definitely NOT for H2O (edf=4)
anova(ARA.gam1)
anova(ARA.loggamF)
anova(ARA.loggam1)
## The main non-linear portion for Cells is between 0-values and the rest: non-zero values are rather linear.
op <- par(mfrow=c(2,2))
plot(ARA.gam)
plot(ARA.gam1)
plot(ARA.loggam1)
## polynomial for moisture (H2O)?
ARA.gam2 <- gam(Y.log ~ s(X.log) + s(poly(H2O, 2)) + Block + Time*Chamber*Frag*Position, data=SECCa)
anova(ARA.gam2)
## what about sqrt-transforming Cell Density (on a whim)
ARA.Xsqrt.gam <- gam(Y ~ s(sqrt(X)) + s(H2O) + Block * Time*Chamber*Position + Frag, data=SECCa)
ARA.logsqrt.gam <- gam(Y.log ~ s(sqrt(X)) + s(H2O) + Block * Time*Chamber*Position + Frag, data=SECCa)
anova(ARA.Xsqrt.gam)
anova(ARA.logsqrt.gam)
## Definitely linear for sqrt(Cell Density)!! wtf?
plot(ARA.logsqrt.gam)

par(op)                                # restore original settings



##==============================================================
## Model Fitting: What Transformations do I need?
##==============================================================
## fixed effects only for assessment (How much variation is explained?)
Y.lmain <- lm( Y.main , data=SECCa)
summary(Y.lmain)
anova(Y.lmain)                          # ORDER MATTERS! (see Zuur et al. 2009 pg. 540))
Y.fmain <- lm( Y.fixed , data=SECCa)   # no replication of all interaction combinations
anova.full <- anova(Y.fmain)
anova.full$PropVar <- anova.full[, "Sum Sq"] / sum(anova.full[, "Sum Sq"])
anova.full <- anova.full[order(anova.full$PropVar, decreasing=TRUE), ] # sort by Prop. Var
print(sum(anova.full$PropVar))
print(anova.full)

## Reasonable model of effects & Interactions for assessment (thanks John Connelly)
Y.lm <- lm( Y.trans ~ X.trans * H2O + Block + Time * Chamber * Frag * Position , data=SECCa)
summary(Y.lm)                          # R2 = 0.59
anova(Y.lm)
## Frag interactions are not significant.  Let's look at interactions with Blocks instead:
Y.lm1  <- lm( Y.trans ~ X.trans * H2O + Block * Time * Chamber * Position + Frag + Block:Time, data=SECCa)
Y.lmH  <- lm( Y.trans ~ X.trans * H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag, data=SECCa)
Y.lmHB <- lm( Y.trans ~ X.trans * H2O * I(H2O^2) * Block + Time * Chamber * Position + Frag + 
              Block:Time + Chamber:X.trans, data=SECCa)
Y.lmHC <- lm( Y.trans ~ X.trans * H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag + 
              Chamber:X.trans, data=SECCa)
summary(Y.lm1)                         # Block * Time ... R2 = 0.75 !!
summary(Y.lmH)                         # Block * Time ... R2 = 0.77 !!
summary(Y.lmHB)                        # Block * Cells... R2 = 0.75
summary(Y.lmHC)                        # Block * Time + Chamber:Cells ... R2 = 0.77 !
anova(Y.lmH, Y.lmHC)                   # Adding Chamber:Cells interaction is not a sig. improvement
anova(Y.lmH)                          # ORDER MATTERS! (see Zuur et al. 2009 pg. 540)
Anova(Y.lmH)
AIC(Y.lm, Y.lm1, Y.lmH, Y.lmHC, Y.lmHB) # Block interactions are a better fit than Frag interactions (more degrees of freedom?)

## Transformation?
Y.rawlm   <- lm( Y ~ X * H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag , data=SECCa)
summary(Y.rawlm)
Y.translm <- lm( Y.trans ~ X * H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag , data=SECCa)
summary(Y.translm)                     # better fit than log-log transformation?
if (F) { ## Transformation via link function in glm? (invalid Y values causes error)
  Y.logglm <- glm( Y ~ X * H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag , 
                  data=SECCa[SECCa$Y>0, ], family=gaussian(link="log"))
  summary(Y.logglm)
}
Y.lmtrans <- lm( Y ~ X.trans * H2O + Block * Time * Chamber * Position + Frag , data=SECCa)
## incorporating sqrt(Cells), quadratic H2O term, and better interactions ****
Y.lmsqrt <- lm( Y.log ~ I(sqrt(X)) * H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag , data=SECCa)

AIC(Y.lmain, Y.lm, Y.lm1, Y.lmH, Y.rawlm, Y.lmtrans, Y.translm, Y.lmsqrt)
summary(Y.lmsqrt)                         # The best model, according to AIC**
## AIC suggests a log-transformation on only the response (ARA/N-fixation) is a good fit.
## However, this implies an *exponential*, rather than saturating log- relationship
## between N-fixation and Cell abundance, which we predicted / hypothesised.
## Granted, the AIC improvement is about 1-2% (745 -> 729; 721 -> 729)
## and the R^2 is only marginally higher (0.769 vs. 0.764)
## Both are substantially better than the model with untransformed (response) variables.
## It is probably more defensible to use the log-log model, which makes more sense theoretically.
## There is so much noise in the data that the marginally better log-Y fit may be spurious or inconsequential.
## The model with more interactions seems to be a *poorer* fit than main effects only, 
##   according to AIC?  But the R^2 is a little higher (0.59 vs. 0.57)
## Correction: depends on which interactions are included ;-) (R^2= 0.769 vs. 0.57)

Y.model  <- Y.lmH



##==============================================================
## MODEL SELECTION: Which TERMS do I need in the model?
##==============================================================
## multi-model averaging and inference with glmulti!
library(glmulti)                       #  v1.0, april 2011; v1.0.3, Nov 2011
## library(MASS)                          # ?
## library(leaps)
## Note: glmulti performs an exhaustive search of all candidate models.  
##       2^7 = 128 candidate models (not including interactions). 256 according to glmulti
## method="l" (leaps) uses a *much* faster algorithm, but not with factors or interactions :(
ARA.glmulti1 <- glmulti(Y.fixed, data=SECCa, crit=aic, level=1, fitfunc=lm, marginality = TRUE,
                        method="h", confsetsize=256, plotty=FALSE, report=FALSE)
print(ARA.glmulti1)

if (T) {
  ## level=2 for 2-way interactions
  ## 2^(7 + choose(7,2) ) = 268,435,456 candidate models with 2-way interactions!!
  ## 416,869,200 according to glmulti (without marginality); 286,192,513 with marginality.
  ## An exhaustive exploration would take ~ 1 month on my computer.  
  ## Try with 'genetic' algorithm to speed things up (method="g"). ~ 30 minutes
  ## Best to run 2-4+ replicate genetic algorithms, and take consensus.
  ## use method="d" to print a summary of candidate models (no fitting)
  ARA.multi <- list()
  for (i in 1:4) {
    cat("\n================ glmulti: Genetic Algorithm Run", i, "================\n\n")
    ARA.multi[[i]] <- glmulti(Y.fixed, data=SECCa, crit=aic, level=2, fitfunc=lm, 
                              marginality = TRUE, method="g", confsetsize=256, 
                              plotty=FALSE, report=TRUE)
  }
  ARA.glmulti2 <- consensus(ARA.multi, confsetsize=256) # more models use more memory!
  rm(ARA.multi)                        # clean-up
  ## save object to speed up loading for future analysis.  This is still a big file: >1 GB!
  save(ARA.glmulti1, ARA.glmulti2, file="./save/ARA-cb.glmulti.R")
} else {
  ## load saved object to speed things up
  load("./save/ARA-cb.glmulti.R")
}
print(ARA.glmulti2)
## summary(ARA.glmulti2)

ARA.best2 <- as.formula(summary(ARA.glmulti2)$bestmodel)
ARA.best2lm <- lm(ARA.best2, data=SECCa)
summary(ARA.best2lm)


getCoef.glmulti <- function(glmObj, minImportance=0) {
  glm.coef <- as.data.frame(coef(glmObj))
  ## from coef: these are weights of model *coefficients*, NOT *model terms*...
  ## I think I want weights of model terms (variables, rather than levels of each factor)
  ## ggplot will order the bars by levels of the explanatory factor
  glm.coef$Term <- factor(row.names(glm.coef), levels=unique(row.names(glm.coef))) 
  glm.order <- order(glm.coef$Importance, decreasing=TRUE)
  ## glm.coef$Term <- factor(glm.coef$Term, levels=levels(glm.coef$Term)[glm.order])
  ## Drop terms below the threshold
  glm.minImp <- which(glm.coef$Importance >= minImportance)
  glm.coef <- glm.coef[glm.minImp, ]
  ## facilitate confidence intervals on Estimates
  glm.coef$Emax <- glm.coef$Estimate + glm.coef[, "+/- (alpha=0.05)"]
  glm.coef$Emin <- glm.coef$Estimate - glm.coef[, "+/- (alpha=0.05)"]
  ## Clean up labels
  levels(glm.coef$Term) <- gsub("X.trans", attr(SECC, "labels")[X.col], levels(glm.coef$Term))
  levels(glm.coef$Term) <- gsub("ChamberFull Chamber", "Chamber", levels(glm.coef$Term))
  levels(glm.coef$Term) <- gsub("PositionOuter", "Position", levels(glm.coef$Term))
  levels(glm.coef$Term) <- gsub("Frag(.*)", "Frag (\\1)", levels(glm.coef$Term))
  levels(glm.coef$Term) <- gsub("Time(.*)", "Time (\\1)", levels(glm.coef$Term))
  glm.coef
}
if (F) {
  names(ARA.glmulti2)
  str(summary(ARA.glmulti2))
  summary(ARA.glmulti1)$modelweights
  weightable(ARA.glmulti1)             # models, IC values, and relative weights (for confset)
}

ARA.coef1 <- getCoef.glmulti(ARA.glmulti1)
ARA.coef2 <- getCoef.glmulti(ARA.glmulti2)
nrow(ARA.coef2)                        # 129 terms (128 + intercept)!
## defaults
ARA.glmulti <- ARA.glmulti1
ARA.coef <- ARA.coef1

## Output graphs
par(mfrow=c(1,1))
plot(ARA.glmulti1, type="p")           # red line at ~2 IC units above best model: should consider at least all models below this line.
plot(ARA.glmulti1, type="w")           # red line where cumulative evidence weight = 95%
## plot(ARA.glmulti1, type="r")           # diagnostic plots (windows only?)
plot(ARA.glmulti1, type="s")           # term weights (v 1.0.3)
plot(ARA.glmulti2, type="s")           # term weights (v 1.0.3)
## sum of relative evidence weights over all models that include each term?
barplot(ARA.coef1[, "Importance"], horiz=TRUE, names.arg=ARA.coef1$Term, las=2) 
plot(ARA.glmulti2, type="p")
plot(ARA.glmulti2, type="w")

## ggplot2: theme settings
bar.import <- function(glmObj) {
  ggplot(glmObj, aes(x=Term, y=Importance), 
                          stat="identity", xlab = "Model Terms") +
         list(geom_bar(), coord_flip(), scale_y_continuous(expand=c(0,0)),
              opts(panel.border=theme_blank(), axis.line=theme_segment(),
                   plot.title=theme_text(size = 16, face = "bold"),
                   title="Model-averaged importance of effects")
         )
}
est.confint <- function(glmObj) {
  conf.wd <- glmObj[, "+/- (alpha=0.05)"]
  conf.int <- aes(ymax = Estimate + conf.wd, ymin = Estimate - conf.wd)
  ggplot(glmObj, aes(x=Term, y=Estimate) ) +
  geom_hline(yintercept=0, colour="grey") + 
  list(geom_point(),
       geom_errorbar(aes(ymax = Emax, ymin = Emin), width=0.2),
       coord_flip())
}

ARA.importance1 <- bar.import(ARA.coef1)
ARA.est1 <-est.confint(ARA.coef1)
print(ARA.importance1)
print(ARA.est1)

ARA.coef2plot <- ARA.coef2[ARA.coef2$Importance>=0.5, ]  #  sharp jump from 0.2->0.3->99
ARA.coef2plot$Term <- factor(ARA.coef2plot$Term, levels=unique(ARA.coef2plot$Term))
ARA.importance2 <- bar.import(ARA.coef2plot) + 
                              opts(axis.text.y=theme_text(size=8, hjust=1))
ARA.est2 <-est.confint(ARA.coef2plot) + 
                              opts(axis.text.y=theme_text(size=8, hjust=1))
print(ARA.importance2)
print(ARA.est2)

## Important 2-way interactions:
## Block:Time
## Chamber:Position
## Chamber:Time
## H2O:Cells
## H2O^2:Cells
## Block:Cells
## Block:H2O
## Block:H2O^2
## Chamber:Cells
## H2O^2:Time
## Position:Time
## Implied Higher-order interactions:
##   Cells * H2O * I(H2O^2) * Block
##   Time * Chamber * Position





##==============================================================
## MODEL FITTING: Final Model Structure with higher interactions?
##==============================================================
## use glmulti output to identify which 2-way interactions are important
## include relevant higher-order interactions as necessary, and feasible
Y.formula <- ARA.best2
Y.model <- Y.lmH
Y.model <- Y.lmHB
Y.model <- ARA.best2lm

## Add mixed effects extensions to avoid violating model assumptions?




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

diagnostics(Y.lmain)
diagnostics(Y.lm1)
diagnostics(Y.lmH)
diagnostics(Y.lmsqrt)
diagnostics(ARA.best2lm)

RE <- diagnostics(Y.model, resType="pearson", more=TRUE) # full diagnostics, just to be sure
op <- par(mfrow=c(2,2))
plot(Y.model)
par(op)
residualPlots(Y.model)                 # car

## Check for spatial patterns in residuals?

## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)



################################################################
## ANALYSIS: GET RESULTS
################################################################
## the anova() function performs sequential (Type I) tests: order matters.
anova(Y.model)                         # ORDER MATTERS! (see Zuur et al. 2009 pg. 540))
Anova(Y.model, type=2)                 # Type II: car package**
summary(Y.model)
## effects of single-term deletions?
drop1(Y.model)
drop1(Y.lmain)

## Partial effects of each variable (Zuur et al. 2009, pg. 400)
## termplot() only works for main effects, not when interactions are present :(
## termplot(Y.model, se=T, rug=T, partial.resid=T)
library(effects)
plot(allEffects(Y.model), ask=FALSE)   # interesting, but messy

intermean <- function (vec) {
  vec1 <- rep(NA, length(vec) -1)
  for (i in 1:length(vec1) ) {
    vec1[i] <- mean(vec[c(i, i+1)])
  }
  vec1
}
H2O.breaks9 <- seq(0, max(SECCa$H2O), length.out=10) # 9 groups
H2O.breaks4 <- seq(0, max(SECCa$H2O), length.out=5 ) # 4 groups
H2O.9lvls   <- intermean(H2O.breaks9)  # 9 groups
H2O.4lvls   <- intermean(H2O.breaks4)  # 4 groups

X.eff     <- effect("Chamber:X.trans", Y.model)
X.TCP.eff <- effect("Time : Chamber : Position : X.trans", Y.model)
X.BTC.eff <- effect("Block : Time : Chamber : X.trans", Y.model)
F.eff     <- effect("Time:Frag", Y.model)
X.H.eff   <- effect("X.trans:H2O:I(H2O^2)", Y.model, xlevels=list(H2O=H2O.9lvls))
X.HB.eff  <- effect("Block:X.trans:H2O:I(H2O^2)", Y.model,
                    xlevels=list(Block=1:8, H2O=H2O.4lvls))
plot(X.eff,     ask=FALSE)
plot(X.TCP.eff, ask=FALSE)
plot(X.BTC.eff, ask=FALSE)
plot(F.eff,     ask=FALSE)
plot(X.H.eff,   ask=FALSE)
plot(X.HB.eff,   ask=FALSE)

plot(effect("X.trans",       Y.model), ask=FALSE) # Warning: averaged over interactions
plot(effect("Block:X.trans", Y.model), ask=FALSE)
plot(effect("Block:Time",    Y.model), ask=FALSE)
plot(effect("Time*Chamber",  Y.model), ask=FALSE)
plot(effect("I(H2O^2)",      Y.model), ask=FALSE) # wtf?
plot(effect("H2O",           Y.model), ask=FALSE) # wtf?

if (F) {
  plot(X.eff$x$X.trans, X.eff$fit, ylim=range(c(X.eff$lower, X.eff$upper)), type="n")
  X.eff.ambient <- which(X.eff$x$Chamber=="Ambient")
  lines(X.eff$x$X.trans[ X.eff.ambient], X.eff$fit[   X.eff.ambient], lty=1, col="black")
  lines(X.eff$x$X.trans[ X.eff.ambient], X.eff$lower[ X.eff.ambient], lty=2, col="black")
  lines(X.eff$x$X.trans[ X.eff.ambient], X.eff$upper[ X.eff.ambient], lty=2, col="black")
  lines(X.eff$x$X.trans[-X.eff.ambient], X.eff$fit[  -X.eff.ambient], lty=1, col="red4")
  lines(X.eff$x$X.trans[-X.eff.ambient], X.eff$lower[-X.eff.ambient], lty=2, col="red4")
  lines(X.eff$x$X.trans[-X.eff.ambient], X.eff$upper[-X.eff.ambient], lty=2, col="red4")
  ## effect method for glmulti objects?
  ## getS3method("effect", "lm")
  ## analyze.model? debug(effects), then call "analyze.model" while debugging
  ARA.X.eff <- effect.glmulti("X.trans", ARA.glmulti1)
  ARA.X.eff <- effect("X.trans", Y.model)
}

##================================================
## Partial regression
##================================================
## http://intersci.ss.uci.edu/wiki/index.php/R_partial_regression_plots ??
avPlots(Y.lmH, terms= ~ X.trans * I(H2O^2), ask=FALSE) # car
avPlots(Y.model, terms= ~ X.trans * I(H2O^2), ask=FALSE) # car

if (F) {
  ## using update() gives same output as avPlots, above
  ARA.part <- update(Y.model, .~. -X.trans)
  cb.part  <- update(Y.model, X.trans ~ . -X.trans)
  ## Removing interaction terms gives a slightly different result (positive relationship)
  ARA.part <- update(Y.model, .~. -X.trans -X.trans:H2O -X.trans:I(H2O^2) -Block:X.trans -Chamber:X.trans)
  cb.part  <- update(Y.model, X.trans ~ . -X.trans -X.trans:H2O -X.trans:I(H2O^2) -Block:X.trans -Chamber:X.trans)
  ## Partial out higher-order interactions (cleaner +ve relationship)
  ARA.part <- lm(Y.trans ~ H2O * I(H2O^2) + Block * Time * Chamber * Position + Frag, data=SECCa)
  cb.part  <- lm(X.trans ~ H2O * I(H2O^2) + Block * Time *Chamber*Position + Frag, data=SECCa)
  ## dynamic (simple interactions)
  terms.drop <- grep("X\\.trans(:.+)?", value=TRUE, 
                     strsplit(as.character(formula(Y.model)[3]), " \\+ ")[[1]] )
  terms.drop <- paste(terms.drop, collapse=" -")
  ARA.part <- paste("update(Y.model, .~. -", terms.drop, ")", sep="")
  cb.part  <- paste("update(Y.model, X.trans ~. -", terms.drop, ")", sep="")
}

## Removing main + interaction terms dynamically
RHS.part <- as.character(formula(Y.model)[3]) 
RHS.part <- gsub("\\s?\\*?\\s?X\\.trans\\s?\\*?\\s?", "", RHS.part) # *-notation
RHS.part <- gsub("\\+\\+", "+", RHS.part) # leftovers
RHS.part <- gsub(" [^ ]+\\:\\+", "", RHS.part) # :-notation
RHS.part <- gsub(" \\+ [^ ]+:$", "", RHS.part) # :-notation
ARA.part <- paste("update(Y.model, .~ ", RHS.part, ")", sep="")
cb.part  <- paste("update(Y.model, X.trans ~ ", RHS.part, ")", sep="")
ARA.part <- eval(parse(text=ARA.part))
cb.part  <- eval(parse(text=cb.part))

ARA.re   <- resid(ARA.part, type="response")
cb.re    <- resid(cb.part, type = "response")
ARA.cb   <- lm(ARA.re ~ cb.re)
x.ord <- order(cb.re)
ARA.cb.pred <- predict(ARA.cb, interval="confidence", level=0.95) # 95% CI bands

plot(cb.re, ARA.re, pch=20)
## points(cb.re[SECCa$Chamber=="Full Chamber"], ARA.re[SECCa$Chamber=="Full Chamber"], pch=19, col="red4")
## abline(ARA.cb, col="red")
## 95% CI?
lines(cb.re[x.ord], ARA.cb.pred[x.ord, 1], col="red", lty=1, lwd=2)
lines(cb.re[x.ord], ARA.cb.pred[x.ord, 2], col="red", lty=2)
lines(cb.re[x.ord], ARA.cb.pred[x.ord, 3], col="red", lty=2)
## the relationship is visually clearer when partialling out higher-order interactions between Block & other experimental treatments, but the r-squared is lower!
## i.e. Y.model <- Y.lmH
summary(ARA.cb)                        # R^2 = 0.04 ! :(
ARA.cb.r2 <- format(summary(ARA.cb)$adj.r.squared, digits=2)



################################################################
## SAVE OUTPUT
################################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  capture.output(cat(Save.header), 
                 print(ARA.glmulti2),
				 print(Y.formula),              # model
				 Anova(Y.model),                # model summary
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
                      X.trans  =seq(0, max(SECCa$X.trans), length.out=10 ),
                      H2O      =seq(0, max(SECCa$H2O), length.out=10 ) 
                      )
Y.pred.response <- predict(Y.model, newdata=Y.pred, type="response", interval="confidence", level=0.95)  # newdata must have same explanatory variable names for predict to work.
Y.pred$predicted <- Y.pred.response[, 1]
Y.pred$pred.lwr  <- Y.pred.response[, 2]
Y.pred$pred.upr  <- Y.pred.response[, 3]
Y.pred.terms <- predict(Y.model, type="terms")
## type=="response" for full predictions (including interactions)
## type=="terms" for partial predictions (?)

## Note: there is a predict() method for glmulti objects...
Y.multipred <- predict(ARA.glmulti2, se.fit=TRUE)
Y.pred$multi.fit <- Y.multipred$averages[1]
Y.pred$multi.lwr <- Y.multipred$averages[1] - Y.multipred$variability[, "+/- (alpha=0.05)"]
Y.pred$multi.upr <- Y.multipred$averages[1] + Y.multipred$variability[, "+/- (alpha=0.05)"]
## Use effect() to get predicted effects for specific terms
## summary(eff.obj$term) : $effect, $lower, $upper for 95% CI

if (F) {                               # average predictions for subsets of terms?
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

## spaghetti plots (too many other terms)
par(mfrow=c(1,1))
pred.Y <- with( Y.pred, 
               aggregate(cbind(predicted), list(Chamber = Chamber, X = X.trans), mean)
)  # I should be getting direct predictions, not means of predictions. *****
	# pred | augpred | ?
	plot(SECCa$X.trans, SECCa$Y.trans, type="p",
		ylab=Y.plotlab, xlab=X.plotlab,
		pch=SECCa$pt, col=SECCa$colr, bg=SECCa$fill
        )
	lines(predicted ~ X.trans, data=subset(Y.pred, Chamber == "Ambient"), 
          col = Chamber.map$col[1], 
          lty = Chamber.map$lty[1]
          )
	lines(predicted ~ X.trans, data=subset(Y.pred, Chamber == "Full Chamber"), 
          col = as.character(Chamber.map$col[2]), 
          lty = Chamber.map$lty[2]
          )
	legend( "topright", legend=Chamber.map$label, pch=point, col=as.character(Chamber.map$col), pt.bg=as.character(Chamber.map$bg) )



##==============================================================
## Plot fitted on observed, by factor?
##==============================================================
Chamber.label <- attr(SECC, "labels")[["Chamber"]]
ChamberPts  <- ggPts.SECC(Chamber.map, Chamber.label) 
TopLegend   <- opts(legend.position = "top", legend.direction = "horizontal")
X.label <- paste("\"", attr(SECCa, "labels")[X.col], ": \"*log[10](", attr(SECCa, "units")[X.col], ")", sep="")
Y.label <- paste("\"", attr(SECCa, "labels")[Y.col], ": \"*log[10](", attr(SECCa, "units")[Y.col], ")", sep="")
X.label <- parse(text=X.label)
Y.label <- parse(text=Y.label)
XY.axlab <- list( xlab(bquote(.(X.label))), ylab(bquote(.(Y.label))) )

df.eff <- function (effObj) {
  eff.df <- data.frame(fit     = effObj$fit,
                       lower   = effObj$lower,
                       upper   = effObj$upper)
  for (v in names(effObj$x)) {
    eff.df[, v] <- effObj$x[, v]
  }
  eff.df
}
eff.layers <- function(effObj, conf.int=TRUE, colour="black", bin.name="H2Obin", bin.lvls=NULL, bin.var="H2O") {
  eff.df <- df.eff(effObj)
  if (!is.null(bin.lvls)) {
    eff.df[, bin.name] <- factor(eff.df[, bin.var])
    levels(eff.df[, bin.name]) <- bin.lvls
  }
  ## no $ in aes() !!!!
  result <- list(geom_line(data=eff.df, colour=colour, aes(x=X.trans, y=fit, group=NULL, shape=NULL)) )
  if (conf.int == TRUE) {
    result <- c(result, 
                geom_line(data=eff.df, colour=colour, aes(x=X.trans, y=lower, lty=2)), 
                geom_line(data=eff.df, colour=colour, aes(x=X.trans, y=upper, lty=2)) 
                )
  }
  result
}
eff.C.layers <- function(effObj, conf.int=TRUE) {
  eff.df <- df.eff(effObj)
  ## no $ in aes() !!!!
  result <- list(geom_line(data=eff.df, aes(x=X.trans, y=fit,  group=Chamber, colour=Chamber)) )
  if (conf.int == TRUE) {
    result <- c(result, 
                geom_line(data=eff.df, aes(x=X.trans, y=lower, group=Chamber, colour=Chamber, lty=2)), 
                geom_line(data=eff.df, aes(x=X.trans, y=upper, group=Chamber, colour=Chamber, lty=2)) 
                )
  }
  result
}
eff.F.layers <- function(effObj, conf.int=TRUE, boxes=TRUE) { # for factors
  eff.df <- df.eff(effObj)
  ## no $ in aes() !!!!
  if (boxes==TRUE) {
    result <- list(geom_crossbar(data=eff.df, size=0.6, 
                                 aes(x=Frag, y=fit, ymin=lower, ymax=upper, group=Time)) 
    )
  } else {
    result <- list(geom_line(data=eff.df, size=1.5, aes(x=Frag, y=fit, group=Time)) )
    if (conf.int == TRUE) {
      result <- c(result, 
                  geom_line(data=eff.df, size=1.5, aes(x=Frag, y=lower, lty=2, group=Time)), 
                  geom_line(data=eff.df, size=1.5, aes(x=Frag, y=upper, lty=2, group=Time)) 
                  )
    }
  }
  result
}

ARA.plot <- qplot(X.trans, Y.trans, data=SECCa, geom = "point", size = I(3),
                  group = Chamber, colour = Chamber, shape = Chamber,
                  xlab = bquote(.(X.label)), ylab = bquote(.(Y.label))) + 
                  jaw.ggplot() + ChamberPts + TopLegend
ARA.plot <- ggplot(data=SECCa, aes(x=X.trans, y=Y.trans)) +
            geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) + 
            XY.axlab + jaw.ggplot() + ChamberPts + TopLegend
ARA.C.plot <- ARA.plot + geom_line(data=X.eff$x, aes(x=X.trans, y=X.eff$fit, group=Chamber)) +
              geom_line(data=X.eff$x, aes(x=X.trans, y=X.eff$lower, group=Chamber, lty=2)) +
              geom_line(data=X.eff$x, aes(x=X.trans, y=X.eff$upper, group=Chamber, lty=2))
ARA.C.plot <- ARA.plot + eff.C.layers(X.eff)
ARA.facet.plot <- ARA.plot + eff.C.layers(X.TCP.eff) + facet_grid(facets=Position~Time)
ARA.Block.plot <- ARA.plot + eff.C.layers(X.BTC.eff, conf=F) + facet_grid(facets=Block~Time)

ARA.Frag.plot <- ARA.plot + aes(x=Frag, y=Y.trans) + XY.axlab + xlab("Fragmentation Treatment")
ARA.Frag.plot <- ARA.Frag.plot + eff.F.layers(F.eff) + facet_grid(facets=.~Time) +
                 geom_jitter(size=3, aes(group=Chamber, colour=Chamber, shape=Chamber))# optional

SECCa$H2Obin9 <- cut(SECCa$H2O, breaks=H2O.breaks9)
SECCa$H2Obin4 <- cut(SECCa$H2O, breaks=H2O.breaks4)
SECCa$H2Obin4 <- factor(SECCa$H2Obin4, levels=rev(levels(SECCa$H2Obin4)) ) # arrange in decreasing order, for top-down rows in plot
ARA.H2O.plot <- ggplot(data=SECCa, aes(x=X.trans, y=Y.trans)) +
                geom_point(size = 3, aes(group = Chamber, 
                               colour = Chamber, shape = Chamber)) + 
                XY.axlab + jaw.ggplot() + ChamberPts + TopLegend
ARA.HB.plot  <- ARA.H2O.plot + 
                eff.layers(X.HB.eff, bin.name="H2Obin4", bin.lvls=levels(SECCa$H2Obin4)) + 
                facet_grid(H2Obin4 ~ Block)
ARA.H2O.plot <- ARA.H2O.plot + 
                eff.layers(X.H.eff, bin.name="H2Obin9", bin.lvls=levels(SECCa$H2Obin9)) + 
                facet_wrap(~ H2Obin9)

print(ARA.C.plot)
print(ARA.facet.plot)
print(ARA.Block.plot)
print(ARA.Frag.plot)
print(ARA.H2O.plot)
print(ARA.HB.plot)



if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
