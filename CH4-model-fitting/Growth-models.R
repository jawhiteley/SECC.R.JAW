################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Regression models for Moss Growth
###   WITHOUT Nfix as a predictor variable (highly correlated with H2O)
### Jonathan Whiteley     R v2.12     2012-07-17
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
source('./CH4-model-fitting/Growth_setup.R')
Save.results   <- FALSE


library(nlme)
################################################################
## MODEL FORMULA
################################################################
Y.main   <- Y.trans ~ Block + Chamber + Frag + H2O + logTAN
Y.fixed  <- Y.trans ~ Block * Chamber * Frag * H2O * logTAN
Y.full   <- Y.fixed                    # not enough replication to test full range of interactions
Y.random <-  ~ 1 | Block/Chamber/Frag 

##==============================================================
cat("- Processing Data (excluding NAs)\n")
## Remove NA rows from data, but only for variables in the model ...
## Easier to do it once here, than in every single fitting function in the script ;)
SECCf <- SECCa
Mod.cols <- unlist( strsplit("SampleID, Block, Time, Chamber, Frag, Position, TempC, H2O, ARA.m, Nfix, logNfix, TAN, logTAN, Growth, Y, Y.trans", 
                             ", ", fixed = TRUE) )
SECCa <- na.exclude( SECCa[, Mod.cols] ) # only exclude NAs from relevant (used) columns
attr(SECCa, "labels") <- attr(SECCf, "labels")
attr(SECCa, "units")  <- attr(SECCf, "units")
levels(SECCa$Chamber)[levels(SECCa$Chamber) == "Full Chamber"] <- "Chamber"



## generate grid to add predicted values to (X-values in all combinations of factors).
## - watch length.out: if it's too long, R will choke.
## - real replication is 1 anyway, so it doesn't need to be so big in this case.
Y.pred <- expand.grid(Block    = levels(SECCa$Block) , 
                      Frag     = levels(SECCa$Frag), 
                      Chamber  = levels(SECCa$Chamber), 
                      logNfix  = seq(0, max(SECCa$logNfix, na.rm = TRUE), length.out=3 ),
                      H2O      = seq(0, max(SECCa$H2O, na.rm = TRUE), length.out=3 ), 
                      logTAN  = seq(0, max(SECCa$logTAN, na.rm = TRUE), length.out=3 ) 
                      )


cat("- Fitting models:", Y.col, "\n")
################################################################
## MODEL FITTING
################################################################
## Transformations & Linearity
##==============================================================
## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)
library(mgcv)
Y.gam     <- gam(Y.trans ~ Block + Chamber + Frag + s(H2O)  + s(TAN), data = SECCa)
Y.log.gam <- gam(Y.trans ~ Block + Chamber + Frag + s(H2O)  + s(logTAN), data = SECCa)
anova(Y.gam)
anova(Y.log.gam)
AIC(Y.gam, Y.log.gam)
## Nfix, TAN, and H2O are linear (no need to transform?)
## I log-transformed Nfix in other analyses: I should probably continue to do so, unless I have a really good reason not to.
## The non-linearity in logNfix may be driven largely by a big gap between small & large values

op <- par(mfrow=c(2,2))
plot(Y.gam, ask = FALSE)
plot(Y.log.gam, ask = FALSE)
par(op)


##==============================================================
## MODEL STRUCTURE
##==============================================================
Y.lmain <- lm(Y.main, data = SECCa)
summary(Y.lmain)
anova(Y.lmain)                         # ORDER MATTERS! (see Zuur et al. 2009 pg. 540))
Anova(Y.lmain, type = 2)               # H2O, Nfixation **
Y.lmfull <- lm( Y.full , data=SECCa)   # n = 2 for all interaction combinations (perfect fit)
anova.full <- anova(Y.lmfull)
anova.full$PropVar <- anova.full[, "Sum Sq"] / sum(anova.full[, "Sum Sq"])
anova.full <- anova.full[order(anova.full$PropVar, decreasing=TRUE), ] # sort by Prop. Var
print(sum(anova.full$PropVar))
print(anova.full)

lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 500, msMaxIter = 100, opt="optim")
if (TRUE)
{ ## Mixed Modelling: may need to establish which interaction terms I need to keep first
  Y.gls  <- gls(Y.main, data = SECCa, method = "REML")
  Y.lme  <- lme(Y.main, random = Y.random, data = SECCa)
  Y.lme2 <- lme(Y.main, random = ~ 1 | Block/Frag, data = SECCa)
  anova(Y.gls, Y.lme)                  # random structure not an improvement
  anova(Y.gls, Y.lme2)                 # random structure not an improvement
  Y.lmB  <- gls(Y.main, weights = varIdent(form = ~ 1 | Block), data = SECCa, method = "REML")
  Y.lmeB <- lme(Y.main, random = Y.random, weights = varIdent(form = ~ 1 | Block), 
                data = SECCa, method = "REML")
  anova(Y.gls, Y.lmB, Y.lmeB)          # allowing heterogeneity among Blocks is an improvement (p < 0.001).
  Y.lmBC <- gls(Y.main, weights = varIdent(form = ~ 1 | Block * Chamber), 
                data = SECCa, method = "REML")
  anova(Y.lmB, Y.lmBC)                 # allowing heterogeneity among Blocks AND Chambers is a sig. improvement (p = 0.026)
  ## Nevertheless, acocunting for heterogeneity in both Blocks and Chambers causes convergence problems in some models :(
}



##==============================================================
## MODEL SELECTION: Which TERMS do I need?
##==============================================================
## Normally, this would come AFTER choosing a random structure.
## Unfortunately, I have so little replication that I cannot include the full range of interaction terms
## Therefore, I first need to figure out which terms are more important to include, and which ones can be dropped.
## multi-model averaging and inference with glmulti!
## WARNING: this uses **A LOT** of memory (~ 2GB!!)
library(glmulti)                       #  v1.0, april 2011; v1.0.3, Nov 2011

if (file.exists(Save.glmulti)) 
{ ## load saved object to speed things up
  cat("- loading glmulti objects from previous run\n")
  load(Save.glmulti)
} else {
  ## Note: glmulti performs an exhaustive search of all candidate models, by default.  
  ## If I can ignore the random (nested) structure, I can use gls() as the fitfunc, and account for heterogeneity.
  ## gls() and nlme() don't have the same coef() methods, which makes it harder to extract coefficients after :(
  ## remember: REML for RANDOM effects, ML for Fixed effects ;)
  ## This is quite a bit slower than using lm, but should only take a few minutes.
  Y.glmultiB <- glmulti(Y.main, data=SECCa, crit=aic, level=1, fitfunc=gls.glmulti, 
                        marginality = TRUE, method="h", confsetsize=256, 
                        plotty=FALSE, report=FALSE,
                        weights = varIdent(form = ~ 1 | Block)
                        )
  print(Y.glmultiB)
  ## using lme to account for nested (random) structure: may not actually improve the model (see above), but does account for the experimental design.
  Y.glmultiBr <- glmulti(Y.main, data=SECCa, crit=aic, level=1, fitfunc=lme.glmulti, 
                        marginality = TRUE, method="h", confsetsize=256, 
                        plotty=FALSE, report=TRUE,
                        weights = varIdent(form = ~ 1 | Block),
                        random = Y.random, control = lmc
                        )
  print(Y.glmultiBr)
  Y.glmulti1 <- Y.glmultiBr

  ## 2-WAY INTERACTIONS (level = 2)
  ## An exhaustive exploration would take ~ 1 month on my computer.  
  ## Try with 'genetic' algorithm to speed things up (method="g"). ~ 30 minutes
  ## Best to run 2-4+ replicate genetic algorithms, and take consensus.
  ## use method="d" to print a summary of candidate models (no fitting)
  ## larger confsetsize -> more memory useage
  Y.multi <- list()
  for (i in 1:6) 
  {
    cat("\n================ glmulti: Genetic Algorithm Run", i, "================\n\n")
    Y.multi[[i]] <- glmulti(Y.main, data=SECCa, crit=aic, level=2, fitfunc=lm, 
                            marginality = TRUE, method="g", confsetsize=256, 
                            plotty=FALSE, report=TRUE
                            )
  }
  Y.glmulti2 <- consensus(Y.multi, confsetsize=256) # more models use more memory!

  rm(Y.multi)                        # clean-up

  ## capture output to save before deleting glmulti objects
  Y.pmulti1 <- capture.output(print(Y.glmulti1))
  Y.pmulti2 <- capture.output(print(Y.glmulti2))
  Y.best2 <- as.formula(summary(Y.glmulti2)$bestmodel)
  Y.best2lm <- lm(Y.best2, data=SECCa)
  Y.mtable1 <- weightable(Y.glmulti1)             # models, IC values, and relative weights (for confset)
  Y.mtable2 <- weightable(Y.glmulti2)             # models, IC values, and relative weights (for confset)
  summary(Y.best2lm)


  if (F) 
  {
    summary(Y.glmulti1)
    summary(Y.glmulti2)
    names(Y.glmulti2)
    str(summary(Y.glmulti2))
    summary(Y.glmulti1)$modelweights
  }

  ##   Y.coef1 <- getCoef.glmulti(Y.glmulti1)   # not for gls or lme
  Y.coef2 <- getCoef.glmulti(Y.glmulti2)
  nrow(Y.coef2)                        # 97 terms (96 + intercept)!
  Y.imp1 <- importance.glmulti(Y.glmulti1)
  Y.imp2 <- importance.glmulti(Y.glmulti2)
  ## Extract model-averaged predictions: needs a lot of memory (>600 MB)
  ## May produce errors, but I don't really understand why :(
  Y.multipred <- predict(Y.glmulti2, newdata=Y.pred, se.fit=TRUE)

  ## Output graphs
  par(mfrow=c(1,1))
  plot(Y.glmulti1, type="p")           # red line at ~2 IC units above best model: should consider at least all models below this line.
  plot(Y.glmulti1, type="w")           # red line where cumulative evidence weight = 95%
  ## plot(Y.glmulti1, type="r")           # diagnostic plots (windows only?)
  plot(Y.glmulti1, type="s")           # term weights (v 1.0.3)
  plot(Y.glmulti2, type="s")           # term weights (v 1.0.3)
  plot(Y.glmulti2, type="p")
  plot(Y.glmulti2, type="w")
  ## sum of relative evidence weights over all models that include each term?
  ##   barplot(Y.coef1[, "Importance"], horiz=TRUE, names.arg=Y.coef1$Term, las=2) 



  ## save derivative objects to speed up loading for future analysis.  
  ## The raw glmulti objects make for a big file (and a lot of memory): >1 GB!
  save(Y.pmulti1, Y.pmulti2, Y.best2, Y.best2lm, 
	   Y.mtable1, Y.mtable2, Y.coef2, Y.imp1, Y.imp2, # Y.coef1, 
	   Y.multipred, file=Save.glmulti)

  rm(Y.glmulti1, Y.glmulti2, Y.glmultiB, Y.glmultiBr) # save memory? not right away, but maybe eventually :(
}



clean.term.labels <- function(tls, coef.labels = FALSE) {
  ## custom function for text replacement in model output: useful for making readable graphs and other output.
  if (coef.labels)
  {  ## Coefficient term labels (interactions, etc.)
    tls <- gsub("ChamberFull Chamber", "Chamber", tls)
    tls <- gsub("PositionOuter", "Position", tls)
    tls <- gsub("Block(.*)", "Block \\1", tls)
    tls <- gsub("Frag([^:]*)", "Frag (\\1)", tls)
    tls <- gsub("Time([^:]*)", "Time (\\1)", tls)
  }
  ## Term labels
  tls <- gsub("Cells.m", "Cyanobacteria", tls, fixed=TRUE)
  tls <- gsub("logCells", "Cyanobacteria", tls, fixed=TRUE)
  tls <- gsub("logNfix", "N-fixation", tls, fixed=TRUE)
  tls <- gsub("logTAN", "Total N", tls, fixed=TRUE)
  tls <- gsub("logTAN", "Total N", tls, fixed=TRUE)
  tls <- gsub("I(H2O^2)", "H2O^2", tls, fixed=TRUE)
  tls <- gsub("H2O", "Moisture", tls, fixed=TRUE)
  tls <- gsub("TempC", "Temperature", tls, fixed=TRUE)
  tls <- gsub("Frag", "Fragmentation", tls)
  if (FALSE) { # math expressions for axis labels?
    tls <- gsub(":", "%*%", tls)       # interactions
    tls <- expression(tls)             # expressions
  } else {
  }
  tls
}

## ggplot2: theme settings
library(ggplot2)
bar.import <- function(glm.imp) 
{
  ##     glm.imp <- importance.glmulti(glmObj)
  ggplot(glm.imp, aes(x=Term, y=Importance), stat="identity", xlab = "Model Terms") +
  list(geom_bar(colour="#333333", fill="#999999"), coord_flip(), 
       scale_y_continuous(expand=c(0,0)), scale_x_discrete(expand=c(0.01,0)),
       geom_hline(aes(yintercept=0.8), colour="#000000", lty=3),
       opts(title = "Model-averaged importance of effects",
            panel.border=theme_blank(), axis.line=theme_segment(),
            plot.title = theme_text(size = 16, lineheight = 1.2, face = "bold")
            )
       )
}
est.confint <- function(glmObj) 
{
  ##   conf.wd <- glmObj[, "+/- (alpha=0.05)"]
  ##   conf.int <- aes(ymax = Estimate + conf.wd, ymin = Estimate - conf.wd)
  ggplot(glmObj, aes(x=Term, y=Estimate) ) +
  geom_hline(yintercept=0, colour="grey") + 
  list(geom_point(),
       geom_errorbar(aes(ymax = Emax, ymin = Emin), width=0.2),
       coord_flip()
       )
}

## MAIN effects
Y.importance1 <- bar.import(Y.imp1) + jaw.ggplot()
## Y.est1 <- est.confint(Y.coef1) + jaw.ggplot()    # nlme objects don't provide compatible coef() output
print(Y.importance1)
## print(Y.est1)

## ALL model terms
Y.imp2$Term <- clean.term.labels(Y.imp2$Term)            # DO NOT run this more than once ;)
Y.imp2$Term <- factor(Y.imp2$Term, levels = Y.imp2$Term) # for ggplot order

Y.importance2 <- bar.import(Y.imp2) + jaw.ggplot()       # + opts(axis.text.y=theme_text(size=8, hjust=1))
Y.coef2plot <- Y.coef2[Y.coef2$Importance>=0.2, ]        # Plot the "most important"
Y.coef2plot$Term <- clean.term.labels(Y.coef2plot$Term, coef = TRUE)
Y.coef2plot$Term <- factor(Y.coef2plot$Term, levels=Y.coef2plot$Term)
Y.est2 <- est.confint(Y.coef2plot) + jaw.ggplot() + 
opts(axis.text.y=theme_text(size=8, hjust=1))
print(Y.importance2)
print(Y.est2)

anova(Y.best2lm)                       # order matters (Type 1)
Anova(Y.best2lm, type=2)               # Type II: car package**

## Including Nfix does improve the R-square (explains more variance), 
##   but is it over-fitting or masking other terms?
summary(lm(Y.trans ~ Block + Frag + Chamber * H2O + logTAN, data = SECCa))
summary(lm(Y.trans ~ Block + Frag + Chamber * H2O * logNfix + logTAN, data = SECCa))


## Important 2-way interactions (no mixed effects):
## H2O:Chamber

## Interactions of interest (could be more or less than Y.best2, depending on how that goes)
##   I'm less interested in Block interactions
Y.fixed <- Y.trans ~ Block + Chamber + Frag + H2O + logTAN + Chamber:H2O

## Implied higher-order interactions
## Actually, the same model...
Y.fixHi <- Y.trans ~ Block + Frag + Chamber * H2O + logTAN


##==============================================================
## MODEL FITTING: Final Fixed & Random Effects?
##==============================================================
##	Model-fitting often fails with many interaction terms + heterogeneity + random nested structure (+ ML)
TryMM <- TRUE                         # try fitting complex Mixed Models with interaction terms?

Y.fxlm <- lm(Y.fixed, data = SECCa)
Y.fHlm <- lm(Y.fixHi, data = SECCa)
anova(Y.best2lm, Y.fxlm, Y.fHlm)
Y.f2   <- gls(Y.best2, data = SECCa, method = "ML")
Y.fx   <- gls(Y.fixed, data = SECCa, method = "ML")
Y.fH   <- gls(Y.fixHi, data = SECCa, method = "ML")
anova(Y.fx, Y.f2, Y.fH)                # Adding higher-order interactions is not better
anova(Y.fx, Y.fH)                      
Y.fx   <- update(Y.fx, method = "REML")
Y.fH   <- update(Y.fH, method = "REML")
## I will likely need to account for heterogeneity among Blocks, maybe Chambers as well (TempC/Warming)
Y.mH <- gls(Y.fixHi, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
Y.mr <- lme(Y.fixHi, random = Y.random, data = SECCa, method ="REML")
anova(Y.fH, Y.mr)                      # Adding random effects
anova(Y.fH, Y.mH)                      # Accounting for heterogeneity
anova(Y.mH, Y.mr)                      # are these really nested?  (might not be a valid comparison)
if (TryMM)
{                                      # higher-order interaction models often crash here: trying to do too much!
  Y.me <- lme(Y.fixHi, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
  anova(Y.fH, Y.mH, Y.me, Y.mr)        # compare random structures (using REML)
  anova(Y.fH, Y.me)
}
## The nested structure of random effects is technically not necessary, but it did improve the AIC in glmulti (using ML), 
## and probably should be done on theoretical grounds, to account for the structure of the experiment.

Y.mHf <- gls(Y.fixed, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
anova(Y.fx, Y.mHf)                     # 2-way interactions might be better without heterogeneity!
if (TryMM)
{  # ML estimation may fail with many interactions.
  Y.mef <- lme(Y.fixed, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method = "REML")
  anova(Y.mHf, Y.mef)                  # adding random effects to Y.fixed
  ## Compare fuller model with best model from glmulti
  Y.m2  <- gls(Y.best2, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
  Y.ml2 <- update(Y.m2, method ="ML")
  Y.mb2 <- lme(Y.best2, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="ML")
  Y.mlf <- lme(Y.fixed, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="ML")
  ## Y.mlh <- update(Y.me, method ="ML")    # X Singularity in backsolve at level 0, block 1 :(
  anova(Y.mlf, Y.mb2)      # adding the extra interaction terms does improve the AIC, but not significantly so.
}


## SO, Which model should I use?
##  All I really need a model for at this point is predictions and graphs.
##  Adding heterogeneity improves the residuals, and the fit.
##  Adding random effects doesn't improve the fit, but might improve the residuals
## Y.fit  <- Y.me
##  Unfortunately, effects() throws an error if I try to include both heterogeneity and random structure :(
Y.fit  <- Y.mH




################################################################
## CHECK ASSUMPTIONS: MODEL VALIDATION
################################################################
cat("- Validating models\n")

## residualPlots(Y.best2lm)               # car: lm only
diagnostics(Y.best2lm, X.cols = c("logNfix", "H2O", "logTAN")) # Best model from glmulti
diagnostics(Y.mHf, X.cols = c("logNfix", "H2O", "logTAN"))     # 2-way interactions and heterogeneity
diagnostics(Y.mH,  X.cols = c("logNfix", "H2O", "logTAN"))     # higher-order interactions and heterogeneity
diagnostics(Y.mr,  X.cols = c("logNfix", "H2O", "logTAN"))     # higher-order interactions and nested structure
if (TryMM) diagnostics(Y.me,  X.cols = c("logNfix", "H2O", "logTAN"))     # higher-order interactions + mixed effects
## I do prefer the distribution and lack of patterns in residuals in the mixed effects models
## There are some disturbing patterns in the residuals for the lm model
## Allowing heterogeneity is a modest improvement: a nested random structure is about the same (but theoretically justified).
## There may also be some heterogeneity with H2O, but I might just have to live with it :(

RE <- diagnostics(Y.fit, resType = "pearson", 
                  X.cols = c("logNfix", "H2O", "TAN"), more = TRUE) # Full diagnostics; requires all ID columns
op <- par(mfrow=c(2,2))
plot(Y.fit)
par(op)




################################################################
## GET RESULTS
################################################################
cat("- Generating results & predictions\n")
summary(Y.fit)
anova(Y.fit)                           # ORDER MATTERS! (see Zuur et al. 2009 pg. 540))
if (inherits(Y.fit, "lm")) {
  ## Why does this work sometimes, and throw errors others? (instability; low memory?)
  Anova(Y.fit, type=2)                 # Type II: car package**
  ## effects of single-term deletions?
  drop1(Y.fit)
  drop1(Y.lmain)
}

## Partial effects of each variable (Zuur et al. 2009, pg. 400)
## Use effect() to get predicted effects for specific terms
library(effects)

## Break continuous predictors into discrete groups for facetting.
H2O.breaks9 <- seq(0, max(SECCa$H2O), length.out=10) # 9 groups
H2O.breaks4 <- seq(0, max(SECCa$H2O), length.out=5 ) # 4 groups
T.breaks    <- seq(min(SECCa$TempC), max(SECCa$TempC), length.out=2 ) # 2 groups
H2O.9lvls   <- intermean(H2O.breaks9)  # 9 groups
H2O.4lvls   <- intermean(H2O.breaks4)  # 4 groups

if (!inherits(Y.fit, "lm"))
{  ## effect() throws an error if the model isn't specified explicitly :(
  Y.call <- deparse(Y.fit$call)
  Y.calf <- regexpr("(?<=formula = )[^,]+", Y.call, perl = TRUE)
  if (all(Y.calf == -1)) Y.calf <- regexpr("(?<=model = )[^,]+", Y.call, perl = TRUE)
  if (all(Y.calf == -1)) Y.calf <- regexpr("(?<=fixed = )[^,]+", Y.call, perl = TRUE)
  Y.form <- substr(Y.call, Y.calf, Y.calf + attr(Y.calf, "match.length") -1 )
  Y.form <- Y.form[Y.form != ""]
  Y.newf <- paste( deparse(get(Y.form)), collapse = "\n")
  Y.call <- sub(Y.form, Y.newf, paste(Y.call, collapse = "\n"), fixed = TRUE)
  Y.fit  <- eval( parse( text = Y.call ) )
}

## extract effects (order of terms in interactions matter)
Y.eff     <- allEffects(Y.fit)
T.eff     <- effect("Chamber", Y.fit)
F.eff     <- effect("Frag", Y.fit)
H.eff     <- effect("H2O", Y.fit)
## A.eff     <- effect("logNfix", Y.fit)
N.eff     <- effect("logTAN", Y.fit)
## AH.eff    <- effect("H2O:logNfix", Y.fit, xlevels=list(H2O=H2O.9lvls))
TH.eff    <- effect("Chamber:H2O", Y.fit)
## AT.eff    <- effect("Chamber:logNfix", Y.fit)
## ATC.eff   <- effect("Chamber:H2O:logNfix", Y.fit, xlevels = list(H2O=H2O.9lvls) )

plot(T.eff, ask = FALSE)
plot(F.eff, ask = FALSE)
plot(H.eff, ask = FALSE)
## plot(A.eff, ask = FALSE)
plot(N.eff, ask = FALSE)
## plot(AH.eff, x.var = "logNfix", ask = FALSE)
## plot(AH.eff, x.var = "H2O", ask = FALSE)
plot(TH.eff, x.var = "H2O", ask = FALSE)
## plot(AT.eff, x.var = "logNfix", ask = FALSE)
## plot(ATC.eff, x.var = "logNfix", ask = FALSE)
## plot(ATC.eff, x.var = "H2O", ask = FALSE)



##==============================================================
## Predictions
##==============================================================
## Note: there is a predict() method for glmulti objects...
Y.pred$multi.fit <- Y.multipred$averages[1]
Y.pred$multi.lwr <- Y.multipred$averages[1] - Y.multipred$variability[, "+/- (alpha=0.05)"]
Y.pred$multi.upr <- Y.multipred$averages[1] + Y.multipred$variability[, "+/- (alpha=0.05)"]


##==============================================================
## Partial regression
##==============================================================
## gls (mixed effects) model causes problems in this section
Y.lfit <- if (inherits(Y.fit, "lm")) Y.fit else Y.fHlm
avPlots(Y.best2lm, terms= ~ H2O, ask=FALSE) # car
avPlots(Y.lfit, terms= ~ H2O + logTAN, ask=FALSE) # car


##______________________________________________________________
## Partial regression on logNfix?
## Removing main + interaction terms dynamically (gls should be ok here)
## This should be the same as Growth~Nfix ...
Parts <- PartialFormula("Y.fit", x.var = "logNfix")
Y.part <- eval(Parts$y) # problems fitting with gls? :(
X.part <- eval(Parts$x)

Y.re     <- resid(Y.part, type = "response")
X.re     <- resid(X.part,  type = "response")
Y.X      <- lm(Y.re ~ X.re)
x.ord    <- order(X.re)
Y.X.pred <- predict(Y.X, interval="confidence", level=0.95) # 95% CI bands

plot(X.re, Y.re, pch=20)
## points(X.re[SECCa$Chamber=="Full Chamber"], Y.re[SECCa$Chamber=="Full Chamber"], pch=19, col="red4")
## abline(Y.X, col="red")
## 95% CI?
lines(X.re[x.ord], Y.X.pred[x.ord, 1], col="red", lty=1, lwd=2)
lines(X.re[x.ord], Y.X.pred[x.ord, 2], col="red", lty=2)
lines(X.re[x.ord], Y.X.pred[x.ord, 3], col="red", lty=2)

residualPlots(Y.X)                 # car

summary(Y.X)                        # R^2 = 0.02 ! :(
Y.X.r2 <- format(summary(Y.X)$adj.r.squared, digits=2)
Y.X.df <- data.frame(X=X.re, Y=Y.re, fit=Y.X.pred[, "fit"], 
                        lower=Y.X.pred[, "lwr"], upper=Y.X.pred[, "upr"])


##______________________________________________________________
## Partial regression on H2O
## This will be much stronger without removing the effect of NFix!
Parts <- PartialFormula("Y.fit", x.var = "H2O")
Y.part <- eval(Parts$y) # problems fitting with gls? :(
H.part <- eval(Parts$x)

Y.re     <- resid(Y.part, type = "response")
H.re     <- resid(H.part,  type = "response")
Y.H      <- lm(Y.re ~ H.re)
x.ord    <- order(H.re)
Y.H.pred <- predict(Y.H, interval="confidence", level=0.95) # 95% CI bands

plot(H.re, Y.re, pch=20)
## points(H.re[SECCa$Chamber=="Full Chamber"], Y.re[SECCa$Chamber=="Full Chamber"], pch=19, col="red4")
## abline(Y.H, col="red")
## 95% CI?
lines(H.re[x.ord], Y.H.pred[x.ord, 1], col="red", lty=1, lwd=2)
lines(H.re[x.ord], Y.H.pred[x.ord, 2], col="red", lty=2)
lines(H.re[x.ord], Y.H.pred[x.ord, 3], col="red", lty=2)

residualPlots(Y.H)                 # car

summary(Y.H)                        # R^2 = 0.028 ! :(
Y.H.r2 <- format(summary(Y.H)$adj.r.squared, digits=2)
Y.H.df <- data.frame(H=H.re, Y=Y.re, fit=Y.H.pred[, "fit"], 
                        lower=Y.H.pred[, "lwr"], upper=Y.H.pred[, "upr"])


##______________________________________________________________
## Partial regression on logTAN?
Parts <- PartialFormula("Y.fit", x.var = "logTAN")
Y.part <- eval(Parts$y) # problems fitting with gls? :(
N.part <- eval(Parts$x)

Y.re     <- resid(Y.part, type = "response")
N.re     <- resid(N.part,  type = "response")
Y.N      <- lm(Y.re ~ N.re)
x.ord    <- order(N.re)
Y.N.pred <- predict(Y.N, interval="confidence", level=0.95) # 95% CI bands

plot(N.re, Y.re, pch=20)
## points(N.re[SECCa$Chamber=="Full Chamber"], Y.re[SECCa$Chamber=="Full Chamber"], pch=19, col="red4")
## abline(Y.N, col="red")
## 95% CI?
lines(N.re[x.ord], Y.N.pred[x.ord, 1], col="red", lty=1, lwd=2)
lines(N.re[x.ord], Y.N.pred[x.ord, 2], col="red", lty=2)
lines(N.re[x.ord], Y.N.pred[x.ord, 3], col="red", lty=2)

residualPlots(Y.N)                 # car

summary(Y.N)                        # R^2 = 0.028 ! :(
Y.N.r2 <- format(summary(Y.N)$adj.r.squared, digits=2)
Y.N.df <- data.frame(N=N.re, Y=Y.re, fit=Y.N.pred[, "fit"], 
                        lower=Y.N.pred[, "lwr"], upper=Y.N.pred[, "upr"])




################################################################
## OUTPUT
################################################################
## Save Text Output
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  cat("- Saving Results\n")
  capture.output(cat(Save.head.txt), 
                 print(anova.full),				   # Full model: variance partitioning
				 cat("\n\n"),                      # for output
                 cat(Y.pmulti2, fill=TRUE),		   # multi-model selection
                 print(Y.imp2),					   # model-averaged estimates & weights
				 cat("\n\n"),                      # for output
				 print(formula(Y.fit)),            # model
                 if (inherits(Y.fit, "lm")) Anova(Y.fit, type=2) else anova(Y.fit),
				 summary(Y.fit),                   # model summary
				 cat("\n\n"),                      # for output
				 Anova(Y.X),                       # partial regression
				 summary(Y.X),                     # model summary
				 cat("\n\n"),                      # for output
				 Anova(Y.N),                       # partial regression
				 summary(Y.N),                     # model summary
				 cat("\n\n"),                      # for output
				 Anova(Y.H),                       # partial regression
				 summary(Y.H),                     # model summary
				 cat("\n\n"),                      # for output
				 cat(Save.end.txt),                # END OUTPUT #
				 file = Save.text
				)
}



##==============================================================
## PUBLICATION GRAPHS
##==============================================================
cat("- Final Graphics\n")
Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
levels(Chamber.map$label)[levels(Chamber.map$label) == "Full Chamber"] <- "Chamber"
point <- 21	# 21 for circles with col & bg ; 16 for solid circles
Chamber.map$pch <- c(21, 16)  # use circles for both treatments

Chamber.label <- attr(SECC, "labels")[["Chamber"]]
ChamberPts  <- ggPts.SECC(Chamber.map, Chamber.label) 
TopLegend   <- opts(legend.position = "top", legend.direction = "horizontal")
## Axis Labels:
Y.label <- SECC.axislab(SECCa, col = Y.col, parens=TRUE)
Y.lim   <- range(SECCa$Y)
Y.lim   <- c(-7, 30)

##______________________________________________________________
## utility functions
## (back)transformation functions: defined in SECC.functions.R
eff.Tlayer <- function(eff.df = NULL, conf.int = TRUE)
{   # assumes certain columns in eff.df (using effect.to.df())
  require(ggplot2)
  result <- list(geom_line(data=eff.df, aes(y=effect, group = Chamber, colour = Chamber, lwd = Chamber )) )
  if (conf.int == TRUE) {
    result <- c(result, 
                geom_line(data=eff.df, aes(y=lower, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )), 
                geom_line(data=eff.df, aes(y=upper, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )) 
                )
  }
  result
}



##______________________________________________________________
## Chamber effects
T.pdata <- effect.to.df(T.eff)
T.plot <- ggplot(SECCa, aes(y = Y, x = Chamber)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber),
					   position = position_jitter(width = 0.1)) +
			eff.layer(eff = T.pdata, conf.int = FALSE) +
			xlab(SECC.axislab(SECCa, col = "Chamber", parens=FALSE)) + ylab(Y.label) +
			jaw.ggplot() + ChamberPts + TopLegend # yes, the order matters :/


##______________________________________________________________
## H2O effects
H.pdata <- effect.to.df(H.eff)
H.plot <- ggplot(SECCa, aes(y = Y, x = H2O)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
			eff.layer(eff = H.pdata, conf.int = TRUE) +
			xlab(SECC.axislab(SECCa, col = "H2O", parens=TRUE)) + ylab(Y.label) +
			jaw.ggplot() + ChamberPts + TopLegend # yes, the order matters :/

TH.pdata <- effect.to.df(TH.eff)
TH.plot  <- ggplot(SECCa, aes(y = Y, x = H2O)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
			geom_line(data=TH.pdata, aes(y=effect, group = Chamber, colour = Chamber, lwd = Chamber )) +
			geom_line(data=TH.pdata, aes(y=lower, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )) +
			geom_line(data=TH.pdata, aes(y=upper, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )) +
			xlab(SECC.axislab(SECCa, col = "H2O", parens=TRUE)) + ylab(Y.label) +
			jaw.ggplot() + ChamberPts + TopLegend

## N-fixation lower in Contiguous chambers on average: likely due to less disturbance (DeLuca pers. comm.)?
## FH.pdata <- effect.to.df(FH.eff)
## FH.plot  <- ggplot(SECCa, aes(y = Y, x = H2O)) + ylim(Y.lim) +
##             geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
##             eff.layer(eff = FH.pdata, conf.int = TRUE) + facet_wrap(~ Frag) +
##             xlab(SECC.axislab(SECCa, col = "H2O", parens=TRUE)) + ylab(Y.label) +
##             jaw.ggplot() + ChamberPts + TopLegend


##______________________________________________________________
## Available N effects
N.pdata <- effect.to.df(N.eff)
N.pdata <- within(N.pdata, TAN <- 10^logTAN)
N.plot  <- ggplot(SECCa, aes(y = Y, x = TAN)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
			eff.layer(eff = N.pdata, conf.int = TRUE) +
			xlab(SECC.axislab(SECCa, col = "TAN", parens=TRUE)) + ylab(Y.label) +
			jaw.ggplot() + ChamberPts + TopLegend # yes, the order matters :/

##______________________________________________________________
## Partial Regression: N-fixation
## adding plotMath to a ggplot graph: https://groups.google.com/forum/?fromgroups#!topic/ggplot2/-Ind8XDqaPQ
Xpart.notes <- RegPlot.annote(Y.X)

X.part.plot <- ggplot(data=Y.X.df, aes(x=X, y=Y)) +
                 geom_point(size=3, pch=20) + jaw.ggplot()   +
				 geom_text(aes(min(X), max(Y), label = Xpart.notes[1] ), 
						   size = 4, hjust = 0, vjust = 0, parse = TRUE) +
				 geom_text(aes(min(X), max(Y), label = Xpart.notes[2] ), 
						   size = 4, hjust = 0, vjust = 1.5, parse = TRUE) +
				 geom_text(aes(min(X), max(Y), label = Xpart.notes[3] ), 
						   size = 4, hjust = 0, vjust = 2.7, parse = TRUE) +
                 xlab("N-fixation | others") + 
                 ylab("Moss Growth | others") 
X.part.plot <- X.part.plot + geom_line(aes(y=fit), size=1, lty=1, colour="#CC0000") +
                 geom_line(aes(y=lower), size=0.5, lty=2, colour="#CC0000") + 
                 geom_line(aes(y=upper), size=0.5, lty=2, colour="#CC0000")


##______________________________________________________________
## Partial Regression: H2O
Hpart.notes <- RegPlot.annote(Y.H)

H.part.plot <- ggplot(data=Y.H.df, aes(x=H, y=Y)) +
                 geom_point(size=3, pch=20) + jaw.ggplot()   +
				 geom_text(aes(min(H), max(Y), label = Hpart.notes[1] ), 
						   size = 4, hjust = 0, vjust = 0, parse = TRUE) +
				 geom_text(aes(min(H), max(Y), label = Hpart.notes[2] ), 
						   size = 4, hjust = 0, vjust = 1.5, parse = TRUE) +
				 geom_text(aes(min(H), max(Y), label = Hpart.notes[3] ), 
						   size = 4, hjust = 0, vjust = 2.7, parse = TRUE) +
                 xlab("Moisture Contents | others") + 
                 ylab("Moss Growth | others") 
H.part.plot <- H.part.plot + geom_line(aes(y=fit), size=1, lty=1, colour="#CC0000") +
                 geom_line(aes(y=lower), size=0.5, lty=2, colour="#CC0000") + 
                 geom_line(aes(y=upper), size=0.5, lty=2, colour="#CC0000")


##______________________________________________________________
## Partial Regression: TAN
## adding plotMath to a ggplot graph: https://groups.google.com/forum/?fromgroups#!topic/ggplot2/-Ind8XDqaPQ
Npart.notes <- RegPlot.annote(Y.N)

N.part.plot <- ggplot(data=Y.N.df, aes(x=N, y=Y)) +
                 geom_point(size=3, pch=20) + jaw.ggplot()   +
				 geom_text(aes(max(N), max(Y), label = Npart.notes[1] ), 
						   size = 4, hjust = 1, vjust = 0, parse = TRUE) +
				 geom_text(aes(max(N), max(Y), label = Npart.notes[2] ), 
						   size = 4, hjust = 1, vjust = 1.5, parse = TRUE) +
				 geom_text(aes(max(N), max(Y), label = Npart.notes[3] ), 
						   size = 4, hjust = 1, vjust = 2.7, parse = TRUE) +
                 xlab("Total N | others") + 
                 ylab("Moss Growth | others") 
N.part.plot <- N.part.plot + geom_line(aes(y=fit), size=1, lty=1, colour="#CC0000") +
                 geom_line(aes(y=lower), size=0.5, lty=2, colour="#CC0000") + 
                 geom_line(aes(y=upper), size=0.5, lty=2, colour="#CC0000")



##==============================================================
## SAVE GRAPHS
##==============================================================
if (Save.results == TRUE) 
{
  ## glmulti plots
  ggsave(filename = sprintf("%sImportance1.eps", Fig.filename), plot = Y.importance1, 
         width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sImportance2.eps", Fig.filename), plot = Y.importance2, 
         width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sEstimates.eps",   Fig.filename), plot = Y.est2, 
         width = 4, height = 4, scale = 1.5)
  ## Effects plots
  ggsave(filename = sprintf("%sTemp.eps",   Fig.filename), plot = T.plot, 
         width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sH2O.eps",    Fig.filename), plot = H.plot, 
         width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sTxH2O.eps",  Fig.filename), plot = TH.plot, 
         width = 4, height = 4, scale = 1.5)
  ##   ggsave(filename = sprintf("%sFxH2O.eps",  Fig.filename), plot = FH.plot, 
  ##          width = 4, height = 4, scale = 1.5)
  ##   ggsave(filename = sprintf("%sNfix.eps",   Fig.filename), plot = A.plot, 
  ##          width = 4, height = 4, scale = 1.5)
  ##   ggsave(filename = sprintf("%sTxNfix.eps", Fig.filename), plot = AT.plot, 
  ##          width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sTAN.eps",    Fig.filename), plot = N.plot, 
         width = 4, height = 4, scale = 1.5)
  ## Partial Regression plots
  ##   ggsave(filename = sprintf("%sNfix-partial.eps", Fig.filename), plot = X.part.plot, 
  ##          width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sTAN-partial.eps",  Fig.filename), plot = N.part.plot, 
         width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sH2O-partial.eps",  Fig.filename), plot = H.part.plot, 
         width = 4, height = 4, scale = 1.5)
} else {
  print(T.plot)
  print(H.plot)
  print(TH.plot)
  print(FH.plot)
  print(A.plot)
  ##   print(AT.plot)
  print(N.plot)
  print(NT.plot)
  ## Partial Regression plots
  ##   print(X.part.plot)
  print(N.part.plot)
  print(H.part.plot)
}
cat("- Finished Model Fitting:", Y.col, "-\n")
