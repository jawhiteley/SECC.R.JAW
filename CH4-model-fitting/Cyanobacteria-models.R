################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Regression models for Cyanobacteria density
### Jonathan Whiteley     R v2.12     2012-08-05
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
source('./CH4-model-fitting/Cyanobacteria_setup.R')


library(nlme)
################################################################
## MODEL FORMULA
################################################################
Y.main   <- Y.trans ~ Block + Chamber + Frag + H2O + I(H2O^2) + log10(TAN)
Y.fixed  <- Y.trans ~ Block + Chamber * Frag * H2O * I(H2O^2) * log10(TAN)
Y.full   <- Y.fixed                    # not enough replication to test full range of interactions
Y.random <-  ~ 1 | Block/Chamber/Frag  # does this still work with Climate pseudo-factor?
## Y.random <-  ~ 1 | Block/Frag/Climate  # Spatially, Climate might be bigger than Frag, but it makes less sense conceptually: patches *within Frag treatments* are subject to different climate conditions, not the other way around.

##==============================================================
cat("- Processing Data (excluding NAs)\n")
## Remove NA rows from data, but only for variables in the model ...
## Easier to do it once here, than in every single fitting function in the script ;)
SECCf <- SECCa
if (TRUE)
{
  Mod.cols <- unlist( strsplit("SampleID, Block, Time, Chamber, Frag, Position, Climate, H2O, Cells, Cells.m, TAN, Y, Y.log, Y.trans", 
                       ", ", fixed = TRUE) ) # include untransformed columns to use their NA values for filtering
  SECCa <- na.exclude( SECCa[, Mod.cols] ) # only exclude NAs from relevant (used) columns
} else {
  ModFrame <- model.frame(Y.main, data = SECCa, na.action = na.exclude) # includes transformations, drops all other attributes :(
  Mod.cols <- attr(ModFrame, "variables")
  SECCa <- ModFrame                    # includes transformations, missing some columns I want to keep
}
attr(SECCa, "labels") <- attr(SECCf, "labels")
attr(SECCa, "units")  <- attr(SECCf, "units")
levels(SECCa$Chamber)[levels(SECCa$Chamber) == "Full Chamber"] <- "Chamber"



## generate grid to add predicted values to (X-values in all combinations of factors).
## - watch length.out: if it's too long, R will choke.
## - real replication is 1 anyway, so it doesn't need to be so big in this case.
Y.pred <- expand.grid(Block    = levels(SECCa$Block) , 
                      Frag     = levels(SECCa$Frag), 
                      Chamber  = levels(SECCa$Chamber), 
                      ##                       Position = levels(SECCa$Position), 
                      ##                       Climate = levels(SECCa$Climate), 
                      H2O      = seq(0, max(SECCa$H2O, na.rm = TRUE), length.out=3 ), 
                      TAN      = seq(0, max(SECCa$TAN, na.rm = TRUE), length.out=3 ) 
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
Y.gam     <- gam(Y.trans ~ Block + Chamber + Frag + s(H2O) + s(TAN), data = SECCa)
Y.log.gam <- gam(Y.trans ~ Block + Chamber + Frag + s(H2O) + s(log10(TAN)), data = SECCa)
anova(Y.gam)
anova(Y.log.gam)
AIC(Y.gam, Y.log.gam)
## I log-transformed Cells in original analysis: I should probably continue to do so, unless I have a really good reason not to.

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
Anova(Y.lmain, type = 2)               # 
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
  anova(Y.gls, Y.lmB, Y.lmeB)          # allowing heterogeneity among Blocks is an improvement (p = 0.01).
  Y.lmBC <- gls(Y.main, weights = varIdent(form = ~ 1 | Block * Chamber), 
                data = SECCa, method = "REML")
  anova(Y.lmB, Y.lmBC)                 # allowing heterogeneity among Blocks AND Chambers IS a sig. improvement (p < 0.0001)
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
## library(MASS)                          # ?
## library(leaps)

if (file.exists(Save.glmulti)) 
{ ## load saved object to speed things up
  cat("- loading glmulti objects from previous run\n")
  load(Save.glmulti)
} else {
  ## Note: glmulti performs an exhaustive search of all candidate models, by default.  
  ## gls() and nlme() don't have the same coef() methods, which makes it harder to extract coefficients after :(
  ## remember: REML for RANDOM effects, ML for Fixed effects ;)
  ## This is quite a bit slower than using lm, but should only take a few minutes.
  ## Allowing heterogeneity among Blocks AND Chambers causes => Error in eigen(val) : infinite or missing values in 'x'
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
	   Y.mtable1, Y.mtable2, Y.coef2, Y.imp1, Y.imp2, 
	   ## Y.coef1, Y.importance1, Y.est1, Y.importance2, Y.est2,
	   Y.multipred, file=Save.glmulti)
  cat("- glmulti objects saved.\n")

  rm(Y.glmulti1, Y.glmulti2, Y.glmultiB, Y.glmultiBr) # save memory? not right away, but maybe eventually :(
  cat("- glmulti objects clean-up.\n")
}



clean.term.labels <- function(tls, coef.labels = FALSE) {
  ## custom function for text replacement in model output: useful for making readable graphs and other output.
  if (coef.labels)
  {  ## Coefficient term labels (interactions, etc.)
    tls <- gsub("ChamberFull Chamber", "Chamber", tls)
    tls <- gsub("ChamberChamber", "Chamber", tls)
    tls <- gsub("PositionOuter", "Position", tls)
    tls <- gsub("Climate(.*)", "\\1", tls)
    tls <- gsub("Block(.*)", "Block \\1", tls)
    tls <- gsub("Frag([^:]*)", "Isolation (\\1)", tls)
    tls <- gsub("Time([^:]*)", "Time (\\1)", tls)
  }
  ## Term labels
  tls <- gsub("logCells", "Cyanobacteria", tls) # attr(SECC, "labels")[X.col]
  tls <- gsub("log10(TAN)", "Total N", tls, fixed=TRUE)
  tls <- gsub("I(H2O^2)", "H2O^2", tls, fixed=TRUE)
  tls <- gsub("TempC", "Temperature", tls, fixed=TRUE)
  tls <- gsub("H2O", "Moisture", tls)
  # tls <- gsub("Frag", "Isolation", tls)  # Fragmentation or Isolation
  tls <- gsub("Frag", "Fragmentation", tls)
  tls <- gsub("Fragmentationmentation", "Fragmentation", tls) # fix double-substitution
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
       coord_flip())
}

## MAIN effects
Y.importance1 <- bar.import(Y.imp1) + jaw.ggplot()
## Y.est1 <- est.confint(Y.coef1) + jaw.ggplot()
print(Y.importance1)
## print(Y.est1)

## ALL model terms
Y.imp2$Term <- clean.term.labels(Y.imp2$Term)
Y.imp2$Term <- factor(Y.imp2$Term, levels = Y.imp2$Term)

Y.importance2 <- bar.import(Y.imp2) + jaw.ggplot() # + opts(axis.text.y=theme_text(size=8, hjust=1))
Y.coef2plot <- Y.coef2[Y.coef2$Importance>=0.5, ] #  sharp jump from 0.25->0.75->0.8
Y.coef2plot$Term <- clean.term.labels(Y.coef2plot$Term, coef = TRUE)
Y.coef2plot$Term <- factor(Y.coef2plot$Term, levels=Y.coef2plot$Term)
Y.est2 <- est.confint(Y.coef2plot) + jaw.ggplot() + 
opts(axis.text.y=theme_text(size=8, hjust=1))
print(Y.importance2)
print(Y.est2)

anova(Y.best2lm)                       # order matters (Type 1)
Anova(Y.best2lm, type=2)               # Type II: car package**


## Important 2-way interactions (no mixed effects):

## Be careful with the H2O:I(H2O^2) interaction: it does screwy things to predicted values
##  e.g. fitting values higher than the data, particularly in dry Ambient patches (where there is NO DATA).
##  It's basically fitting a cubic polynomial function, which isn't really what I want :P
##  Better to avoid it in the model fitting (but I will need it for plotting effects).

##  I should probably avoid Chamber interactions (especially Chamber:H20 AND H2O^2), given that I have no data on dry Ambient patches :(
##  See Nfix-models.R for more notes about this.

Y.fixed <- Y.trans ~ Block + Chamber + Frag + H2O + I(H2O^2) + log10(TAN) + 
            ##  Not interested in Block interactions; Position is out of the model.
            ## Block:Chamber + Block:Frag + Block:Position + Chamber:Position + Frag:Position +
			Block:H2O # + Block:I(H2O^2) + # + Block:log10(TAN) + 
            ## Chamber:I(H2O^2) + Chamber:H2O # + 
			## Position:H2O + Position:I(H2O^2) + Position:log10(TAN)

## Implied higher-order interactions
## Add 3-way interactions to Y.fixed
Y.fixHi <- update(Y.fixed, .~. + Chamber:Frag + Frag:H2O + Frag:I(H2O^2) # + Chamber:Frag:H2O + Chamber:Frag:I(H2O^2)
                  ##+  Chamber:Frag:Position + Chamber:Frag:log10(TAN) +
                  ## Chamber:Position:H2O + Chamber:Position:I(H2O^2) # NO! avoid H2O:I(H2O^2)
)
## Y.fixHi <- Y.trans ~ Block + Chamber * Frag * H2O * log10(TAN) + Chamber * Frag * I(H2O^2) * log10(TAN) + Block:H2O


##==============================================================
## MODEL FITTING: Final Fixed & Random Effects?
##==============================================================
TryMM <- TRUE                         # try fitting complex Mixed Models with interaction terms?

Y.fxlm <- lm(Y.fixed, data = SECCa)
Y.fHlm <- lm(Y.fixHi, data = SECCa)
anova(Y.best2lm, Y.fxlm, Y.fHlm)
AIC(Y.best2lm, Y.fxlm, Y.fHlm)
Y.f2   <- gls(Y.best2, data = SECCa, method = "ML")
Y.fx   <- gls(Y.fixed, data = SECCa, method = "ML")
Y.fH   <- gls(Y.fixHi, data = SECCa, method = "ML")
anova(Y.fx, Y.f2, Y.fH)                # Adding higher-order interactions
Y.f2   <- update(Y.f2, method = "REML")
Y.fx   <- update(Y.fx, method = "REML")
Y.fH   <- update(Y.fH, method = "REML")
## I will likely need to account for heterogeneity among Blocks, maybe Chambers as well (TempC/Warming)
Y.mH2 <- gls(Y.best2, weights = varIdent(form = ~ 1 | Block * Chamber), data = SECCa, method ="REML")
Y.mHf <- gls(Y.fixed, weights = varIdent(form = ~ 1 | Block * Chamber), data = SECCa, method ="REML")
Y.mH  <- gls(Y.fixHi, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
Y.mr  <- lme(Y.fixHi, random = Y.random, data = SECCa, method ="REML")
anova(Y.f2, Y.mH2)                     # Adding heterogeneity to best model from glmulti
anova(Y.fx, Y.mHf)                     # Accounting for heterogeneity: 2-way interactions
anova(Y.fH, Y.mH)                      # Accounting for heterogeneity: higher-order interactions
anova(Y.fH, Y.mr)                      # Adding random effects
anova(Y.mH, Y.mr)                      # are these really nested?  (might not be a valid comparison)
if (TryMM)
{                                      # higher-order interaction models often crash here: trying to do too much!
  Y.me <- lme(Y.fixHi, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
  anova(Y.fH, Y.mH, Y.me, Y.mr)        # compare random structures (using REML)
  anova(Y.fH, Y.me)
}
## The nested structure of random effects is technically not necessary, 
## but probably should be done on theoretical grounds, to account for the structure of the experiment.

if (TryMM)
{  # ML estimation may fail with many interactions.
  Y.mef <- lme(Y.fixed, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method = "REML")
  anova(Y.mHf, Y.mef)
  ## Compare fuller model with best model from glmulti
  Y.mb2 <- lme(Y.best2, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="ML")
  Y.mlf <- lme(Y.fixed, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="ML")
  ## Y.mlh <- update(Y.me, method ="ML")    # X Singularity in backsolve at level 0, block 1 :(
  anova(Y.mlf, Y.mb2)      # adding the extra interaction term does increase the AIC (worse), but not significantly so.
}

## SO, Which model should I use?
##  All I really need a model for at this point is predictions and graphs.
##  Ideally, I would prefer the full model with higher order interactions and all random effects (Y.me),
##    however, this model may be too complex: several coefficients have 0 df & p = NaN (see summary(Y.me))
##	Model-fitting often fails with many interaction terms + heterogeneity + random nested structure (+ ML)
##    Effects also can't handle nesting structure? X can't handle higher-order interactions :(
##  Adding heterogeneity to some models makes the residuals worse (less normal, more patterns?)
##  Higher-order interactions seem to be a better fit, but the residuals may be worse, and the predicted values get really weird.
## Y.fit  <- Y.me
Y.fit  <- Y.mHf




################################################################
## CHECK ASSUMPTIONS: MODEL VALIDATION
################################################################
cat("- Validating models\n")

## residualPlots(Y.best2lm)               # car: lm only
diagnostics(Y.best2lm, X.cols = c("H2O", "TAN")) # Best model from glmulti
diagnostics(Y.mH2, X.cols = c("H2O", "TAN"))     # Best model from glmulti + heterogeneity
diagnostics(Y.mHf, X.cols = c("H2O", "TAN"))     # 2-way interactions + heterogeneity
diagnostics(Y.mH,  X.cols = c("H2O", "TAN"))     # higher-order interactions + heterogeneity
diagnostics(Y.mr,  X.cols = c("H2O", "TAN"))     # higher-order interactions + nested structure
if (TryMM) diagnostics(Y.me,  X.cols = c("H2O", "TAN"))     # higher-order interactions + mixed effects
## I do prefer the distribution and lack of patterns in residuals in the mixed effects models
## There are some disturbing patterns in the residuals for the lm model
## Allowing heterogeneity is an improvement: a nested random structure is not.
## There may also be some heterogeneity with H2O, but I might just have to live with it :(

RE <- diagnostics(Y.fit, resType = "pearson", 
                  X.cols = c("H2O", "TAN"), more = TRUE) # Full diagnostics; requires all ID columns
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
  Anova(Y.fit, type=2)                 # Type II: car package**
  ## effects of single-term deletions?
  drop1(Y.fit)
  drop1(Y.lmain)
}

## Partial effects of each variable (Zuur et al. 2009, pg. 400)
## Use effect() to get predicted effects for specific terms
## Effects to plot: 
##   Block * logCells * H2O * H2O^2
##   Frag * H2O(^2)
##   TempC * TAN
##   logCells * H2O(^2)
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
H.eff     <- effect("H2O:I(H2O^2)", Y.fit)
N.eff     <- effect("log10(TAN)", Y.fit)
TH.eff    <- effect("Chamber:H2O:I(H2O^2)", Y.fit)
## TPH.eff    <- effect("Chamber:Position:H2O:I(H2O^2)", Y.fit)
##   FP.eff    <- effect("Frag:Position", Y.fit)
##   NT.eff    <- effect("Chamber:log10(TAN)", Y.fit)
##   NP.eff    <- effect("Position:log10(TAN)", Y.fit)
##   NTP.eff    <- effect("Chamber:Position:log10(TAN)", Y.fit)

plot(F.eff, ask = FALSE)
plot(H.eff, ask = FALSE)
plot(T.eff, ask = FALSE)
plot(N.eff, ask = FALSE)
plot(TH.eff, x.var = "H2O", ask = FALSE)
## plot(TPH.eff, x.var = "H2O", ask = FALSE)
##   plot(FP.eff, x.var = "Frag", ask = FALSE)
##   plot(FP.eff, x.var = "Position", ask = FALSE)
##   plot(NT.eff, x.var = "TAN", ask = FALSE)
##   plot(NP.eff, x.var = "TAN", ask = FALSE)
##   plot(NTP.eff, x.var = "TAN", ask = FALSE)


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
Y.lfit <- if (inherits(Y.fit, "lm")) Y.fit else Y.fHlm
avPlots(Y.best2lm, terms= ~ H2O * I(H2O^2) + log10(TAN), ask=FALSE) # car
avPlots(Y.lfit,    terms= ~ H2O * I(H2O^2) + log10(TAN), ask=FALSE) # car


##______________________________________________________________
## Partial regression on H2O
Parts <- PartialFormula("Y.fit", x.var = "H2O")
Y.part <- eval(Parts$y)
H.part <- eval(Parts$x)

Y.re     <- resid(Y.part, type = "response")
H.re     <- resid(H.part,  type = "response")
x.ord    <- order(H.re)
Y.H      <- lm(Y.re ~ H.re)
Y.H.pred <- predict(Y.H, interval="confidence", level=0.95) # 95% CI bands

plot(H.re, Y.re, pch=20)
## 95% CI
lines(H.re[x.ord], Y.H.pred[x.ord, 1], col="red", lty=1, lwd=2)
lines(H.re[x.ord], Y.H.pred[x.ord, 2], col="red", lty=2)
lines(H.re[x.ord], Y.H.pred[x.ord, 3], col="red", lty=2)

residualPlots(Y.H)                 # car

summary(Y.H)                        # R^2 = 0.028 ! :(
Y.H.r2 <- format(summary(Y.H)$r.squared, digits=2)
Y.H.df <- data.frame(H=H.re, Y=Y.re, fit=Y.H.pred[, "fit"], 
                        lower=Y.H.pred[, "lwr"], upper=Y.H.pred[, "upr"])



##______________________________________________________________
## Partial regression on TAN?
Parts <- PartialFormula("Y.fit", x.var = "log10(TAN)")
Y.part <- eval(Parts$y) # problems fitting with gls? :(
N.part <- eval(Parts$x)

Y.re     <- resid(Y.part, type = "response")
N.re     <- resid(N.part,  type = "response")
Y.N      <- lm(Y.re ~ N.re)
x.ord    <- order(N.re)
Y.N.pred <- predict(Y.N, interval="confidence", level=0.95) # 95% CI bands

plot(N.re, Y.re, pch=20)
## 95% CI
lines(N.re[x.ord], Y.N.pred[x.ord, 1], col="red", lty=1, lwd=2)
lines(N.re[x.ord], Y.N.pred[x.ord, 2], col="red", lty=2)
lines(N.re[x.ord], Y.N.pred[x.ord, 3], col="red", lty=2)

residualPlots(Y.N)                 # car

summary(Y.N)                        # R^2 = 0.02 ! :(
Y.N.r2 <- format(summary(Y.N)$r.squared, digits=2)
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
				 Anova(Y.H),                       # partial regression
				 summary(Y.H),                     # model summary
				 cat("\n\n"),                      # for output
				 Anova(Y.N),                       # partial regression
				 summary(Y.N),                     # model summary
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
Y.lim   <- c(0, 100000)

##______________________________________________________________
## utility functions
## (back)transformation functions: defined in SECC.functions.R
scale_y_Cells <- function()
{
  require(ggplot2)
  scale_y_continuous(trans = "log10", 
                     breaks = c(0 , 100, 200, 400, 600, 800, 1000, 
                                2000, 4000, 6000, 8000, 10000, 
                                20000, 40000, 60000, 80000, 100000), 
                     labels = c(0 , 100 , ""  , ""  , ""  , "" , 1000 , 
                                ""  , ""  , ""  , "" , 10000,
                                ""  , ""  , ""  , "" , 100000),
                     )
}

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
T.pdata <- effect.to.df(T.eff, fun = alog0)
T.plot <- ggplot(SECCa, aes(y = Y, x = Chamber)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber),
					   position = position_jitter(width = 0.1)) +
			geom_crossbar(aes(y = effect, ymin = lower, ymax = upper), width = 0.3, data = T.pdata) +
			## eff.layer(eff = T.pdata, conf.int = FALSE) +
			xlab(SECC.axislab(SECCa, col = "Chamber", parens=FALSE)) + ylab(Y.label) +
			scale_y_Cells() + jaw.ggplot() + ChamberPts + TopLegend # yes, the order matters :/


##______________________________________________________________
## Frag effects
F.pdata <- effect.to.df(F.eff, fun = alog0)
F.plot <- ggplot(SECCa, aes(y = Y, x = Frag)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber),
					   position = position_jitter(width = 0.1)) +
			geom_crossbar(aes(y = effect, ymin = lower, ymax = upper), width = 0.3, data = F.pdata) +
			## eff.layer(eff = F.pdata, conf.int = FALSE) +
			xlab(SECC.axislab(SECCa, col = "Frag", parens=FALSE)) + ylab(Y.label) +
			scale_y_Cells() + jaw.ggplot() + ChamberPts + TopLegend # yes, the order matters :/


##______________________________________________________________
## H2O effects
H.pdata <- effect.to.df(H.eff, fun = alog0)
H.plot <- ggplot(SECCa, aes(y = Y, x = H2O)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
			eff.layer(eff = H.pdata, conf.int = TRUE) +
			xlab(SECC.axislab(SECCa, col = "H2O", parens=TRUE)) + ylab(Y.label) +
			scale_y_Cells() + jaw.ggplot() + ChamberPts + TopLegend # yes, the order matters :/

TH.pdata <- effect.to.df(TH.eff, fun = alog0)
TH.plot  <- ggplot(SECCa, aes(y = Y, x = H2O)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
			eff.Tlayer(eff = TH.pdata, conf.int = TRUE) + # facet_wrap(~ Climate) +
			xlab(SECC.axislab(SECCa, col = "H2O", parens=TRUE)) + ylab(Y.label) +
			scale_y_Cells() + jaw.ggplot() + ChamberPts + TopLegend

            ## TPH.pdata <- effect.to.df(TPH.eff, fun = alog0)
            ## TPH.plot  <- ggplot(SECCa, aes(y = Y, x = H2O)) + ylim(Y.lim) +
            ##             geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
            ##             geom_line(data=TPH.pdata, aes(y=effect, group = Chamber, colour = Chamber, lwd = Chamber )) +
            ##             geom_line(data=TPH.pdata, aes(y=lower, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )) +
            ##             geom_line(data=TPH.pdata, aes(y=upper, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )) +
            ##             scale_y_Cells() + facet_wrap(~ Position) +
            ##             xlab(SECC.axislab(SECCa, col = "H2O", parens=TRUE)) + ylab(Y.label) +
            ##             jaw.ggplot() + ChamberPts + TopLegend


##______________________________________________________________
## Available N effects
N.pdata <- effect.to.df(N.eff, fun = alog0)
N.plot  <- ggplot(SECCa, aes(y = Y, x = TAN)) + ylim(Y.lim) +
			geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
			eff.layer(eff = N.pdata, conf.int = TRUE) +
			xlab(SECC.axislab(SECCa, col = "TAN", parens=TRUE)) + ylab(Y.label) +
			scale_x_log10() + scale_y_Cells() + jaw.ggplot() + ChamberPts + TopLegend

            ## NTP.pdata <- effect.to.df(NTP.eff, fun = alog0)
            ## NTP.plot  <- ggplot(SECCa, aes(y = Y, x = TAN)) + ylim(Y.lim) +
            ##             geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) +
            ##             geom_line(data=NTP.pdata, aes(y=effect, group = Chamber, colour = Chamber, lwd = Chamber )) +
            ##             geom_line(data=NTP.pdata, aes(y=lower, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )) +
            ##             geom_line(data=NTP.pdata, aes(y=upper, lty=2, group = Chamber, colour = Chamber, lwd = Chamber )) +
            ##             facet_wrap(~ Position) +
            ##             xlab(SECC.axislab(SECCa, col = "TAN", parens=TRUE)) + ylab(Y.label) +
            ##             scale_x_log10() + scale_y_Cells() + jaw.ggplot() + ChamberPts + TopLegend


##______________________________________________________________
## Partial Regression: H2O
Y.plim <- c(10^-5, 10^4)
Hpart.notes <- RegPlot.annote(Y.H)
## Back-transform: this is simple 10^x, despite the original log10(x+1) transformation.  Subtracting the 1 leads to -ve values, which won't be plotted on log-axes :(
Y.H.df <- within(Y.H.df, {
                 Y <- 10^(Y)
                 fit <- 10^(fit)
                 lower <- 10^(lower)
                 upper <- 10^(upper)
})

H.part.plot <- ggplot(data=Y.H.df, aes( x=H, y=Y )) +
                 geom_point(size=3, pch=20) + jaw.ggplot()   +
				 geom_text(aes(max(H), min(Y.plim), label = Hpart.notes[1] ), 
						   size = 3, hjust = 1, vjust = -1.3, parse = TRUE) +
				 geom_text(aes(max(H), min(Y.plim), label = Hpart.notes[2] ), 
						   size = 3, hjust = 1, vjust = 0, parse = TRUE) +
				 geom_text(aes(max(H), min(Y.plim), label = Hpart.notes[3] ), 
						   size = 3, hjust = 1, vjust = -2.5, parse = TRUE) +
                 xlab(H2O.labpart) + 
                 ylab(Cells.labpart) + 
                 scale_y_log10(limits = Y.plim)
H.part.plot <- H.part.plot + geom_line(aes(y=fit), size=1, lty=1, colour="#CC0000") +
                 geom_line(aes(y=lower), size=0.5, lty=2, colour="#CC0000") + 
                 geom_line(aes(y=upper), size=0.5, lty=2, colour="#CC0000")


##______________________________________________________________
## Partial Regression: TAN
## adding plotMath to a ggplot graph: https://groups.google.com/forum/?fromgroups#!topic/ggplot2/-Ind8XDqaPQ
Npart.notes <- RegPlot.annote(Y.N)
## Back-transform: this is simple 10^x, despite the original log10(x+1) transformation.  Subtracting the 1 leads to -ve values, which won't be plotted on log-axes :(
Y.N.df <- within(Y.N.df, {
                 N <- 10^(N)
                 Y <- 10^(Y)
                 fit <- 10^(fit)
                 lower <- 10^(lower)
                 upper <- 10^(upper)
})

N.part.plot <- ggplot(data=Y.N.df, aes(x=N, y=Y)) +
                 geom_point(size=3, pch=20) + jaw.ggplot()   +
				 geom_text(aes(max(N), min(Y.plim), label = Npart.notes[1] ), 
						   size = 3, hjust = 1, vjust = -1.3, parse = TRUE) +
				 geom_text(aes(max(N), min(Y.plim), label = Npart.notes[2] ), 
						   size = 3, hjust = 1, vjust = 0, parse = TRUE) +
				 geom_text(aes(max(N), min(Y.plim), label = Npart.notes[3] ), 
						   size = 3, hjust = 1, vjust = -2.5, parse = TRUE) +
                 xlab(TAN.labpart) + 
                 ylab(Cells.labpart) + 
                 scale_y_log10(limits = Y.plim) + scale_x_log10() 
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
         width = 4, height = 4, scale = 1)
  ggsave(filename = sprintf("%sImportance2.eps", Fig.filename), plot = Y.importance2, 
         width = 4, height = 4, scale = 1)
  ggsave(filename = sprintf("%sEstimates.eps",	 Fig.filename), plot = Y.est2, 
         width = 4, height = 4, scale = 1)
  ## Effects plots
  ggsave(filename = sprintf("%sChamber.eps",	 Fig.filename), plot = T.plot, 
         width = 4, height = 4, scale = 1)
  ggsave(filename = sprintf("%sFrag.eps",	 Fig.filename), plot = F.plot, 
         width = 4, height = 4, scale = 1)
  ggsave(filename = sprintf("%sH2O.eps",	 Fig.filename), plot = H.plot, 
         width = 4, height = 4, scale = 1)
  ggsave(filename = sprintf("%sCxH2O.eps",	 Fig.filename), plot = TH.plot, 
         width = 4, height = 4, scale = 1)
  ##   ggsave(filename = sprintf("%sCPxH2O.eps",	 Fig.filename), plot = TPH.plot, 
  ##          width = 4, height = 4, scale = 1.5)
  ggsave(filename = sprintf("%sTAN.eps",	 Fig.filename), plot = N.plot, 
         width = 4, height = 4, scale = 1)
  ##   ggsave(filename = sprintf("%sCPxTAN.eps",	 Fig.filename), plot = NTP.plot, 
  ##          width = 4, height = 4, scale = 1.5)
  ## Partial Regression plots
  ggsave(filename = sprintf("%sTAN-partial.eps",   Fig.filename), plot = N.part.plot, 
         width = 4, height = 4, scale = 1)
  ggsave(filename = sprintf("%sH2O-partial.eps",   Fig.filename), plot = H.part.plot, 
         width = 4, height = 4, scale = 1)
} else {
  print(T.plot)
  print(H.plot)
  print(TH.plot)
  ##   print(TPH.plot)
  print(N.plot)
  ##   print(NTP.plot)
  ## Partial Regression plots
  print(H.part.plot)
  print(N.part.plot)
}
cat("- Finished Model Fitting:", Y.col, "-\n")
