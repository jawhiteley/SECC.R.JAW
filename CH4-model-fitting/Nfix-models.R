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


library(nlme)
################################################################
## MODEL FORMULA
################################################################
Y.main   <- Y.trans ~ Block + Frag + TempC + H2O + I(H2O^2) + logCells + log10(TAN)
Y.fixed  <- Y.trans ~ Block + Frag * TempC * H2O * I(H2O^2) * logCells * log10(TAN)
Y.full   <- Y.fixed                    # not enough replication to test full range of interactions
Y.random <-  ~ 1 | Block/Chamber/Frag 

## generate grid to add predicted values to (X-values in all combinations of factors).
## - watch length.out: if it's too long, R will choke.
## - real replication is 1 anyway, so it doesn't need to be so big in this case.
Y.pred <- expand.grid(Block    = levels(SECCa$Block) , 
                      Frag     = levels(SECCa$Frag), 
                      logCells = seq(0, max(SECCa$logCells), length.out=3 ),
                      TempC    = seq(min(SECCa$TempC), max(SECCa$TempC), length.out=2 ),
                      H2O      = seq(0, max(SECCa$H2O), length.out=3 ), 
                      TAN      = seq(0, max(SECCa$TAN), length.out=3 ) 
                      )


cat("- Fitting models\n")
################################################################
## MODEL FITTING
################################################################
## Transformations & Linearity
##==============================================================
## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)
library(mgcv)
Y.gam <- gam(Y.trans ~ TempC + s(H2O) + s(Cells.m) + s(TAN) + Frag + Block, data = SECCa)
Y.log.gam <- gam(Y.trans ~ TempC + s(H2O) + s(logCells) + s(log10(TAN)) + Frag + Block, data = SECCa)
anova(Y.gam)
anova(Y.log.gam)
AIC(Y.gam, Y.log.gam)
## Cells.m is linear (no need to transform?), TAN might be, H2O is NOT
## I log-transformed Cells in original analysis: I should probably continue to do so, unless I have a really good reason not to.
## The non-linearity in H2O is driven largely by a big gap between very dry and other patches.
## log(TAN) looks a little smoother, but I'm not sure I can justify it biologically? 
##   (regardless, the distribution is highly skewed)

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
Y.lmfix <- lm( Y.fixed , data=SECCa)   # n = 2 for all interaction combinations (perfect fit)
anova.full <- anova(Y.lmfix)
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
  Y.lmeB <- lme(Y.main, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method = "REML")
  anova(Y.gls, Y.lmB, Y.lmeB)          # allowing heterogeneity among Blocks is an improvement (p = 0.004).
  Y.lmBC <- gls(Y.main, weights = varIdent(form = ~ 1 | Block * Chamber), data = SECCa, method = "REML")
  anova(Y.lmB, Y.lmBC)                 # allowing heterogeneity among Blocks AND Chambers is not a sig. improvement (p = 0.088)
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

if (file.exists(Save.glmulti)) { ## load saved object to speed things up
  cat("- loading glmulti objects from previous run\n")
  load(Save.glmulti)
} else {
  ## Note: glmulti performs an exhaustive search of all candidate models, by default.  
  ##       2^7 = 128 candidate models (not including interactions). 256 according to glmulti
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
                        weights = varIdent(form = ~ 1 | Block), control = lmc,
                        random = Y.random
                        )
  print(Y.glmultiBr)
  Y.glmulti1 <- Y.glmultiBr

  ## 2-WAY INTERACTIONS (level = 2)
  ## An exhaustive exploration would take ~ 1 month on my computer.  
  ## Try with 'genetic' algorithm to speed things up (method="g"). ~ 30 minutes
  ## Best to run 2-4+ replicate genetic algorithms, and take consensus.
  ## use method="d" to print a summary of candidate models (no fitting)
  ## larger confsetsize -> more memory useage
  if (FALSE)
  { ## ML estimation crashes with weights (heterogeneity) and all terms: false convergence / Singularity at backsolve
    Y.mm <- gls(Y.trans ~ 1 + Block + Frag + TempC + H2O + I(H2O^2) + logCells + 
                log10(TAN) + Frag:Block + log10(TAN):TempC + Block:TempC + 
                Block:I(H2O^2) + Block:log10(TAN) + Frag:TempC + Frag:H2O + 
                Frag:I(H2O^2) + Frag:logCells, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="ML")
  }    

  Y.multi <- list()
  for (i in 1:4) {
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
  summary(Y.best2lm)


  if (F) {
    summary(Y.glmulti1)
    summary(Y.glmulti2)
    names(Y.glmulti2)
    str(summary(Y.glmulti2))
    summary(Y.glmulti1)$modelweights
    weightable(Y.glmulti1)             # models, IC values, and relative weights (for confset)
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
              Y.coef2, Y.imp1, Y.imp2, 
              ## Y.coef1, Y.importance1, Y.est1, Y.importance2, Y.est2,
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
    tls <- gsub("Frag(.*)", "Frag (\\1)", tls)
    tls <- gsub("Time(.*)", "Time (\\1)", tls)
  }
  ## Term labels
  tls <- gsub("logCells", "Cyanobacteria", tls) # attr(SECC, "labels")[X.col]
  tls <- gsub("log10(TAN)", "Total N", tls, fixed=TRUE)
  tls <- gsub("I(H2O^2)", "H2O^2", tls, fixed=TRUE)
  tls <- gsub("TempC", "Temperature", tls, fixed=TRUE)
  tls <- gsub("H2O", "Moisture", tls)
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
## Frag:H2O
## Frag:H2O^2
## Block:logCells
## H2O:H2O^2    (WTF?)
## logTAN:TempC
## H2O:logCells
## Block:H2O
## Block:H2O^2
## H2O^2:logCells

## Implied higher-order interactions?
## Block * logCells * H2O * I(H2O^2)
## except that ML estimation will fail with these interactions :(
Y.fixHi <- Y.trans ~ Block * logCells * H2O * I(H2O^2) + TempC + Frag + log10(TAN) +
                     Frag:H2O + Frag:I(H2O^2) + TempC:log10(TAN)

## based on Y.best2, plus extra terms that have similar "importance" & are related (Block:H2O)
Y.fixed <- Y.trans ~ Block + Frag + TempC + H2O + I(H2O^2) + logCells + log10(TAN) +
                     Frag:H2O + Frag:I(H2O^2) + Block:logCells + H2O:I(H2O^2) +
                     log10(TAN):TempC + H2O:logCells + I(H2O^2):logCells + Block:H2O + Block:I(H2O^2)


##==============================================================
## MODEL FITTING: Final Fixed & Random Effects?
##==============================================================
Y.fxlm <- lm(Y.fixed, data = SECCa)
Y.fHlm <- lm(Y.fixHi, data = SECCa)
anova(Y.best2lm, Y.fxlm, Y.fHlm)
Y.f2   <- gls(Y.best2, data = SECCa, method = "ML")
Y.fx   <- gls(Y.fixed, data = SECCa, method = "ML")
Y.fH   <- gls(Y.fixHi, data = SECCa, method = "ML")
anova(Y.fx, Y.f2, Y.fH)                # Adding higher-order interactions is better (without mixed effects)
Y.fx   <- update(Y.fx, method = "REML")
Y.fH   <- update(Y.fH, method = "REML")
## I will likely need to account for heterogeneity among Blocks, maybe Chambers as well (TempC/Warming)
Y.mH <- gls(Y.fixHi, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
Y.me <- lme(Y.fixHi, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
Y.mr <- lme(Y.fixHi, random = Y.random, data = SECCa, method ="REML")
anova(Y.fH, Y.mH, Y.me, Y.mr)          # compare random structures (using REML)
anova(Y.fH, Y.me)
anova(Y.fH, Y.mr)
anova(Y.mH, Y.mr)                      # are these really nested?  (might not be a valid comparison)
## The nested structure of random effects is technically not necessary, but it did improve the AIC in glmulti (using ML), 
## and probably should be done on theoretical grounds, to account for the structure of the experiment.

## Compare fuller model with best model from glmulti
Y.mb2 <- lme(Y.best2, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="ML")
Y.mlf <- lme(Y.fixed, random = Y.random, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="ML")
## Y.mlh <- update(Y.me, method ="ML")    # X Singularity in backsolve at level 0, block 1 :(
anova(Y.mlf, Y.mb2)      # adding the extra interaction term does increase the AIC (worse), but not significantly so.
Y.mHf <- gls(Y.fixed, weights = varIdent(form = ~ 1 | Block), data = SECCa, method ="REML")
Y.mef <- update(Y.mlf, method = "REML")
anova(Y.mHf, Y.mef)

## SO, Which model should I use?
##  All I really need a model for at this point is predictions and graphs.
##  Ideally, I would prefer the full model with higher order interactions and all random effects (Y.me),
##    however, this model may be too complex: several coefficients have 0 df & p = NaN (see summary(Y.me))
##    Y.mel is a little better, and the gls versions don't seem to have this problem :/
##  Sadly, effect() throws an error with the mixed models (lme including nested effects & heterogeneity)
##  Basically, I can have either nested effects OR heterogeneity, but not both (to work with `effects`)
Y.fit  <- Y.me
Y.fit  <- Y.mH




################################################################
## CHECK ASSUMPTIONS: MODEL VALIDATION
################################################################
cat("- Validating models\n")

## residualPlots(Y.best2lm)               # car: lm only
diagnostics(Y.best2lm, X.cols = c("logCells", "H2O", "TAN")) # Best model from glmulti
diagnostics(Y.mHf, X.cols = c("logCells", "H2O", "TAN"))     # 2-way interactions and heterogeneity
diagnostics(Y.mH,  X.cols = c("logCells", "H2O", "TAN"))     # higher-order interactions and heterogeneity
diagnostics(Y.mr,  X.cols = c("logCells", "H2O", "TAN"))     # higher-order interactions and nested structure
diagnostics(Y.me,  X.cols = c("logCells", "H2O", "TAN"))     # higher-order interactions + mixed effects
## I do prefer the distribution and lack of patterns in residuals in the mixed effects models
## There are some disturbing patterns in the residuals for the lm model
## Allowing heterogeneity is a modest improvement: a nested random structure is about the same (but theoretically justified).
## There may also be some heterogeneity with H2O, but I might just have to live with it :(

RE <- diagnostics(Y.fit, resType = "pearson", 
                  X.cols = c("logCells", "H2O", "TAN"), more = TRUE) # Full diagnostics
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
BH.eff    <- effect("Block:H2O:I(H2O^2)", Y.fit)
FH.eff    <- effect("H2O:I(H2O^2):Frag", Y.fit)
CH.eff    <- effect("logCells:H2O:I(H2O^2)", Y.fit, xlevels=list(H2O=H2O.9lvls))
C.BH.eff  <- effect("Block:logCells:H2O:I(H2O^2)", Y.fit, 
                    xlevels=list(Block=1:8, H2O=H2O.4lvls))

plot(BH.eff, x.var = "H2O", ask = FALSE)
plot(FH.eff, x.var = "H2O", ask = FALSE)
plot(CH.eff, x.var = "logCells", ask = FALSE)
plot(C.BH.eff, x.var = "logCells", ask = FALSE)


##==============================================================
## Predictions
##==============================================================
## Note: there is a predict() method for glmulti objects...
Y.pred$multi.fit <- Y.multipred$averages[1]
Y.pred$multi.lwr <- Y.multipred$averages[1] - Y.multipred$variability[, "+/- (alpha=0.05)"]
Y.pred$multi.upr <- Y.multipred$averages[1] + Y.multipred$variability[, "+/- (alpha=0.05)"]


################################################################
## OUTPUT
################################################################
## Save Text Output


##==============================================================
## PUBLICATION GRAPHS
##==============================================================


