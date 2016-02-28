################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Modelling: GLMM (regression)
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2012-11-04
################################################################
## INITIALISE
################################################################
## Working Directory: see lib/init.R below [\rd in Vim]
if (FALSE) {  # do not run automatically
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("..")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place

  ## Load data, functions, etc.  Process data & setup config. values.  
  ## Includes rm(list=ls()) to clear memory
  source('./Nfix-Cyanobacteria/3_ARA-cb_setup-Models.R')
  source('./Nfix-Cyanobacteria/3-x0_ARA-cb_setup-Models.R')
  Save.results  <- FALSE
}

library(car)
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme



################################################################
## PROCESS DATA: planned
################################################################
## generate grid to add predicted values to (X-values in all combinations of factors).
## - watch length.out: if it's too long, R will choke.
## - real replication is 1 anyway, so it doesn't need to be so big in this case.
Y.pred <- expand.grid(Block    = levels(SECCa$Block) , 
                      Time     = levels(SECCa$Time) , 
                      Chamber  = levels(SECCa$Chamber) , 
                      Frag     = levels(SECCa$Frag), 
                      Position = levels(SECCa$Position), 
                      X.trans  =seq(0, max(SECCa$X.trans), length.out=3 ),
                      H2O      =seq(0, max(SECCa$H2O), length.out=3 ) 
                      )
## levels(Y.pred$Frag)[1] <- "Continuous" # glmulti models fitted using old level names - won't correct for old factor names in saved models :(



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

## Reasonable model of effects & Interactions for assessment (thanks John Connolly)
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
  ## method="l" (leaps) uses a *much* faster algorithm, but not with factors or interactions :(
  ARA.glmulti1 <- glmulti(Y.fixed, data=SECCa, crit=aic, level=1, fitfunc=lm, 
                          marginality = TRUE, method="h", confsetsize=256, 
                          plotty=FALSE, report=FALSE)
  ## level=2 for 2-way interactions
  ## 2^(7 + choose(7,2) ) = 268,435,456 candidate models with 2-way interactions!!
  ## 416,869,200 according to glmulti (without marginality); 286,192,513 with marginality.
  ## An exhaustive exploration would take ~ 1 month on my computer.  
  ## Try with 'genetic' algorithm to speed things up (method="g"). ~ 30 minutes
  ## Best to run 2-4+ replicate genetic algorithms, and take consensus.
  ## use method="d" to print a summary of candidate models (no fitting)
  ## larger confsetsize -> more memory useage?
  ARA.multi <- list()
  for (i in 1:4) {
    cat("\n================ glmulti: Genetic Algorithm Run", i, "================\n\n")
    ARA.multi[[i]] <- glmulti(Y.fixed, data=SECCa, crit=aic, level=2, fitfunc=lm, 
                              marginality = TRUE, method="g", confsetsize=256, 
                              plotty=FALSE, report=TRUE)
  }
  ARA.glmulti2 <- consensus(ARA.multi, confsetsize=256) # more models use more memory!

  rm(ARA.multi)                        # clean-up

  if (F) save(ARA.glmulti1, ARA.glmulti2, file=Save.glmulti)

  ## requires glmulti objects
  ## capture output to save before deleting glmulti objects
  ARA.pmulti1 <- capture.output(print(ARA.glmulti1))
  ARA.pmulti2 <- capture.output(print(ARA.glmulti2))
  ARA.best2 <- as.formula(summary(ARA.glmulti2)$bestmodel)
  ARA.best2lm <- lm(ARA.best2, data=SECCa)
  summary(ARA.best2lm)


  term.labels <- function(tls) {
    tls <- gsub("X.trans", "Cyanobacteria", tls) # attr(SECC, "labels")[X.col]
    tls <- gsub("I(H2O^2)", "H2O^2", tls, fixed=TRUE)
    tls <- gsub("H2O", "Moisture", tls)
    tls <- gsub("Frag", "Fragmentation", tls)
    if (FALSE) { # ggplot2 doesn't support math expressions (yet)?
      tls <- gsub(":", "%*%", tls)       # interactions
      tls <- expression(tls)             # expressions
    } else {
    }
    tls
  }
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
    levels(glm.coef$Term) <- gsub("ChamberFull Chamber", "Chamber", levels(glm.coef$Term))
    levels(glm.coef$Term) <- gsub("PositionOuter", "Position", levels(glm.coef$Term))
    levels(glm.coef$Term) <- gsub("Block(.*)", "Block \\1", levels(glm.coef$Term))
    levels(glm.coef$Term) <- gsub("Frag(.*)", "Frag (\\1)", levels(glm.coef$Term))
    levels(glm.coef$Term) <- gsub("Time(.*)", "Time (\\1)", levels(glm.coef$Term))
    levels(glm.coef$Term) <- term.labels(levels(glm.coef$Term))
    glm.coef
  }
  importance.glmulti <- function(x) {
    ## collect importances (see plot.glmulti( type="s" )
    ww = exp(-(x@crits - x@crits[1])/2)
    ww = ww/sum(ww)
    clartou = function(x) {
      pieces <- sort(strsplit(x, ":")[[1]])
      if (length(pieces) > 1) 
        paste(pieces[1], ":", pieces[2], sep = "")
      else x
    }
    tet = lapply(x@formulas, function(x) sapply(attr(delete.response(terms(x)), 
                                                     "term.labels"), clartou))
    allt <- unique(unlist(tet))
    imp <- sapply(allt, function(x) sum(ww[sapply(tet, function(t) x %in% t)]))
    allt <- term.labels(allt)
    order.imp <- order(imp)
    glm.imp <- data.frame(Term=allt[order.imp], Importance=imp[order.imp], stringsAsFactors=F)
    glm.imp$Term <- factor(glm.imp$Term, levels=glm.imp$Term) # sorted levels?
    glm.imp
  }

  if (F) {
    summary(ARA.glmulti1)
    summary(ARA.glmulti2)
    names(ARA.glmulti2)
    str(summary(ARA.glmulti2))
    summary(ARA.glmulti1)$modelweights
    weightable(ARA.glmulti1)             # models, IC values, and relative weights (for confset)
  }

  ARA.coef1 <- getCoef.glmulti(ARA.glmulti1)
  ARA.coef2 <- getCoef.glmulti(ARA.glmulti2)
  nrow(ARA.coef2)                        # 129 terms (128 + intercept)!
  ARA.imp1 <- importance.glmulti(ARA.glmulti1)
  ARA.imp2 <- importance.glmulti(ARA.glmulti2)
  ## Extract model-averaged predictions: needs a lot of memory (>600 MB)
  ## May produce errors, but I don't really understand why :(
  Y.multipred <- predict(ARA.glmulti2, newdata=Y.pred, se.fit=TRUE)

  ## Output graphs
  par(mfrow=c(1,1))
  plot(ARA.glmulti1, type="p")           # red line at ~2 IC units above best model: should consider at least all models below this line.
  plot(ARA.glmulti1, type="w")           # red line where cumulative evidence weight = 95%
  ## plot(ARA.glmulti1, type="r")           # diagnostic plots (windows only?)
  plot(ARA.glmulti1, type="s")           # term weights (v 1.0.3)
  plot(ARA.glmulti2, type="s")           # term weights (v 1.0.3)
  plot(ARA.glmulti2, type="p")
  plot(ARA.glmulti2, type="w")
  ## sum of relative evidence weights over all models that include each term?
  barplot(ARA.coef1[, "Importance"], horiz=TRUE, names.arg=ARA.coef1$Term, las=2) 



  ## save derivative objects to speed up loading for future analysis.  
  ## The raw glmulti objects make for a big file (and a lot of memory): >1 GB!
  save(ARA.pmulti1, ARA.pmulti2, ARA.best2, ARA.best2lm, 
              ARA.coef1, ARA.coef2, ARA.imp1, ARA.imp2, 
              ## ARA.importance1, ARA.est1, ARA.importance2, ARA.est2,
              Y.multipred, file=Save.glmulti)

  rm(ARA.glmulti1, ARA.glmulti2)         # save memory? not right away, but maybe eventually :(
}

## ggplot2: theme settings - in jaw.graph_functions.R
bar.import <- function(glm.imp, imp.line=0.8) 
{  # horizontal bar graph of importance of model terms from glmulti
  Terms <- levels(glm.imp$Term)
  Terms <- gsub("\\s?:\\s?", " %*% ", Terms)
  ##     glm.imp <- importance.glmulti(glmObj)
  ## include passed parameter, so the actual value is included in the ggplot call, not the variable name 
  ## (ggplot functions build structures, they aren't evaluated until a print() call to actually draw the graph) :P
  geom_impline <- eval(substitute( geom_hline(aes(yintercept = imp.line), colour="#000000", lty=3) )) 
  ggplot(glm.imp, aes(x=Term, y=Importance), stat="identity", xlab = "Model Terms") +
  list(geom_bar(colour="#333333", fill="#999999"), coord_flip(), 
       scale_y_continuous(expand=c(0, 0)), 
       eval(substitute( scale_x_discrete(expand=c(0.01, 0),
                                         breaks=levels(glm.imp$Term), 
                                         labels=parse(text=Terms)) )),
       geom_impline,
       labs(x=NULL, y="Importance"),
       opts(title = "Model-averaged importance of effects",
            plot.title = theme_text(size = 16, lineheight = 1.2, face = "bold")#,
            ##panel.border=theme_blank(), axis.line=theme_segment(),
            ##axis.text.x=theme_text(size=10),
            ##axis.title.y=theme_text(hjust=0.6)
            )
       )
}
est.confint <- function(glmObj) 
{
  Terms <- levels(glmObj$Term)
  Terms <- gsub("-", "*'-'*", Terms)
  Terms <- gsub("\\s", "~", Terms)
  Terms <- gsub("\\s?:\\s?", " %*% ", Terms)
  ##   conf.wd <- glmObj[, "+/- (alpha=0.05)"]
  ##   conf.int <- aes(ymax = Estimate + conf.wd, ymin = Estimate - conf.wd)
  ggplot(glmObj, aes(x=Term, y=Estimate) ) +
  geom_hline(yintercept=0, colour="grey") + 
  list(geom_point(),
       geom_errorbar(aes(ymax = Emax, ymin = Emin), width=0.2),
       eval(substitute( scale_x_discrete(expand=c(0.01, 0),
                                         breaks=levels(glmObj$Term), 
                                         labels=parse(text=Terms)) )),
       coord_flip())
}

## MAIN effects
ARA.importance1 <- bar.import(ARA.imp1) + jaw.ggplot()
ARA.est1 <- est.confint(ARA.coef1) + jaw.ggplot()
print(ARA.importance1)
print(ARA.est1)

## ALL model terms
ARA.importance2 <- bar.import(ARA.imp2, imp.line = 0.5) + jaw.ggplot() 
## remove grid lines for Oecologia
ARA.importance2 <- ARA.importance2 + 
  opts(panel.grid.major=theme_blank(), panel.grid.minor=theme_blank(),
       axis.text.y=theme_text(size=8, hjust=1),
       axis.title.x=theme_text(size=12, hjust=0.75))
ARA.coef2plot <- ARA.coef2[ARA.coef2$Importance>=0.5, ] #  sharp jump from 0.2->0.3->0.99
ARA.coef2plot$Term <- factor(ARA.coef2plot$Term, levels=unique(ARA.coef2plot$Term))
ARA.est2 <- est.confint(ARA.coef2plot) + jaw.ggplot() + 
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




################################################################
## ADD MIXED EFFECTS?
################################################################
## Add mixed effects extensions to avoid violating model assumptions?
library(nlme)
lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 500, msMaxIter = 100, opt="optim")

Y.gls <- gls(Y.formula, data=SECCa, method="REML")
Y.lmeBT <- gls(formula(Y.formula), data=SECCa, weights=varIdent(form = ~ 1 | Block * Time), 
             control=lmc, method="REML")
Y.lme <- gls(Y.main, data=SECCa, weights=varIdent(form = ~ 1 | Block * Time), 
             control=lmc, method="REML")

anova(Y.gls, Y.lmeBT)                  # Do random effects improve the model?
AIC(Y.gls, Y.lmeBT, Y.lme)             # By how much do random effects improve the model?

## effect() produces errors if model is not specified explicitly :(
## Error in x$formula : object of type 'symbol' is not subsettable
## Wrapping the formula reference in formula() seems to help, 
## but then it won't work with avPlots :(

if (UseMM) Y.model <- Y.lmeBT




################################################################
## CHECK ASSUMPTIONS: MODEL VALIDATION
################################################################
cat("- Validating model\n")
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
if (inherits(Y.model, "lm")) residualPlots(Y.model)                 # car: lm only

## Check for spatial patterns in residuals?

## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)



################################################################
## ANALYSIS: GET RESULTS
################################################################
cat("- Generating results & predictions\n")
## the anova() function performs sequential (Type I) tests: order matters.
summary(Y.model)
anova(Y.model)                         # ORDER MATTERS! (see Zuur et al. 2009 pg. 540))
if (inherits(Y.model, "lm")) {
  Anova(Y.model, type=2)                 # Type II: car package**
  ## effects of single-term deletions?
  drop1(Y.model)
  drop1(Y.lmain)
}

## Partial effects of each variable (Zuur et al. 2009, pg. 400)
## termplot() only works for main effects, not when interactions are present :(
## termplot(Y.model, se=T, rug=T, partial.resid=T)
## Use effect() to get predicted effects for specific terms
## summary(eff.obj$term) : $effect, $lower, $upper for 95% CI
library(effects)
if (F) plot(allEffects(Y.model), ask=FALSE)   # interesting, but messy (and slow for gls)

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


##==============================================================
## Predictions
##==============================================================
if (inherits(Y.model, "lm")) {
  Y.pred.response <- predict(Y.model, newdata=Y.pred, type="response", interval="confidence", level=0.95) # newdata must have same explanatory variable names for predict to work.
  Y.pred$predicted <- Y.pred.response[, 1]
  Y.pred$pred.lwr  <- Y.pred.response[, 2]
  Y.pred$pred.upr  <- Y.pred.response[, 3]
  Y.pred.terms <- predict(Y.model, type="terms")
  ## type=="response" for full predictions (including interactions)
  ## type=="terms" for partial predictions (?)
}

## Note: there is a predict() method for glmulti objects...
Y.pred$multi.fit <- Y.multipred$averages[1]
Y.pred$multi.lwr <- Y.multipred$averages[1] - Y.multipred$variability[, "+/- (alpha=0.05)"]
Y.pred$multi.upr <- Y.multipred$averages[1] + Y.multipred$variability[, "+/- (alpha=0.05)"]



##==============================================================
## Partial regression
##==============================================================
## gls (Mixed Effects) model causes problems in this section
Y.lmodel <- if (inherits(Y.model, "lm")) Y.model else ARA.best2lm
avPlots(Y.lmH, terms= ~ X.trans * I(H2O^2), ask=FALSE) # car
avPlots(Y.lmodel, terms= ~ X.trans * I(H2O^2), ask=FALSE) # car

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
ARA.part <- paste("update(Y.lmodel, .~ ", RHS.part, ")", sep="")
cb.part  <- paste("update(Y.lmodel, X.trans ~ ", RHS.part, ")", sep="")
ARA.part <- eval(parse(text=ARA.part)) # problems fitting with gls :(
cb.part  <- eval(parse(text=cb.part))

ARA.re   <- resid(ARA.part, type = "response")
cb.re    <- resid(cb.part,  type = "response")
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

residualPlots(ARA.cb)                 # car

summary(ARA.cb)                        # R^2 = 0.04 ! :(
ARA.cb.r2 <- format(summary(ARA.cb)$adj.r.squared, digits=2)
ARA.cb.df <- data.frame(Cells=cb.re, ARA=ARA.re, fit=ARA.cb.pred[, "fit"], 
                        lower=ARA.cb.pred[, "lwr"], upper=ARA.cb.pred[, "upr"])



################################################################
## SAVE OUTPUT
################################################################
cat("- Saving Results (if requested)\n")
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  capture.output(cat(Save.head.txt), 
                 print(anova.full),    # Full model: variance partitioning
				 cat("\n\n"),                      # for output
                 cat(ARA.pmulti2, fill=TRUE), # multi-model selection
                 print(ARA.imp2),      # model-averaged estimates & weights
				 cat("\n\n"),                      # for output
				 print(Y.formula),                 # model
                 if (inherits(Y.model, "lm")) Anova(Y.model) else anova(Y.model),
				 summary(Y.model),                 # model summary
				 cat("\n\n"),                      # for output
				 Anova(ARA.cb),                    # partial regression
				 summary(ARA.cb),                  # model summary
				 cat("\n\n"),                      # for output
				 cat(Save.end.txt),                # END OUTPUT #
				 file = Save.text
				)
}


################################################################
## FINAL GRAPHICS
################################################################
cat("- Final Graphics\n")
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




##==============================================================
## Plot fitted on observed, by factor?
##==============================================================
Chamber.label <- attr(SECC, "labels")[["Chamber"]]
ChamberPts  <- ggPts.SECC(Chamber.map, Chamber.label) 
TopLegend   <- opts(legend.position = "top", legend.direction = "horizontal")
## Axis Labels: could also use X.plotlab, and Y.plotlab, but this is older code - and I'm back-transforming :P
## X.label <- paste("\"", attr(SECCa, "labels")[X.col], " \"*log[10](", attr(SECCa, "units")[X.col], ")", sep="")
## Y.label <- paste("\"", attr(SECCa, "labels")[Y.col], " \"*log[10](", attr(SECCa, "units")[Y.col], ")", sep="")
X.label <- paste("\"", attr(SECCa, "labels")[X.col], " \"*(", attr(SECCa, "units")[X.col], ")", sep="")
Y.label <- paste("\"", attr(SECCa, "labels")[Y.col], " \"*(", attr(SECCa, "units")[Y.col], ")", sep="")
X.label <- parse(text=X.label)
Y.label <- parse(text=Y.label)
XY.axlab <- list( xlab(bquote(.(X.label))), ylab(bquote(.(Y.label))) )

## (back)transformation functions
log0 <- function(x, base = exp(1))
{   # log-transform, but 0 stays as 0
  y <- log(x, base)
  y[x <= 0] <- 0
  y
}
log1p <- function(x, base = exp(1))
{   # OVERRIDING INTERNAL DEFAULT, JUST FOR ggplot2 GRAPHS
  log0(x, base)
}
alog0 <- function(y)
{   # back-transform log(X), where 0s stay as 0s
  x <- 10^y
  ##   x[y == 0] <- NA                      # MISSING (these cause glitches in transformed axes in ggplot)
  x
}

## functions for extracting data, and adding predictions to ggplot2 graphs
df.eff <- function (effObj, fun.trans = NULL, cols.trans = "numeric") {
  eff.df <- data.frame(fit     = effObj$fit,
                       lower   = effObj$lower,
                       upper   = effObj$upper)
  for (v in names(effObj$x)) {
    eff.df[, v] <- effObj$x[, v]       # by name, not index ;)
  }
  ## apply (back) transformation function
  if (!is.null(fun.trans)) 
  {
    if (all(cols.trans == "numeric")) cols.trans <- which( sapply(eff.df, class) %in% cols.trans )
    eff.df[, cols.trans] <- do.call(fun.trans, list(eff.df[, cols.trans]))
  }
  eff.df
}
eff.H.layers <- function(effObj, conf.int=TRUE, colour="black", bin.name="H2Obin", bin.lvls=NULL, bin.var="H2O") {
  eff.df <- df.eff(effObj)
  eff.df <- if (class(effObj) == "eff") df.eff(effObj, fun = alog0, cols.trans = c("fit", "lower", "upper", "X.trans")) else effObj
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
  eff.df <- if (class(effObj) == "eff") df.eff(effObj, fun = alog0, cols = c("fit", "lower", "upper", "X.trans")) else effObj
  ## no $ in aes() !!!!
  result <- list(geom_line(data=eff.df, aes(x=X.trans, y=fit,  group=Chamber, colour=Chamber, size = Chamber)),
                 scale_size_manual(name = Chamber.label,
                                   values = Chamber.map$lwd*0.5, 
                                   breaks = levels(Chamber.map$label)
                                   )
                 )
  if (conf.int == TRUE) {
    result <- c(result, 
                geom_line(data=eff.df, aes(x=X.trans, y=lower, group=Chamber, colour=Chamber, size = Chamber), lty=2), 
                geom_line(data=eff.df, aes(x=X.trans, y=upper, group=Chamber, colour=Chamber, size = Chamber), lty=2) 
                )
  }
  result
}
eff.F.layers <- function(effObj, conf.int=TRUE, boxes=TRUE) { # for factors
  eff.df <- if (class(effObj) == "eff") df.eff(effObj, fun = alog0, cols = c("fit", "lower", "upper")) else effObj
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

scale_y_ARA <- function()
{
  require(ggplot2)
  scale_y_continuous(trans = "log10", 
                     breaks = c(0 , 2  , 4  , 6  , 8  , 10 , 20 , 40 , 60 , 80 , 100 , 200 , 400 , 600 , 800 , 1000), 
                     labels = c(0 , "" , "" , "" , "" , 10 , "" , "" , "" , "" , 100 , ""  , ""  , ""  , ""  , 1000),
                     minor_breaks = c(2, 4, 6, 8, 20, 40, 60, 80, 200, 400, 600, 800) # minor_breaks having no effect :(  simulated with blank labels above
                     )
}

## extract, and backtransform data (moot)
## Plot.data <- SECCa[, strsplit("SampleID, Block, Time, Chamber, Frag, Pos, Position, X, X.trans, Y, Y.trans", ", ")[[1]]]

### with ggplot(), and back-transformed?
ARA.plot <- ggplot(data=SECCa, aes(x=X, y=Y)) +
            geom_point(size = 3, aes(group = Chamber, colour = Chamber, shape = Chamber)) + 
            ##             coord_trans(x = "log10", y = "log10") + 
            ##             coord_trans(x = "log1p", y = "log1p") + 
            ##             scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
            scale_x_log10() + scale_y_ARA() +
            XY.axlab + jaw.ggplot() + ChamberPts + TopLegend
ARA.C.plot     <- ARA.plot + eff.C.layers(X.eff)
ARA.facet.plot <- ARA.plot + eff.C.layers(X.TCP.eff) + facet_grid(facets=Position~Time)
ARA.Block.plot <- ARA.plot + eff.C.layers(X.BTC.eff, conf=F) + facet_grid(facets=Block~Time)

ARA.Frag.plot <- ARA.plot + aes(x=Frag, y=Y) + XY.axlab + xlab("Fragmentation Treatment")
ARA.Frag.plot <- ARA.Frag.plot + 
##   coord_trans(y = "log10") + scale_x_continuous() + scale_y_continuous() + 
  scale_x_continuous() + scale_y_ARA() +
  geom_jitter(size=3, aes(group=Chamber, colour=Chamber, shape=Chamber)) + # optional
  eff.F.layers(F.eff) + facet_grid(facets=.~Time)

## Load eps graphics for plot labels
library(grImport)
FragIcons <- SECCicons()[1:4]       # Frag 1:4

## Add imported graphics as x-axis tick labels :D
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
ARA.Frag.plot <- ARA.Frag.plot + scale_x_discrete(labels = names(FragIcons),
                                                  breaks = levels(SECCa$Frag)) +
opts(axis.ticks.margin = unit(0.2, "lines"),
     axis.text.x = picture_axis(FragIcons, icon.size = unit(1.4, "lines")) 
)


SECCa$H2Obin9 <- cut(SECCa$H2O, breaks=H2O.breaks9)
SECCa$H2Obin4 <- cut(SECCa$H2O, breaks=H2O.breaks4)
SECCa$H2Obin4 <- factor(SECCa$H2Obin4, levels=rev(levels(SECCa$H2Obin4)) ) # arrange in decreasing order, for top-down rows in plot
ARA.H2O.plot <- ggplot(data=SECCa, aes(x=X, y=Y)) +
                geom_point(size = 3, aes(group = Chamber, 
                               colour = Chamber, shape = Chamber)) + 
                scale_x_log10() + scale_y_ARA() +
                XY.axlab + jaw.ggplot() + ChamberPts + TopLegend
ARA.HB.plot  <- ARA.H2O.plot + 
                eff.H.layers(X.HB.eff, bin.name="H2Obin4", bin.lvls=levels(SECCa$H2Obin4)) + 
                facet_grid(H2Obin4 ~ Block)
ARA.H2O.plot <- ARA.H2O.plot + 
                eff.H.layers(X.H.eff, bin.name="H2Obin9", bin.lvls=levels(SECCa$H2Obin9)) + 
                facet_wrap(~ H2Obin9)

## Partial Regression graph
##library(scales)
ARA.part.plot <- ggplot(data=as.data.frame(alog0(ARA.cb.df)), aes(x=Cells, y=ARA)) +
                 geom_point(size=3, pch=20, colour="black", fill="#999999") + jaw.ggplot()   +
                 ##geom_point(size=2, pch=21, colour="black", fill="#999999") + jaw.ggplot()   +
                 xlab( bquote(paste( italic("Residual "), .(attr(SECC, "labels")[[X.col]]) ))#, 
                                    ##" ", (.(attr(SECC, "units")[[X.col]])), "" ))   # different parens in text vs. literal :/
                 ) + 
                 ylab( bquote(paste( italic("Residual "), .(attr(SECC, "labels")[[Y.col]]) ))#, 
                                    ##" ", (.(attr(SECC, "units")[[Y.col]])), "" )) # residuals of a log-linear regression = ratio to geometric mean, not original units?
                 ) 
ARA.part.plot <- ARA.part.plot + geom_line(aes(y=fit), size=1, lty=1, colour="#990000") +
                 geom_line(aes(y=lower), size=0.5, lty=2, colour="#990000") + 
                 geom_line(aes(y=upper), size=0.5, lty=2, colour="#990000")
ARA.part.plot <- ARA.part.plot + scale_y_log10(limits=10^c(-2, 1.5)) + scale_x_log10()
## I wanted to format the numbers in the exponent, but could not figure out how. mathformat() (library(scales)) keeps returning an error.
##ARA.part.plot <- ARA.part.plot + scale_y_log10(limits=10^c(-2, 1.5), breaks=10^seq(-2, 2, by=0.5), labels=parse(text=paste("10^", format(seq(-2, 2, by=0.5), digits=2), sep="")) ) + scale_x_log10()
## Remove grid for Oecologia
ARA.part.plot <- ARA.part.plot +
opts(panel.grid.major = theme_blank(),
     panel.grid.minor = theme_blank()
)


Save.plot.dir <- "./graphs/"               # for output
if (Save.results == TRUE) {
  ggsave(filename=paste(Fig.filename, "-Importance.eps", sep=""), 
         plot = ARA.importance2, width=4, height=5, scale=1.2)
  ggsave(filename=paste(Fig.filename, "-Estimates.eps", sep=""), 
         plot = ARA.est2, width=4, height=6, scale=1.2)
  ggsave(filename=paste(Fig.filename, "*Time*Position*Chamber.eps", sep=""), 
         plot = ARA.facet.plot, width=5, height=4, scale=1.5)
  ggsave(filename=sprintf("%sFigure-%s~%s", SaveDir.plots(), Y.col, "Frag.eps"),
         plot = ARA.Frag.plot, width=4, height=2, scale=2)
  ggsave(filename=paste(Fig.filename, "*Chamber.eps", sep=""), 
         plot = ARA.C.plot, width=4, height=4, scale=1.5)
  ggsave(filename=paste(Suppl.filename, "*Block*Time.eps", sep=""), 
         plot = ARA.Block.plot, width=4, height=6, scale=1.5)
  ggsave(filename=paste(Suppl.filename, "*H2O.eps", sep=""), 
         plot = ARA.H2O.plot, width=4, height=4, scale=1.5)
  ggsave(filename=paste(Suppl.filename, "*H2O*Block.eps", sep=""), 
         plot = ARA.HB.plot, width=8, height=6, scale=1.5)
  ggsave(filename=paste(Fig.filename, "-partial.eps", sep=""), 
         plot = ARA.part.plot, width=4, height=4, scale=1.5)
} else {
  print(ARA.importance2)
  print(ARA.est2)
  print(ARA.C.plot)
  print(ARA.facet.plot)
  print(ARA.Block.plot)
  print(ARA.Frag.plot)
  print(ARA.H2O.plot)
  print(ARA.HB.plot)
  print(ARA.part.plot)
}
cat("- Finished Model Fitting -\n")
