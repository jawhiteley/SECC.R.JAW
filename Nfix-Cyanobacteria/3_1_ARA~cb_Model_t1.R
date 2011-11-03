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
                Climate <- factor( paste(Position, Chamber) ) # psuedo-factor to simplify modelling: fewer interactions to deal with.
})
## use only data from t1
SECCt <- SECCa[SECCa$Time==levels(SECCa$Time)[1], ]

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
##  - Partial regression to remove effect of moisture FIRST, then model the resulting residuals
## H2O
## - does it affect N-fixation directly, or mediate the effect of cyanobacteria?

## Subsume Chamber & Position into "Climate" pseudo-treatment?
## Time as a fixed factor, or separate analysis on each Time?

##==============================================================
## Model Formula
##==============================================================
### Fixed effects
## Include H2O^2 to test for unimodal relationship?
## Including Time as a factor?
## - It might be better to analyse each time period separately, and avoid the higher-order interactions that I know exist.
  Y.fixed <- Y.trans ~ X.trans * H2O * Chamber * Frag * Position
  Y.fixCl <- Y.trans ~ X.trans * H2O * Climate * Frag
  ## Main effects only
  Y.main  <- Y.trans ~ X.trans + H2O + Chamber + Frag + Position
  Y.mainCl<- Y.trans ~ X.trans + H2O + Climate + Frag
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




################################################################
## GLMM - Hierarchical / Multilevel Mixed Model               **
################################################################
##==============================================================
## Model Fitting
##==============================================================
## Using a mixed-effects model to account for treatments of different sizes? 
##  Block / Chamber / Frag / Position
Y.fm <- gls(Y.fixed, data=SECCt, method="REML") # fixed effects only for comparison
lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000) # takes a while to converge...
Y.rim  <- lme(Y.fixed, random=Y.Ri,  data=SECCt, control=lmd, method="REML")
## SLOW! :(
## more complex models take a few minutes (or hours) to fit.  Use with discretion.
Y.rism <- lme(Y.fixed, random=Y.Ris, data=SECCt, control=lmc, method="REML")
Y.rie  <- lme(Y.fixed, random=Y.Ri,  
              weights=varIdent(form=~ 1 | Block),
              data=SECCt, control=lmc, method="REML")
if (UseClimateFac) {
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block),
                 correlation=corAR1(form = ~ 1 | Block/Frag/Climate),
                 data=SECCt, control=lmc, method="REML")
} else {
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 ##                  weights=varIdent(form=~ 1 | Block),
                 correlation=corAR1(form = ~ 1 | Block/Chamber/Frag),
                 data=SECCt, control=lmc, method="REML")
}

if (FALSE) {                           # these attempts are not working :(
  Y.rise <- lme(Y.fixed, random=Y.Ris,  
                weights=varIdent(form=~ 1 | Block),
                data=SECCt, control=lmc, method="REML")
  ## Error in getGroups.data.frame(data, form, level = length(splitFormula(grpForm,  : Invalid formula for groups
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block/Time/Chamber/Frag),
                 data=SECCt, control=lmc, method="REML")
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block/Time/Frag/Climate),
                 data=SECCt, control=lmc, method="REML")

  ## I want to use varConstPower (variance covariate with values of 0), but that produces an error (false convergence) :(
  Y.rieXP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varConstPower(form = ~ X.trans),
                 data=SECCt, control=lmd, method="REML")

  ## varFixed gives no errors, but takes >12 hours to fit (if at all)
  Y.rieXP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varFixed(~ X.trans),
                 data=SECCt, control=lmd, method="REML")
  ## varPower works, but not for X.trans, which has 0's!   :/
  Y.rieXP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varPower(form = ~ X.trans),
                 data=SECCt, control=lmd, method="REML")
  Y.rieHP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varPower(form = ~ H2O),
                 data=SECCt, control=lmd, method="REML")

  Y.riceX <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ X.trans)),
                 data=SECCt, control=lmc, method="REML")
  Y.riceH <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ H2O)),
                 data=SECCt, control=lmc, method="REML")
  Y.rice  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ X.trans),
                                 varConstPower(form=~ H2O)),
                 data=SECCt, control=lmc, method="REML")
}

## Main effects only: ignore interactions?
Y.mainMM <- lme(Y.main, data = SECCt, random = Y.Ri, method="REML")




##==============================================================
## MODEL SELECTION
##==============================================================
## RANDOM structure
anova(Y.fm, Y.rim)                     # do random effects improve the model?
anova(Y.fm, Y.rism)                    # do random effects improve the model?
anova(Y.fm, Y.rie)                     # do random effects improve the model?
anova(Y.rieN, Y.rie, Y.rim, Y.rism)    # do we need random slopes or error terms?
anova(Y.rieN, Y.rim)                   # do we need nested error terms?

Y.mm <- Y.rim                          # Optimal random structure
## The biggest improvement seems to come from random intercepts across Blocks.
## However, that model shows major heterogeneity in the residuals.
## I may need other random factors to maintain a valid model fit.

## optimize FIXED factors
if (FALSE) {
  Y.ml  <- lme(Y.fixed, data=SECCt, random=Y.Ri, method="ML") # re-fit with ML
  Y.rieML <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block),
                 data=SECCt, control=lmc, method="ML")
}
Y.ml  <- update(Y.mm, method="ML")     # re-fit with ML; some models can't be :(
Y.mainML <- update(Y.mainMM, method="ML")

drop1(Y.ml)
## Y.step <- step(Y.ml)                # stepwise back & forward model selection?  Not for lme
if (!UseClimateFac) {                  # All factors, or Climate pseudo-factor
  Y.ml1 <- update(Y.ml, .~. - X.trans:H2O:Chamber:Frag:Position) # <-
  anova(Y.ml, Y.ml1)
  ## 4-way interactions
  drop1(Y.ml1)
  Y.ml1.1 <- update(Y.ml1, .~. - X.trans:H2O:Chamber:Frag) # <-
  Y.ml1.2 <- update(Y.ml1, .~. - X.trans:H2O:Chamber:Position)
  Y.ml1.3 <- update(Y.ml1, .~. - X.trans:H2O:Frag:Position)
  Y.ml1.4 <- update(Y.ml1, .~. - X.trans:Chamber:Frag:Position)
  Y.ml1.5 <- update(Y.ml1, .~. - H2O:Chamber:Frag:Position)
  anova(Y.ml1, Y.ml1.1)
  anova(Y.ml1, Y.ml1.2)
  anova(Y.ml1, Y.ml1.3)
  anova(Y.ml1, Y.ml1.4)
  anova(Y.ml1, Y.ml1.5)
  Y.ml2 <- Y.ml1.1
  drop1(Y.ml2)
  Y.ml2.1 <- update(Y.ml2, .~. - X.trans:H2O:Chamber:Position)
  Y.ml2.2 <- update(Y.ml2, .~. - X.trans:H2O:Frag:Position)
  Y.ml2.3 <- update(Y.ml2, .~. - X.trans:Chamber:Frag:Position)
  Y.ml2.4 <- update(Y.ml2, .~. - H2O:Chamber:Frag:Position) # <-
  anova(Y.ml2, Y.ml2.1)
  anova(Y.ml2, Y.ml2.2)
  anova(Y.ml2, Y.ml2.3)
  anova(Y.ml2, Y.ml2.4)
  Y.ml3 <- Y.ml2.4
  drop1(Y.ml3)
  Y.ml3.1 <- update(Y.ml3, .~. - X.trans:H2O:Chamber:Position)  # *
  Y.ml3.2 <- update(Y.ml3, .~. - X.trans:H2O:Frag:Position)     # <-
  Y.ml3.3 <- update(Y.ml3, .~. - X.trans:Chamber:Frag:Position)
  anova(Y.ml3, Y.ml3.1)
  anova(Y.ml3, Y.ml3.2)
  anova(Y.ml3, Y.ml3.3)
  Y.ml4 <- Y.ml3.2
  drop1(Y.ml4)
  Y.ml4.1 <- update(Y.ml4, .~. - X.trans:H2O:Chamber:Position)  # *
  Y.ml4.2 <- update(Y.ml4, .~. - X.trans:Chamber:Frag:Position) # ~
  anova(Y.ml4, Y.ml4.1)
  anova(Y.ml4, Y.ml4.2)
  ## 3-way interactions?
  Y.ml4.3 <- update(Y.ml4, .~. - X.trans:H2O:Frag)
  Y.ml4.4 <- update(Y.ml4, .~. - H2O:Chamber:Frag)
  Y.ml4.5 <- update(Y.ml4, .~. - H2O:Frag:Position) # <-
  anova(Y.ml4, Y.ml4.3)
  anova(Y.ml4, Y.ml4.4)
  anova(Y.ml4, Y.ml4.5)
  Y.ml5 <- Y.ml4.5
  drop1(Y.ml5)
  Y.ml5.1 <- update(Y.ml5, .~. - X.trans:H2O:Frag) # <-
  Y.ml5.2 <- update(Y.ml5, .~. - H2O:Chamber:Frag)
  anova(Y.ml5, Y.ml5.1)
  anova(Y.ml5, Y.ml5.2)
  Y.ml6 <- Y.ml5.1
  drop1(Y.ml6)
  Y.ml6.1 <- update(Y.ml6, .~. - H2O:Chamber:Frag) # X
  anova(Y.ml6, Y.ml6.1)
  Y.mm <- update(Y.ml6, method="REML")


  ## No interactions?
  drop1(Y.mainML)
  Y.main1 <- update(Y.mainML, .~. - X.trans)  # *
  Y.main2 <- update(Y.mainML, .~. - H2O)      # *
  Y.main3 <- update(Y.mainML, .~. - Chamber)  # <-
  Y.main4 <- update(Y.mainML, .~. - Frag)     # *
  Y.main5 <- update(Y.mainML, .~. - Position)
  anova(Y.mainML, Y.main1)
  anova(Y.mainML, Y.main2)
  anova(Y.mainML, Y.main3)
  anova(Y.mainML, Y.main4)
  anova(Y.mainML, Y.main5)
  Y.m1 <- Y.main3
  drop1(Y.m1)
  Y.m1.1 <- update(Y.m1, .~. - X.trans)  # *
  Y.m1.2 <- update(Y.m1, .~. - H2O)      # *
  Y.m1.3 <- update(Y.m1, .~. - Frag)     # *
  Y.m1.4 <- update(Y.m1, .~. - Position) # <-
  anova(Y.m1, Y.m1.1)
  anova(Y.m1, Y.m1.2)
  anova(Y.m1, Y.m1.3)
  anova(Y.m1, Y.m1.4)
  Y.m2 <- Y.m1.4
  drop1(Y.m2)
  Y.m2.1 <- update(Y.m2, .~. - X.trans)  # *
  Y.m2.2 <- update(Y.m2, .~. - H2O)      # *
  Y.m2.3 <- update(Y.m2, .~. - Frag)     # *
  anova(Y.m2, Y.m2.1)
  anova(Y.m2, Y.m2.2)
  anova(Y.m2, Y.m2.3)
  Y.mainM <- Y.m2
  Y.mainMM <- update(Y.mainM, method="REML")

  ## Forward model selection to add interactions?
  Y.m1     <- update(Y.mainM, .~. + X.trans:H2O)
  Y.m2     <- update(Y.mainM, .~. + X.trans:Frag)
  Y.m3     <- update(Y.mainM, .~. + H2O:Frag)
  Y.m1.1 <- update(Y.mainM, .~. + X.trans:H2O
                     + X.trans:Frag + H2O:Frag
                     + X.trans:H2O:Frag )
  anova(Y.mainM, Y.m1)
  anova(Y.mainM, Y.m2)
  anova(Y.mainM, Y.m3)
  anova(Y.mainM, Y.m1.1)

  Y.mainMM <- update(Y.mainM, method="REML")

  AIC(Y.ml, Y.mainM)
  AIC(Y.mm, Y.mainMM)

} else {                               # Chamber & position lumped into Climate pseudo-factor

  AIC(Y.ml, Y.mainML)
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

diagnostics(Y.mm)
diagnostics(Y.rie)                     # Not Normal :(
diagnostics(Y.rieN)
diagnostics(Y.mainMM)                  # Fewest assumptions violated ?!
diagnostics(Y.m1.1)                    # main effects +interactions?

Y.model <- Y.mm
diagnostics(Y.model)



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
Y.pred <- expand.grid(Chamber  = levels(SECCt$Chamber) , 
                      Frag     = levels(SECCt$Frag), 
                      Position = levels(SECCt$Position), 
                      X=seq(0, max(SECCt$X), length.out=100 ) 
                      )
Y.pred$predicted <- predict(Y.model, newdata=Y.pred, type="response" )  # newdata must have same explanatory variable name for predict to work.

if (FALSE) {
  pred.Chamber <- expand.grid(Chamber = levels(SECCt$Chamber) , 
                              X=seq(0, max(SECCt$X), length.out=100 ) 
  )
  pred.Chamber$predicted <- predict(Y.model, newdata=pred.Chamber, type="response" )
  pred.Frag <- expand.grid(Frag = levels(SECCt$Frag) , 
                           X=seq(0, max(SECCt$X), length.out=100 )
  )
  pred.Frag$predicted <- predict(Y.model, newdata=pred.Frag, type="response" )
  pred.Position <- expand.grid(Position = levels(SECCt$Position) , 
                               X=seq(0, max(SECCt$X), length.out=100 ) 
  )
  pred.Position$predicted <- predict(Y.model, newdata=pred.Position, type="response" )
  pred.FxP <- expand.grid(Frag     = levels(SECCt$Frag), 
                          Position = levels(SECCt$Position), 
                          X=seq(0, max(SECCt$X), length.out=100 ) 
                          )
  pred.FxP$predicted <- predict(Y.model, newdata=pred.Position, type="response" )
}


Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
point <- 21	# 21 for circles with col & bg ; 16 for solid circles
Chamber.map$pch <- c(21, 16)  # use circles for both treatments

SECCt <- within( SECCt,{
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
               aggregate(cbind(predicted), list(Chamber = Chamber, X = X), mean)
)  # I should be getting direct predictions, not means of predictions. *****
with( SECCt,{
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

## lattice panels
print( xyplot( Y ~ X | Frag * Position , data=SECCt, 
              pch = SECCt$pt, col = SECCt$colr, 
              xlab = quote(X.plotlab), ylab = quote(Y.plotlab), 
              panel = function(..., data, subscripts) {
                panel.xyplot(...)  # regular plot of data points
                Frag.lvl <- unique(SECCt$Frag[subscripts]) # get current factor levels
                Pos.lvl  <- unique(SECCt$Position[subscripts])
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




if (FALSE) {
  ## Coplots with linear fits (from Zuur et al. 2007 Chapter 22 R code)
  ## individual lm's within each panel.  Not exactly what I want.
  coplot( Y ~ X | Frag * Position, data=SECCt, 
         pch=SECCt$pt, col=SECCt$colr, # , bg=Chamber.map$bg
         panel = panel.lines2
         )

  ## Plotting: Observed and Fitted from GLMM - from Richard & Zofia's GLMM workshop
  df <- coef( lmList(Y ~ X | Chamber * Position, data=SECCt) )
  cc1 <- as.data.frame(coef(Y.model)$Y)
  names(cc1) <- c("A", "B")
  df <- cbind(df, cc1)
  ff <- fixef(Y.model)

  print( xyplot( Y ~ X | Chamber * Position, data = SECCt, 
                aspect = "xy", layout = c(4,3),
                type = c("g", "p", "r"), coef.list = df[,3:4],
                panel = function(..., coef.list) {
                  panel.xyplot(...)
                  panel.abline(as.numeric( coef.list[packet.number(),] ), 
                               col.line = trellis.par.get("superpose.line")$col[2],
                               lty = trellis.par.get("superpose.line")$lty[2]
                               )
                  panel.abline(fixef(Y.model), 
                               col.line = trellis.par.get("superpose.line")$col[4],
                               lty = trellis.par.get("superpose.line")$lty[4]
                               )
                },
                index.cond = function(x,y) coef(lm(y ~ x))[1],
                xlab = X.plotlab,
                ylab = Y.plotlab,
                key = list(space = "top", columns = 3,
                  text = list(c("Within-subject", "Mixed model", "Population")),
                  lines = list(col = trellis.par.get("superpose.line")$col[c(2:1,4)],
                  lty = trellis.par.get("superpose.line")$lty[c(2:1,4)]
                  ))
                )
  )
}


if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()



