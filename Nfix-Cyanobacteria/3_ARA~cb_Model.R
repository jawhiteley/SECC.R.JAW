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
  source('1_ARA-cb_setup.R')
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

## Subsume Chamber & Position into "Climate" pseudo-treatment?
## Time as a fixed factor, or separtae analysis on each Time?

##==============================================================
## Model Formula
##==============================================================
UseClimateFac <- TRUE

### Fixed effects
### Include H2O^2 to test for unimodal relationship?
## Including Time as a factor?
## It might be better to analyse each time period separately, and avoid the higher-order interactions that I know exist.
if ( length(Time.use) > 1 ) {
  Y.fixed <- Y.log ~ X.log * H2O * Time * Chamber * Frag * Position
  Y.fixCl <- Y.log ~ X.log * H2O * Time * Climate * Frag
  ## Main effects only
  Y.main  <- Y.log ~ X.log + H2O + Time + Chamber + Frag + Position
  Y.mainCl<- Y.log ~ X.log + H2O + Time + Climate + Frag
} else {
  Y.fixed <- Y.log ~ X.log * H2O * Chamber * Frag * Position
  Y.fixCl <- Y.log ~ X.log * H2O * Climate * Frag
  ## Main effects only
  Y.main  <- Y.log ~ X.log + H2O + Chamber + Frag + Position
  Y.mainCl<- Y.log ~ X.log + H2O + Climate + Frag
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
  Y.RisN <- ~ 1 + X.log * H2O | Block/Time/Chamber/Frag # Random intercept + slopes
} else {
  Y.RiN  <- ~ 1|Block/Time/Climate
  Y.RisN <- ~ 1 + X.log * H2O | Block/Time/Climate # Random intercept + slopes
}

## This might be the wrong approach: Block is the only real "random" factor.
## Fixed factors happen to be nested within it.
## Maybe what I should really be doing is accounting for heterogeneity across nested fixed factors?
Y.Ri  <- ~ 1 | Block
Y.Ris <- ~ 1 + X.log * H2O | Block     # Random intercept + slopes




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
Y.model.main <- glm(Y.log ~ X.log + Time + Chamber + Frag + Position + H2O + I(H2O^2), data = SECCa)  # Main factors only; no interaction terms.

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
## GLMM - Hierarchical / Multilevel Mixed Model               **
################################################################
##==============================================================
## Model Fitting
##==============================================================
## Using a mixed-effects model to account for treatments of different sizes? 
##  Block / Chamber / Frag / Position
Y.fm <- gls( Y.fixed, data=SECCa, method="REML") # fixed effects only for comparison
lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000) # takes a while to converge...
Y.rim  <- lme(Y.fixed, random=Y.Ri,  data=SECCa, control=lmd, method="REML")
## SLOW! :(
## more complex models take a few minutes (or hours) to fit.  Use with caution.
Y.rism <- lme(Y.fixed, random=Y.Ris, data=SECCa, control=lmc, method="REML")
Y.rie  <- lme(Y.fixed, random=Y.Ri,  
              weights=varIdent(form=~ 1 | Block),
              data=SECCa, control=lmc, method="REML")
Y.rise <- lme(Y.fixed, random=Y.Ris,  
              weights=varIdent(form=~ 1 | Block),
              data=SECCa, control=lmc, method="REML")
if (false) {                           # not working :(
  Y.rieN  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block/Time/Chamber/Frag),
                 data=SECCa, control=lmc, method="REML")
  ## I want to use varConstPower (variance covariate with values of 0), but that produces an error :(
  ## varFixed gives no errors, but takes >12 hours to fit (if at all)
  Y.rieXP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varFixed(~ X.log),
                 data=SECCa, control=lmd, method="REML")
  Y.rieHP <- lme(Y.fixed, random=Y.Ri,  
                 weights=varPower(form=~ H2O),
                 data=SECCa, control=lmd, method="REML")
  Y.riceX <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ X.log)),
                 data=SECCa, control=lmc, method="REML")
  Y.riceH <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ H2O)),
                 data=SECCa, control=lmc, method="REML")
  Y.rice  <- lme(Y.fixed, random=Y.Ri,  
                 weights=varComb(varIdent(form=~ 1 | Block),
                                 varConstPower(form=~ X.log),
                                 varConstPower(form=~ H2O)),
                 data=SECCa, control=lmc, method="REML")
}

## Main effects only: ignore interactions?
Y.mainML <- lme(Y.main, data = SECCa, random = Y.Ri, method="ML")




##==============================================================
## MODEL SELECTION
##==============================================================
## RANDOM structure
anova(Y.fm, Y.rism)                    # do random effects improve the model?
anova(Y.fm, Y.rim)                     # do random effects improve the model?
anova(Y.fm, Y.rie)                     # do random effects improve the model?
anova(Y.rie, Y.rim, Y.rism, Y.rise)    # do we need random slopes or error terms?

Y.mm <- Y.rim                          # Optimal random structure
## The biggest improvement seems to come from random intercepts across Blocks.
## However, that model shows major heterogeneity in the residuals.
## I need other random factors to maintain a valid model fit.

## optimize FIXED factors
if (FALSE) {
  Y.ml  <- lme(Y.fixed, data=SECCa, random=Y.Ri, method="ML") # re-fit with ML
  Y.rieML <- lme(Y.fixed, random=Y.Ri,  
                 weights=varIdent(form=~ 1 | Block),
                 data=SECCa, control=lmc, method="ML")
}
Y.ml  <- update(Y.mm, method="ML")     # re-fit with ML; some models can't be :(
drop1(Y.ml)                            # not encouraging
## Y.step <- step(Y.ml)                # stepwise back & forward model selection?  Not for lme
if (!UseClimateFac) {                  # All factors, or Climate pseudo-factor
  Y.ml1 <- update(Y.ml, ~. - X.log:H2O:Time:Chamber:Frag:Position)

  ## No interactions?
  anova(Y.ml, Y.mlMain)

} else {                               # Chamber & position lumped into Climate pseudo-factor
  Y.ml1 <- update(Y.ml, ~. - X.log:H2O:Time:Climate:Frag)
  anova(Y.ml, Y.ml1)
  ## drop 4-way interactions?
  Y.ml2 <- update(Y.ml1, ~. - X.log:H2O:Time:Climate) # *
  anova(Y.ml1, Y.ml2)
  Y.ml3 <- update(Y.ml1, ~. - X.log:H2O:Time:Frag)
  anova(Y.ml1, Y.ml3)
  Y.ml4 <- update(Y.ml1, ~. - X.log:H2O:Climate:Frag) # <-
  anova(Y.ml1, Y.ml4)
  Y.ml5 <- update(Y.ml1, ~. - X.log:Time:Climate:Frag) # **
  anova(Y.ml1, Y.ml5)
  Y.ml6 <- update(Y.ml1, ~. - H2O:Time:Climate:Frag)
  anova(Y.ml1, Y.ml6)
  ## dropped least significant 4-way interaction: next
  Y.ml2 <- update(Y.ml4, ~. - X.log:H2O:Time:Climate) # *
  anova(Y.ml4, Y.ml2)
  Y.ml3 <- update(Y.ml4, ~. - X.log:H2O:Time:Frag) # <-
  anova(Y.ml4, Y.ml3)
  Y.ml5 <- update(Y.ml4, ~. - X.log:Time:Climate:Frag) # **
  anova(Y.ml4, Y.ml5)
  Y.ml6 <- update(Y.ml4, ~. - H2O:Time:Climate:Frag)
  anova(Y.ml4, Y.ml6)
  ## dropped least significant 4-way interaction: next
  Y.ml2 <- update(Y.ml3, ~. - X.log:H2O:Time:Climate) # *
  anova(Y.ml3, Y.ml2)
  Y.ml5 <- update(Y.ml3, ~. - X.log:Time:Climate:Frag) # *
  anova(Y.ml3, Y.ml5)
  Y.ml6 <- update(Y.ml3, ~. - H2O:Time:Climate:Frag) # <-
  anova(Y.ml3, Y.ml6)
  ## dropped least significant 4-way interaction: next
  Y.ml2 <- update(Y.ml6, ~. - X.log:H2O:Time:Climate) # *
  anova(Y.ml6, Y.ml2)
  Y.ml5 <- update(Y.ml6, ~. - X.log:Time:Climate:Frag) # *
  anova(Y.ml6, Y.ml5)
  Y.ml4 <- Y.ml6                       # optimal model with 4-way interactions
  anova(Y.ml, Y.ml4)
  Y.mm <- update(Y.ml4, method="REML") # re-fit with REML

  ## 3-way interactions?

  ## Ignore interactions?
  Y.main1 <- update(Y.mainML, ~. - X.log)
  Y.main2 <- update(Y.mainML, ~. - H2O)
  Y.main3 <- update(Y.mainML, ~. - Time)
  Y.main4 <- update(Y.mainML, ~. - Climate)
  Y.main5 <- update(Y.mainML, ~. - Frag) # <-
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


## Standard diagnostic plots 
op <- par(mfrow=c(2,2))	 # panel of figures
plot(Y.mm)
par(op)

## Residuals ##
RE.lab <- "Standardized Residuals"
RE <- resid(Y.mm, type="p")            # "pearson" (standardised) residuals
if (FALSE) {                           # type="p" or type="normalized"?
  ## type = "normalized" residuals 
  ##  Zuur et al. 2009 use this for 'standardized' residuals, but it actually does something more complicated.  see ?residuals.lme
  ##  - in the examples in Zuur et al. 2009, there is no difference between normalized and standardised ("pearson"), but it might matter with more complicated correlation structures.
  ## type = "pearson" for 'standardised' residuals (see ?residuals.lme)
  ##  - this is what Zuur et al. 2009 are actually referring to, and is what Zuur et al. 2007 use in their R-code (type="p")
  ## see also rstudent (studentized) and rstandard (standardized) residuals.
  ## Zuur et al. use stdres() & studres() from the MASS library - what's the difference?
}

## Plot REsiduals: see Zuur et al. 2007, pg. 131-133
## REsiduals: Normal distribution?
op <- par(mfrow=c(2,2))
hist(RE)
hist(RE, freq=FALSE)
Xnorm <- seq(min(RE), max(RE), length=40)
lines(Xnorm, dnorm(Xnorm), col="grey70", lwd=2) # reference Normal distribution
qqnorm(RE)
qqPlot(RE)

## Homogeneity, Independence: Pattern in residuals indicates violation
Fit  <- fitted(Y.mm)
Fit0 <- fitted(Y.mm, level=0)
Fit1 <- fitted(Y.mm, level=1)
plot(Fit, RE, xlab="Fitted values", ylab=RE.lab)
plot(SECCa$X.log, RE, ylab=RE.lab)             
plot(SECCa$H2O,   RE, ylab=RE.lab)
## spreadLevelPlot(Y.mm)                  # library(car)
par(op)

qplot(Block, RE, data = SECCa, facets = . ~ Time, geom="boxplot" ) + jaw.ggplot()
qplot(Chamber, RE, data = SECCa, facets = Frag * Position ~ Time, geom="boxplot" ) + jaw.ggplot()

qplot(X.log, RE, data = SECCa, facets = Block ~ Time) + jaw.ggplot()
qplot(H2O, RE, data = SECCa, facets = Block ~ Time) + jaw.ggplot()


## Compare model to GA(M)M to check that a linear fit is the most appropriate
## see Zuur et al. (2007, 2009: Appendix)


## Global Validation of Linear Model Assumptions (gvlma package)
if ("lm" %in% class(Y.mm)) {
  validation <- gvlma(Y.mm)
  plot(validation)
  summary(validation)
}


################################################################
## ANALYSIS: GET RESULTS
################################################################
anova(Y.mm)
summary(Y.mm)

## intervals(Fitted.Model.Object)   # Approx. Confidence intervals.  see ?intervals


##==============================================================
## Variance Components Analysis (variance decomposition)?








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
Y.pred <- expand.grid(Chamber  = levels(SECCa$Chamber) , 
                      Frag     = levels(SECCa$Frag), 
                      Position = levels(SECCa$Position), 
                      X=seq(0, max(SECCa$X), length.out=100 ) 
                      )
Y.pred$predicted <- predict(Y.model, newdata=Y.pred, type="response" )  # newdata must have same explanatory variable name for predict to work.

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
               aggregate(cbind(predicted), list(Chamber = Chamber, X = X), mean)
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




if (FALSE) {
  ## Coplots with linear fits (from Zuur et al. 2007 Chapter 22 R code)
  ## individual lm's within each panel.  Not exactly what I want.
  coplot( Y ~ X | Frag * Position, data=SECCa, 
         pch=SECCa$pt, col=SECCa$colr, # , bg=Chamber.map$bg
         panel = panel.lines2
         )

  ## Plotting: Observed and Fitted from GLMM - from Richard & Zofia's GLMM workshop
  df <- coef( lmList(Y ~ X | Chamber * Position, data=SECCa) )
  cc1 <- as.data.frame(coef(Y.model)$Y)
  names(cc1) <- c("A", "B")
  df <- cbind(df, cc1)
  ff <- fixef(Y.model)

  print( xyplot( Y ~ X | Chamber * Position, data = SECCa, 
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



