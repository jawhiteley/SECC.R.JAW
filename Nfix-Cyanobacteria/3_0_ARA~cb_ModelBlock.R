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
  getwd()  # Check that we're in the right place

  ## Load data, functions, etc.  Process data & setup config. values.  
  ## Includes rm(list=ls()) to clear memory
  source('./Nfix-Cyanobacteria/1_ARA-cb_setup.R')
}


## library(lattice)    # ggplot2 with faceting is easier!
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme

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
## Start simple and add complexity (i.e. forward model-selection


### Fixed effects: All?  No Replication for all combinations!
Y.fixed <- Y.trans ~ X.trans * H2O * I(H2O^2) * Block * Time * Chamber * Frag * Position
Y.fixCl <- Y.trans ~ X.trans * H2O * I(H2O^2) * Block * Time * Climate * Frag
## Main effects only
Y.main  <- Y.trans ~ X.trans + H2O + I(H2O^2) + Block + Time + Chamber + Frag + Position
Y.mainCl<- Y.trans ~ X.trans + H2O + I(H2O^2) + Block + Time + Climate + Frag




##==============================================================
## Model Fitting
##==============================================================
## fixed effects only for assessment (How much variation is explained?)
Y.lmain <- lm( Y.main , data=SECCa)
summary(Y.lmain)
Y.fmain <- lm( Y.fixed , data=SECCa)   # no replication of all interaction combinations
anova(Y.fmain)
## Reasonable model of effects & Interactions for assessment (thanks John Connelly)
Y.lm <- lm( Y.trans ~ X.trans * H2O + Block + Time * Chamber * Frag * Position , data=SECCa)
summary(Y.lm)
## Transformation?
Y.rawlm <- lm( Y ~ X * H2O + Block + Time * Chamber * Frag * Position , data=SECCa)
summary(Y.rawlm)
Y.translm <- lm( Y.trans ~ X * H2O + Block + Time * Chamber * Frag * Position , data=SECCa)
summary(Y.translm)
Y.lmtrans <- lm( Y ~ X.trans * H2O + Block + Time * Chamber * Frag * Position , data=SECCa)
Y.lmsqrt <- lm( Y.log ~ I(sqrt(X)) * H2O + I(H2O^2) + Block + Time * Chamber * Frag * Position , data=SECCa)

AIC(Y.lmain, Y.lm, Y.rawlm, Y.lmtrans, Y.translm, Y.lmsqrt)
summary(Y.lmain)                       # The best model, according to AIC**
## AIC suggests a log-transformation on only the response (ARA/N-fixation) is a good fit.
## However, this implies an *exponential*, rather than saturating log- relationship
## between N-fixation and Cell abundance, which we predicted / hypothesised.
## Granted, the AIC improvement is about 1.7% (905 -> 888), 
## and the R^2 is only marginally higher (0.609 vs. 0.591)
## Both are substantially better than the model with untransformed (response) variables.
## It is probably more defensible to use the log-log model, which makes more sense theoretically.
## There is so much noise in the data that the marginally better log-Y fit may be spurious or inconsequential.
## The model with more interactions seems to be a *poorer* fit than main effects only, 
##   according to AIC?  But the R^2 is a little higher (0.59 vs. 0.57)

Y.model  <- Y.lmain



##==============================================================
## MODEL SELECTION
##==============================================================
## multi-model averaging and inference with glmulti!
library(glmulti)                       #  v1.0, april 2011
library(MASS)
## library(leaps)
## Note: glmulti performs an exhaustive search of all candidate models.  
##       2^7 = 128 candidate models (not including interactions). 256 according to glmulti
## method="l" (leaps) uses a *much* faster algorithm, but not with factors or interactions :(
ARA.glmulti1 <- glmulti(Y.fixed, data=SECCa, crit=aic, level=1, fitfunc=lm, 
                        confsetsize=256, plotty=FALSE, report=FALSE)
print(ARA.glmulti1)
ARA.fit <- ARA.glmulti1@objects

## level=2 for 2-way interactions
## 2^(7 + choose(7,2) ) = 268,435,456 candidate models with 2-way interactions!!
## 416,869,200 according to glmulti
## An exhaustive exploration would take ~ 1 month on my computer.  
## Try with 'genetic' algorithm to speed things up (method="g").
## method="d" to print a summary of candidate models (no fitting)
ARA.glmulti2 <- glmulti(Y.fixed, data=SECCa, crit=aic, level=2, fitfunc=lm, 
                        method="g", confsetsize=1024, plotty=FALSE, report=TRUE)
print(ARA.glmulti2)


getCoef.glmulti <- function(glmObj, minImportance=0) {
  glm.coef <- as.data.frame(coef(glmObj))
  ## ggplot will order the bars by levels of the explanatory factor
  glm.coef$Variable <- factor(row.names(glm.coef), levels=unique(row.names(glm.coef))) 
  glm.order <- order(glm.coef$Importance, decreasing=TRUE)
  ## glm.coef$Variable <- factor(glm.coef$Variable, levels=levels(glm.coef$Variable)[glm.order])
  ## Drop terms below the threshold
  glm.minImp <- which(glm.coef$Importance >= minImportance)
  glm.coef <- glm.coef[glm.minImp, ]
  ## facilitate confidence intervals on Estimates
  glm.coef$Emax <- glm.coef$Estimate + glm.coef[, "+/- (alpha=0.05)"]
  glm.coef$Emin <- glm.coef$Estimate - glm.coef[, "+/- (alpha=0.05)"]
  ## Clean up labels
  levels(glm.coef$Variable) <- gsub("X.trans", "Cells", levels(glm.coef$Variable))
  levels(glm.coef$Variable) <- gsub("Chamber.*", "Chamber", levels(glm.coef$Variable))
  levels(glm.coef$Variable) <- gsub("Position.*", "Position", levels(glm.coef$Variable))
  levels(glm.coef$Variable) <- gsub("Frag(.*)", "Frag(\\1)", levels(glm.coef$Variable))
  levels(glm.coef$Variable) <- gsub("Time(.*)", "Time(\\1)", levels(glm.coef$Variable))
  glm.coef
}

ARA.coef1 <- getCoef.glmulti(ARA.glmulti1)
ARA.coef2 <- getCoef.glmulti(ARA.glmulti2)
nrow(ARA.coef2)                        # 129 terms (128 + intercept)!
## defaults
ARA.glmulti <- ARA.glmulti1
ARA.coef <- ARA.coef1

## Output graphs
par(mfrow=c(1,1))
plot(ARA.glmulti1, type="w")
barplot(ARA.coef1[, "Importance"], horiz=TRUE, names.arg=ARA.coef1$Variable, las=2) 

## ggplot2 theme settings
bar.import <- function(glmObj) {
  ggplot(glmObj, aes(x=Variable, y=Importance), 
                          stat="identity", xlab = "Explanatory Variables") +
         list(geom_bar(), coord_flip(), scale_y_continuous(expand=c(0,0)),
              opts(panel.border=theme_blank(), axis.line=theme_segment(),
                   plot.title=theme_text(size = 16, face = "bold"),
                   title="Model-averaged importance of effects")
         )
}
est.confint <- function(glmObj) {
  conf.wd <- glmObj[, "+/- (alpha=0.05)"]
  conf.int <- aes(ymax = Estimate + conf.wd, ymin = Estimate - conf.wd)
  ggplot(glmObj, aes(x=Variable, y=Estimate) ) +
  list(geom_point(), 
       geom_errorbar(aes(ymax = Emax, ymin = Emin), width=0.1),
       coord_flip())
}

ARA.importance1 <- bar.import(ARA.coef1)
ARA.est1 <-est.confint(ARA.coef1)
print(ARA.importance1)
print(ARA.est1)

ARA.importance2 <- bar.import(ARA.coef2[ARA.coef2$Importance>=1, 
                              c("Variable", "Importance")]) + 
                              opts(axis.text.y=theme_text(size=5))
ARA.est2 <-est.confint(ARA.coef2)
print(ARA.importance2)
print(ARA.est2)



##==============================================================
## GAM: Is the relationship linear or not?
##==============================================================
## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)
library(mgcv)
ARA.gam <- gam(Y.log ~ s(X.log) + s(H2O) + Block + Time*Chamber*Frag*Position, data=SECCa)
anova(ARA.gam)
ARA.gam1 <- gam(Y.log ~ s(X) + s(H2O) + Block + Time*Chamber*Frag*Position, data=SECCa)
anova(ARA.gam1)
## It *might* be linear for log(Cells) (X.log; edf = 2.18), but definitely NOT for H2O (edf=4)
## polynomial for moisture (H2O)?
ARA.gam2 <- gam(Y.log ~ s(X.log) + s(poly(H2O, 2)) + Block + Time*Chamber*Frag*Position, data=SECCa)
anova(ARA.gam2)
## what about sqrt-transforming Cell Density (on a whim)
ARA.Xsqrt.gam <- gam(Y ~ s(sqrt(X)) + s(H2O) + Block + Time*Chamber*Frag*Position, data=SECCa)
anova(ARA.Xsqrt.gam)
ARA.logsqrt.gam <- gam(Y.log ~ s(sqrt(X)) + s(H2O) + Block + Time*Chamber*Frag*Position, data=SECCa)
anova(ARA.logsqrt.gam)
## Definitely linear for sqrt(Cell Density)!!




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


RE <- diagnostics(Y.model, resType="pearson", more=TRUE) # full diagnostics, just to be sure

## Check for spatial patterns in residuals?

## Compare model to GA(M)M to check that a linear fit is the most appropriate?
## see Zuur et al. (2007, 2009: Appendix)



################################################################
## ANALYSIS: GET RESULTS
################################################################
anova(Y.model)
summary(Y.model)

## intervals(Fitted.Model.Object)   # Approx. Confidence intervals.  see ?intervals
## see multcomp package for multiple comparisons (Tukey's HSD on mixed effects model?)





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
