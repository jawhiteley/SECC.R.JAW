################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 3:
### Balance between Productivity and Decomposition in moss layer
### Jonathan Whiteley		R v2.12		2012-07-03
################################################################
## INITIALIZE
if (FALSE) {  # do not run automatically
  ## Set Working Directory: path in quotes "".
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("./")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place

  source("./lib/load.R")  # (re-)load data
}   

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source("./lib/init.R")  # Initialize - all analysis scripts should start with this.


library(ggplot2)

##==============================================================
## Convert Production and Decomposition to common units (t4 only)
##==============================================================
## Decomposition is recorded as % mass loss,
## BUT, it seems unlikely and unrealistic that 100% of of a patch's dry weight is "available" for decomposition
## I am assuming that a portion of the dry weight of each patch is 'decomposing'
## But, what is a reasonable proportion?  
## 1/2?  At that rate, totall mass loss by decomposition may appear to exceed productivity, even in Ambient patches?
##       (actually, this may have been the result of an error in my math)
## 1/3?  Seems reasonable?
## Analysis of the Moss biomass data (`Moss_biomass.R`) suggests a proportion of 0.115 - 0.118 (Brown + Dead tissue / Total dwt).
## - This might be an underestimate, given that there's likely a bunch of moss litter not attached to the shoot.
## By length: 0.15 - 0.16, but that's based on total segments for each shoot.
## - this may also be biased, because of the variation in total shoot length.
## A more conservative approach might just be to go by depth:
## - 75% of "green" segments are above 6 cm
## - most "brown" or "dead" segments are at 6cm, or below.
## - assuming a total depth of 10 cm, 1/3 is actually a reasonable proportion, assuming mass doesn't change with depth
## - realistically, biomass probably does change with depth, with more litter deeper down, and less total biomass at the top.
## Assuming 1/3 seems like a reasonable balance between the various estimates based on Moss biomass data

## Conversion factors - defined in SECC.functions

SECC <- within(SECC, 
               {
                 Productivity <- Prod12 + Prod23 # Cumulative productivity over 2nd year (same time period as litter bags for Decomposition)
                 if (FALSE)
                 {
                   Productivity <- Productivity * 385.8 # shoots / patch (approximately, on average)
                 } else {
                   ## multiply by fraction of patch weight in ~1 shoot (half of sample used for cyanobacteria Cells)
                   Productivity <- Productivity * Patch.dwt / (Cells.dwt / 1000 * 0.5) 
                 }
                 Productivity <- Productivity / 1000 # mg -> g (yes, I still need to do this)
                 ## Analysis of the Moss biomass data (`Moss_biomass.R`) suggests a proportion of 0.115 - 0.118.
                 Decompositng <- Decomposition * Patch.dwt * (1/3) # only dead tissue is really decomposing?
                 ## Convert each g / patch / yr -> g / m2 / yr
                 Productivity <- Productivity / patch.m2
                 Decompositng <- Decompositng / patch.m2
                 ## remove outliers; 44C-1.I, 64C-1.I, if extrapolation based on Patch.dwt
                 Productivity[which(Productivity > 500)] <- NA
                 # Net Production
                 PD.bal <- Productivity - Decompositng
               }
)

attr(SECC, "labels")[["PD.bal"]] <- "Net Moss Production" # "Net Biomass Production - Decomposition" ?
attr(SECC, "units" )[["PD.bal"]] <- quote(g %.% m^-2 %.% yr^-1)



##==============================================================
## Data Exploration & Checking
##==============================================================
if (FALSE)
{
  with(SECC, hist(Productivity) )
  with(SECC, hist(Decompositng) )
  with(SECC, hist(PD.bal) )

  SECCsub <- subset(SECC, Position %in% c("Inner", "Outer") ) 
  SECCsub <- subset(SECCsub, Chamber %in% c("Ambient", "Full Chamber") ) 
  SECCsub$Position <- factor(SECCsub$Position)
  SECCsub$Chamber  <- factor(SECCsub$Chamber)

  library(ggplot2)
  PD.plot <- ggplot(SECCsub, 
                    aes(x = Chamber, y = PD.bal, colour = Position)
  ) +
               stat_summary(fun.data = "mean_cl_boot", geom = "errorbar" ) +
               stat_summary(fun.y = mean, geom = "line" )

               print(PD.plot)

  library(mgcv)
  PxD.plot <- ggplot(SECCsub,
                    aes(x = Decompositng, y = Productivity, colour = Chamber, shape = Position)
  ) +
  geom_point() + stat_smooth(method = "gam") +
  geom_abline(intercept = 0, slope = 1, colour = "#999999") + coord_equal()
  print(PxD.plot)

  ## Possible Productivity Outliers
  Psort <- order(SECC$Productivity, decreasing = TRUE)
  SECC[Psort[1:5], ]

  ## quick regression of Net Biomass Production vs. Patch.dwt?
  Patchreg.plot <- ggplot(SECCsub, aes(x = Patch.dwt, y = PD.bal, colour = Chamber, shape = Position) ) +
                    geom_point() + stat_smooth(method = "gam")
  print(Patchreg.plot)

}





################################################################
## RUN STANDARD NESTED ANALYSIS (ANOVA)
################################################################
##==============================================================
## CONFIGURE BASIC ANALYSIS
##==============================================================

### Response Variable *****
Y.col <- 'PD.bal'     # Column to analyze as response variable           *****
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****

### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[3]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include (all available, but we'll exclude partial for simplicity?)
Position.use <- levels(SECC$Position)[c(1, 3)]  # Patch Positions to include: Inner, Outer

## Define Labels
Y.units <- bquote( .(Y.units) )     # store as quote(expression)  *****

## Output Results?
Save.results  <- TRUE


### Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

##================================================
## CUSTOM LABELS
##================================================



##################################################
### RUN STANDARD nested ANOVA
##################################################
source("./SECCanova/SECC - nested lm.R", echo = FALSE)

##   print( effect("Chamber:Position", Yp.fit) )
print( effect("Chamber:Position", Yp.fit) ) # works fine here


## ALERT: Tukey HSD family-wise CI are quite a bit wider than CIs in basic effects graphs!
##        Treatment differences may not be as strong as suggested by CIs


##================================================
## check for Block effects: I can't include it in the full model, because of lack of replication
Yp.b   <- gls(Y.trans ~ Block*Chamber*Position, data = SECCp,
              method="REML", na.action = na.omit
)
anova(Yp.b)
Yp.b0 <- update(Yp.b , method = "ML")
Yp.b1 <- update(Yp.b0, . ~ . -Block:Chamber:Position)
anova(Yp.b0, Yp.b1)                # Maybe a main effect, but no sig. interaction.  *pfew*

Block.plot <- ggplot(SECCp, aes(y = Y.trans, x = Chamber, group = Position)) + 
                geom_point(aes(colour = Position, shape = Position)) +
                stat_summary(fun.y = mean, geom = "line", aes(colour = Position, size = Position)) +
                facet_wrap(~ Block)
print(Block.plot)


##================================================
## MODEL SELECTION?
##================================================
## control for heterogeneity between blocks & years (see data exploration graphs)
Yp.hB  <- update(Yp.fit, weights = varIdent(form= ~1 | Block))
Yp.hF  <- update(Yp.fit, weights = varIdent(form= ~1 | Frag))
Yp.hBF <- update(Yp.fit, weights = varIdent(form= ~1 | Block * Frag))

Yp.gls <- gls(Yp.fixed, data = SECCp, method="REML", na.action = na.omit)

## Check Random Structure
anova(Yp.hBF, Yp.hB)
anova(Yp.hBF, Yp.hF)
anova(Yp.hBF, Yp.hF, Yp.fit, Yp.gls)              # random structure might make a difference ...

## controlling for heterogeneity between fragmentation treatments improves the model fit by AIC *slightly*, but not BIC.
## However, doing so means that Chamber:Frag:Position is no longer significant ... :/
## It also prevents effect() from working, making plotting WAY more difficult! :(
## Yp.fit  <- Yp.hF
anova(Yp.hF)

if (FALSE)
{   ## Backwards-selection (using likelihood-ratio tests)
  ML <- update(Yp.fit, method = "ML")
  ML1 <- update(ML, . ~ . - Chamber:Frag:Position)
  anova(ML, ML1)                       # * p = 0.0497 (without allowing heterogeneity across Frag trts) !!
  ML2 <- update(ML1, . ~ . - Chamber:Frag)
  ML3 <- update(ML1, . ~ . - Chamber:Position)
  ML4 <- update(ML1, . ~ . - Frag:Position)
  anova(ML1, ML2)
  anova(ML1, ML3)                      # ***
  anova(ML1, ML4)
  ML5 <- update(ML2, . ~ . - Frag:Position)
  anova(ML2, ML5)
  ML6 <- update(ML5, . ~ . - Frag)
  anova(ML5, ML6)
  ## Basically, the anova output produces largely the same results (slightly different p-values, but same terms are significant, or not).
}



###===============================================
### Prepare results for graphing
###===============================================
## effects() processing functions in SECC.functions.R

CxP.data <- effect.to.df(effect("Chamber:Position", Yp.fit)) # CxP.eff)
CxP.data$Chamber <- factor(CxP.data$Chamber, labels = c("Ambient", "Chamber"))
CFP.data <- effect.to.df(effect("Chamber:Frag:Position", Yp.fit))
CFP.data$Chamber <- factor(CFP.data$Chamber, labels = c("Ambient", "Chamber"))


##################################################
### PUBLICATION GRAPHS
##################################################
Y.lim <- range(c(CxP.data$lower, CxP.data$upper))
Plot.Title <- bquote("Patch means " %+-% "95% Confidence Intervals")
Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = c("Inner", "other", "Outer") )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
PositionPts  <- ggPts.SECC(Position.map, Position.label) 


CxP.plot <- ggplot(data = CxP.data, 
                   aes(x = Chamber, y = effect, 
                       group = Position, colour = Position, fill = Position, shape = Position)
                   ) + ylim(Y.lim) +
                    labs(x = attr(SECC, "labels")[["Chamber"]], y = Y.plotlab) +
                    opts(title = Plot.Title) +
                    geom_hline(yintercept = 0, size = 0.5, colour = "#999999") +
                    geom_line(aes(group = Position, size = Position)) +
                    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
                    geom_point(aes(group = Position), size = 3) 

CxP.plot <- CxP.plot + PositionPts +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label) +
scale_size_manual(name = Position.label,
                  values = Position.map$lwd, 
                  breaks = Position.map$label)
CxP.plot <- CxP.plot + jaw.ggplot()
print(CxP.plot)


## Chamber-Frag-Position interaction?
FragIconList <- list(FragIcon1 = FragIcon1,
                     FragIcon2 = FragIcon2,
                     FragIcon3 = FragIcon3,
                     FragIcon4 = FragIcon4
                     )
Y.lim <- range(c(CFP.data$lower, CFP.data$upper))
Plot.Title <- bquote("Patch means " %+-% "95% Confidence Intervals")
Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = c("Inner", "other", "Outer") )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
PositionPts  <- ggPts.SECC(Position.map, Position.label) 


CFP.plot <- ggplot(data = CFP.data, 
                   aes(x = Frag, y = effect, 
                       group = Position, colour = Position, fill = Position, shape = Position)
                   ) + ylim(Y.lim) +
                    labs(x = attr(SECC, "labels")[["Chamber"]], y = Y.plotlab) +
                    opts(title = Plot.Title) +
                    geom_hline(yintercept = 0, size = 0.5, colour = "#999999") +
                    geom_line(aes(group = Position, size = Position)) +
                    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
                    geom_point(aes(group = Position), size = 3) + facet_grid(.~ Chamber)

CFP.plot <- CFP.plot + PositionPts +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label) +
scale_size_manual(name = Position.label,
                  values = Position.map$lwd, 
                  breaks = Position.map$label) +
scale_x_discrete(labels = names(FragIconList),
                 breaks = levels(CFP.data$Frag))
CFP.plot <- CFP.plot + jaw.ggplot() + 
opts(axis.ticks.margin = unit(0.2, "lines"),
     axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
)

print(CFP.plot)


## stacked bar graph of Total Production - Decomposition?
## I currently don't have a graphical representation of the relative magnitudes of these two processes: only the combined net effect.
## ...


if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  Save.final <- paste(SaveDir.plots(), "Figure - ", "PD-balance", sep="")
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 3, height = 3, scale = 1.5)
  ggsave(file = paste(Save.final, "- CFP.eps"), plot = CFP.plot, width = 3, height = 3, scale = 1.5)
}

