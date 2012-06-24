################################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Main control Script - Chapter 3:
### Balance between Productivity and Decomposition in moss layer
### Jonathan Whiteley		R v2.12		2012-06-24
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
} else {
  ## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
  source("./lib/init.R")  # Initialize - all analysis scripts should start with this.
}


library(ggplot2)

##==============================================================
## Convert Production and Decomposition to common units (t4 only)
##==============================================================
## Decomposition is recorded as % mass loss,
## BUT, it seems unlikely and unrealistic that 100% of the moss is "eligible" for decomposition
## I am assuming that a portion of the dry weight of each patch is 'decomposing'
## But, what is a reasonable proportion?  
## 1/2?  At that rate, totall mass loss by decomposition may appear to exceed productivity, even in Ambient patches?
##       (actually, this may have been the result of an error in my math)
## 1/3?  Seems reasonable?
## I should really use Heather's data (Moss_biomass) to estimate what proportion, or at what depth moss tissue is "brown", and use that.

## Conversion factors - defined in SECC.functions

SECC <- within(SECC, 
               {
                 Productivity <- Prod12 + Prod23
                 if (FALSE)
                 {
                   Productivity <- Productivity * 385.8 # shoots / patch (appoximately, on average)
                 } else {
                   Productivity <- Productivity * Patch.dwt / (Cells.dwt / 1000 * 0.5) # multiply by fraction of patch weight in ~1 shoot (half of sample used for cyanobacteria Cells)
                 }
                 Productivity <- Productivity / 1000 # mg -> g (yes, I still need to do this)
                 Decompositng <- Decomposition * Patch.dwt/2 # only dead tissue is really decomposing?
                 PD.bal <- Productivity - Decompositng  # g / patch over 1 year
                 PD.bal <- PD.bal / (patch.m2 * 1)            # patch area -> square m ; 1000 g / kg X
               }
)

attr(SECC, "labels")[["PD.bal"]] <- "Net Biomass Production (- Decomposition)"
attr(SECC, "units" )[["PD.bal"]] <- quote(g %.% m^-2 %.% y^-1)



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
}





################################################################
## PRODUCTION - DECOMPOSITION MASS BALANCE (ANOVA)
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

###===============================================
### Prepare results for graphing
###===============================================
## effects() processing functions in SECC.functions.R

CxP.data <- effect.to.df(CxP.eff)
CxP.data$Chamber <- factor(CxP.data$Chamber, labels = c("Ambient", "Chamber"))

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


if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  Save.final <- paste(SaveDir.plots(), "Figure - ", "PD-balance", sep="")
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 3, height = 3, scale = 1.5)
}

