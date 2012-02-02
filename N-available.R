##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Template for basic analyses of experimental data
### Available Nitrogen: NH4, NO3, TAN  @ time 4
### Jonathan Whiteley     R v2.12     2012-01-23
##################################################
## INITIALISE
##################################################
## This script is used in a generic way for most univariate analyses
## Working Directory: see lib/init.R below
if (FALSE) {  # do not run automatically
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("..")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')


##################################################
## CONFIGURE BASIC ANALYSIS
##################################################

### Response Variable *****
##  - NH4, NO3, TAN
Y.col <- 'TAN'      # Column to analyze as response variable           *****
Y.use <- 'Y.log'        # Which transformation is being used (for labels)? *****

### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[3]            # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)               # Frag treatments to include
Position.use <- levels(SECC$Position)           # Patch Positions to include
## ion resin capsules were placed in "other" patches several months after "Inner" & "Outer" 
##   patches, so they were deployed for a much shorter period
##   (see Ndays column).

## Define Labels
Y.units <- bquote( log(.(Y.units)) )  # sqrt(.(Y.units), 4)  # store as quote(expression) *****

## Output Results to File?
Save.results  <- FALSE  

## File names are determined automatically (in labels script).  Specify custom filenames *after* default labels are calculated below.


##================================================
## CUSTOM CALCULATIONS 
##================================================
SECC.prime <- SECC    # save a copy of the original for reference.

## str(SECC)
sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample
ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * (12.5^2)      # patch area
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2
Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.

SECC <- within( SECC, { 
			   ## change negative ARA values to 0 - should I wait until after aggregation?
			   ARA.ml[ARA.ml < 0] <- 0
			   ARA.m[ ARA.m  < 0] <- 0
               ARA.g[ ARA.g  < 0] <- 0
			   Nfix <- ARA.m * Nfix.ARA.ratio
               Nyrs <- Ndays / 365
})


### Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

##================================================
## CUSTOM LABELS
##================================================




##################################################
### RUN STANDARD nested ANOVA
##################################################
## No data from t1-2
cat("\n\n\nProcessing Time:", Time.use, "\n")

## RUN STANDARD nested ANOVA
source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)


###===============================================
### Include Time as a factor in nested ANOVA
###===============================================
## Note that Samples at different times are actually independent
## in this design, due to destructive sampling.
## For available N data, time is somewhat confounded with Position
## - capsules in "other" patches were deployed only for the summer before
##   being collected with the rest.
## - see Ndays column for actual number of days each capsule was deployed.
## - It does appear that most of the adsorption happens in the summer 
##   (other patches have much higher rates, over a shorter period of time)

if (FALSE)
{
  Time.use     <- levels(SECC$Time)      # Include *ALL* Times (as a Treatment)
  source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) # Load default Labels. *****
  source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE) # RUN STANDARD nested ANOVA
}


###===============================================
### Remove control for time deployed
###===============================================
## Should be no difference at all: therefore most of the adsorption happens in the summer, since the extra 10 months made no difference to the capsules 
## (or the method got seriously messed up somehow)

SECCt <- SECC                          # save a copy
SECC <- within( SECC, 
               { 
                 Nyrs <- Ndays / 365
                 NH4 <- NH4 * Nyrs
                 NO3 <- NO3 * Nyrs
                 TAN <- TAN * Nyrs
               }
)
attr(SECC, "units")[[Y.col]]  <- bquote(g %.% m^-2)

### Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

## RUN STANDARD nested ANOVA
cat("\n\n\nProcessing without Time\n")
source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)




##################################################
### Further Exploration
##################################################
if (F)
{
  ARA.N <- qplot(y = log1p(ARA.m), x = log1p(TAN), data = subset(SECCp, log1p(TAN) < 0.2) )
  print(ARA.N)
  mean(SECCp[SECCp$Position=="other", "Ndays"]) #  57.6 days =  58 days
  mean(SECCp[SECCp$Position!="other", "Ndays"]) # 395.4 days = 395 days
}



##################################################
### PUBLICATION GRAPHS
##################################################

Y.lim <- c(-0.01, 0.1)
Plot.Title <- bquote("Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
Position.map$label <- c("395d Inner", " 58d other", "395d Outer")
Position.map$label <- c("Inner - 1 year", "other  - summer", "Outer - 1 year")
Position.map <- Position.map[c(1, 3, 2), ]


## data frame of plot values (for ggplot2).
## might be able to accomplish much the same effect with stat_summary using means in ggplot2?
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Position=SECCp$Position, Time=SECCp$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means$error <- as.numeric(msd["Chamber:Position"]/2)
levels(plot.means$Chamber)[2] <- "Chamber"
plot.means$Position <- factor(plot.means$Position, 
                              levels = levels(plot.means$Position)[c(1, 3, 2)], 
                              labels = Position.map$label
                              )

CxP.plot <- qplot(Chamber, x, data = plot.means, group = Position, 
                    geom = "point", ylim = Y.lim, size = I(3), 
                    colour = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = .~Time)
## CxP.plot <- CxP.plot + geom_point(aes(Chamber, x), size = 2)
CxP.plot <- CxP.plot + geom_line(aes(group = Position), size = 0.8)
CxP.plot <- CxP.plot + geom_errorbar(aes(ymin = x - error, ymax = x + error), 
                                         width = 0.2, size = 0.5)
CxP.plot <- CxP.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_shape_manual(name = Position.label,
                                           values = Position.map$pch, 
                                           breaks = Position.map$label)
CxP.plot <- CxP.plot + jaw.ggplot()
print(CxP.plot)

if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 3, scale = 1.5)
}
