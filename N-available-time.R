##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Template for basic analyses of experimental data
### Available Nitrogen: NH4, NO3, TAN  @ time 4
### Jonathan Whiteley     R v2.12     2012-07-18
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
Y.col <- 'logTAN'      # Column to analyze as response variable           *****
Y.use <- 'Y'        # Which transformation is being used (for labels)? *****

##================================================
## CUSTOM CALCULATIONS 
##================================================
SECC.prime <- SECC    # save a copy of the original for reference.

## Conversion factors defined in SECC.functions

SECC <- within( SECC, { 
			   ## change negative ARA values to 0 - should I wait until after aggregation?
			   ARA.ml[ARA.ml < 0] <- 0
			   ARA.m[ ARA.m  < 0] <- 0
               ARA.g[ ARA.g  < 0] <- 0
			   Nfix <- ARA.m * Nfix.ARA.ratio
               Nyrs <- Ndays / 365     # 365 days / year
               NH4 <- NH4 * Nyrs
               NO3 <- NO3 * Nyrs
               TAN <- TAN * Nyrs
               logTAN <- log10(TAN)    # No log(x +1) or any BS
})
attr(SECC, "labels")[["logTAN"]] <- attr(SECC, "labels")[["TAN"]]
attr(SECC, "units")[["logTAN"]]  <- quote(g %.% m^-2) # I'm going to back-transform it for plotting :P

## Define Labels
## Y.units <- attr(SECC, "units")[[Y.col]]
## Y.units <- bquote( log(.(Y.units)) )  # sqrt(.(Y.units), 4)  # store as quote(expression) *****

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

## Output Results to File?
Save.results  <- TRUE  

## File names are determined automatically (in labels script).  Specify custom filenames *after* default labels are calculated below.


### Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 


##================================================
## CUSTOM LABELS
##================================================
Save.filename <- gsub("TAN", "TAN-time", Save.filename, fixed=TRUE)
Save.text     <- gsub("TAN", "TAN-time", Save.text,     fixed=TRUE)
Save.plots    <- gsub("TAN", "TAN-time", Save.plots,    fixed=TRUE)
Save.final    <- gsub("TAN", "TAN-time", Save.final,    fixed=TRUE)



##################################################
### RUN STANDARD nested ANOVA
##################################################
## No data from t1-2
cat("\n\n\nProcessing Time:", Time.use, "\n")
cat("* without controlling for duration of sampling\n")

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




##################################################
### Further Exploration
##################################################
if (F)
{
  ## are TAN or logTAN colinear with moisture? (given apparent inverse of CHabmer x Position interaction, albeit NS)
  pairplot(SECCp[, c("logTAN", "TAN", "H2O")]) # NO!  Thank God.  Still interesting, though.
}



##################################################
### PUBLICATION GRAPHS
##################################################

Y.lim <- c(0.0, 0.04)
Plot.Title <- bquote("Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
## Position.map$label <- c("395d Inner", " 58d other", "395d Outer")
Position.map$label <- c("Inner - 1 year", "other  - summer", "Outer - 1 year")
Position.map$label <- c("Inner (395d)", "other  \n(58d summer)", "Outer (395d)")
Position.map <- Position.map[c(3, 1, 2), ]


## data frame of plot values (for ggplot2).
plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Chamber=SECCp$Chamber, 
                                     Position=SECCp$Position, 
                                     Time=SECCp$Time), 
                                error = msd["Chamber:Position"]
                                )
## back-transform
plot.means[, c("x", "lower", "upper")] <- sapply(plot.means[, c("x", "lower", "upper")], function (x) 10^x)
levels(plot.means$Chamber)[2] <- "Chamber"
plot.means$Position <- factor(plot.means$Position, 
                              levels = levels(plot.means$Position)[c(3, 1, 2)], 
                              labels = Position.map$label
                              )

CxP.plot <- qplot(Chamber, x, data = plot.means, group = Position, 
                    geom = "point", size = I(3), ylim = Y.lim, 
                    colour = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE)
## CxP.plot <- CxP.plot + geom_point(aes(Chamber, x), size = 2)
CxP.plot <- CxP.plot + geom_line(aes(group = Position, lwd = Position) ) +
                        geom_errorbar(aes(ymin = lower, ymax = upper), 
                                      width = 0.2, size = 0.5)
CxP.plot <- CxP.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label) +
                        scale_fill_manual(name = Position.label,
                                          values = Position.map$bg, 
                                          breaks = Position.map$label) +
                        scale_shape_manual(name = Position.label,
                                           values = Position.map$pch, 
                                           breaks = Position.map$label) +
                        scale_size_manual(name = Position.label,
                                          values = Position.map$lwd/2, 
                                          breaks = Position.map$label)
CxP.plot <- CxP.plot + jaw.ggplot()
print(CxP.plot)


## Position main effect (this is what is signficant, according to ANOVA)
plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Position=SECCp$Position, 
                                     Time=SECCp$Time), 
                                error = msd["Position"]
                                )
## back-transform
plot.means[, c("x", "lower", "upper")] <- sapply(plot.means[, c("x", "lower", "upper")], function (x) 10^x)
plot.means$Position <- factor(plot.means$Position, 
                              levels = levels(plot.means$Position)[c(3, 1, 2)], 
                              labels = Position.map$label
                              )

P.plot <- qplot(x = Position, y = x, data = plot.means, 
                    geom = "bar", stat="identity", ylim = Y.lim, 
                    fill = I("#999999"),
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Pos"]],
                    ylab = Y.plotlab,
                    legend = FALSE)
## P.plot <- P.plot + geom_point(aes(Chamber, x), size = 2)
P.plot <- P.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                         width = 0.2, size = 0.5)
P.plot <- P.plot + jaw.ggplot()
print(P.plot)


if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 4, height = 4, scale = 1.5)
  ggsave(file = paste(Save.final, "- P.eps"),   plot =   P.plot, width = 4, height = 4, scale = 1.5)
}
