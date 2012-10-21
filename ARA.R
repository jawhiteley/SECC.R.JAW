##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Acetylene Reduction Assay (ARA: N-fixation)  @ time #s
### Jonathan Whiteley     R v2.12     2012-07-03
##################################################
## INITIALISE
##################################################
## This script is used in a generic way for most univariate analyses
## Working Directory: see lib/init.R below
if (FALSE) {  # do not run automatically
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (Mac).
  getwd()  # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')


##################################################
## CONFIGURE BASIC ANALYSIS
##################################################

### Response Variable *****
Y.col <- 'ARA.m'     # Column to analyze as response variable           *****
Y.use <- 'Y.sqrt'    # Which transformation is being used (for labels)? *****

### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[1]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include

## Define Labels
## Y.units <- bquote( sqrt(.(Y.units)) )     # store as quote(expression)  *****

## Output Results?
Save.results  <- TRUE  


##================================================
## CUSTOM CALCULATIONS 
##================================================
## !is.na(SECC$Time) ; NAs in factors are annoying
SECC.prime <- SECC    # save a copy of the original for reference.

if (FALSE)
{                                      # Conversion factors defined in SECC.functions
  ## str(SECC)
  sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
                  #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
                  #  or 6.4 for 20 shoots, based on density survey.
  sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
  sample_ml    <- 50  # 50 ml sample
  ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
  patchA   <- pi * ((12.5 / 2)^2)      # patch area, in cm^2 (12.5 cm diameter patch)
  patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2
  Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.
}

SECC <- within( SECC, { 
  ## change negative ARA values to 0 - should I wait until after aggregation?
  ARA.ml[ARA.ml < 0] <- 0
  ARA.m[ ARA.m  < 0] <- 0
  ARA.g[ ARA.g  < 0] <- 0
  Nfix <- ARA.m * Nfix.ARA.ratio
})


### Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

##================================================
## CUSTOM LABELS
##================================================

Y.lim <- c(-2, 18) # consistent Y-axis (transformed)



##================================================
## Data Exploration
##================================================
if (FALSE) {
  SECC.t1 <- SECC[which(SECC$Time==levels(SECC$Time)[1] & 
                  SECC$Chamber==levels(SECC$Chamber)[3]), ]
  boxplot(ARA.ml~Block, data =SECC.t1)
  SECC.t1[which(SECC.t1$ARA.ml == max(SECC.t1$ARA.ml, na.rm = TRUE)), c('SampleID', 'ARA.ml')]
  SECC.t1[which(SECC.t1$ARA.ml == min(SECC.t1$ARA.ml, na.rm = TRUE)), c('SampleID', 'ARA.ml')]
}



##################################################
### RUN STANDARD nested ANOVA
##################################################
### Run analysis on each Time point in sequence.
for ( Time.i in 1:length(levels(SECC$Time)) ) {
  ## Specify which treatment levels to include (by index is probably easiest)
  Time.use     <- levels(SECC$Time)[Time.i]      # Time (index: 1-3) to include in this run
  cat("\n\n\nProcessing Time:", Time.use, "\n")

  ## Load default Labels - dependent on above settings. *****
  source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

  ## RUN STANDARD nested ANOVA
  source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)

}


###===============================================
### Include Time as a factor in nested ANOVA
###===============================================
## Note that Samples at different times are actually independent
## in this design, due to destructive sampling.

Time.use     <- levels(SECC$Time)      # Include *ALL* Times (as a Treatment)
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) # Load default Labels. *****
source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE) # RUN STANDARD nested ANOVA


##################################################
### PUBLICATION GRAPHS
##################################################

Y.lim <- c(-10, 400)
Plot.Title <- bquote(.(Time.label) * "Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 
Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]

## data frame of plot values (for ggplot2).
## might be able to accomplish much the same effect with stat_summary using means in ggplot2?
plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Chamber=SECCp$Chamber, 
                                     Position=SECCp$Position, Time=SECCp$Time), 
                                error = msd["Time:Chamber:Position"]
                                )
levels(plot.means$Chamber)[2] <- "Chamber"

CxP.plot <- qplot(Chamber, x, data = plot.means, group = Position, 
                    geom = "line", ylim = Y.lim, size = Position,
                    colour = Position, shape = Position, fill = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = .~Time)
CxP.plot <- CxP.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                         width = 0.2, size = 0.5)
CxP.plot <- CxP.plot + geom_point(aes(group = Position), size = 3)
CxP.plot <- CxP.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_shape_manual(name = Position.label,
                                          values = Position.map$pch, 
                                          breaks = Position.map$label)
CxP.plot <- CxP.plot + scale_size_manual(name = Position.label,
                                         values = Position.map$lwd*0.5, 
                                         breaks = Position.map$label)
CxP.plot <- CxP.plot + jaw.ggplot()
print(CxP.plot)

## Same plot with internal legend, for Oecologia
CxP.legend <- CxP.plot + 
                opts(legend.position = c(0.54, 0.6),  # position legend inside main graph (for export dimensions)
                     legend.text = theme_text(size = 10),
                     legend.key.size = unit(1.5, "lines"),
                     legend.title = theme_blank()
                )
                # opts(legend.position = "none")
print(CxP.legend)


if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 4, scale = 1)
  ggsave(file = paste(Save.final, "- CPL.eps"), plot = CxP.legend, width = 6, height = 4, scale = 1) # for Oecologia
}
