##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Aggregate Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2012-04-01
###===============================================
### ** Typically called from 'Fauna.R'  **
###    - all the pre-processing & filtering happens there.
### Species identified to morphospecies (usually family or genus level)
### Counts / sample converted to # / g dwt of moss substrate (using 'Patch.dwt' column)
### This script aggregates data over Fragmentation treatments ("Meta-communities"),
### and performs basic univariate analyses

##################################################
## INITIALISE
##################################################

if (FALSE) {  # do not run automatically
  ## Working Directory: see lib/init.R below
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (Mac).
  getwd()  # Check that we're in the right place

  ## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
  source('./lib/init.R')

}

## Should already have the following in memory:
## SECC.fauna      raw data (counts / g dwt)
## SECC.fauna.meta metadata about species
## SECC.fauna.sum  summary data by major taxonomic groups
## SECC.sp         filtered & aggregated to 'real' species
## SECC.sp.sum     summary data for filtered data


## Most scripts assume analysis on 'SECC' data frame, 
## but relevant fauna data is stored in 'SECC.fauna'
SECC.df <- SECC     # temporary storage for this script (restored at the end)
if (FALSE) {
  SECC <- SECC.df   # restore
}

### ANOVA Response Variable *****
if (!exists('Y.col')) Y.col <- 'Richness' # Column to analyze as response variable           *****
if (!exists('Y.use')) Y.use <- 'Y'        # Which transformation is being used (for labels)? *****

##==============================================================
## Process Data
##==============================================================
## Aggregate Fauna counts across regions / meta-communities (Frag treatments)
MCID <- substr(rownames(Fauna), 1, 5)
Factors.mc <- data.frame(SampleID = rownames(Fauna), # for recoding & checking
                         Block   = substr(MCID, 1, 1),
                         Time    = substr(MCID, 2, 2),
                         Chamber = substr(MCID, 3, 3),
                         Frag    = substr(MCID, 5, 5),
                         Pos     = substr(rownames(Fauna), 7, 7) # for recoding only.
                         )
Factors.mc <- checkSECCdata( Factors.mc, "Fauna.factors" )
Factors.mc <- recodeSECC( Factors.mc )
## the lack of some Pos values will trigger warnings.  That's ok, it doesn't really matter.
if ( any( MCID != substr(Factors.mc$SampleID, 1, 5) ) )
  ## check that things were not rearranged
  stop("Fauna factors were rearranged during processing!")
## Keep only necessary columns, and every other row (aggregated)
Factors.mc <- Factors.mc[seq(1, NROW(Factors.mc), by = 2), c("Block", "Time", "Chamber", "Frag")]

Fauna.mc <- aggregate(Fauna, by = list(Sample = MCID), FUN = sum)
rownames(Fauna.mc) <- Fauna.mc$Sample
Fauna.mc <- Fauna.mc[, -1] ## drop 'Sample' column
## assuming both data frames are still in the same order...
Fauna.mc.sum <- within(Factors.mc, {
        Richness <- apply(Fauna.mc, 1, function(x) length(which(x>0)) )  # observed # spp.
        Evenness <- diversity(Fauna.mc, index = "invsimpson")
                         })

## Calculate Group totals
##  as in lib/load_fauna.R
Taxa.groups <- unique(SECC.fauna.meta$Taxonomic.Group)

Fauna.mc.sum <- within(Fauna.mc.sum, {
### Mesostigs
### Collembola
### Prostigs
### Other
    for (taxa in Taxa.groups) {
      taxa.cols <- SECC.fauna.meta$ID[which(SECC.fauna.meta$Taxonomic.Group == taxa)]
      taxa.cols <- intersect(colnames(Fauna.mc), taxa.cols)
      assign( taxa, apply(Fauna.mc[, taxa.cols], 1, sum) )
    }

### Uropodina (non-predatory Mesostigs)
    Uropodina <- apply(Fauna.mc[, SECC.fauna.meta$ID[which(SECC.fauna.meta$Major.Taxa == "Uropodina")] ], 1, sum)
### Mesostig.preds = Mesostigs - Uropodina
    Mesostig.preds <- Mesostigmata - Uropodina
### Predators = Mesostigs.preds + Prostigs
    Predators <- Mesostig.preds + Prostigmata
### Grazers = Collembola + Uropodina (+ Oribatids)
    Grazers <- Collembola + Uropodina
### fauna.jaw = Mesostigs + Collembola (+ Prostigs?)
    fauna.jaw <- Mesostigmata + Collembola
### fauna = all (-Other) * including ZL data?
    taxa.cols <- intersect(colnames(Fauna.mc), SECC.fauna.meta$ID)
    taxa.cols <- setdiff(taxa.cols, SECC.fauna.meta$ID[which(SECC.fauna.meta$Taxonomic.Group == "Other")])
    fauna <- apply(Fauna.mc[, taxa.cols], 1, sum)
    rm(taxa, taxa.cols)
})
Fauna.mc.sum <- Fauna.mc.sum[, c('Block', 'Time', 'Chamber', 'Frag',
                              'Mesostigmata', 'Collembola', 'Prostigmata', 'Other',
                              'Uropodina', 'Mesostig.preds', 'Predators', 'Grazers', 
                              'fauna.jaw', 'fauna',
                              'Richness', 'Evenness')]  # manually reorder columns 


##################################################
## DATA EXPLORATION
##################################################

if (FALSE) {

  plot(Fauna.mc.sum[, which( sapply(Fauna.mc.sum, is.numeric) )])

  MC.Rich.Frag <- ggplot(SECC.mc.sum, aes(x = Frag, y = Richness, colour = Chamber)) +
    stat_summary(fun.data = "mean_cl_boot")
  print(MC.Rich.Frag)


}





##################################################
## CONFIGURE BASIC ANALYSIS
##################################################



### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

##================================================
## CUSTOM SETTINGS 
##================================================
## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[3]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include: A, C
Frag.use     <- levels(SECC$Frag)[c(1, 2, 4)]   # Fragmentation treatments to include: 1, 2, 4
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

attr(SECC, "labels")[["Frag"]] <- "Habitat Isolation" # different interpretation, particularly as far as the fauna is concerned.


##################################################
### RUN STANDARD nested ANOVA
##################################################
source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)

if (FALSE) 
{

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

}


##################################################
### PUBLICATION GRAPHS
##################################################
Plot.Title <- bquote(.(Time.label) * "\nPatch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
FragIconList <- list(FragIcon1 = FragIcon1,
                     FragIcon2 = FragIcon2,
                     FragIcon4 = FragIcon4
                     )

## 3-way interaction: NS
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, Position=SECCp$Position, Time=SECCp$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means <- within(plot.means, 
                     {
                       error <- as.numeric(msd["Chamber:Frag:Position"]/2)
                       upper <- x + error
                       lower <- x - error
                       levels(Chamber)[2] <- "Chamber"
                     })

## if (!exists('Y.lim')) Y.lim <- c(-1, 20)
Y.lim <- with(plot.means, range(lower, upper))
Y.lim <- c(floor(Y.lim[1]/10), ceiling(Y.lim[2]/10) ) *10

CFP.plot <- qplot(Frag, x, data = plot.means, group = Position, 
                    geom = "line", ylim = Y.lim, size = I(0.8), 
                    colour = Position, fill = Position, 
                    shape = Position, # lty = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Frag"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = .~Chamber)
CFP.plot <- CFP.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                         width = 0.2, size = 0.5)
CFP.plot <- CFP.plot + geom_point(aes(group = Position), size = 3)
CFP.plot <- CFP.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
CFP.plot <- CFP.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
CFP.plot <- CFP.plot + scale_shape_manual(name = Position.label,
                                           values = Position.map$pch, 
                                           breaks = Position.map$label)
## CFP.plot <- CFP.plot + scale_linetype_manual(name = Position.label,
##                                              values = Position.map$lty, 
##                                              breaks = Position.map$label)
CFP.plot <- CFP.plot + jaw.ggplot()
## Add imported graphics as x-axis tick labels :D
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
CFP.plot <- CFP.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 4), 
                                        breaks = levels(plot.means$Frag)) +
opts(axis.ticks.margin = unit(0.2, "lines"),
     axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
)
print(CFP.plot)


## Frag x Position
plot.means <- aggregate(SECCp$Y.trans, list(Frag=SECCp$Frag, Position=SECCp$Position, Time=SECCp$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means <- within(plot.means, 
                     {
                       error <- as.numeric(msd["Frag:Position"]/2)
                       upper <- x + error
                       lower <- x - error
                     })

if (exists('Y.lim1'))
{
  Y.lim <- Y.lim1
} else {
  Y.lim <- with(plot.means, range(lower, upper))
  Y.lim <- c(floor(Y.lim[1]/5), ceiling(Y.lim[2]/5) ) *5
}

FxP.plot <- qplot(Frag, x, data = plot.means, group = Position, 
                  geom = "line", ylim = Y.lim, size = I(0.8), 
                  colour = Position, fill = Position, 
                  shape = Position, # lty = Position,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Frag"]],
                  ylab = Y.plotlab,
                  legend = FALSE)
FxP.plot <- FxP.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                     width = 0.2, size = 0.5)
FxP.plot <- FxP.plot + geom_point(aes(group = Position), size = 3)
FxP.plot <- FxP.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
FxP.plot <- FxP.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
FxP.plot <- FxP.plot + scale_shape_manual(name = Position.label,
                                          values = Position.map$pch, 
                                          breaks = Position.map$label)
## FxP.plot <- FxP.plot + scale_linetype_manual(name = Position.label,
##                                              values = Position.map$lty, 
##                                              breaks = Position.map$label)
FxP.plot <- FxP.plot + jaw.ggplot()
## Add imported graphics as x-axis tick labels :D
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
FxP.plot <- FxP.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 4), 
                                        breaks = levels(plot.means$Frag)) +
     opts(axis.ticks.margin = unit(0.2, "lines"),
          axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
     )
print(FxP.plot)


## Chamber x Position
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Position=SECCp$Position, Time=SECCp$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means <- within(plot.means, 
                     {
                       error <- as.numeric(msd["Chamber:Position"]/2)
                       upper <- x + error
                       lower <- x - error
                       levels(Chamber)[2] <- "Chamber\n"
                     })

if (exists('Y.lim1'))
{
  Y.lim <- Y.lim1
} else {
  Y.lim <- with(plot.means, range(lower, upper))
  Y.lim <- c(floor(Y.lim[1]/5), ceiling(Y.lim[2]/5) ) *5
}

CxP.plot <- qplot(Chamber, x, data = plot.means, group = Position, 
                  geom = "line", ylim = Y.lim, size = I(0.8), 
                  colour = Position, fill = Position, 
                  shape = Position, # lty = Position,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Chamber"]],
                  ylab = Y.plotlab,
                  legend = FALSE)
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
## CxP.plot <- CxP.plot + scale_linetype_manual(name = Position.label,
##                                              values = Position.map$lty, 
##                                              breaks = Position.map$label)
CxP.plot <- CxP.plot + jaw.ggplot() +
  ## adjust x-axis labelling to allow the plot to line up better next to corresponding FxP plot
  opts(axis.ticks.margin = unit(0.2, "lines"),
        axis.text.x = theme_text(lineheight = 1.1, vjust = 1) # trial & error 
  ## plot.margin = unit(c(1, 1, 1.4, 0.5), "lines")  # this lines up exactly, but moves the axis title up as well: that's not what I want.
  )
print(CxP.plot)



if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- FxP.eps"), plot = FxP.plot, width = 3, height = 4, scale = 1.5)
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 3, height = 4, scale = 1.5)
  ##   ggsave(file = paste(Save.final, "- CFP.eps"), plot = CFP.plot, width = 6, height = 4, scale = 1.5)
}



##################################################
### CLEAN-UP / HOUSEKEEPING
##################################################
SECC <- SECC.df # restore original
