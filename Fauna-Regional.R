##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Aggregate Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2012-04-12
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


##==============================================================
## Process Data
##==============================================================
## Aggregate Fauna counts across regions / meta-communities (Frag treatments)
SECC.fauna.coded <- checkSECCdata(SECC.fauna[Samples.fauna, ], "SECC.fauna")
SECC.fauna.coded <- recodeSECC(SECC.fauna.coded)
## the lack of some Pos values may trigger warnings.  That's ok, it doesn't really matter.

MCID <- substr(rownames(SECC.fauna.coded), 1, 5)
SECC.fauna.mc <- aggregate(SECC.fauna.coded[, Spp.cols], 
                      by = list(Block   = SECC.fauna.coded$Block,
                                Time    = SECC.fauna.coded$Time,
                                Chamber = SECC.fauna.coded$Chamber,
                                Frag    = SECC.fauna.coded$Frag, 
                                Sample  = MCID # @ end for sorting
                                ), FUN = sum)
rownames(SECC.fauna.mc) <- SECC.fauna.mc$Sample
Factors.mc    <- SECC.fauna.mc[,  (1:4)] # keep trt columns
SECC.fauna.mc <- SECC.fauna.mc[, -(1:5)] # drop trt columns
Fauna.mc      <- SECC.fauna.mc[, Spp.fauna] # exclude Prostigs + Other
## assuming both data frames are still in the same order...
Fauna.mc.sum <- within(Factors.mc, 
                       {
                         ## Note that this is now GAMMA (regional) richness
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
      taxa.cols <- intersect(colnames(SECC.fauna.mc), taxa.cols)
      assign( taxa, apply(SECC.fauna.mc[, taxa.cols], 1, sum) )
    }

### Uropodina (non-predatory Mesostigs)
    Uropodina <- apply(SECC.fauna.mc[, SECC.fauna.meta$ID[which(SECC.fauna.meta$Major.Taxa == "Uropodina")] ], 1, sum)
### Mesostig.preds = Mesostigs - Uropodina
    Mesostig.preds <- Mesostigmata - Uropodina
### Predators = Mesostigs.preds + Prostigs
    Predators <- Mesostig.preds + Prostigmata
### Grazers = Collembola + Uropodina (+ Oribatids)
    Grazers <- Collembola + Uropodina
### fauna.jaw = Mesostigs + Collembola (+ Prostigs?)
    fauna.jaw <- Mesostigmata + Collembola
### fauna = all (-Other) * including ZL data?
    taxa.cols <- intersect(colnames(SECC.fauna.mc), SECC.fauna.meta$ID)
    taxa.cols <- setdiff(taxa.cols, SECC.fauna.meta$ID[which(SECC.fauna.meta$Taxonomic.Group == "Other")])
    fauna <- apply(SECC.fauna.mc[, taxa.cols], 1, sum)
    rm(taxa, taxa.cols)
})
Fauna.mc.sum <- Fauna.mc.sum[, c('Block', 'Time', 'Chamber', 'Frag',
                              'Mesostigmata', 'Collembola', 'Prostigmata', 'Other',
                              'Uropodina', 'Mesostig.preds', 'Predators', 'Grazers', 
                              'fauna.jaw', 'fauna',
                              'Richness', 'Evenness')]  # manually reorder columns 

attr(Fauna.mc.sum, "labels") <- attr(SECC.sp.sum, "labels") # :)
attr(Fauna.mc.sum,  "units") <- attr(SECC.sp.sum,  "units") # :)
attr(SECC.fauna.coded, "labels") <- attr(SECC.sp.sum, "labels") # :)
attr(SECC.fauna.coded,  "units") <- attr(SECC.sp.sum,  "units") # :)

##################################################
## DATA EXPLORATION
##################################################

if (FALSE) {

  plot(Fauna.mc.sum[, which( sapply(Fauna.mc.sum, is.numeric) )])

  MC.Rich.Frag <- ggplot(Fauna.mc.sum, aes(x = Frag, y = Richness, colour = Chamber)) +
    stat_summary(fun.data = "mean_cl_boot")
  print(MC.Rich.Frag)
  ## Looks like there might be an interesting interaction, but too much noise (sample size too small) :(
  ## Ambient treatment apppears to lose gamma richness with isolation; Chamber *increases*?
  ## suggests lack of nestedness...

}





##################################################
## CONFIGURE BASIC ANALYSIS
##################################################
## Most scripts assume analysis on 'SECC' data frame, 
## but relevant fauna data is stored in 'SECC.fauna'
SECC.df <- SECC     # temporary storage for this script (restored at the end)
if (FALSE) {
  SECC <- SECC.df   # restore
}

SECC  <- Fauna.mc.sum                  # for attributes
SECCa <- Fauna.mc.sum                  # Meta-Community (Frag) scale data

### ANOVA Response Variable *****
if (!exists('Y.col')) Y.col <- 'Richness' # Column to analyze as response variable           *****
if (!exists('Y.use')) Y.use <- 'Y'        # Which transformation is being used (for labels)? *****




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
## Position.use <- levels(SECC$Position)[c(1, 3)]  # Patch Positions to include: Inner, Outer

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
source("./SECCanova/SECC - nested ANOVA Frag.R", echo = FALSE)

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

Chamber.label <- "Chamber\nTreatment" # attr(SECC, "labels")[["Chamber"]]
Chamber.map <- plotMap( factor = "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label[2] <- "Chamber"
FragIconList <- list(FragIcon1 = FragIcon1,
                     FragIcon2 = FragIcon2,
                     FragIcon4 = FragIcon4
                     )

## Chamber x Frag
plot.means <- aggregate(SECCmc$Y.trans, list(Chamber=SECCmc$Chamber, Frag=SECCmc$Frag, Time=SECCmc$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means <- within(plot.means, 
                     {
                       error <- as.numeric(msd.mc["Chamber:Frag"]/2)
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

CxF.plot <- qplot(Frag, x, data = plot.means, group = Chamber, 
                  geom = "line", ylim = Y.lim, size = I(0.8), 
                  colour = Chamber, fill = Chamber, 
                  shape = Chamber, # lty = Chamber,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Frag"]],
                  ylab = Y.plotlab,
                  legend = FALSE)
CxF.plot <- CxF.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                     width = 0.2, size = 0.5)
CxF.plot <- CxF.plot + geom_point(aes(group = Chamber), size = 3)
CxF.plot <- CxF.plot + scale_colour_manual(name = Chamber.label,
                                           values = Chamber.map$col, 
                                           breaks = Chamber.map$label)
CxF.plot <- CxF.plot + scale_fill_manual(name = Chamber.label,
                                         values = Chamber.map$bg, 
                                         breaks = Chamber.map$label)
CxF.plot <- CxF.plot + scale_shape_manual(name = Chamber.label,
                                          values = Chamber.map$pch, 
                                          breaks = Chamber.map$label)
## CxF.plot <- CxF.plot + scale_linetype_manual(name = Chamber.label,
##                                              values = Chamber.map$lty, 
##                                              breaks = Chamber.map$label)
CxF.plot <- CxF.plot + jaw.ggplot()
## Add imported graphics as x-axis tick labels :D
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
CxF.plot <- CxF.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 4), 
                                        breaks = levels(plot.means$Frag)) +
     opts(axis.ticks.margin = unit(0.2, "lines"),
          axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
     )
print(CxF.plot)



if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxF.eps"), plot = CxF.plot, width = 6, height = 4, scale = 1.5)
}



##################################################
### CLEAN-UP / HOUSEKEEPING
##################################################
SECC <- SECC.df # restore original
