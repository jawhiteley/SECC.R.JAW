##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Moisture Content (H2O % of moss dry wt)
### Jonathan Whiteley     R v2.12     2012-07-21
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
Y.col <- 'H2O'  # Column to analyze as response variable           *****
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****


##================================================
## CUSTOM CALCULATIONS 
##================================================

SECC <- within( SECC, { 
    H2O.asq <- asin(sqrt(H2O.wwt))
    H2O     <- H2O     * 100  # convert to %
    H2O.wwt <- H2O.wwt * 100
})

attr(SECC, "labels")[["H2O.asq"]] <- "Moisture"
attr(SECC, "units" )[["H2O.asq"]] <- quote(asin(sqrt("% "* H[2]*O)))
attr(SECC, "labels")[["H2O"]] <- "Moisture"
attr(SECC, "units" )[["H2O"]] <- quote(""* H[2]*O *" as % of dry wt. moss")

##==============================================================
## Data Exploration & Checking
##==============================================================
if (FALSE)
{
  with(SECC, hist(H2O) )

  library(ggplot2)
  library(mgcv)

  ## I originally wanted to try a quick regression of H2O on temp, but then I realized I don't have that kind of temperature data.
  H2O_temp.plot <- ggplot(SECC,
                    aes(y = Prod12 + Prod23, x = H2O, colour = Chamber, shape = Position)
  ) +
  geom_point() + stat_smooth(method = "gam")
  print(H2O_temp.plot)

}


##================================================
## CUSTOM SETTINGS 
##================================================
### Load default settings (based on response variable) *****
source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

## delete lines to use the defaults.

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)               # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include

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
### sub-analysis of final inner & outer only (for comparison with other analyses)
###===============================================
Position.Main <- Position.use                   # store temporarilly
Position.use  <- levels(SECC$Position)[c(1, 3)] # Inner & Outer only
Time.use     <- levels(SECC$Time)[3]      # Time (index: 1-3) to include in this run
## Load default Labels - dependent on above settings. *****
source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 
Save.text  <- gsub(".txt", "-IO.txt", Save.text,  fixed=TRUE)
Save.plots <- gsub(".pdf", "-IO.txt", Save.plots, fixed=TRUE)

  ## RUN STANDARD nested ANOVA
  source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)

msd.4IO <- msd                                  # store for later?
Position.Main -> Position.use                   # restore original value


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
library(ggplot2)

Y.lim <- c(0, 800)
Plot.Title <- bquote(.(Time.label) * "Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

if (FALSE) {
  postscript(file = Save.final, width = 6, height = 2)


Chamber.map <- plotMap( factor = "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Frag.map <- plotMap( factor = "Frag", labels = levels(SECC$Frag) )
Frag.map <- Frag.map[ levels(SECC$Frag) %in% Frag.use, ]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
Position.col <- Position.map$col
names(Position.col) <- Position.map$label

old.par <- par(mfrow=c(1,3), las = 1, oma = c(3, 2, 3, 1), mar = c(3, 3, 2, 0) +0.1 )

for(Time.for in 1:length(Time.use)) {
with( SECCp[SECCp$Time == Time.use[Time.for], ], {
  plot.error <- matrix( as.numeric(msd["Chamber:Position"]/2),
                       nrow = length(levels(Chamber)),  # rows: x-axis
                       ncol = length(levels(Position))  # cols: y-axis
                       )
  ## custom plotMeans function, with custom error bars (LSD)
  plotMeans( Y.trans, Chamber, Position, 
            error.bars = "custom", level = plot.error, ylim = Y.lim,
            cex = 2, lwd = 2, lty = Position.map$lty, pch = Position.map$pch,
            yaxt = ifelse(Time.for > 1, "n", "s"),
            col = as.character(Position.map$col),
            bg  = as.character(Position.map$bg),
            main = Time.use[Time.for],
            sub  = "",
            xlab = "",
            ylab = ""
            )   
  })
}
mtext(Plot.Title, side = 3, outer = TRUE)
mtext(attr(SECC, "labels")[["Chamber"]], side = 1, outer = TRUE)
mtext(Sub.msd, side = 1, padj = 1.5, outer = TRUE)
mtext(Y.plotlab, side = 2, outer = TRUE, las = 0)

par(old.par)
dev.off()
}

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
plot.means$label <- as.character(plot.means$Position)
plot.means$label[plot.means$Chamber == "Ambient"] <- NA
plot.means$label[plot.means$Time != "August\n12 months"] <- NA
plot.means$label[is.na(plot.means$label)] <- ""
plot.means$label[plot.means$label == "Inner"] <- "Inner\n\n"

CxP.plot <- qplot(Chamber, x, data = plot.means, group = Position, 
                    geom = "line", ylim = Y.lim, size = Position,
                    colour = Position, shape = Position, fill = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Chamber"]],
                    ylab = Y.plotlab,
                    legend = FALSE) +
            facet_grid(facets = .~Time)
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
                ## geom_text(aes(label = label), size = 4, hjust = 1.5, vjust = 0.5) +
                ## facet_grid(facets = .~Time, scales = "free_x", space = "free") +  # adjust panel sizes to show text annotations?
                opts(legend.position = c(0.24, 0.69),  # position legend inside main graph (for export dimensions)
                     legend.title = theme_blank(),
                     legend.text  = theme_text(size = 10, lineheight=1),
                     legend.key.width  = unit(1.5, "lines"),
                     legend.key.height = unit(1.0, "lines"),
                     panel.grid.major = theme_blank(),
                     panel.grid.minor = theme_blank()
                     )
                # opts(legend.position = "none")
print(CxP.legend)


## t4-only for follow-up manuscripts
CP4.plot <- qplot(Chamber, x, data = subset(plot.means, Time == "August\n24 months"), 
                  group = Position,  size = Position,
                  colour = Position, shape = Position, fill = Position,
                  geom = "line", ylim = Y.lim,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Chamber"]],
                  ylab = Y.plotlab, legend = FALSE)
CP4.plot <- CP4.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                         width = 0.2, size = 0.5)
CP4.plot <- CP4.plot + geom_point(aes(group = Position), size = 3)
CP4.plot <- CP4.plot + scale_colour_manual(name = Position.label,
                                           values = Position.map$col, 
                                           breaks = Position.map$label)
CP4.plot <- CP4.plot + scale_fill_manual(name = Position.label,
                                         values = Position.map$bg, 
                                         breaks = Position.map$label)
CP4.plot <- CP4.plot + scale_shape_manual(name = Position.label,
                                          values = Position.map$pch, 
                                          breaks = Position.map$label)
CP4.plot <- CP4.plot + scale_size_manual(name = Position.label,
                                         values = Position.map$lwd*0.5, 
                                         breaks = Position.map$label)
CP4.plot <- CP4.plot + jaw.ggplot()
print(CP4.plot)


# Inner & Outer only
IO.map <- subset(Position.map, label!="intermediate")

CPIO.plot <- qplot(Chamber, x, data = droplevels( subset(plot.means, Position!="intermediate") ), 
                  group = Position,  size = Position,
                  colour = Position, shape = Position, fill = Position,
                  geom = "line", ylim = Y.lim,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Chamber"]],
                  ylab = Y.plotlab, legend = FALSE) +
            facet_grid(facets = .~Time)
CPIO.plot <- CPIO.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                         width = 0.2, size = 0.5)
CPIO.plot <- CPIO.plot + geom_point(aes(group = Position), size = 3)
CPIO.plot <- CPIO.plot + scale_colour_manual(name = Position.label,
                                               values = IO.map$col, 
                                               breaks = IO.map$label)
CPIO.plot <- CPIO.plot + scale_fill_manual(name = Position.label,
                                             values = IO.map$bg, 
                                             breaks = IO.map$label)
CPIO.plot <- CPIO.plot + scale_shape_manual(name = Position.label,
                                              values = IO.map$pch, 
                                              breaks = IO.map$label)
CPIO.plot <- CPIO.plot + scale_size_manual(name = Position.label,
                                             values = IO.map$lwd*0.5, 
                                             breaks = IO.map$label)
CPIO.plot <- CPIO.plot + jaw.ggplot()
print(CPIO.plot)

plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Chamber=SECCp$Chamber, 
                                     Position=SECCp$Position, Time=SECCp$Time), 
                                error = msd.4IO["Chamber:Position"]
                                )
levels(plot.means$Chamber)[2] <- "Chamber"
CP4IO.plot <- qplot(Chamber, x, data = droplevels( subset(plot.means, 
                                                          Time == "August\n24 months" & Position!="intermediate") ), 
                  group = Position,  size = Position,
                  colour = Position, shape = Position, fill = Position,
                  geom = "line", ylim = Y.lim,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Chamber"]],
                  ylab = Y.plotlab, legend = FALSE)
CP4IO.plot <- CP4IO.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                         width = 0.2, size = 0.5)
CP4IO.plot <- CP4IO.plot + geom_point(aes(group = Position), size = 3)
CP4IO.plot <- CP4IO.plot + scale_colour_manual(name = Position.label,
                                               values = IO.map$col, 
                                               breaks = IO.map$label)
CP4IO.plot <- CP4IO.plot + scale_fill_manual(name = Position.label,
                                             values = IO.map$bg, 
                                             breaks = IO.map$label)
CP4IO.plot <- CP4IO.plot + scale_shape_manual(name = Position.label,
                                              values = IO.map$pch, 
                                              breaks = IO.map$label)
CP4IO.plot <- CP4IO.plot + scale_size_manual(name = Position.label,
                                             values = IO.map$lwd*0.5, 
                                             breaks = IO.map$label)
CP4IO.plot <- CP4IO.plot + jaw.ggplot()
print(CP4IO.plot)



## Frag x Pos Interaction

Y.lim <- c(-100, 1000)
plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, Position=SECCp$Position, Time=SECCp$Time), mean)
levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
plot.means$error <- as.numeric(msd["Time:Chamber:Frag:Position"]/2)
plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, 
                                     Position=SECCp$Position, Time=SECCp$Time), 
                                error = msd["Time:Chamber:Frag:Position"]
                                )
levels(plot.means$Chamber)[2] <- "Chamber"

FragIconList <- list(FragIcon1 = FragIcon1,
                     FragIcon2 = FragIcon2,
                     FragIcon3 = FragIcon3,
                     FragIcon4 = FragIcon4
                     )

FxP.plot <- qplot(Frag, x, data = plot.means, group = Position, 
                    geom = "line", ylim = Y.lim, size = Position, 
                    colour = Position, fill = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Frag"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = Chamber~Time)
FxP.plot <- FxP.plot + 
geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
geom_point(aes(shape = Position), size = 3) +
scale_colour_manual(name = Position.label,
                    values = Position.map$col, 
                    breaks = Position.map$label) +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label) +
scale_shape_manual(name = Position.label,
                   values = Position.map$pch, 
                   breaks = Position.map$label) +
scale_size_manual(name = Position.label,
                  values = Position.map$lwd*0.5, 
                  breaks = Position.map$label) +
scale_x_discrete(labels = c(1, 2, 3, 4), 
                 breaks = levels(plot.means$Frag)) + jaw.ggplot() 
## Add imported graphics as x-axis tick labels :D
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
FxP.plot <- FxP.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 3, 4), 
                                        breaks = levels(plot.means$Frag)) +
opts(axis.ticks.margin = unit(0.2, "lines"),
     axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
)
print(FxP.plot)


## t4-only for follow-up manuscripts
FP4.plot <- qplot(Frag, x, data = subset(plot.means, Time == "August\n24 months"), 
                  group = Position,  size = Position, 
                  colour = Position, fill = Position, shape = Position,
                  geom = "line", ylim = Y.lim,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Frag"]],
                  ylab = Y.plotlab,
                  legend = FALSE,
                  facets = .~Chamber)
FP4.plot <- FP4.plot + 
geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
geom_point(aes(shape = Position), size = 3) +
scale_colour_manual(name = Position.label,
                    values = Position.map$col, 
                    breaks = Position.map$label) +
     scale_fill_manual(name = Position.label,
                       values = Position.map$bg, 
                       breaks = Position.map$label) +
     scale_shape_manual(name = Position.label,
                        values = Position.map$pch, 
                        breaks = Position.map$label) +
     scale_size_manual(name = Position.label,
                       values = Position.map$lwd*0.5, 
                       breaks = Position.map$label) +
     scale_x_discrete(labels = c(1, 2, 3, 4), 
                      breaks = levels(plot.means$Frag)) + jaw.ggplot() 
     ## Add imported graphics as x-axis tick labels :D
     ## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
     FP4.plot <- FP4.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 3, 4), 
                                             breaks = levels(plot.means$Frag)) +
     opts(axis.ticks.margin = unit(0.2, "lines"),
          axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
     )
print(FP4.plot)



if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 3, scale = 1.2)
  ggsave(file = paste(Save.final, "- CxP-IO.eps"), plot = CPIO.plot, width = 6, height = 3, scale = 1.2)
  ggsave(file = paste(Save.final, "- FxP.eps"), plot = FxP.plot, width = 6, height = 4, scale = 1.2)
  ggsave(file = paste(Save.final, "- FP4.eps"), plot = FP4.plot, width = 6, height = 4, scale = 1)
  ggsave(file = paste(Save.final, "- CPL.eps"), plot = CxP.legend, width = 6, height = 4, scale = 1) # for Oecologia
  ggsave(file = paste(Save.final, "- CP4.eps"), plot = CP4.plot, width = 4, height = 4, scale = 1)
  ggsave(file = paste(Save.final, "- CP4-IO.eps"), plot = CP4IO.plot, width = 4, height = 4, scale = 1)
}
