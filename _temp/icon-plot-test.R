rm(list=ls())

library(grImport)                      # import eps files
setwd("./save/")
fragIcons <- PostScriptTrace("Frag-Black-icons.eps")
fragIcons <- readPicture("Frag-Black-icons.eps.xml")
Hex <- PostScriptTrace("hexagon.eps")
Hex <- readPicture("hexagon.eps.xml")
setwd("..")
FragIcon1 <- fragIcons[49:50]          # 1. Continuous 
FragIcon2 <- fragIcons[29:48]          # 2. Corridors
FragIcon3 <- fragIcons[ 9:28]          # 3. Pseudo-corridors
FragIcon4 <- fragIcons[ 1:8 ]          # 4. Isolated

## str(fragIcons)
picturePaths(fragIcons)
if (FALSE)
{
  plot.new()
  grid.picture(fragIcons)
  plot.new()
  grid.picture(fragIcons[1])             #  (first path)
  grid.picture(fragIcons[1:8])           # 4. Isolated
  grid.picture(fragIcons[9:28])          # 3. Pseudo-corridors
  grid.picture(fragIcons[29:48])         # 2. Corridors
  grid.picture(fragIcons[49:50])         # 1. Continuous
}

plot.new()
grid.picture(FragIcon1)
grid.picture(FragIcon2)
grid.picture(FragIcon3)
grid.picture(FragIcon4)

plot.new()
picturePaths(Hex)
plot.new()
grid.picture(Hex, width = 0.5, height = 0.5, just = "bottom", 
             use.gc = FALSE, gp = gpar(col = "#555555", fill = NA, 
                                       lwd = 4, lty = 3)
)
grid.picture(Hex, width = 0.5, height = 0.5, just = "top", 
             use.gc = FALSE, gp = gpar(col = "#000000", fill = "#CC9999", 
                                       lwd = 4, lty = 1)
)

vp1 <- viewport(width = 0.1, height = 0.1, x = 1/6, y = -0.2)




## including frag labels in H2O frag plot
source("H2O.R", echo=FALSE)

if (TRUE)                              # build basic plot, without x-axis labels (blank)
{
  Y.lim <- c(-100, 900)
  plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, Position=SECCp$Position, Time=SECCp$Time), mean)
  levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
  ## levels(plot.means$Frag) <- c(FragIcon1, FragIcon2, FragIcon3, FragIcon4)
  levels(plot.means$Frag) <- c(1, 2, 3, 4)
  plot.means$error <- as.numeric(msd["Time:Chamber:Frag:Position"]/2)
  levels(plot.means$Chamber)[2] <- "Chamber"

  FxP.plot <- qplot(Frag, x, data = plot.means, group = Position, 
                    geom = "point", ylim = Y.lim, size = I(3), 
                    colour = Position, shape = Position,
                    main = Plot.Title, sub = Sub.msd,
                    xlab = attr(SECC, "labels")[["Frag"]],
                    ylab = Y.plotlab,
                    legend = FALSE,
                    facets = Chamber~Time)
  FxP.plot <- FxP.plot + geom_line(aes(group = Position), size = 0.8)
  FxP.plot <- FxP.plot + geom_errorbar(aes(ymin = x - error, ymax = x + error), 
                                       width = 0.2, size = 0.5)
  FxP.plot <- FxP.plot + scale_colour_manual(name = Position.label,
                                             values = Position.map$col, 
                                             breaks = Position.map$label)
  FxP.plot <- FxP.plot + scale_fill_manual(name = Position.label,
                                           values = Position.map$bg, 
                                           breaks = Position.map$label)
  FxP.plot <- FxP.plot + scale_shape_manual(name = Position.label,
                                            values = Position.map$pch, 
                                            breaks = Position.map$label)
  FxP.plot <- FxP.plot + scale_x_discrete(labels = c(" ", "  ", "   ", "    "), 
                                          breaks = levels(plot.means$Frag))
  FxP.plot <- FxP.plot + jaw.ggplot()
  ##  + opts(axis.text.x = theme_blank())  # removes axis text, but does not leave space ...
}
print(FxP.plot)

current.vpTree(all=TRUE)
grid.ls(grobs = FALSE, viewports = TRUE)
grid.ls()

## move current viewport to axis text, and draw inside?
downViewport("axis_h-7-3")
iconSize  <- 0.9
iconHt    <- 0.3 
grid.symbols(FragIcon1, x = 0.13, y = iconHt, size = iconSize)
# how do I find the axis labels to line the symbols up with them auto-magically?

FxP.g <- ggplotGrob(FxP.plot)
axis.x.1 <- getGrob(FxP.g, gPath("axis_h-7-3"), grep=TRUE)

## Add Frag icons to axis labels: HACK!!
##  units relative to main viewport
iconSize  <- 0.035
iconHt    <- 0.054 
tickSp    <- 0.057
grid.symbols(FragIcon1, x = 0.137             , y = iconHt, size = iconSize)
grid.symbols(FragIcon2, x = 0.137 + (tickSp*1), y = iconHt, size = iconSize)
grid.symbols(FragIcon3, x = 0.137 + (tickSp*2), y = iconHt, size = iconSize)
grid.symbols(FragIcon4, x = 0.137 + (tickSp*3), y = iconHt, size = iconSize)
grid.symbols(FragIcon1, x = 0.378             , y = iconHt, size = iconSize)
grid.symbols(FragIcon2, x = 0.378 + (tickSp*1), y = iconHt, size = iconSize)
grid.symbols(FragIcon3, x = 0.378 + (tickSp*2), y = iconHt, size = iconSize)
grid.symbols(FragIcon4, x = 0.378 + (tickSp*3), y = iconHt, size = iconSize)
grid.symbols(FragIcon1, x = 0.620             , y = iconHt, size = iconSize)
grid.symbols(FragIcon2, x = 0.620 + (tickSp*1), y = iconHt, size = iconSize)
grid.symbols(FragIcon3, x = 0.620 + (tickSp*2), y = iconHt, size = iconSize)
grid.symbols(FragIcon4, x = 0.620 + (tickSp*3), y = iconHt, size = iconSize)



## Save it
postscript(file = paste(Save.final, "- FxP.eps"), width = 6, height = 4)
print(FxP.plot)
## draw grid.symbols ...
dev.off()
## dev.copy2eps() ?dev.copy2eps

