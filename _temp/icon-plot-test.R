rm(list=ls())

library(grImport)                      # import eps files
setwd("./save/")
fragIcons <- PostScriptTrace("Frag-Black-icons.eps")
fragIcons <- readPicture("Frag-Black-icons.eps.xml")
Hex <- PostScriptTrace("hexagon.eps")
Hex <- readPicture("hexagon.eps.xml")
setwd("..")
FragIcon1 <- fragIcons[49:50]          # 1. Contiguous 
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
  grid.picture(fragIcons[49:50])         # 1. Contiguous
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



## simple plot
plot(1:25, rep(1, 25))
grid.symbols(FragIcon2, 1:25, rep(1.2, 25), units = "native", size = unit(5, "mm"))



################################################################
## including frag labels in H2O frag plot
source("H2O.R", echo=FALSE)

if (TRUE)                              # build basic plot, without x-axis labels (blank)
{
  Y.lim <- c(-100, 1000)
  plot.means <- aggregate(SECCp$Y.trans, list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, Position=SECCp$Position, Time=SECCp$Time), mean)
  levels(plot.means$Time) <- paste(c("August", "June", "August"), levels(plot.means$Time), sep="\n")
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

## Add Frag Icons: still Hack, but more accurate.
iconSize  <- 0.9
iconHt    <- 0.3 
expand    <- 0.025                     # ggplot2 $options$panel.margin?  expansion factor of axis in ggplot2?
tickSp    <- (1-expand)/8
for (vp in c("axis_h-7-3", "axis_h-7-5", "axis_h-7-7")) 
{
  downViewport(vp)
  grid.symbols(FragIcon1, x = (expand/2 + tickSp*1), y = iconHt, size = iconSize)
  grid.symbols(FragIcon2, x = (expand/2 + tickSp*3), y = iconHt, size = iconSize)
  grid.symbols(FragIcon3, x = (expand/2 + tickSp*5), y = iconHt, size = iconSize)
  grid.symbols(FragIcon4, x = (expand/2 + tickSp*7), y = iconHt, size = iconSize)
  upViewport()
}



current.vpTree(all=TRUE)
grid.ls(grobs = FALSE, viewports = TRUE)
grid.ls()

## move current viewport to graph panel, and draw inside?
## downViewport("panel.background.rect.5927")
downViewport("panel-3-3")
iconSize  <- 1
iconHt    <- 0 
grid.symbols(FragIcon2, x = unit(1, "native"), y = iconHt, size = iconSize, units="native")
## why are units not "native" to viewport?
convertX(unit(0:1, "npc"), "native", valueOnly=TRUE)
## rather: why are there no "native" units (xscale, yscale) in the graphing panel?!
## the scales seem to be associated with geoms, not the viewports themselves?

## move current viewport to axis text, and draw inside?
downViewport("axis_h-7-3")
iconSize  <- 0.9
iconHt    <- 0.3 
expand    <- 0.03                      # expansion factor of axis in ggplot2?
tickSp    <- (1-expand)/8
grid.symbols(FragIcon1, x = 1, y = iconHt, size = iconSize, units = "native")
# how do I find the axis labels to line the symbols up with them auto-magically?
grid.text(1:4, x = unit(1:4, "native"), default.units = "native")


FxP.g <- ggplotGrob(FxP.plot)
axis.x.1 <- getGrob(FxP.g, gPath("axis_h-7-3"), grep=TRUE)

FxP <- ggplot_build(FxP.plot)
str(FxP$plot)

FxP <- plot_clone(FxP.plot)
FxP$options$axis.text.x
FxP$facet$clone

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
grid.symbols(FragIcon1, x = 0.619             , y = iconHt, size = iconSize)
grid.symbols(FragIcon2, x = 0.619 + (tickSp*1), y = iconHt, size = iconSize)
grid.symbols(FragIcon3, x = 0.619 + (tickSp*2), y = iconHt, size = iconSize)
grid.symbols(FragIcon4, x = 0.619 + (tickSp*3), y = iconHt, size = iconSize)



## Save it
dev.copy2eps(file = paste(Save.final, "- FxP_icon.eps"), width = 9, height = 6) # equivalent relative output to ggsave( 6x4, scale=1.5), except ... it's 1.5x bigger in absolute size :P

dev.copy2eps(file = paste(Save.final, "- FxP_icon.eps"))

postscript(file = paste(Save.final, "- FxP_icon.eps"), width = 6, height = 4)
print(FxP.plot)
## draw grid.symbols ...
downViewport("axis_h-7-3")
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
grid.symbols(FragIcon1, x = 0.619             , y = iconHt, size = iconSize)
grid.symbols(FragIcon2, x = 0.619 + (tickSp*1), y = iconHt, size = iconSize)
grid.symbols(FragIcon3, x = 0.619 + (tickSp*2), y = iconHt, size = iconSize)
grid.symbols(FragIcon4, x = 0.619 + (tickSp*3), y = iconHt, size = iconSize)
dev.off()
## dev.copy2eps() ?dev.copy2eps





##==============================================================
## example for Stack Overflow question
## How can I use a graphic imported with grImport as axis tick labels in ggplot2 (using grid functions)?
## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
## http://github.com/hadley/ggplot2/wiki/Editing-raw-grid-objects-from-a-ggplot

word1 <- regexpr("\\s", paste(row.names(mtcars), "")) 
cars <- within(mtcars, Company <- substr(row.names(mtcars), 1, word1-1) )
cars <- within(mtcars, cyl <- factor(cyl) )
numLvls <- length(levels(cars$cyl))


library(ggplot2)
## library(grImport)  # not needed for this example, but would be for grid.symbols()


p <- ggplot(mtcars, aes(cyl, mpg)) + stat_summary(fun.data = "mean_cl_boot")
print(p)

## Replace (in this case, overlay) x-axis tick labels with a graphic / grob
iconSize  <- 0.05
iconHt    <- 0.2 
padding   <- 0.09                      # horizontal padding around axis: trial & error
tickSp    <- (1-padding)/(4*2)
downViewport("axis_h-5-3")
## I would use grid.symbols() with an imported Picture in place of grid.circle(),
## but the idea is the same: draw a shape at the ticks along the axis.
for (i in 0:(max(mtcars$cyl) - min(mtcars$cyl)) )
{
  grid.circle(x = (padding/2 + tickSp*(i*2)), y = iconHt, 
              r = iconSize*(min(mtcars$cyl)+i), gp = gpar(fill="black"))
}

upViewport()




################################################################
## kohske: custom grob axis labels (improved).  GENIUS     *****
## http://stackoverflow.com/questions/8905101/how-can-i-use-a-graphic-imported-with-grimport-as-axis-tick-labels-in-ggplot2-u

# convert ps to RGML
PostScriptTrace(file.path(system.file(package = "grImport"), "doc", "GNU.ps"), "GNU.xml")
PostScriptTrace(file.path(system.file(package = "grImport"), "doc", "tiger.ps"), "tiger.xml")
# read xml
pictures <- list(a = readPicture("GNU.xml"), b = readPicture("tiger.xml"))

# custom function for x axis label.
picture_axis <- function (pics, icon.size = unit(1, "lines"), ...) {
  structure(
      function(label, x = 0.5, y = 0.5, ...) {
         absoluteGrob(
           do.call("gList", mapply(symbolsGrob, pics[label], x, y, SIMPLIFY = FALSE)),
           height = icon.size)
    }
)}

p <- qplot(factor(c("a", "b")), 1:2) + opts( axis.text.x = picture_axis(pictures) )

pictures <- list(a = FragIcon2, b = Hex, c = FragIcon4, d = Hex)
dataf <- data.frame(fac = factor(c("a", "b", "a", "b")),
                    y   = 1:4,
                    group = c(1, 1, 2, 2) 
                    )
p <- qplot(fac, y, group = group, data = dataf, facets = .~ group) + 
scale_x_discrete(labels = names(pictures)[1:2], 
                 breaks = levels(dataf$fac)) +
opts( axis.text.x = picture_axis(pictures, icon.size = unit(3, "lines")) ) 
print(p)


## My graph
pictures  <- list(FragIcon1 = FragIcon1,
                  FragIcon2 = FragIcon2,
                  FragIcon3 = FragIcon3,
                  FragIcon4 = FragIcon4
                  )
picture_text_axis <- function (pics, icon.size = unit(1, "lines"), labels, lab.gp = gpar(), ...) {
  structure(
      function(label, x = 0.5, y = 0.5, ...) {
         absoluteGrob(
           do.call("gList", mapply(symbolsGrob, pics[label], x, y, SIMPLIFY = FALSE)),
           height = icon.size)
  ##          textGrob(labels, x, y, gp = lab.gp)  # this replaces pictures above instead of overlaying.
    }
)}

FxP.icons <- FxP.plot + 
scale_x_discrete(labels = names(pictures), 
                 breaks = levels(plot.means$Frag)) +
opts(axis.ticks.margin = unit(0.2, "lines"),
     ##      axis.text.x = picture_axis(pictures, icon.size = unit(1.4, "lines")) 
     axis.text.x = picture_text_axis(pictures, icon.size = unit(1.4, "lines"), labels = 1:4, lab.gp = gpar(col = c("white", "black", "black", "black")) ) +
)
## Overlay text on top of icons?  Or add a 'legend' for symbols?

print(FxP.icons)

################################################################


##==============================================================
## baptiste answer #1: custom grob axis labels
## http://stackoverflow.com/questions/8905101/how-can-i-use-a-graphic-imported-with-grimport-as-axis-tick-labels-in-ggplot2-u

library(grid)
library(ggplot2)

## convert the labels to some parameter to be used in the custom grob
## here simply an index that will be interpreted as color
mapping <- function(x, ...){
  seq_along(x) ## but could be something fancier ...
}

library(grImport)

hourglass <- new("Picture",
paths= list(new("PictureFill",
x=c(0, 1, 0, 1),
y=c(0, 0, 1, 1))),
summary= new("PictureSummary",
numPaths=1,
xscale=c(0, 1),
yscale=c(0, 1)))

## bare bones edit of theme_text()
my_axis <- function () 
{
    structure(function(label, x = 0.5, y = 0.5, default.units = "npc", ...) {
      cols <-mapping(label)

      symbolsGrob(hourglass, x, 0*x + unit(0.5, "npc"),
                  use.gc=FALSE,size=unit(5,"mm"), gp=gpar(fill=cols))

    }, class = "theme", type = "custom", call = match.call())
}

icon_axis1 <- function (picture, ...) 
{
  structure(function(label, x = 0.5, y = 0.5, default.units = "npc", ...) 
            {
              symbolsGrob(picture, x, 0*x + unit(0.5, "npc"),
                          size=unit(1,"lines"), ... )
            }, class = "theme", type = "custom", call = match.call()
            )
}

## not easy to auto-draw multiple pictures ...
icon_axis <- function (pictureList, ...) 
{
  if (class(pictureList) != "list") pictureList <- list(pictureList)

  structure(function(label, x = 0.5, y = 0.5, default.units = "npc", ...) 
            {
              cols <- seq_along(label)
              pictureIndex <- rep( 1:length(pictureList), length.out = length(x) )
              for (i in 1:length(x))
              {
                icon <- pictureList[[ pictureIndex[i] ]]
                xi <- x[i]
                symbolsGrob(icon, xi, 0*xi + unit(0.5, "npc"),
                            size=unit(2,"lines"), ... )
              }
            }, class = "theme", type = "custom", call = match.call()
            )
}


p <- qplot(1:12, rnorm(12)) +
opts(axis.text.x = icon_axis1(FragIcon2), 
     axis.ticks.margin = unit(0.3, "cm"),
     axis.title.x = theme_text(size = 12, vjust=-0.2)
     )

print(p)

## works well with multiple panels, though!
print(FxP.plot +
      opts(axis.text.x = icon_axis1(FragIcon2), 
           axis.ticks.margin = unit(0.3, "cm"),
           axis.title.x = theme_text(size = 12, vjust=-0.2)
           )
)




##==============================================================
## baptiste answer #2: edit axis grob directly
## http://stackoverflow.com/questions/8905101/how-can-i-use-a-graphic-imported-with-grimport-as-axis-tick-labels-in-ggplot2-u
library(grid)
library(ggplot2)

p <- qplot(1:12, rnorm(12)) +
opts(axis.ticks.margin = unit(0.3, "lines"),
     axis.title.x = theme_text(size = 12, vjust=-0.2)
     )

grid.newpage()
g <- ggplotGrob(p)
grid.draw(g)
g0 <- getGrob(g, gPath("axis.text.x"), grep=TRUE)
grid.set(gPath("axis.text.x"),
         pointsGrob(x = g0$x, y=0*g0$x + unit(0.5,"npc"),
                    pch=19, gp=gpar(col=seq_along(g0$x)),
                    name = g0$name), grep = TRUE)

## multiple pictures?



##==============================================================
## Hybrid solution: extract x-axis tick locations, then use a loop to draw each picture
## multiple panels?
iconNames <- c("FragIcon1", "FragIcon2", "FragIcon3", "FragIcon4")  # g0$label

FxP.icons <- FxP.plot + 
scale_x_discrete(labels = iconNames, 
                 breaks = levels(plot.means$Frag)) +
opts(axis.ticks.margin = unit(0.3, "cm"),
     axis.title.x = theme_text(size = 12, vjust=-0.2)
     )

grid.newpage()
g <- ggplotGrob(FxP.icons)
grid.draw(g)
axis.grobs <- grep("axis.text.x", grid.ls(grobs = TRUE, viewports = FALSE, print = FALSE)$name, value=TRUE)
g0 <- getGrob(g, gPath("axis.text.x"), grep=TRUE)
gx    <- g0$x                          # extract x-axis tick locations :D
icons <- g0$label                      # extract labels to use as picture names :P

# + multiple pictures?
for (axis.text in axis.grobs)
{ ## remove existing x-axis labels (From 1 axis at a time )
  g0 <- getGrob(g, gPath(axis.text))
  grid.set(gPath(axis.text),
           textGrob("", x = g0$x, y=0*g0$x + unit(0.5,"npc"),
                    name = axis.text))
}
## Add Frag Icons: still Hack, but more accurate.
iconSize   <- unit(1.3, "lines")
iconMargin <- unit(0.35, "lines") 
vports <- grep("axis_h", grid.ls(grobs = FALSE, viewports = TRUE, print = FALSE)$name, value=TRUE)
for (vp in unique(vports)) 
{
  downViewport(vp)
  iconHt    <- unit(1, "npc") - (iconSize*0.5) - iconMargin
  for (i in 1:length(icons))
  {
    icon <- icons[i]
    grid.symbols(get(icon), x = gx[i], y = iconHt, size = iconSize)
  }
  upViewport()
}





##==============================================================
## tikzDevice?
library(tikzDevice)
options(tikzLatexPackages = 
        c(getOption('tikzLatexPackages'),"\\usepackage{graphicx}\n")) 

d = data.frame(x=1:10, y=1:10, f=factor(sample(letters[1:2], 10, repl=TRUE))) 

p <- qplot(x,y,data=d) + theme_bw() + 
  opts(plot.margin = unit(c(1, 1, 5, 1), "lines"), 
       axis.text.x = theme_text(size = 12 * 
        0.8, lineheight = 0.9, vjust = 10)) + 
  scale_x_continuous(breaks = c(2, 8), labels=c("\\phe{15}", "\\leu{15}")) 

tikz("annotation.tex",standAlone=T,width=4,height=4) 
print(p) 
dev.off()                              # close tikz Device



FxP.icons <- FxP.plot + 
scale_x_discrete(labels = c("\\includegraphics[trim=5 5 0 0, clip]{Frag-Black-icons}", 
                            "", "", ""), 
                 breaks = levels(plot.means$Frag)) +
opts(axis.ticks.margin = unit(0.3, "cm"),
     axis.title.x = theme_text(size = 12, vjust=-0.2)
     )

tikz("./save/icons.tex", standAlone=T, width=4, height=4) 
print(FxP.icons)
dev.off()
## produces LaTeX input file, which still needs to be processed?
