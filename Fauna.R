##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Analyses of Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2011-08-22
###===============================================
### Species identified to morphospecies (usually family-level)
### Counts / sample converted to # / g dwt of moss substrate (using 'Patch.dwt' column)
##################################################
## INITIALISE
##################################################
if (FALSE) {  # do not run automatically
  ## Working Directory: see lib/init.R below
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (Mac).
  getwd()  # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')

## Load packages
library(vegan)	# load external package
library(ggplot2)	# load external package


##################################################
## CHECK & PROCESS DATA
##################################################
SECC.fauna.full <- SECC.fauna   # save a copy, just in case
## Filter by Time pt.  
## I do this here to maintain correspondence in row numbers between the various data frames.  
## I'm not really using the other data yet anyway.
SECC.fauna <- SECC.fauna[which(SECC.fauna$Time == 4), ]

## Aggregate count data by morphospecies with "confidence" 
## (i.e. lumping things together that are probably the same)
SECC.fauna.sp <- data.frame(SampleID = SECC.fauna[['SampleID']])

Taxa.groups <- rev( unique(SECC.fauna.meta$sp_alias) )
SECC.fauna.sp <- within(SECC.fauna.sp, {
    for (taxa in Taxa.groups) {
      taxa.sp <- SECC.fauna.meta$ID[which(SECC.fauna.meta$sp_alias == taxa)]
      assign( taxa, 
             if (length(taxa.sp) > 1) 
               apply(SECC.fauna[, taxa.sp], 1, sum) 
             else 
               SECC.fauna[, taxa.sp] 
      )
    }
    rm(taxa, taxa.sp)
})
rownames(SECC.fauna.sp) <- SECC.fauna.sp$SampleID
  
## Data matrix only: more convenient name, minus ID column (rownames still contain SampleIDs) :P
Fauna <- SECC.fauna.sp[, -1]  



##################################################
## DATA EXPLORATION
##################################################

## check variation & transformations
boxplot( Fauna            , main = "raw data")
boxplot( sqrt( Fauna )    , main = "sqrt")
boxplot( Fauna ^0.25      , main = "4th rt")
boxplot( log( Fauna +1 )  , main = "log_e(X +1)")
boxplot( log( Fauna +1 , base = 2)  , main = "log_2(X +1)")
boxplot( decostand( Fauna, method = 'log' )         , main = "log (decostand)") # ?? wtf?
boxplot( decostand( Fauna, method = 'normalize' )   , main = "normalized") 


##################################################
## NMDS
##################################################
## Transformation
Fauna.trans <- log( Fauna +1 ) 

## Distance/Similarity matrix (included in MDS call)
## Compute MDS
Fauna.mds  <- metaMDS( Fauna.trans, distance = "bray", autotransform = FALSE, k = 2)
Fauna.mds3 <- metaMDS( Fauna.trans, distance = "bray", autotransform = FALSE, k = 3)

## Check MDS quality (Stress)
stress.label <- paste( "Stress = ", round( Fauna.mds$stress, 1), "%", sep="" )
stressplot(Fauna.mds, main=stress.label)

stress.label <- paste( "Stress = ", round( Fauna.mds3$stress, 1), "%", sep="" )
stressplot(Fauna.mds3, main=stress.label)

## Look at the pretty graph
Fauna.plot <- ordiplot( Fauna.mds, type="text" )







##################################################
### PUBLICATION GRAPHS
##################################################
## Prepare data to plot
stress.label <- paste( "Stress = ", round( Fauna.mds$stress, 1)/100 )	# Make text label of Stress value, rounded to 1 decimal place, and converted from % to a decimal number.
MDS.pts <- as.data.frame(Fauna.mds$points)
colnames(MDS.pts) <- c('x', 'y')
MDS.pts <- cbind(MDS.pts, Time = SECC.fauna$Time, Chamber = SECC.fauna$Chamber, 
                 Frag = SECC.fauna$Frag, Pos = SECC.fauna$Pos
)

## plot symbol map
## Inner = full
## Outer = open
## Chamber = Red
## Ambient = Black
## Pos: pch (shape); Chamber: colour; Frag: shape?
Position.map <- plotMap( factor = "Position" )
Position.map <- Position.map[ Position.map$label %in% levels(SECC.fauna$Pos), ]
Chamber.map  <- plotMap( factor = "Chamber"  )
Chamber.map  <- Chamber.map[ Chamber.map$label %in% levels(SECC.fauna$Chamber), ]

Position.label <- attr(SECC, "labels")[["Pos"]] # "Patch\nPosition"
Chamber.label  <- "Chamber"         # attr(SECC, "labels")[["Chamber"]]

## Build plot in ggplot
Fauna.plot <- qplot(x, y, data = MDS.pts, group = Chamber,
                    shape = Pos, colour = Chamber, size = I(3), lwd = I(1.5),
                    xlab = NULL, ylab = NULL,
                    )   # facets=Frag~Time
Fauna.plot <- Fauna.plot + scale_shape_manual(name = Position.label,
                                              values = Position.map$pch, 
                                              breaks = Position.map$label,
                                              labels = c("Wet", "Dry")
                                              )
Fauna.plot <- Fauna.plot + scale_colour_manual(name = Chamber.label,
                                               values = Chamber.map$col, 
                                               breaks = Chamber.map$label,
                                               labels = c("Ambient", "Chamber")
                                               )
## options: see theme_get() for available theme options.
Fauna.plot <- Fauna.plot + theme_bw() + coord_equal()
Fauna.plot <- Fauna.plot + opts(axis.ticks = theme_blank(), 
                                axis.text.x = theme_blank(), 
                                axis.text.y = theme_blank(), 
                                legend.key = theme_rect(colour = NA),
                                legend.position = "right",
                                legend.direction = "vertical"
                                )
print(Fauna.plot)
mtext( stress.label , side=3, line=0, adj=0 )

Fauna.plot.frag <- Fauna.plot + coord_cartesian() + facet_grid(facets = Frag~.)
print(Fauna.plot.frag)

##################################################
### CLEAN-UP / HOUSEKEEPING
##################################################
