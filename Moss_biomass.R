################################################################
## Schefferville Experiment on Climate Change (SEC-C)
## Exploration of Moss Biomass, for use in converting
## Moss Growth to Moss Productivity
##   Moss weight data collected by Heather McIntosh: 
##   Volunteer & work study in the Gonzalez lab.  Thanks Heather!
## Jonathan Whiteley		R v2.12		2011-12-12
################################################################
## INITIALIZE
################################################################
if (FALSE) {  # do not run automatically
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("./")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')

library(MASS)                          # truehist
library(lattice)                       # densityplot
library(ggplot2)
theme_set(theme_bw())                  # change global ggplot2 theme
library(mgcv)                          # gam()
library(nlme)                          # Mixed effects
library(car)                           # Anova()

## load raw biomass data
Moss.biomass <- read.csv("./data/Moss_biomass.csv", na.strings=c("NA", "", "."),
						 stringsAsFactors=FALSE)


##==============================================================
## Check Data
##==============================================================
if (F) {
  head(Moss.biomass)
  str(Moss.biomass)
  SECCstr(Moss.biomass)
}

##==============================================================
## Process Data
##==============================================================
DropRows <- grep("^[^1-8]", Moss.biomass$Sample) # illegal Sample ID values ("MISSING...")
Moss.biomass <- Moss.biomass[-DropRows, ]
Moss.biomass <- strip_empty_dims(Moss.biomass, dim=1, cols=c(3, 4, 5, 6)) # cols 1--2 default to value of previous row if empty in raw values.
Moss.biomass <- within( Moss.biomass, 
					   {
						 dwt.mg <- gsub("<0\\.1", "0.05", dwt.mg) # values of "<0.1" ; [^0-9.]([0-9.]+)
						 dwt.mg <- as.numeric(dwt.mg)
						 Width  <- factor(Width,  levels=c("N", "W"), labels=c("<1cm", ">=1cm"))
						 Colour <- factor(Colour, levels=c("A", "B", "D"), labels=c("Green", "Brown", "Dead"))
						 ## extract treatment columns
						 Block   <- factor(substr(Sample, 1, 1))
						 Time    <- factor(substr(Sample, 2, 2))
						 Chamber <- factor(substr(Sample, 3, 3))
						 Frag    <- factor(substr(Sample, 5, 5))
						 Pos     <- factor(substr(Sample, 7, 7))
					   })
Moss.biomass <- Moss.biomass[, c("Sample", "Block", "Time", "Chamber", "Frag", "Pos", 
								 "Total.length", "Segment", "Width", "Colour", "dwt.mg")]
##  Moss.biomass <- checkSECCdata(Moss.biomass) # will fail due to duplicates



################################################################
## CALCULATIONS
################################################################
## + Time 2,4 ; Block 3, 5, 6, 7, 8 ; Chamber A, C ; Frag 1, 4 ; Position I, O
## + Some preliminary data from: Time 4, Chamber B ; Block 1, 2, 3, 6, 7, 8 ; Frag 1, 2, 3, 4 ; Position I, O, S
Moss.biomass <- within( Moss.biomass, 
					   {
					   })



################################################################
## EXPLORATION
################################################################
summary(Moss.biomass)

pairplot(Moss.biomass[, names(Moss.biomass)!="Sample"])

## Boxplots
op <- par(mfrow=c(3, 3))
boxplot(dwt.mg ~ Block   , data = Moss.biomass , xlab = "Block")
boxplot(dwt.mg ~ Time    , data = Moss.biomass , xlab = "Time")
boxplot(dwt.mg ~ Chamber , data = Moss.biomass , xlab = "Chamber")
boxplot(dwt.mg ~ Frag    , data = Moss.biomass , xlab = "Frag")
boxplot(dwt.mg ~ Pos     , data = Moss.biomass , xlab = "Pos")
boxplot(dwt.mg ~ Width   , data = Moss.biomass , xlab = "Width")
boxplot(dwt.mg ~ Colour  , data = Moss.biomass , xlab = "Colour")
boxplot(dwt.mg ~ Segment , data = Moss.biomass , xlab = "Segment")
boxplot(Segment ~ Colour , data = Moss.biomass , xlab = "Colour")
boxplot(Segment ~ Width  , data = Moss.biomass , xlab = "Width")
ar(op)

Point.smooth <- list(geom_point(), stat_smooth(method="gam"))
Biomass.plot <- ggplot(Moss.biomass, aes(x=Segment, y=dwt.mg)) + Point.smooth + facet_grid(facets = Time*Width ~ Block)
print(Biomass.plot)


################################################################
## ANALYSIS
################################################################
## Models with many variables take a while to fit ...

## Nested ANOVA?
biomass.aov <- aov(dwt.mg ~ Block * Time * Chamber * Frag * Pos * Width * Colour * Segment + Error(Time/Block/Chamber/Frag), data = Moss.biomass)

biomass.lm  <- lm( dwt.mg ~ Block * Time * Chamber * Frag * Pos * Width * Colour * Segment, data = Moss.biomass)

## Mixed model, mostly for variance decomposition: where is most of the variance?
lmc <- lmeControl(niterEM = 500, msMaxIter = 100, opt="optim")
biomass.mm  <- gls(dwt.mg ~ Segment * Width * Colour, weights = varIdent(form = ~ 1 | Block * Time * Chamber * Frag * Pos), 
				   data = Moss.biomass, method="REML", control = lmc)

biomass.gam <- gam(dwt.mg ~ Block * Time * Chamber * Frag * Pos + Width + Colour + s(Segment), data = Moss.biomass)
plot(biomass.gam)

AIC(biomass.lm, biomass.mm, biomass.gam)

summary(biomass.aov)
summary(biomass.lm)
summary(biomass.mm)
summary(biomass.gam)

## Anova(biomass.aov)
Anova(biomass.lm)
Anova(biomass.mm)

##==============================================================
## Categorical Analyses
##==============================================================
Segments <- 2:3

Width.contingency <- with(subset(Moss.biomass, Segment %in% Segments), table(Width, Block, Time, Chamber, Frag, Pos))
Width.Block <- with(subset(Moss.biomass, Segment %in% Segments), table(Width, Block))
Width.BC    <- with(subset(Moss.biomass, Segment %in% Segments), table(Width, Block, Chamber))

chisq.test(Width.Block)
chisq.test(Width.BC)
chisq.test(Width.contingency)

## Where are the Narrow shoots?
Narrow.segments <- Moss.biomass[Moss.biomass$Width=="<1cm" & 
								Moss.biomass$Segment %in% Segments & 
								Moss.biomass$Chamber != "B", 
								]
Wide.segments   <- Moss.biomass[Moss.biomass$Width==">=1cm" & 
								Moss.biomass$Segment %in% Segments & 
								Moss.biomass$Chamber != "B", 
								]
op <- par(mfcol=c(2,2))
hist(Moss.biomass[, "dwt.mg"], xlim=c(0, 20))
hist(Narrow.segments[, "dwt.mg"], xlim=c(0, 20))
hist(Moss.biomass[Moss.biomass$Width==">=1cm", "dwt.mg"], xlim=c(0, 20))
hist(Wide.segments[, "dwt.mg"], xlim=c(0, 20))
par(op)



################################################################
## RESULTS / CONCLUSIONS
################################################################
## segments 2-4 seem roughly similar (gam, exploration graphs)
## The shortest stem is 4 cm, so this would capture every shoot
## Not surprisingly, Width & Colour make a big difference, although:
## - Colour is pretty closely related to segment: as expected, most Green segments are near the top (1-9), most Dead ones near the bottom (segment 3+)
## - How are these variables correlated with Block, etc.?
##   - Most segments 2--4 are >=1cm
##   - Block 2, Chamber B was mostly <1 cm
##     - as were Blocks 7 & 8, Chamber B

## Grand Means: not as useful as estimates for each patch (or group of patches), but better than nothing :/
Moss.segments <- subset(Moss.biomass, subset = Segment %in% Segments & Chamber!="B" & Pos %in% c("0", "1") )
Segment.mean  <- mean(Moss.segments$dwt.mg, na.rm=TRUE)
Wide.mean     <- mean(Wide.segments$dwt.mg, na.rm=TRUE) 
Narrow.mean   <- mean(Narrow.segments$dwt.mg, na.rm=TRUE) 

## Moss growth data does contain notes indicating which shoots had "skinny" tips (<<1cm width)
## - I should ensure these use a value for "skinny" shoots
## - This may also explain why liner growth is roughly the same over summer vs. winter
##   - Summer values may reflect skinny extension, and overall lower biomass production
##   - Most productivity may be occuring during the winter!

### ALGORITHM to provide biomass for 1 cm of linear extension for a given patch:
## Basic value: average of mass of Segments 2--3
## Order of preference:
## - Values from that patch (if present in this data)
## - Average of values from 4 nearest patches (need spatial data for this)
##   - Matched by experimental treatments?
##     Chamber, Pos, Frag? (Time?) ; Block is implicit in the nearest neighbor criteria, and therefore unnecessary
## - if "skinny" in notes:
##   - only use values for Narrow shoots (<1cm)
##   - use half the average value (value / 2)?  Even this might be generous!
##   - what I really need is an approximate weight of the main shoot, without lateral branches for these :(
