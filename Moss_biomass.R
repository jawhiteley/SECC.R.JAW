################################################################
## Schefferville Experiment on Climate Change (SEC-C)
## Exploration of Moss Biomass, 
## and calculations for converting Moss Growth to Productivity
##   Moss weight data collected by Heather McIntosh: 
##   Volunteer & work study in the Gonzalez lab. Thanks Heather!
## Jonathan Whiteley		R v2.12		2011-12-14
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
cat("Processing moss biomass data\n")

DropRows <- grep("^[^1-8]", Moss.biomass$Sample) # illegal Sample ID values ("MISSING...")
Moss.biomass <- Moss.biomass[-DropRows, ]
Moss.biomass <- strip_empty_dims(Moss.biomass, dim=1, cols=c(3, 4, 5, 6)) # cols 1--2 default to value of previous row if empty in raw values.
Moss.biomass <- within( Moss.biomass, 
					   {
						 dwt.mg <- gsub("<0\\.1", "0.05", dwt.mg) # values of "<0.1" ; [^0-9.]([0-9.]+)
						 dwt.mg <- as.numeric(dwt.mg)
						 Width  <- factor(Width,  levels=c("N", "W"), labels=c("<1cm", ">=1cm"))
						 Colour <- factor(Colour, levels=c("A", "B", "D"), labels=c("Green", "Brown", "Dead"))
						 ## manually fix ID codes
						 Sample <- gsub("\\.0", ".O", Sample)
						 Sample <- gsub("\\.1", ".I", Sample)
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
## levels(Moss.biomass$Pos)[levels(Moss.biomass$Pos)=="0"] <- "O"
## levels(Moss.biomass$Pos)[levels(Moss.biomass$Pos)=="1"] <- "I"


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
cat("Exploring moss biomass data\n")

summary(Moss.biomass)

pairplot(Moss.biomass[, names(Moss.biomass)!="Sample"])

## Boxplots
op <- par(mfrow=c(3, 3))
boxplot(dwt.mg ~ Block   , data = Moss.biomass , xlab = "Block"   , ylab = "mg / cm")
boxplot(dwt.mg ~ Time    , data = Moss.biomass , xlab = "Time"    , ylab = "mg / cm")
boxplot(dwt.mg ~ Chamber , data = Moss.biomass , xlab = "Chamber" , ylab = "mg / cm")
boxplot(dwt.mg ~ Frag    , data = Moss.biomass , xlab = "Frag"    , ylab = "mg / cm")
boxplot(dwt.mg ~ Pos     , data = Moss.biomass , xlab = "Pos"     , ylab = "mg / cm")
boxplot(dwt.mg ~ Width   , data = Moss.biomass , xlab = "Width"   , ylab = "mg / cm")
boxplot(dwt.mg ~ Colour  , data = Moss.biomass , xlab = "Colour"  , ylab = "mg / cm")
boxplot(Segment ~ Colour , data = Moss.biomass , xlab = "Colour"  , ylab = "Segment")
boxplot(Segment ~ Width  , data = Moss.biomass , xlab = "Width"   , ylab = "Segment")
par(op)
boxplot(dwt.mg ~ Segment , data = Moss.biomass , xlab = "Segment" , ylab = "mg / cm")

Point.smooth <- list(geom_point(), stat_smooth(method="gam"))
Biomass.plot <- ggplot(Moss.biomass, aes(x=Segment, y=dwt.mg)) + Point.smooth + facet_grid(facets = Time*Width ~ Block)
print(Biomass.plot)



################################################################
## ANALYSIS
################################################################
Segments <- 2:3
Moss.segments <- subset(Moss.biomass, subset = Segment %in% Segments & Pos %in% c("O", "I") )  #  & Chamber!="B"
Moss.segmeans <-  aggregate(Moss.segments$dwt.mg, 
							by=list(SampleID = Moss.segments$Sample,
									Block    = Moss.segments$Block,
									Time     = Moss.segments$Time,
									Chamber  = Moss.segments$Chamber,
									Frag     = Moss.segments$Frag,
									Pos      = Moss.segments$Pos), 
							FUN=mean)
names(Moss.segmeans)[names(Moss.segmeans)=="x"] <- "dwt"
Moss.segmeans <- checkSECCdata(Moss.segmeans)

## Models with many variables take a while to fit ...
cat("Fitting models to moss biomass data.  This may take a few minutes...\n")
## Nested ANOVA?
biomass.aov <- aov(dwt.mg ~ Block * Time * Chamber * Frag * Pos * Width * Colour * Segment + Error(Time/Block/Chamber/Frag), data = Moss.biomass)

biomass.lm  <- lm( dwt.mg ~ Block * Time * Chamber * Frag * Pos * Width * Colour * Segment, data = Moss.biomass)
seg.xpt.lm  <- lm( dwt.mg ~ Block * Time * Chamber * Frag * Pos, data = Moss.segments)

## Mixed model, mostly for variance decomposition: where is most of the variance?
lmc <- lmeControl(niterEM = 500, msMaxIter = 100, opt="optim")
biomass.mm  <- gls(dwt.mg ~ Segment * Width * Colour, weights = varIdent(form = ~ 1 | Block * Time * Chamber * Frag * Pos), 
				   data = Moss.biomass, method="REML", control = lmc)

biomass.gam <- gam(dwt.mg ~ Block * Time * Chamber * Frag * Pos + Width + Colour + s(Segment), data = Moss.biomass)
plot(biomass.gam)

AIC(biomass.lm, biomass.mm, biomass.gam)

summary(biomass.aov)
summary(biomass.lm)
summary(seg.xpt.lm)
summary(biomass.mm)
summary(biomass.gam)

## Anova(biomass.aov)
Anova(biomass.lm)
anova(seg.xpt.lm)                      # Anova throws Error
## Anova(biomass.mm)

cat("Drawing interaction plots\n")
## Effects plots: SLOW & buggy
library(effects)
if (F) {
CxP.eff <- effect("Chamber:Pos" , seg.xpt.lm) # ***
CxP.eff <- effect("Chamber:Pos" , biomass.lm) # ***
plot(effect("Chamber:Pos"                     , biomass.lm)) # ***
plot(effect("Time:Frag"                       , biomass.lm)) # ***
plot(effect("Time:Frag:Pos"                   , biomass.lm)) # ***
plot(effect("Block:Time:Chamber:Pos"          , biomass.lm)) # ***
plot(effect("Block:Chamber:Frag:Segment"      , biomass.lm)) # *
plot(effect("Block:Pos:Colour:Segment"        , biomass.lm))
plot(effect("Block:Time:Frag:Segment"         , biomass.lm))
plot(effect("Block:Time:Chamber:Frag:Segment" , biomass.lm))
plot(effect("Frag:Width"                      , biomass.lm)) # *
plot(effect("Time:Frag:Pos:Width"             , biomass.lm)) # *
plot(effect("Time:Chamber:Frag:Colour"        , biomass.lm)) # *
plot(effect("Time:Chamber:Pos:Segment"        , biomass.lm))
plot(effect("Time:Frag:Pos:Segment"           , biomass.lm)) #
plot(effect("Block:Chamber:Colour"            , biomass.lm)) # *
plot(effect("Block:Chamber:Segment"           , biomass.lm))
}

## Interaction plots with ggplot
Pts.Chamber <- list(geom_point(shape=1),
			   stat_summary(aes(colour=Chamber), fun.data="mean_cl_normal", geom="line"), # colour="#888888" 
			   stat_summary(aes(colour=Chamber), fun.data="mean_cl_normal", geom="errorbar"), # colour="#888888" 
			   scale_colour_manual(name="Chamber Treatment", 
								   values=c("A"="#990000", "B"="#000099", "C"="#000000"))
			   )
Pts.CI  <- list(geom_point(shape=1),
			   stat_summary(fun.data="mean_cl_normal", geom="line", colour="#888888"),
			   stat_summary(fun.data="mean_cl_normal", geom="errorbar", colour="#888888"),
			   scale_colour_manual(name="Chamber Treatment", 
								   values=c("A"="#990000", "B"="#000099", "C"="#000000"))
			   )
Pos.plot <- ggplot(Moss.segments, aes(x=Pos, y=dwt.mg))
CP.plot  <- ggplot(Moss.segments, aes(x=Pos, y=dwt.mg, colour=Chamber, group=Chamber))
CF.plot  <- ggplot(Moss.segments, aes(x=Frag, y=dwt.mg, colour=Chamber, group=Chamber))
CxP.plot  <- CP.plot + Pts.Chamber + facet_grid(facets= ~ Chamber)
BTCP.plot <- CP.plot + Pts.Chamber + jaw.ggplot() + facet_grid(facets= Block ~ Time)
BTCF.plot <- CF.plot + Pts.Chamber + jaw.ggplot() + facet_grid(facets= Block ~ Time)
TPCF.plot <- CF.plot + Pts.Chamber + jaw.ggplot() + facet_grid(facets= Pos ~ Time)
Frag.plot <- CF.plot + aes(group=NULL) + Pts.CI + jaw.ggplot()
TFP.plot  <- ggplot(Moss.segments, aes(x=Frag, y=dwt.mg, colour=Pos, group=Pos)) + 
			 stat_summary(aes(colour=Pos), fun.data="mean_cl_normal", geom="line") +
			 stat_summary(aes(colour=Pos), fun.data="mean_cl_normal", geom="errorbar") + 
			 geom_point(shape=21, fill="white") + 
			 scale_colour_manual(name="Position", 
								 values=c("O"="#990000", "I"="#000000"), 
								 breaks=c("O", "I")) +
			 jaw.ggplot() + facet_grid(facets= ~ Time)


print(CxP.plot)
print(BTCP.plot)
print(BTCF.plot)
print(TPCF.plot)
print(Frag.plot)
print(TFP.plot)

cat("Follow-up analyses\n")
##==============================================================
## Categorical Analyses
##==============================================================
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

##==============================================================
## Spatial Analysis
##==============================================================
## Bubble Plot of raw values / residuals in space?
Moss.segmeans <- merge(Moss.segmeans[, c("SampleID", "dwt")], SECC.xy, by="SampleID", all.x = TRUE, all.y = FALSE)
Moss.bubble <- ggplot(Moss.segmeans, aes(x=xE, y=yN)) + geom_point(shape=1, aes(size=dwt))
print(Moss.bubble)

## Variogram of segment biomass by distance
Moss.vary <- gls(dwt ~ xE + yN, data = Moss.segmeans)
Moss.variogram <- Variogram(Moss.vary, form = ~xE + yN)
plot(Moss.variogram)
## No spatial autocorrelation evident


################################################################
## RESULTS / CONCLUSIONS
################################################################
## segments 2-4 seem roughly similar (gam, exploration graphs)
## The shortest stem is 4 cm, so this would capture every shoot
## Not surprisingly, Width & Colour make a big difference, although:
## - Colour is pretty closely related to segment: 
##   as expected, most Green segments are near the top (1-9), 
##   most Dead ones near the bottom (segment 3+)
## - How are these variables correlated with Block, etc.?
##   - Most segments 2--4 are >=1 cm wide
##   - Block 2, Chamber B was mostly <1 cm wide
##     - as were Blocks 7 & 8, Chamber B
## Width is perhaps the single variable with the most effect on average biomass for a given 1cm segment

## There is A LOT of variation, even within a single patch.  e.g.
Moss.segments[Moss.segments$Frag==3, ]
summary(Moss.segments[Moss.segments$Frag==3, ])
summary(Moss.segments)
summary(Moss.biomass[Moss.biomass$Pos=="S", ]) # all <1cm (except for 1)

## Grand Means: not as ideal as estimates for each patch (or group of patches), but better than nothing :/
Segment.mean  <- mean(Moss.segments$dwt.mg, na.rm=TRUE)
Wide.mean     <- mean(Wide.segments$dwt.mg, na.rm=TRUE) 
Narrow.mean   <- mean(Narrow.segments$dwt.mg, na.rm=TRUE) 

## Moss growth data does contain notes indicating which shoots had "skinny" tips (<<1cm width)
## - I should ensure these use an appropriate biomass value for "skinny" shoots
## - This may also explain why liner growth is roughly the same over summer vs. winter
##   - Summer values may reflect skinny extension, and overall lower biomass production
##   - Most productivity may be occuring during the winter!

##==============================================================
### ALGORITHM to provide biomass for 1 cm of linear extension for a given patch:
## Basic value: average of mass of Segments 2--3
## Order of preference:
## - Values from that patch (if present in this data)
## - Average of values from 4? nearest patches, up to max. 10m(?) (need spatial data for this)
##   - Matched by experimental treatments?
##     Chamber, Pos, Frag? (Time?) ; 
##     - Block is implicit in the nearest neighbor criteria, and therefore unnecessary
##     - Chamber & Pos are important in some Blocks...
##     - Not sure I really have enough to judge if there is a Frag effect
##       - When in doubt, prefer Continuous patches over Isolated (more crowded)?
## - if "skinny" in notes:
##   - only use values for Narrow shoots (<1cm); 'skinny' shoots are often even skinnier!
##   - use half the average value (value / 2)?  Even this might be generous!
##   - what I really need is an approximate weight of the main shoot, without lateral branches for these :(
cat("Calculating Conversion Values\n")

NumNeighbs <- 3                        # Number of Nearest Neighbour values to use to interpolate patch estimates.
Moss.mg <- within( subset(SECC.xy, Time %in% 3:4 & Chamber!="B" ), 
				  mg.cm <- NA )
SECC.dist <- dist(SECC.xy[, c("xE", "yN")])
Moss.dist <- dist(Moss.segmeans[, c("xE", "yN")]) # distance matrix of patches with biomass estimates
Dist.matrix <- as.matrix(SECC.dist)
ValuesDs    <- which( colnames(Dist.matrix) %in% Moss.segmeans$SampleID )

for (i in 1:nrow(Moss.mg))
{
  IDi <- Moss.mg[i, "SampleID"]
  iChamber <- Moss.mg[i, "Chamber"]
  iPos     <- if (iChamber=="A") "." else Moss.mg[i, "Pos"] # any Ambient patch will do
  if (!(iPos %in% c("I", "O"))) iPos <- "." # No useable data for 'other' patches
  ## currently match any closest patch for 'other' patches (not Inner or Outer)
  ## - Ideally, I would average an equal number of Inner AND Outer patches,
  ##   But hopefully this is 'good enough'
  ## try for an exact match?
  iHaveThat <- IDi %in% Moss.segmeans$SampleID
  if (iHaveThat) 
  {                                    # Extract patch values if exact match
	PatchValues  <- subset(Moss.segmeans, SampleID == IDi, select="dwt")
	PatchWeights <- rep(1, length(PatchValues)) 
  } else {                             # Extract other values to use for patch estimate
	## find nearest neighbours, with biomass values AND same Chamber & Position values
	Distances   <- Dist.matrix[IDi, ValuesDs]
	MatchDs     <- grep(sprintf("..%s-.\\.%s", iChamber, iPos), names(Distances))
	Distances   <- sort( Distances[MatchDs] ) #     Neighbours  <- order(Distances)
	NeiIDs      <- names(Distances)[1:NumNeighbs]
	## extract nearest neighbor values
	PatchValues  <- subset(Moss.segmeans, SampleID %in% NeiIDs, select="dwt")[, "dwt"]
	PatchWeights <- 1/Distances[1:NumNeighbs] 
  }
  ## Average of Patch Values (exact or neighbours)
  if (length(PatchValues) < 1)
  {
	PatchEstimate <- NA
  } else {
	PatchEstimate <- weighted.mean(PatchValues, w=PatchWeights, na.rm=TRUE)
  }
  Moss.mg[i, "mg.cm"] <- PatchEstimate
}
## Moss.mg$mg.cm[is.nan(Moss.mg$mg.cm)] <- NA

hist(Moss.mg$mg.cm)
summary(Moss.mg)

## Attributes (optional)
attr(Moss.mg, "SECC columns") <- "mg.cm"
attr(Moss.mg, "labels") <- "Moss biomass"
attr(Moss.mg, "units")  <- quote(mg %.% cm^-1)

## Save Conversion values
cat("Saving Conversion Data\n")
Filename <- sprintf( "%sMoss_mg.%%s", SaveDir.obj() )
save(Moss.mg, file = sprintf(Filename, "R") )
write.csv(Moss.mg, file = sprintf(Filename, "csv"), row.names = FALSE) 

