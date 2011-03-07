##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Template for basic analyses of experimental data
# Response Variable(s) @ time #s
# R v2.10.1		;	2010-05-04
##################################################
## INITIALISE
##################################################
rm(list=ls())	# clear memory
# setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis")	# Set Working Directory: replace with a path in quotes "".	# iMac@McGill
# setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis")	# JAW-MBP
library(car)		# load external package 'car', for recode()
library(lattice)	# mostly for xyplot
library(nlme)		# mixed-effect models
#=================================================
## DEFINE FUNCTIONS
#=================================================
source("SECC.functions.R")	#2010-04-14

##################################################
## LOAD DATA
##################################################
SECC.full <- read.csv("SECC.csv", na.string=".")	# specify missing values with na.string=""

## Define Labels
Y.label <- "Response Variable"	# response variable label for this script.
Chamber.label = "Chamber"
Frag.label = "Fragmentation Treatment"
pos.label = "patch position"
dataset.list <- c("SECC", "SECCr")
dataset.labels <- c( "patch scale", "Regional scale" )

##################################################
## PROCESS DATA: planned
##################################################
# str(SECC.full)
## Generate factor columns, re-order factor levels, etc.  Creating new columns does not work inside with().
SECC.full <- within( SECC.full, {	## instead of attach, this places focus within the indicated data frame / object, and returns object with changes applied (in reverse order...)
	# Generic Response Variable (used in remainder of script)
	Y <- as.numeric(Y)
	## Rename columns /convert to standard informative names for the rest of this script
	Block 	<- as.factor(Block)
	TimePt	<- TimePt
	Warming <- as.factor(Warming)
	Frag 	<- as.factor(Frag)
	pos 	<- as.factor(pos)
	# new factor column with meaningful values
	Chamber <- as.factor(Warming)
	levels(Chamber) <- list( 'Ambient'='A', 'Partial Chamber'='B', 'Full Chamber'='C' )	# rename and reorder factor levels to informative names.  Unspecified values converted to 'NA'.
	## rename and reorder factor levels the easy way - maintains empty values if empty factor specified, otherwise converts to 'NA'.  Requires package 'car'
	Chamber <- recode( Warming, 
		"'A'='Ambient'; 'B'='Partial Chamber'; 'C'='Full Chamber'; else=''", 
		as.factor.result=TRUE, 
		levels=c( "Ambient", "Partial Chamber", "Full Chamber" ) 
	)
	levels(Frag) <- list( 'Continuous'=1, 'Full Corridors'=2, 'Pseudo-Corridors'=3, 'Isolated'=4 )
	pos <- factor(pos, levels=c('1', 'S', 'W', 'E', 'N', '0'))	# safely reorder factor levels
	levels(pos) <- list( 'I'='1', 'S'='S', 'W'='W', 'E'='E', 'N'='N', 'O'='0' )	# rename some factor levels (omitted levels are dropped and replaced with empty strings)
	# new factor column with simplified values for patch position
	pos1 <- recode( pos, 
		"'I'='Inner'; 'O'='Outer'; else='other'", 
		as.factor.result=TRUE, 
		levels=c( "Inner", "other", "Outer" ) 
	)
})

## Filter data for analysis.  If no filter is desired, retain assignment to new object for consistent naming & reference for the rest of this script.
SECC.full <- SECC.full[SECC.full$TimePt==1, ]	# only one time point at a time (not repeated measures)
SECC <- SECC.full[SECC.full$SampleControl=="Sample",]	# only rows for samples, exclude controls

## Summarize data by means across different (or all) positions to prevent unbalanced effects?
SECCr <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, TimePt=TimePt, Chamber=Chamber, Frag=Frag), mean ) )	# for regional-level analyses (ignoring position)
SECC <- with( SECC, aggregate( cbind(Y) , by=list(Block=Block, TimePt=TimePt, Chamber=Chamber, Frag=Frag, pos1=pos1), mean ) )	# using cbind() on the response variables allows multiple columns to be summarized, and also preserves column names.

#=================================================
# PROCESS DATA: unplanned?
# Transformations, calculated columns, etc.  May depend on the results of exploratory graphs (& assumptions), below.
for ( dataset in dataset.list ){
	assign( dataset,
		within( get(dataset), {
			Y <- replace( Y, which( Y<0 ), 0 )	# negative values are problematic for log and other transformations :-(
			Y.ln  <- log( Y +1 )	# defaults to base e=exp(1).
			Y.log <- log( Y +1 , 10 )	# base 10. (stdev prop to mean).
			Y.sqrt <- sqrt( Y )	# useful for Poisson-distributed data (mean prop. to variance).
			Y.4rt <-Y^(0.25)	# fourth-root
			Y.use <- Y	# assign which column to work with for analyses****
		})	
	)
}
# Which response variable is being used (for labels)?****
Y.used <- "Y"

##################################################
## CHECK DATA
##################################################
head(SECC)		# have a peek at the first 6 rows & columns: is this what you expected?
str(SECC)		# check structure: are the appropriate variables factors, numeric, etc.?
summary(SECC)	# summary statistics
## Regional analyses
head(SECCr)		# have a peek at the first 6 rows & columns: is this what you expected?
str(SECCr)		# check structure: are the appropriate variables factors, numeric, etc.?
summary(SECCr)	# summary statistics


##################################################
## EXPLORE: PLOTS
##################################################
## make some meaningful plots of data to check for predicted (expected) patterns.
par( mfrow=c(2,2), cex=0.8)	# panel of figures: 2 rows & 2 columns
## Patch analyses
with( SECC, plot(Y.log ~ Chamber*Frag*pos1, 
	main="Patch") )	# fixed effects only, no nesting
plot.new()	# empty panel
plot.new()	# empty panel
par( mfrow=c(2,2), cex=0.8)	# panel of figures: 2 rows & 2 columns
## Is the response variable normally-distributed? Check skew in histograms, eyeball linearity on QQ-plots
for( i in 1:2 ){
	dataset <- dataset.list[i]
	with( get(dataset), {
		hist( Y )
		hist( Y.log )
		hist( Y.sqrt )
		hist( Y.4rt )
		qqnorm( Y, main=paste("untransformed: ", dataset.labels[i]) )
		qqline( Y, col="grey50" )	# draw line through 1st & 3rd quartiles.
		qqnorm( Y.log, main=paste("log10-transformed: ", dataset.labels[i]) )
		qqline( Y.log, col="grey50" )
		qqnorm( Y.sqrt, main=paste("sqrt-transformed: ", dataset.labels[i]) )
		qqline( Y.sqrt, col="grey50" )
		qqnorm( Y.4rt, main=paste("4th-root-transformed: ", dataset.labels[i]) )
		qqline( Y.4rt, col="grey50" )
	})
}


##################################################
## DEFINE MODEL FORMULA
##################################################
# Nested Fixed Effects, with error term for ANOVA using aov() 
Yp.model <- Y.use ~ Chamber*Frag*pos1 +Error(Block/Chamber/Frag)
# ignoring effect of position: 'regional' effects only
Yr.model <- Y.use ~ Chamber*Frag +Error(Block/Chamber/Frag)
# Mixed Effects model using lme() (package 'nlme')
Yp.fixed <- Y.use ~ Chamber*Frag*pos1
Yp.random <- ~1|Block/Chamber/Frag

##################################################
## ANALYSIS: design
##################################################
Yp.aov <- aov( Yp.model, data=SECC )	# simple ANOVA with nested fixed effects.  Nesting means an object of type "aovlist" is produced* (therefore regular TukeyHSD & other functions won't work )
Yr.aov <- aov( Yr.model, data=SECCr )	# regional effects only.
Yp.lme <- lme( fixed=Yp.fixed, random=Yp.random, data=SECC )

##################################################
## CHECK ASSUMPTIONS: analyse residuals, standard diagnostic plots
##################################################
## anova
# independence?
	# experimental design: random position of Frag within Chambers.
	# possibility of Block gradient (7&8 SW -> 1-6 E:N), which also corresponds roughly to order of samples.
## Patch analyses
	## trellis plots: any pattern across blocks, within frag & chambers?
	xyplot( Y.use ~ Block | Frag + Chamber, data=SECC, 
		pch=21, col="black", bg="grey", cex=0.8,
		main = dataset.labels[1]
	)
par( mfrow=c(2,2), cex=0.8)	# panel of figures: 2 rows & 2 columns
# homogeneity of variances?
with( SECC, plot(Y.use ~ Chamber*Frag*pos1) )	# fixed effects only, no nesting
plot.new()
# normal distribution?
with(SECC, qqnorm(resid(Yp.aov$Within), main="Residuals - patch" ) )	# are residuals normally distributed?
par( mfrow=c(1,1) )
hist(resid(Yp.aov$Within))	# plot residuals
# with(SECC, shapiro.test( Y.use ) )	# normality?

## REGIONAL analyses
	xyplot( Y.use ~ Block | Frag + Chamber, data=SECCr, 
		pch=21, col="black", bg="grey", cex=0.8,
		main = dataset.labels[2]
	)
par( mfrow=c(2,2), cex=0.8)	# panel of figures: 2 rows & 2 columns
# homogeneity of variances?
with( SECCr, plot(Y.use ~ Chamber*Frag) )	# fixed effects only, no nesting
# normal distribution?
with(SECCr, qqnorm(resid(Yp.aov$Within), main="Residuals - Regional" ) )	# are residuals normally distributed?
par( mfrow=c(1,1) )
hist(resid(Yr.aov$"Block:Chamber:Frag"))	# plot residuals
# with(SECCr, shapiro.test( Y.use ) )	# normality?


##################################################
## ANALYSIS: GET RESULTS
##################################################
## Patch analyses
# names(Yp.aov)
summary(Yp.aov)		# summary statistics
model.tables(Yp.aov, "means")	# effect sizes
# Interaction Plots
par(mfrow=c(2,2))	# panel of figures: 2 rows & 2 columns
with( SECC, interaction.plot( Frag, Chamber, Y.use, ylab=paste("mean of ", Y.used) ) )
with( SECC, interaction.plot( pos1, Chamber, Y.use, ylab=paste("mean of ", Y.used) ) )
with( SECC, interaction.plot( pos1, Frag, Y.use, ylab=paste("mean of ", Y.used) ) )

## Planned Multiple Comparisons using Least Significant Differences (LSD) -> comparison intervals for graphical display.
# Chamber x pos Interaction
lsd <- LSD( Yp.aov$Within, Yp.model, data=SECC, alpha=0.05, mode="pairwise" )	# compute LSDs based on a 5% error rate (alpha), 2-tailed.  Manual, due to unbalanced data (this is an estimate).
lsd.Cxp <- lsd["Chamber:pos1"]

## Regional analyses
# names(Yr.aov)
summary(Yr.aov)		# summary statistics
model.tables(Yr.aov, "means")	# effect sizes
# Interaction Plots
par(mfrow=c(1,1))	# panel of figures: 1 rows & 1 columns
with( SECCr, interaction.plot( Frag, Chamber, Y.use, ylab=paste("mean of ", Y.used) ) )

## Planned Multiple Comparisons using Least Significant Differences (LSD) -> comparison intervals for graphical display.
lsd.r <- LSD( Yr.aov$"Block:Chamber:Frag", Yr.model, data=SECCr, alpha=0.05, mode="pairwise" )	# compute LSDs based on a 5% error rate (alpha), 2-tailed.
lsd.r.FxC <- lsd.r["Chamber:Frag"]
lsd.r <- LSD( Yr.aov$"Block:Chamber", Yr.model, data=SECCr, alpha=0.05, mode="pairwise" )	# compute LSDs based on a 5% error rate (alpha), 2-tailed.
lsd.r.C <- lsd.r["Chamber"]

## Mixed Effects Model analysis Results
# names(Yp.lme)
# summary(Yp.lme)		# summary statistics

##################################################
## FINAL GRAPHICS
##################################################
Chamber.map <-	data.frame( label=levels(SECC$Chamber), col=c("#000000","#000099","#990000"), bg=c("#FFFFFF","#FFFFFF","#FFFFFF"), pch=c(21,23,18), lty=c(3,2,1) )
	# Ambient = black, open circles with dotted line ; 
	# Partial = blue, open diamonds with dashed line ; 
	# Full	  = red, solid diamond with solid line.

## Patch results
plot.means <- with( SECC, 
	aggregate( cbind( Y.use ), 
		list(pos=pos1, Chamber=Chamber), 
		mean 
	) 
)
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
	# using custom plotMeans function, with custom error bars (LSD)
	plot.error <- matrix( as.numeric(lsd.Cxp/2), nrow = length(levels(pos)), ncol = length(levels(Chamber)) )
	plotMeans( Y.use , pos , Chamber, 
		error.bars="custom", level=plot.error, cex=2, lwd=2,
		lty=Chamber.map$lty, pch=Chamber.map$pch,
		col=as.character(Chamber.map$col),
		bg=as.character(Chamber.map$bg),
		main="Means of patch values ± LSD (95% comparison intervals)",
		xlab=pos.label, ylab=Y.label
	)	# as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## REGIONAL results
plot.means <- with( SECCr, 
	aggregate( cbind( Y.use ), 
		list(Frag=Frag, Chamber=Chamber), 
		mean 
	) 
)
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
	# using custom plotMeans function, with custom error bars (LSD)
	plot.error <- matrix( as.numeric(lsd.r.FxC/2), nrow = length(levels(Frag)), ncol = length(levels(Chamber)) )
	plotMeans( Y.use , Frag , Chamber, 
		error.bars="custom", level=plot.error, cex=2, lwd=2,
		lty=Chamber.map$lty, pch=Chamber.map$pch,
		col=as.character(Chamber.map$col),
		bg=as.character(Chamber.map$bg),
		main="Means of Regional values ± LSD (95% comparison intervals)",
		xlab=Frag.label, ylab=Y.label
	)	# as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## Chamber Main Effects
plot.means <- with( SECCr, 
	aggregate( cbind( Y.use ), 
		list(Chamber=Chamber), 
		mean 
	) 
)
plot.error <- rep( as.numeric(lsd.r.C/2), length(levels(plot.means$Chamber)) )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
	# using custom plotMeans function, with custom error bars (LSD)

	plotMeans( Y.use , Chamber, 
		error.bars="custom", level=plot.error, cex=2, lwd=2,
		lty=Chamber.map$lty[3], pch=Chamber.map$pch,
		col=as.character(Chamber.map$col),
		bg=as.character(Chamber.map$bg),
		main="Means of Regional values ± LSD (95% comparison intervals)",
		xlab=Chamber.label, ylab=Y.label
	)	# as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})
