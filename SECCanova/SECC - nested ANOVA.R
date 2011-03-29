##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Standard Nested ANOVA of experimental data
### Response Variable(s)  @ time #s
### Jonathan Whiteley     R v2.12     2011-03-26
##################################################
##
##  This script is used to carry out standard
##  nested analysis of variance, according to
##  the structure of the SECC Experiment.
##  Calling scripts should configure settings
##  before source()'ing this file.
##
##################################################
## INITIALISE
##################################################
library(lattice)    # mostly for xyplot
library(ggplot2)    # grammar of graphics?
## library(effects)    # attractive interaction plots.  Does not work with aovlist.

if ( exists('SECC') == FALSE ) stop(
	"No SECC data found to analyze!  Don't forget to source(./lib/init.R)"
  )
if ( exists('Y.col') == FALSE ) stop(
	"Please specify a data column to analyze (Y.col)."
  )

##================================================
## Standard Labels
##================================================

Header.msd  <- "\n95% Minimum Significant Ranges (MSR):\n"


##################################################
## PROCESS DATA: planned
##################################################
## strip empty rows (rows with only NAs)
SECC.use <- strip_empty_dims( SECC, dim = 1, col.class = "numeric" )  

## Filter data for analysis, according to settings above.
SECC.use <- SECC.use[SECC.use$Time     %in% Time.use     &
                     SECC.use$Chamber  %in% Chamber.use  &
                     SECC.use$Frag     %in% Frag.use     &
                     SECC.use$Position %in% Position.use 
                     , ]
## drop unused factor levels (for plotting)
SECC.use <- within( SECC.use, {
  Time     <- factor(Time,     levels = Time.use)
  Chamber  <- factor(Chamber,  levels = Chamber.use)
  Frag     <- factor(Frag,     levels = Frag.use)
  Position <- factor(Position, levels = Position.use)
})

## Summarize data by means across different (or all) positions to prevent unbalanced effects?
## Meta-community (regional) -level analyses (ignoring position):
SECCmc <- SECC_aggregate( SECC.use, trt = 'Frag' )

## aggregate 'other' patch Positions
SECCp  <- SECC_aggregate( SECC.use, trt = 'Position' )


##=================================================
## PROCESS DATA: unplanned?
##=================================================
## plot label for transformed data.
Y.tlabel <- paste( Y.use, ": ", Y.label, sep=""  )  # label for transformation plots

## Transformations, calculated columns, etc.  May depend on the results of exploratory graphs (& assumptions), below.
for ( dataset in Dataset.list ){
  DataObject <- get(dataset)
  DataObject <- within( DataObject, {
	## Generic Response Variable (used in remainder of script)
	Y <- as.numeric( get(Y.col) )

    Y.trans  <- Y                # temporary, for transformations.
    Y.trans[Y.trans < 0] <- 0    # negative values are problematic for log and other transformations :-(
    Y.sqrt   <- sqrt( Y.trans )  # useful for Poisson-distributed data (mean prop. to variance).
    Y.4rt    <- Y.trans^(0.25)   # fourth-root
    Y.ln     <- log(Y.trans +1)  # defaults to base e=exp(1).
    Y.log    <- log(Y.trans +1, 10 )    # base 10. (stdev prop to mean).
    
	Y.trans  <- get(Y.use)      # assign which column to work with for analyses ****
  })    
  ## attach relevant info to object attributes.
  attr(DataObject, "response variable") <- Y.col
  attr(DataObject, "transformation")    <- Y.use
  assign(dataset, DataObject)
}


##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run if source()d
  head(SECCp)     # have a peek at the first 6 rows & columns: is this what you expected?
  str( SECCp)     # check structure: are the appropriate variables factors, numeric, etc.?
  summary(SECCp)  # summary statistics
  ## Regional analyses
  head(SECCmc)    # have a peek at the first 6 rows & columns: is this what you expected?
  str( SECCmc)    # check structure: are the appropriate variables factors, numeric, etc.?
  summary(SECCmc) # summary statistics
}

##################################################
## EXPLORE: PLOTS
##################################################
## make some meaningful plots of data to check for predicted (expected) patterns.

if (Save.results == TRUE && is.null(Save.plots) == FALSE) pdf( file = Save.plots )

par( mfrow=c(2,2), cex=0.8) # panel of figures: 2 rows & 2 columns
## Patch analyses
with( SECCp, plot(Y.trans ~ Chamber * Frag * Position, 
                  main=Dataset.labels[1], ylab = Y.tlabel
                  )
     )  # fixed effects only, no nesting
##* PROMPT *##
## plot.new()   # empty panel
plot.new()  # empty panel
par( mfrow=c(2,2), cex=0.8) # panel of figures: 2 rows & 2 columns
## Is the response variable normally-distributed? Check skew in histograms, eyeball linearity on QQ-plots
for( i in 1:length(Dataset.list) ){
    dataset <- Dataset.list[i]
    with( get(dataset), {
        hist( Y,      sub = Dataset.labels[i] )
        hist( Y.log,  sub = Dataset.labels[i] )
        hist( Y.sqrt, sub = Dataset.labels[i] )
        hist( Y.4rt,  sub = Dataset.labels[i])
        qqnorm( Y,      main = "untransformed",        sub = Dataset.labels[i] )
        qqline( Y,      col="grey50" )  # draw line through 1st & 3rd quartiles.
        qqnorm( Y.log,  main = "log10-transformed",    sub = Dataset.labels[i] )
        qqline( Y.log,  col="grey50" )
        qqnorm( Y.sqrt, main = "sqrt-transformed",     sub = Dataset.labels[i] )
        qqline( Y.sqrt, col="grey50" )
        qqnorm( Y.4rt,  main = "4th-root-transformed", sub = Dataset.labels[i] )
        qqline( Y.4rt,  col="grey50" )
    })
}



##################################################
## DEFINE MODEL FORMULA
##################################################
## Nested ANOVA: specify Error terms, with largest unit first
## in descending order of size (least to most replication).
## Do not include smallest unit/factor in Error.
## see Crawley, pp 469-470

## Including Time as a factor?
if ( length(Time.use) > 1 ) {
  Yp.model <- Y.trans ~ Time*Chamber*Frag*Position +Error(Block/Time/Chamber/Frag)
  Ymc.model <- Y.trans ~ Time*Chamber*Frag +Error(Block/Time/Chamber)
} else {
  ## Nested Fixed Effects, with error term for ANOVA using aov() 
  Yp.model <- Y.trans ~ Chamber*Frag*Position +Error(Block/Chamber/Frag)
  ## ignoring effect of position: 'regional' effects only
  Ymc.model <- Y.trans ~ Chamber*Frag +Error(Block/Chamber)
}


##################################################
## ANALYSIS: design
##################################################
Yp.aov <- aov( Yp.model, data=SECCp )   # simple ANOVA with nested fixed effects.  Nesting means an object of type "aovlist" is produced* (therefore regular TukeyHSD & other functions won't work )
Ymc.aov <- aov( Ymc.model, data=SECCmc )    # regional effects only.


##################################################
## CHECK ASSUMPTIONS: analyse residuals, standard diagnostic plots
##################################################
## anova
## independence?
    ## experimental design: random position of Frag within Chambers.
    ## possibility of Block gradient (7&8 SW -> 1-6 E:N), which also corresponds roughly to order of samples.
## Patch analyses
## trellis plots: any pattern across blocks, within frag & chambers?
xyplot(Y.trans ~ Block | Frag + Chamber, data=SECCp, 
       pch=21, col="black", bg="grey", cex=0.8,
       main = Dataset.labels[1]
       )
par( mfrow=c(2,2), cex=0.8) # panel of figures: 2 rows & 2 columns
## homogeneity of variances?
with( SECCp, plot(Y.trans ~ Chamber*Frag*Position) )    # fixed effects only, no nesting
##* PROMPT *##
## plot.new()  # put next plot in empty panel?
## normal distribution?
Yp.residuals <- resid(Yp.aov$Within)
with(SECCp, qqnorm( Yp.residuals, main="Residuals", sub=Dataset.labels[1] ) )   # are residuals normally distributed?
qqline(Yp.residuals,  col="grey50")
par( mfrow=c(1,1) )
hist(Yp.residuals)  # plot residuals
# with(SECCp, shapiro.test( Yp.residuals ) )    # normality?


##================================================
## REGIONAL analyses
xyplot( Y.trans ~ Block | Frag + Chamber, data=SECCmc, 
       pch=21, col="black", bg="grey", cex=0.8,
       main = Dataset.labels[2]
       )
par( mfrow=c(2,2), cex=0.8) # panel of figures: 2 rows & 2 columns
# homogeneity of variances?
with( SECCmc, plot(Y.trans ~ Chamber*Frag) )    # fixed effects only, no nesting
# normal distribution?
##* PROMPT *##
## plot.new()  # for stupid "Hit <Return> ..." prompt when going line-by-line :(
Ymc.residuals <- resid(Ymc.aov$Within)
with(SECCmc, qqnorm( Ymc.residuals, main="Residuals", sub=Dataset.labels[2] ) ) # are residuals normally distributed?
qqline(Ymc.residuals,  col="grey50")
par( mfrow=c(1,1) )
hist(Ymc.residuals) # plot residuals
## with(SECCmc, shapiro.test( Ymc.residuals ) )  # normality?



##################################################
## ANALYSIS: GET RESULTS
##################################################
## Patch analyses
# names(Yp.aov)
cat("\n\n")
print(Yp.model)                 # for output
print( summary(Yp.aov) )        # summary statistics
model.tables(Yp.aov, "means")   # effect sizes
# Interaction Plots
par(mfrow=c(2,2))   # panel of figures: 2 rows & 2 columns
with( SECCp, interaction.plot(Frag, Chamber, Y.trans,
                              ylab = paste("mean of", Y.use)
                              )
     )
with( SECCp, interaction.plot(Position, Chamber, Y.trans,
                              ylab=paste("mean of", Y.use)
                              )
     )
with( SECCp, interaction.plot(Position, Frag, Y.trans,
                              ylab=paste("mean of", Y.use)
                              )
     )

if (FALSE) {    # effects package doesn't work with aovlist either :(
  Yp.effects <- allEffects(Yp.aov)
  plot(Yp.effect, "Chamber:Frag")
  plot(Yp.effect, "Chamber:Position")
  plot(Yp.effect, "Chamber:Frag:Position")
}


##________________________________________________
## (un)planned Multiple Comparisons using Least Significant Differences (LSD) -> comparison intervals for graphical display.
# Chamber x Pos Interaction
msd <- MSD( Yp.aov, alpha=0.05, method="unplanned" )  # compute MSRs based on a 5% error rate (alpha), 2-tailed.  mode.df="manual" if unbalanced data ( provide n as an estimate).
msd.CxP <- msd["Chamber:Position"]
msd.FxP <- msd["Frag:Position"]
msd.CxFxP <- msd["Chamber:Frag:Position"]


##================================================
## Regional analyses
## names(Ymc.aov)
cat("\n\n")
print(Ymc.model)                # for output
print( summary(Ymc.aov) )       # summary statistics
model.tables(Ymc.aov, "means")  # effect sizes
# Interaction Plots
par(mfrow=c(1,1))   # panel of figures: 1 rows & 1 columns
with( SECCmc, interaction.plot( Frag, Chamber, Y.trans,
                               ylab=paste("mean of ", Y.use)
                               )
     )

##________________________________________________
## (un)planned Multiple Comparisons using Least Significant Differences (LSD) -> comparison intervals for graphical display.
msd.mc <- MSD( Ymc.aov, alpha=0.05, method="unplanned" )   # compute LSDs based on a 5% error rate (alpha), 2-tailed.
msd.mc.FxC <- msd.mc["Chamber:Frag"]
msd.mc.C <- msd.mc["Chamber"]



##################################################
## SAVE OUTPUT
##################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  capture.output(cat(Save.header, Save.patch.header, sep=""),
				 print(Yp.model),                 # model
				 summary(Yp.aov),                 # model summary
				 cat("\n\n"),                     # for output
				 model.tables(Yp.aov, "means"),   # effect sizes
				 cat(Header.msd),
				 print(msd),
				 cat(Save.mc.header),             # Meta-Community Results #
				 print(Ymc.model),                # model
				 summary(Ymc.aov),                # model summary
				 cat("\n\n"),                     # for output
				 model.tables(Ymc.aov, "means"),  # effect sizes
				 cat(Header.msd),
				 print(msd.mc),
				 cat(Save.end),                   # END OUTPUT #
				 file = Save.text
				)
}

if (Save.results == TRUE && is.null(Save.plots) == FALSE && Save.plots != Save.final) dev.off()



##################################################
## FINAL GRAPHICS
##################################################
if (Save.results == TRUE && is.null(Save.final) == FALSE && Save.plots != Save.final) pdf( file = Save.final )


Chamber.map <- plotMap( factor = "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Frag.map <- plotMap( factor = "Frag", labels = levels(SECC$Frag) )
Frag.map <- Frag.map[ levels(SECC$Frag) %in% Frag.use, ]
Position.map <- plotMap( factor = "Position", labels = levels(SECC$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]

if (length(Time.use) > 1) {
  Time.label <- ""
} else {
  Time.label <- paste(Time.use, ": ", sep="")
}
Plot.Title <- bquote(.(Time.label) * "Patch means " %+-% "95% Comparison Intervals")
Sub.msd <- "95% comparison intervals (MSR)" 

par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )


## Patch results: Chamber x Position
plot.means <- with( SECCp, 
                   aggregate( cbind( Y.trans ), 
                             list(Position = Position, Chamber = Chamber), 
                             mean 
                             ) 
                   )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plot.error <- matrix( as.numeric(msd["Chamber:Position"]/2),
                       nrow = length(levels(Chamber)),
                       ncol = length(levels(Position))
                       )
  plotMeans( Y.trans, Chamber, Position, 
            error.bars = "custom", level = plot.error, cex = 2, lwd = 2,
            lty = Position.map$lty, pch = Position.map$pch,
            col = as.character(Position.map$col),
            bg  = as.character(Position.map$bg),
            main = Plot.Title,
            sub  = Sub.msd,
            xlab = attr(SECC, "labels")[["Chamber"]],
            ylab = Y.plotlab
            )   
  ## as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## Patch results: Frag x Position (significant in t4)
plot.means <- with( SECCp, 
                   aggregate( cbind( Y.trans ), 
                             list(Pos = Position, Frag = Frag), 
                             mean 
                             ) 
                   )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plot.error <- matrix( as.numeric(msd.FxP/2),
                       nrow = length(levels(Frag)),
                       ncol = length(levels(Pos))
                       )
  plotMeans( Y.trans, Frag, Pos, 
            error.bars = "custom", level = plot.error, cex = 2, lwd = 2,
            lty = Position.map$lty, pch = Position.map$pch,
            col = as.character(Position.map$col),
            bg  = as.character(Position.map$bg),
            main = Plot.Title,
            sub  = Sub.msd,
            xlab = attr(SECC, "labels")[["Frag"]],
            ylab = Y.plotlab
            )
})

## Patch results: Chamber x Frag
plot.means <- with( SECCp, 
                   aggregate( cbind( Y.trans ), 
                             list(Frag=Frag, Chamber=Chamber), 
                             mean 
                             ) 
                   )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plot.error <- matrix( as.numeric(msd["Chamber:Frag"]/2),
                       nrow = length(levels(Frag)),
                       ncol = length(levels(Chamber))
                       )
  plotMeans( Y.trans , Frag , Chamber, 
            error.bars="custom", level=plot.error, cex=2, lwd=2,
            lty=Chamber.map$lty, pch=Chamber.map$pch,
            col=as.character(Chamber.map$col),
            bg=as.character(Chamber.map$bg),
            main = Plot.Title,
            sub  = Sub.msd,
            xlab = attr(SECC, "labels")[["Frag"]],
            ylab = Y.plotlab
            )   
})



##================================================
## META-COMMUNITY results
Plot.Title <- bquote(.(Time.label) * "Meta-Community means " %+-% "95% Comparison Intervals")

## Chamber Main Effects
plot.means <- with( SECCmc, 
                   aggregate( cbind( Y.trans ), 
                             list(Chamber=Chamber), 
                             mean 
                             ) 
                   )
plot.error <- rep( as.numeric(msd.mc.C/2), length(levels(plot.means$Chamber)) )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plotMeans( Y.trans , Chamber, 
            error.bars="custom", level=plot.error, cex=2, lwd=2,
            lty=1, pch=Chamber.map$pch,
            col=as.character(Chamber.map$col),
            bg=as.character(Chamber.map$bg),
            main = Plot.Title,
            sub  = Sub.msd,
            xlab = attr(SECC, "labels")[["Chamber"]],
            ylab = Y.plotlab
            )   # as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## Fragmentation Main Effects
plot.means <- with( SECCmc, 
                   aggregate( cbind( Y.trans ), 
                             list(Frag=Frag), 
                             mean 
                             ) 
                   )
plot.error <- rep( as.numeric(msd.mc.C/2), length(levels(plot.means$Frag)) )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plotMeans( Y.trans , Frag, 
            error.bars="custom", level=plot.error, cex=2, lwd=2,
            lty=1, pch=Frag.map$pch,
            col=as.character(Frag.map$col),
            bg=as.character(Frag.map$bg),
            main = Plot.Title,
            sub  = Sub.msd,
            xlab = attr(SECC, "labels")[["Frag"]],
            ylab = Y.plotlab
            )   # as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## Interaction: Chamber x Frag
plot.means <- with( SECCmc, 
                   aggregate( cbind( Y.trans ), 
                             list(Frag=Frag, Chamber=Chamber), 
                             mean 
                             ) 
                   )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plot.error <- matrix( as.numeric(msd.mc.FxC/2),
                       nrow = length(levels(Frag)),
                       ncol = length(levels(Chamber))
                       )
  plotMeans( Y.trans , Frag , Chamber, 
            error.bars="custom", level=plot.error, cex=2, lwd=2,
            lty=Chamber.map$lty, pch=Chamber.map$pch,
            col=as.character(Chamber.map$col),
            bg=as.character(Chamber.map$bg),
            main = Plot.Title,
            sub  = Sub.msd,
            xlab = attr(SECC, "labels")[["Frag"]],
            ylab = Y.plotlab
            )
})

if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
