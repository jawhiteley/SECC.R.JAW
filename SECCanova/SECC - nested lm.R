##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Standard Nested modelling of experimental data (unbalanced data)
### Response Variable(s)  @ time #s
### Jonathan Whiteley     R v2.12     2012-02-14
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
## library(nlme)       # non-linear (& nested) linear models - convergence problems
library(lme4)       # linear models with S4 classes
library(lattice)    # mostly for xyplot
library(ggplot2)    # grammar of graphics?

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
## drop rows with NA in column for analysis
##  this also affects factor levels, and therefore model structure...
SECC.use <- SECC.use[-which(is.na(SECC.use[[Y.col]])), ]  

## Filter data for analysis, according to settings above.
SECC.use <- SECC.use[SECC.use$Time     %in% Time.use     &
                     SECC.use$Chamber  %in% Chamber.use  &
                     SECC.use$Frag     %in% Frag.use     &
                     SECC.use$Position %in% Position.use 
                     , ]
## drop unused factor levels (for plotting)
SECC.use <- within( SECC.use, {
  Time     <- factor(Time,     exclude=NULL)
  Chamber  <- factor(Chamber,  exclude=NULL)
  Frag     <- factor(Frag,     exclude=NULL)
  Position <- factor(Position, exclude=NULL)
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
## Nested: specify Error terms, with largest unit first
## in descending order of size (least to most replication).
## Do not include smallest unit/factor in Error.
## see Crawley, pp 469-470

## Including Time as a factor?
if ( length(unique(SECCp$Time)) > 1 ) {
  Yp.model   <- Y.trans ~ Time*Chamber*Frag*Position
  Yp.random  <- ~ 1 | Block/Time/Chamber/Frag
} else {
  Yp.model  <- Y.trans ~ Chamber*Frag*Position
  Yp.random <- ~ 1 | Block/Chamber/Frag 
}


##################################################
## ANALYSIS: design
##################################################
## neither nlme nor lmer seem to be able to converge (with unbalanced data)!!

lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 500, msMaxIter = 100, opt="optim")

Yp.lme  <- lme( Yp.model, random = Yp.random, data=SECCp, method="REML", 
              control = lmc, na.action = na.omit)  # Linear model with nested error terms??

Yp.lmer <- lmer( Y.trans ~ Chamber*Frag*Position + (1 | Block/Chamber/Frag), data=SECCp, REML = TRUE)  # Linear model with nested error terms??
## Frag*Position is causing the problem!!

Yp.lm <- lm( Y.trans ~ Chamber*Frag*Position, data=SECCp)  # Linear model with nested error terms??

Yp.aov <- aov( Y.trans ~ Chamber*Frag*Position + Error(Block/Chamber/Frag), data=SECCp)

Yp.fit <- Yp.lme

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
if (inherits(Yp.fit, "aovlist"))
{
  Yp.residuals <- resid(Yp.fit$Within)
} else {
  Yp.residuals <- resid(Yp.fit)
}
with(SECCp, qqnorm( Yp.residuals, main="Residuals", sub=Dataset.labels[1] ) )   # are residuals normally distributed?
qqline(Yp.residuals,  col="grey50")
par( mfrow=c(1,1) )
hist(Yp.residuals)  # plot residuals
## with(SECCp, shapiro.test( Yp.residuals ) )    # normality?




##################################################
## ANALYSIS: GET RESULTS
##################################################
## Patch analyses
# names(Yp.aov)
cat("\n\n")
print(Yp.model)                 # for output
print( summary(Yp.fit) )        # summary statistics
if( isTRUE( paste("anova", class(Yp.fit), sep=".") %in% methods(anova) ) ) print( anova(Yp.fit) )
Yp.mtab <- try( model.tables(Yp.fit, "means")   # effect sizes
               , silent = TRUE)        # wrapped in try() statement, because unbalanced designs throw errors :(
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


##________________________________________________
## unplanned Multiple Comparisons using multcomp package
##   -> comparison intervals for graphical display.
library(multcomp)   # multiple comparisons
op <- par(mfrow=c(1, 1), mai = par("mai") * c(1, 3, 1, 1))

Yp.glht <- glht(Yp.fit)
confint(Yp.glht)
plot(Yp.glht)
## Interaction terms?
X.grid <- expand.grid(Chamber = unique(SECCp$Chamber), Position = unique(SECCp$Position))
X <- model.matrix(~ Chamber * Position, data = X.grid)
## X.pred <- predict(Yp.fit, newdata = X.grid) # only for the full set of interactions, not subset :(
Tukey <- contrMat(table(SECCp[, c("Chamber", "Position")]), "Tukey") # is this useful?
Tukey <- contrMat(table(SECCp[, "Chamber"]), "Tukey")
Tukey0 <- matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)) 
K1 <- cbind( Tukey, Tukey0 )
rownames(K1) <- paste(levels(SECCp$Chamber)[1], rownames(K1), sep=":")
K2 <- cbind( Tukey0, Tukey )
rownames(K2) <- paste(levels(SECCp$Chamber)[2], rownames(K2), sep=":")
K <- rbind(K1, K2)
colnames(K) <- c(colnames(Tukey), colnames(Tukey))

## Yp.mcp <- glht(Yp.fit, linfct = K)                 # not all terms included?

## fake factor for Chamber x Position interaction
SECCp$CxP <- with(SECCp, interaction(Chamber, Position) )
CxP.fit <- lme(Y.trans ~ CxP, random = Yp.random, data = SECCp) #  with or without intercept?
Yp.mcp <- glht(CxP.fit, linfct = mcp(CxP="Tukey"))
K <- mcp(CxP="Tukey")                  # differences between each combination of factors

## custom contrasts to get estimates & CI of factor levels, not differences
CxP.fit <- lme(Y.trans ~ Chamber * Position, random = Yp.random, data = SECCp) #  with or without intercept?
rownames(X.grid) <- paste(X.grid[, 1], X.grid[, 2], sep=":")
K <- model.matrix(~ Chamber * Position, data = X.grid)
Yp.mcp <- glht(CxP.fit, linfct = K)
## this produces the right type of information, but there's something wrong, because the estimate for Ambient:Inner is way too high.

summary(Yp.mcp)
Yp.ci <- confint(Yp.mcp)                        # estimates of DIFFERENCES and confidence intervals, adjusted for multiple comparisons.  This is close to what I want. (depends on model / contrast matrix K)
plot(Yp.ci)                           # This is the data I want, but not quite the graph layout I want?
Yp.ci$confint[, 3] - Yp.ci$confint[, 1] # uneven confidence widths, due to uneven sample sizes

par(op)


##________________________________________________
## Effects display
if (!inherits(Yp.fit, "aovlist"))
{
  library(effects)    # attractive interaction plots.  Does not work with aovlist.

  Yp.effects <- allEffects(Yp.fit)
  plot(Yp.effects, "Chamber:Frag:Position")
  CxP.eff <- effect("Chamber:Position", Yp.fit)
  CxPs   <- summary(CxP.eff)
  plot(CxP.eff)
  if (F)
  {
    CxP.eff2 <- effect("Chamber:Position", CxP.fit) # model with only interaction terms...
    plot(CxP.eff2)
  }
}


##________________________________________________
## (un)planned Multiple Comparisons using Least Significant Differences (LSD) -> comparison intervals for graphical display.
## Requires balanced data :/
if (F)
{
  ## Chamber x Pos Interaction
  msd <- MSD( Yp.aov, alpha=0.05, method="unplanned" )  # compute MSRs based on a 5% error rate (alpha), 2-tailed.  mode.df="manual" if unbalanced data ( provide n as an estimate).
  msd.CxP <- msd["Chamber:Position"]
  msd.FxP <- msd["Frag:Position"]
  msd.CxFxP <- msd["Chamber:Frag:Position"]
}



##################################################
## SAVE OUTPUT
##################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  capture.output(cat(Save.header, Save.patch.header, sep=""),
				 print(Yp.model),                 # model
				 summary(Yp.fit),                 # model summary
				 cat("\n\n"),                     # for output
				 Yp.mtab,                         # effect sizes
				 cat(Save.end),                   # END OUTPUT #
				 file = Save.text
				)
}



##################################################
## POLISHED GRAPHICS (almost FINAL)
##################################################

if (exists("Y.lim") == FALSE) Y.lim <- range(SECCp$Y.trans) # consistent scale of Y axis

Chamber.map <- plotMap( factor = "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Frag.map <- plotMap( factor = "Frag", labels = levels(SECC$Frag) )
Frag.map <- Frag.map[ levels(SECC$Frag) %in% Frag.use, ]
Position.map <- plotMap( factor = "Position", labels = c("Inner", "other", "Outer") )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]

if (length(Time.use) > 1) {
  Time.label <- ""
} else {
  Time.label <- paste(Time.use, ": ", sep="")
}
## xyplot takes expression()s for plotmath, but not bquote()!! :(
Plot.Title <- bquote(.(Time.label) * "Patch means " %+-% "95% Confidence Intervals")
Sub.msd <- ""                          # "95% comparison intervals (MSR)"
## Y.plotlab <- paste(attr(SECC, "labels")[[Y.col]], attr(SECC, "units")[[Y.col]])

par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )


## Patch results: Chamber x Position
print( plot(CxP.eff, alternating = FALSE, multiline = TRUE, confint = TRUE,
            colors = Position.map$col, ylim = Y.lim,
            key = list( pch = Position.map$pch, lty = Position.map$lty),
            xlab = attr(SECC, "labels")[["Chamber"]],
            )
)


with( SECCp, 
     {
       ## using custom plotMeans function, with custom error bars (LSD)
       plot.error <- CxPs$effect - CxPs$lower
       plotMeans( Y.trans, Chamber, Position, 
                 error.bars = "custom", level = plot.error, ylim = Y.lim,
                 cex = 2, lwd = 2, lty = Position.map$lty, pch = Position.map$pch,
                 col = as.character(Position.map$col),
                 bg  = as.character(Position.map$bg),
                 main = Plot.Title,
                 sub  = Sub.msd,
                 xlab = attr(SECC, "labels")[["Chamber"]],
                 ylab = Y.plotlab
                 )   
       ## as.character() is needed for string arguments (color hex strings), if they are stored as factors().
     }
)

## Patch results: Frag x Position
FxP.eff <- effect("Frag:Position", Yp.fit)
FxPs <- summary(FxP.eff)
plot.error <- FxPs$effect - FxPs$lower
with( SECCp, 
     {
       ## using custom plotMeans function, with custom error bars (LSD)
       plotMeans( Y.trans, Frag, Position, 
                 error.bars = "custom", level = plot.error, ylim = Y.lim,
                 cex = 2, lwd = 2, lty = Position.map$lty, pch = Position.map$pch,
                 col = as.character(Position.map$col),
                 bg  = as.character(Position.map$bg),
                 main = Plot.Title,
                 sub  = Sub.msd,
                 xlab = attr(SECC, "labels")[["Frag"]],
                 ylab = Y.plotlab
                 )
     }
)

## Patch results: Chamber x Frag
CxF.eff <- effect("Chamber:Frag", Yp.fit)
CxFs <- summary(CxF.eff)
plot.error <- t(CxFs$effect - CxFs$lower)
with( SECCp, 
     {
       ## using custom plotMeans function, with custom error bars (LSD)
       plotMeans( Y.trans , Frag , Chamber, 
                 error.bars="custom", level=plot.error, ylim = Y.lim,
                 cex = 2, lwd = 2, lty=Chamber.map$lty, pch=Chamber.map$pch,
                 col=as.character(Chamber.map$col),
                 bg=as.character(Chamber.map$bg),
                 main = Plot.Title,
                 sub  = Sub.msd,
                 xlab = attr(SECC, "labels")[["Frag"]],
                 ylab = Y.plotlab
                 )   
     }
)



if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
