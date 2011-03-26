##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Template for basic analyses of experimental data
### Response Variable(s)  @ time #s
### Jonathan Whiteley     R v2.12     2011-03-25
##################################################
## INITIALISE
##################################################
## Set Working Directory: path in quotes "".
# setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")    # iMac@McGill
# setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
getwd()  # Check that we're in the right place

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')
library(lattice)    # mostly for xyplot
library(ggplot2)    # grammar of graphics
par(ask = FALSE)    # Stop asking me to hit <Return> to see next plot!
options(device.ask.default = FALSE) # same as above: does this work?

##################################################
## CONFIGURE BASIC ANALYSIS
##################################################
## Can this script be used in a generic way for
## most univariate analyses?

##================================================
## SETTINGS - edit
##================================================

## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[3]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]      # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)         # Frag treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)]     # Patch Positions to include
Y.col        <- 'Cells.m'                    # Column to use for response variable.

## Define Labels
## Which response variable is being used (for labels)? ****
Y.use <- "Y.4rt"
Y.label <- attr(SECC, "labels")[[Y.col]]  # response variable label for this script.
Y.units <- attr(SECC, "units" )[[Y.col]]  # response variable units for this script.
Y.units <- bquote( paste( .(Y.label), " (", 
                         sqrt("cells" %.% m^-2, 4),  # manually (for now)
                         ")",
                         sep=""
                         )
                  )
Dataset.list <- c("SECCp", "SECCmc")
Dataset.labels <- c( "Patch scale data", "Meta-Community scale data" )
ID.cols <- c('SampleID', 'Time', 'Block', 'Chamber', 'Frag', 'Pos', 'Position')
Trt.nested <- c('Time', 'Block', 'Chamber', 'Frag', 'Position')

## Output Files - set to NULL to prevent output.
Out.filename  <- paste("Results - ", Y.col, " - ",
                   which(levels(SECC$Time) == Time.use), sep = ""
                   )
Out.text  <- paste("./output/", Out.filename, ".txt", sep = "")
Out.plots <- paste("./graphs/", Out.filename, ".pdf", sep = "")
Out.final <- Out.plots              # Destination for final plots.
Out.header <- paste(  "Nested ANOVA Results for:", Y.label, "(", Y.col, ")",
                    "\nTransformation used:     ", Y.use,
                    "\nSample Time:  ", paste(Time.use,     collapse = ", "),
                    "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                    "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                    "\nPatches:      ", paste(Position.use, collapse = ", "),
                    "\n\n",
                    date(),
                    "\n\n================================================================\n"
                    )


##================================================
## PREPARE DATA - do not edit
##================================================
## strip empty rows (rows with only NAs)?
empty.rows <- which( apply( SECC[, lapply(SECC, class) == "numeric"], 1,
                           function(x) all(is.na(x))
                           )
                    )
SECC.use <- SECC[-empty.rows, ]  # !is.na(SECC$Time) ; NAs in factors are annoying
## SECC.use <- SECC                 # make a copy, for further processing (save the original for reference).


##================================================
## CALCULATIONS - edit
##================================================
# str(SECC.use)
sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample
ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * (12.5^2)      # patch area
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2
Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.

SECC.use <- within( SECC.use, { 
  ## Generic Response Variable (used in remainder of script)
  Nfix  <- ARA.m * Nfix.ARA.ratio
  Y <- as.numeric( get(Y.col) )  # do not need to edit this
  Y[Y < 0] <- 0    ## negative values cause problems for transformations: replace with 0
  ## change negative ARA values to 0 - should I wait until after aggregation?
})



##################################################
## PROCESS DATA: planned
##################################################

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
Y.tlabel <-      paste( Y.use, ": ", Y.label, sep=""  )  # label for transformation plots
Y.plabel <- expression( paste(Y.label) * " ("* paste(Y.units) * ")" )  # plot label?

## Transformations, calculated columns, etc.  May depend on the results of exploratory graphs (& assumptions), below.
for ( dataset in Dataset.list ){
  DataObject <- get(dataset)
  DataObject <- within( DataObject, {
    Y.trans  <- Y                # temporary, for transformations.
    Y.trans[Y.trans < 0] <- 0        # negative values are problematic for log and other transformations :-(
    Y.sqrt   <- sqrt( Y.trans ) # useful for Poisson-distributed data (mean prop. to variance).
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
head(SECCp)     # have a peek at the first 6 rows & columns: is this what you expected?
str( SECCp)     # check structure: are the appropriate variables factors, numeric, etc.?
summary(SECCp)  # summary statistics
## Regional analyses
head(SECCmc)    # have a peek at the first 6 rows & columns: is this what you expected?
str( SECCmc)    # check structure: are the appropriate variables factors, numeric, etc.?
summary(SECCmc) # summary statistics


##################################################
## EXPLORE: PLOTS
##################################################
## make some meaningful plots of data to check for predicted (expected) patterns.

if (is.null(Out.plots) == FALSE) pdf( file = Out.plots )

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
# Nested Fixed Effects, with error term for ANOVA using aov() 
Yp.model <- Y.trans ~ Chamber*Frag*Position +Error(Block/Chamber/Frag)
# ignoring effect of position: 'regional' effects only
Ymc.model <- Y.trans ~ Chamber*Frag +Error(Block/Chamber/Frag)


##################################################
## ANALYSIS: design
##################################################
Yp.aov <- aov( Yp.model, data=SECCp )   # simple ANOVA with nested fixed effects.  Nesting means an object of type "aovlist" is produced* (therefore regular TukeyHSD & other functions won't work )
Ymc.aov <- aov( Ymc.model, data=SECCmc )    # regional effects only.


##################################################
## CHECK ASSUMPTIONS: analyse residuals, standard diagnostic plots
##################################################
## anova
# independence?
    # experimental design: random position of Frag within Chambers.
    # possibility of Block gradient (7&8 SW -> 1-6 E:N), which also corresponds roughly to order of samples.
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
Ymc.residuals <- resid(Ymc.aov$"Block:Chamber:Frag")
with(SECCmc, qqnorm( Ymc.residuals, main="Residuals", sub=Dataset.labels[2] ) ) # are residuals normally distributed?
qqline(Ymc.residuals,  col="grey50")
par( mfrow=c(1,1) )
hist(Ymc.residuals) # plot residuals
# with(SECCmc, shapiro.test( Ymc.residuals ) )  # normality?


##################################################
## ANALYSIS: GET RESULTS
##################################################
if (is.null(Out.text) == FALSE) {
  sink( file = Out.text, split = TRUE, type = "output" )
  cat(Out.header,
      "================  Patch scale Results  ================\n\n",
      sep=""
      )
}

## Patch analyses
# names(Yp.aov)
Yp.model                        # for output
summary(Yp.aov)                 # summary statistics
cat("\n\n")                     # for output
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
                              ylab=paste("mean of ", Y.use)
                              )
     )

##________________________________________________
## Planned Multiple Comparisons using Least Significant Differences (LSD) -> comparison intervals for graphical display.
# Chamber x Pos Interaction
lsd <- LSD( Yp.aov$Within, Yp.model, data=SECCp, alpha=0.05, mode="pairwise" )  # compute LSDs based on a 5% error rate (alpha), 2-tailed.  Manual, due to unbalanced data (this is an estimate).
lsd.CxP <- lsd["Chamber:Position"]
lsd.FxP <- lsd["Frag:Position"]
lsd.CxFxP <- lsd["Chamber:Frag:Position"]

cat("\n95% Least Significant Differences (LSD):\n")
lsd


##================================================
## Regional analyses
cat("\n",
    "================================================================", 
    "================  Meta-Community scale Results  ================",
    "================================================================\n",
    sep = "\n"
    )

## names(Ymc.aov)
Ymc.model
summary(Ymc.aov)        # summary statistics
cat("\n\n")
model.tables(Ymc.aov, "means")  # effect sizes
# Interaction Plots
par(mfrow=c(1,1))   # panel of figures: 1 rows & 1 columns
with( SECCmc, interaction.plot( Frag, Chamber, Y.trans,
                               ylab=paste("mean of ", Y.use)
                               )
     )

##________________________________________________
## Planned Multiple Comparisons using Least Significant Differences (LSD) -> comparison intervals for graphical display.
lsd.mc <- LSD( Ymc.aov$"Block:Chamber:Frag", Ymc.model, data=SECCmc, alpha=0.05, mode="pairwise" )   # compute LSDs based on a 5% error rate (alpha), 2-tailed.
lsd.mc.FxC <- lsd.mc["Chamber:Frag"]
lsd.mc <- LSD( Ymc.aov$"Block:Chamber", Ymc.model, data=SECCmc, alpha=0.05, mode="pairwise" )    # compute LSDs based on a 5% error rate (alpha), 2-tailed.
lsd.mc.C <- lsd.mc["Chamber"]

cat("\n95% Least Significant Differences (LSD):\n")
lsd.mc


if (is.null(Out.text) == FALSE) {
  cat("\n", "<============================= END ============================>",
      sep = "\n"
      )
  sink()
}

if (is.null(Out.plots) == FALSE && Out.plots != Out.final) dev.off()


##################################################
## FINAL GRAPHICS
##################################################
if (is.null(Out.final) == FALSE && Out.plots != Out.final) pdf( file = Out.final )


Chamber.map <- data.frame( label = levels(SECC$Chamber),
                          col = c("#000000", "#000099", "#990000"),
                          bg  = c("#FFFFFF", "#FFFFFF", "#FFFFFF"),
                          pch = c(21, 23, 18),
                          lty = c(3, 2, 1)
                          )
    ## Ambient = black, open circles with dotted line ; 
    ## Partial = blue, open diamonds with dashed line ; 
    ## Full   = red, solid diamond with solid line.
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Frag.map    <- data.frame( label = levels(SECCp$Frag),
                          col = c("#000000", "#666666", "#000099", "#990000"), 
                          bg  = c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF"), 
                          pch = c(19, 15, 22, 21),
                          lty = c(1, 2, 3, 3)
                          )
    ## Continuous         = black, filled circles with solid line ; 
    ## Full Corridors     =  grey, filled squares with dashed line ; 
    ## Pseudo-Corridors   = blue, open squares with dotted line.
    ## Isolated           =  red, open circles with dotted line.

## Patch results: Chamber x Position
plot.means <- with( SECCp, 
                   aggregate( cbind( Y.trans ), 
                             list(Pos = Position, Chamber = Chamber), 
                             mean 
                             ) 
                   )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plot.error <- matrix( as.numeric(lsd.CxP/2),
                       nrow = length(levels(Pos)),
                       ncol = length(levels(Chamber))
                       )
  plotMeans( Y.trans , Pos , Chamber, 
            error.bars = "custom", level = plot.error, cex = 2, lwd = 2,
            lty = Chamber.map$lty, pch = Chamber.map$pch,
            col = as.character(Chamber.map$col),
            bg  = as.character(Chamber.map$bg),
            main = bquote(paste(.(Time.use), ": Patch means " %+-% "95% LSD", sep="")),
            sub  = "95% comparison intervals (LSD)",
            xlab = attr(SECC, "labels")[["Pos"]],
            ylab = Y.units
            )   # as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## Patch results: Frag x Position (significant in t4)
plot.means <- with( SECCp, 
                   aggregate( cbind( Y.trans ), 
                             list(Pos = Position, Frag = Frag), 
                             mean 
                             ) 
                   )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plot.error <- matrix( as.numeric(lsd.FxP/2),
                       nrow = length(levels(Pos)),
                       ncol = length(levels(Frag))
                       )
  plotMeans( Y.trans , Pos , Frag, 
            error.bars = "custom", level = plot.error, cex = 2, lwd = 2,
            lty = Frag.map$lty, pch = Frag.map$pch,
            col = as.character(Frag.map$col),
            bg  = as.character(Frag.map$bg),
            main = bquote(paste(.(Time.use), ": Patch means " %+-% "95% LSD", sep="")),
            sub  = "95% comparison intervals (LSD)",
            xlab = attr(SECC, "labels")[["Pos"]],
            ylab = Y.units
            )   # as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})


##================================================
## META-COMMUNITY results
plot.means <- with( SECCmc, 
                   aggregate( cbind( Y.trans ), 
                             list(Frag=Frag, Chamber=Chamber), 
                             mean 
                             ) 
                   )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plot.error <- matrix( as.numeric(lsd.mc.FxC/2),
                       nrow = length(levels(Frag)),
                       ncol = length(levels(Chamber))
                       )
  plotMeans( Y.trans , Frag , Chamber, 
            error.bars="custom", level=plot.error, cex=2, lwd=2,
            lty=Chamber.map$lty, pch=Chamber.map$pch,
            col=as.character(Chamber.map$col),
            bg=as.character(Chamber.map$bg),
            main = bquote(paste(.(Time.use),
              ": Meta-Community means " %+-% "95% LSD",
              sep=""
              ) ),
            sub  = "95% comparison intervals (LSD)",
            xlab = attr(SECC, "labels")[["Frag"]],
            ylab = Y.units
            )   # as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## Chamber Main Effects
plot.means <- with( SECCmc, 
                   aggregate( cbind( Y.trans ), 
                             list(Chamber=Chamber), 
                             mean 
                             ) 
                   )
plot.error <- rep( as.numeric(lsd.mc.C/2), length(levels(plot.means$Chamber)) )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plotMeans( Y.trans , Chamber, 
            error.bars="custom", level=plot.error, cex=2, lwd=2,
            lty=1, pch=Chamber.map$pch,
            col=as.character(Chamber.map$col),
            bg=as.character(Chamber.map$bg),
            main = bquote(paste(.(Time.use),
              ": Meta-Community means " %+-% "95% LSD",
              sep=""
              ) ),
            sub  = "95% comparison intervals (LSD)",
            xlab = attr(SECC, "labels")[["Chamber"]],
            ylab = Y.units
            )   # as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

## Fragmentation Main Effects
plot.means <- with( SECCmc, 
                   aggregate( cbind( Y.trans ), 
                             list(Frag=Frag), 
                             mean 
                             ) 
                   )
plot.error <- rep( as.numeric(lsd.mc.C/2), length(levels(plot.means$Frag)) )
par( mfrow=c(1,1), lty=1, cex=1, lwd=1 )
with( plot.means, {
  ## using custom plotMeans function, with custom error bars (LSD)
  plotMeans( Y.trans , Frag, 
            error.bars="custom", level=plot.error, cex=2, lwd=2,
            lty=1, pch=Frag.map$pch,
            col=as.character(Frag.map$col),
            bg=as.character(Frag.map$bg),
            main = bquote(paste(.(Time.use),
              ": Meta-Community means " %+-% "95% LSD",
              sep=""
              ) ),
            sub  = "95% comparison intervals (LSD)",
            xlab = attr(SECC, "labels")[["Frag"]],
            ylab = Y.units
            )   # as.character() is needed for string arguments (color hex strings), but I'm still not entirely sure why.  If it is not used, that argument is essentially ignored, and (ugly) defaults are used instead.
})

if (is.null(Out.plots) == FALSE) dev.off()
