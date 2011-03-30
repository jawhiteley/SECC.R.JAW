##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### GLMM (regression)
### Acetylene Reduction Assay (ARA: N-fixation)
### vs. cyanobacteria density
### Jonathan Whiteley     R v2.12     2011-03-30
##################################################
## INITIALISE
##################################################
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
X.col <- 'Cells.m'   # Column to analyze as explanatory variable        *****
Y.col <- 'ARA.m'     # Column to analyze as response variable           *****

vars.ls   <- c("ARA.m", "Cells.m", "Hcells.m")  # for data exploration?

##================================================
## SETTINGS 
##================================================
## Specify which treatment levels to include (by index is probably easiest)
Time.use     <- levels(SECC$Time)[1]      # Time (index: 1-3) to include in this run
Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]      # Chamber treatments to include
Frag.use     <- levels(SECC$Frag)         # Frag treatments to include
Position.use <- levels(SECC$Position)[c(1, 3)]     # Patch Positions to include


## Output Results?
Save.results  <- FALSE  


##================================================
## CALCULATIONS 
##================================================
SECC.prime <- SECC    # save a copy of the original for reference.

# str(SECC)
sampleA  <- 6   # sample Area, in cm^2:  pi * (2.75/2)^2 ; pi * (2.8 / 2)^2
      #     6 for rough estimate of inner tube diameter (2.8 cm): pi*(2.8/2)^2,
      #  or 6.4 for 20 shoots, based on density survey.
sample.to.m2 <- (100*100)/sampleA   # scale sample area, in cm^2 to m^2
sample_ml    <- 50  # 50 ml sample
ARA.m2   <- sampleA/(100*100)  # ARA sample area,   in (cm^2 to) m^2
patchA   <- pi * (12.5^2)      # patch area
patch.m2 <- patchA/(100*100)   # patch sample area, in (cm^2 to) m^2
Nfix.ARA.ratio <- 1/3  # ratio of N-fixation : ARA.

SECC <- within( SECC, { 
  ## change negative ARA values to 0 - should I wait until after aggregation?
  ARA.ml[ARA.ml < 0] <- 0
  ARA.m[ ARA.m  < 0] <- 0
  ARA.g[ ARA.g  < 0] <- 0
  Nfix <- ARA.m * Nfix.ARA.ratio
})


##================================================
## LABELS
##================================================

X.label <- attr(SECC, "labels")[[X.col]]  # explanatory variable label
X.units <- attr(SECC, "units" )[[X.col]]  # explanatory variable units

Y.label <- attr(SECC, "labels")[[Y.col]]  # response variable label
Y.units <- attr(SECC, "units" )[[Y.col]]  # response variable units

X.plotlab <- bquote( .(X.label) * "  " * .(X.units) *  "" )
Y.plotlab <- bquote( .(Y.label) * "  " * .(Y.units) *  "" )

labels.ls <- c()
for(i in 1:length(vars.ls) ){
  var <- vars.ls[i]
  labels.ls[i] <- attr(SECC, "labels")[[var]]
}


## Save Output to Files - set to NULL to prevent output.
Save.filename <- paste("Results - ", Y.col, "~" , X.col, " - ",
                       paste(which(levels(SECC$Time) == Time.use), collapse=""),
                       sep = ""
                   )
Save.text  <- paste("./output/", Save.filename, ".txt", sep = "")
Save.plots <- paste("./graphs/", Save.filename, ".pdf", sep = "")
Save.final <- Save.plots              # Destination for final plots.

## Output text
Save.divider <-        "================================================================\n"
Save.header  <- paste( "GLM Results for:", Y.label, "(", Y.col, ")",
                     "\n               ~", X.label, "(", X.col, ")",
                     "\nExpt. Time:   ", paste(Time.use,     collapse = ", "),
                     "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                     "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                     "\nPatches:      ", paste(Position.use, collapse = ", "),
                     paste("\n\n", date(), "\n\n", Save.divider, sep = "")
                     )
Save.end      <- paste("\n", 
					   "<============================= END ============================>",
						sep = "\n"
					  )

###===============================================
### Include Time as a factor 
###===============================================
## Note that Samples at different times are actually independent
## in this design, due to destructive sampling.

## Time.use     <- levels(SECC$Time)      # Include *ALL* Times (as a Treatment)
## source("./R", echo = FALSE) # Load default Labels. *****
## source("./R", echo = FALSE) # RUN STANDARD nested ANOVA



##################################################
## PROCESS DATA: planned
##################################################
## filter & clean SECC data frame.
SECC.use <- SECCclean(SECC, Time.use, Chamber.use, Frag.use, Position.use)

## Summarize data by means across different (or all) positions to prevent unbalanced effects?
## aggregate 'other' patch Positions
SECCp  <- SECC.use  # SECC_aggregate( SECC.use, trt = 'Position' )
## Meta-community (regional) -level analyses (ignoring position):
SECCmc <- SECC_aggregate( SECC.use, trt = 'Frag' )

SECC.scale  <- "patch"  # c("patch", "mc")
SECCa <- if(SECC.scale == "patch") SECCp else if (SECC.scale == "mc") SECCmc else SECC


SECCa <- within( SECCa, {
                X <- as.numeric( get(X.col) )
                Y <- as.numeric( get(Y.col) )
                Y.use <- Y  # compatibility with older code
})

##################################################
## CHECK DATA
##################################################
if (FALSE) {  # do not run if source()d
  head(SECCa)     # have a peek at the first 6 rows & columns: is this what you expected?
  str( SECCa)     # check structure: are the appropriate variables factors, numeric, etc.?
  summary(SECCa)  # summary statistics
}


##################################################
## EXPLORE: PLOTS
##################################################
## make some meaningful plots of data to check for predicted (expected) patterns.
if (Save.results == TRUE && is.null(Save.plots) == FALSE) pdf( file = Save.plots )

Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
point <- 21	# 21 for circles with col & bg ; 16 for solid circles
Chamber.map$pch <- c(21, 16)

SECCa <- within( SECCa,{
	colr = ifelse( Chamber == Chamber.map$label[1], 
			Chamber.map$col[1], 
			Chamber.map$col[2] 
		)
	fill = ifelse( Chamber == Chamber.map$label[1], 
			Chamber.map$bg[1], 
			Chamber.map$bg[2] 
		)
	pt = ifelse( Chamber == Chamber.map$label[1], 
			Chamber.map$pch[1], 
			Chamber.map$pch[2]
		)
})

with( SECCa,{
	par(mfrow=c(1,1))

	## scatterplot
	plot(X, Y, type="p",
		ylab=Y.plotlab, xlab=X.plotlab, 
		pch=pt, col=colr
	)
	legend("topright", legend=Chamber.map$label, 
           pch=Chamber.map$pch, col=Chamber.map$col 
    )
})

par(mfcol=c(2,2))
## Check distributions
for(i in 1:length(vars.ls) ){
  var <- vars.ls[i]
  label <- labels.ls[i]
  with( SECCa,{
        cat(var, " ")
        X.var <- get(var)
        X.max  <- max( X.var )
        cat(X.max, "\n")
        freq.max <- length(X.var)/2
        X.maxD <- max( density( X.var )$y )*1.5
        for(Ch.trt in levels(Chamber)){
          X.trt <- X.var[Chamber==Ch.trt]
          X.density <- density( X.trt )
          hist( X.trt,
               main=Ch.trt, xlab=label,
               xlim=c(0, X.max), 
               ylim=c(0, freq.max),
               breaks=seq( 0, X.max, length.out=16 ),
               col="#CCCCCC"
               )
          abline( 5, 0, lty=3, col="#666666" ) # reference line
          plot( X.density,
               main=Ch.trt, xlab=label,
               xlim=c(0, X.max),
               ylim=c(0, X.maxD)
               )
          densityplot( X.trt,
                      main=Ch.trt, xlab=label,
                      xlim=c(0, max( X.var ) )
                      )	# ** TRELLIS plot
        }
    })
	# qqplot(x, y) to compare distributions.
}

coplot( Y ~ X | Frag * Position, data=SECCa, 
       pch=SECCa$pt, col=SECCa$colr	# , bg=Chamber.map$bg
)	# why does recycling Chamber.map work for bg, but not col?
## I think it's because of the way the arguments are passed: coplot has specific arguments for col & pch, but not bg: bg is simply passed on directly to the plotting function (points) within each panel.
coplot( Y ~ X | Chamber * Frag , data=SECCa, 
       pch=SECCa$pt, col=SECCa$colr	# , bg= Chamber.map$bg
)
coplot( Y ~ X | Chamber * Frag , data=SECCa, 
       pch=SECCa$pt, col=SECCa$colr	# , bg=SECC$fill
)

## much the same breakdown, using TRELLIS xyplot over 3 factors.
xyplot( Y ~ X | Chamber * Frag * Position, data=SECCa, 
       pch=point, col=SECCa$colr, bg = SECCa$fill 
)
xyplot( Y ~ X | Chamber * Frag , data=SECCa, col=1 
       ## , key=simpleKey(levels(SECCa$Chamber))
)



##################################################
## ANALYSIS: formula
##################################################
## GLMM

## Including Time as a factor?
if ( length(Time.use) > 1 ) {
  Y.formula <- Y ~ X*Time*Chamber*Frag*Position
} else {
  Y.formula <- Y ~ X*Chamber*Frag*Position
}

##################################################
## ANALYSIS: design
##################################################
# Should I not also be using a mixed-effects model to account for treatments of different sizes? Chamber / Frag / Position
# YES
# Variance Components Analysis (variance decomposition)

Y.model <- glm( Y.formula, data=SECCa, family="quasipoisson" )



##################################################
## CHECK ASSUMPTIONS
##################################################
## analyse residuals, standard diagnostic plots
par(mfrow=c(3,2))	    # panel of figures: 3 rows & 2 columns
plot(Y.model)
hist(resid(Y.model))    # plot residuals



##################################################
## ANALYSIS: GET RESULTS
##################################################
anova(Y.model)
summary(Y.model)


  
  
##################################################
## SAVE OUTPUT
##################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) {
  capture.output(cat(Save.header), 
				 print(Y.formula),              # model
				 anova(Y.model),                # model summary
				 summary(Y.model),              # model summary
				 cat("\n\n"),                   # for output
				 model.tables(Y.model, "means"),# effect sizes
				 cat(Save.end),                 # END OUTPUT #
				 file = Save.text
				)
}

if (Save.results == TRUE && is.null(Save.plots) == FALSE && Save.plots != Save.final) dev.off()



##################################################
## FINAL GRAPHICS
##################################################
if (Save.results == TRUE && is.null(Save.final) == FALSE && Save.plots != Save.final) pdf( file = Save.final )

## generate grid to add predicted values to (X-values in all combinations of factors).
Y.pred 	<- expand.grid( Chamber=levels(SECCa$Chamber) , 
                        Frag=levels(SECCa$Frag), 
                        Position=levels(SECCa$Position), 
                        X=seq(0, max(SECCa$X), length.out=100 ) 
                        )
Y.pred$preds  <- predict(Y.model, newdata=Y.pred, type="response" )	# newdata must have same explanatory variable name for predict to work.


Chamber.map <- plotMap( "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map$label <- factor(Chamber.map$label)
point <- 21	# 21 for circles with col & bg ; 16 for solid circles


par(mfrow=c(1,1))
with( SECCa,{
	colr = ifelse(Chamber=="Ambient", 
                  as.character(Chamber.map$col[1]), 
                  as.character(Chamber.map$col[2]) 
                  )
	bg.colr = ifelse(Chamber=="Ambient", Chamber.map$bg[1], Chamber.map$bg[2] )
	# pred | augpred | ?
	plot(X, Y, type="p",
		ylab=Y.plotlab, xlab=X.plotlab,
		pch=point, col=colr, bg=bg.colr
        )
	lines(preds ~ X, data=subset(Y.pred, Chamber=="Ambient"), 
          col = as.character(Chamber.map$col[1]), 
          lty = Chamber.map$lty[1]
          )
	lines(preds ~ X, data=subset(Y.pred, Chamber=="Full Chamber"), 
          col = as.character(Chamber.map$col[2]), 
          lty = Chamber.map$lty[2]
          )
	legend( "topright", legend=Chamber.map$label, pch=point, col=as.character(Chamber.map$col), pt.bg=as.character(Chamber.map$bg) )
})


if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()
