##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Moss Growth / Productivity
### Jonathan Whiteley     R v2.12     2012-02-05
##################################################
## INITIALISE
##################################################
## Working Directory: see lib/init.R below
if (FALSE) {  # do not run automatically
  setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis/ SECC/")  # iMac@McGill
  setwd("/Users/jaw/Documents/ My Documents/ Academic/McGill/PhD/Analysis/ SECC/")  # JAW-MBP
  setwd("./ SECC/") # relative to my usual default wd in R GUI (MBP).
  setwd("./")       # relative to this file (\rd in Vim-R)
  getwd()           # Check that we're in the right place
}


## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')
library(car)

## Moss growth is a little different than other measure variables:
## - Measured from the same t3 & t4 patches throughout.
##   - True repeated measures, unlike destructive sampling
##   - t3 data not included in SECC data frame
## - Unbalanced data (missing values)
##   - this causes serious problems for model.tables,
##     which is used to calculate the MSD error bars :(
##   - A much older version of the analysis using LSD bars seemed to work, however :/
## All told, it might be easier to try a more sophisticated approach
##  using lme/lmer with a nested error structure
##  and multcomp for multiple comparisons.
## If it works, I might even be able to replace the nested ANOVA script
##  with this approach...


##================================================
## CUSTOM CALCULATIONS 
##================================================
SECC.full <- SECC                     # save a copy for reference
## using full moss growth data (including t3 & t4)
attr(SECC.moss.full, "labels")[["Block"]] <- attr(SECC, "labels")[["Block"]]
attr(SECC.moss.full, "labels")[["Time"]] <- attr(SECC, "labels")[["Time"]]
attr(SECC.moss.full, "labels")[["Chamber"]] <- attr(SECC, "labels")[["Chamber"]]
attr(SECC.moss.full, "labels")[["Frag"]] <- attr(SECC, "labels")[["Frag"]]
attr(SECC.moss.full, "labels")[["Position"]] <- attr(SECC, "labels")[["Position"]]

SECC <- within( SECC.moss.full, 
               {
                 grow02 <- grow01 + grow12
                 grow03 <- grow01 + grow12 + grow23
                 grow13 <- grow12 + grow23
                 Prod02 <- Prod01 + Prod12
                 Prod03 <- Prod01 + Prod12 + Prod23
                 Prod13 <- Prod12 + Prod23
               }
)

attr(SECC, "labels")[["grow02"]] <- "Moss growth 18 months"
attr(SECC, "labels")[["grow03"]] <- "Moss growth 24 months"
attr(SECC, "labels")[["grow13"]] <- "Moss growth 12-24 months"
attr(SECC, "units" )[["grow02"]] <- attr(SECC, "units" )[["grow01"]]
attr(SECC, "units" )[["grow03"]] <- attr(SECC, "units" )[["grow01"]]
attr(SECC, "units" )[["grow13"]] <- attr(SECC, "units" )[["grow01"]]
attr(SECC, "labels")[["Prod02"]] <- "Moss Productivity 18 months"
attr(SECC, "labels")[["Prod03"]] <- "Moss Productivity 24 months"
attr(SECC, "labels")[["Prod13"]] <- "Moss Productivity 12-24 months"
attr(SECC, "units" )[["Prod02"]] <- attr(SECC, "units" )[["Prod01"]]
attr(SECC, "units" )[["Prod03"]] <- attr(SECC, "units" )[["Prod01"]]
attr(SECC, "units" )[["Prod13"]] <- attr(SECC, "units" )[["Prod01"]]

## SECC <- checkSECCdata(SECC)            # too many non-standard IDs
SECC <- recodeSECC(SECC)               # drops t3
SECC <- within( SECC, 
               {
                 Time <- factor(Time, 
                                levels = c(3, levels(Time)[3]),
                                labels = c("t3", "t4")
                                )
                 Time[is.na(Time)] <- levels(Time)[1]
                 Position <- recode(Pos, 
                                    "1='Inner'; 0='Outer'; 'c'='corridor'; 'c2'='corridor'; else='other'", 
                                    levels=c( "Inner", "other", "Outer", "corridor"),
                                    as.factor.result=TRUE
                                    )
               }
)



##################################################
## CONFIGURE BASIC ANALYSIS
##################################################
## Moss growth was measured in the same t4 plots
##  throughout the experiment, so unlike destructive samples,
##  these are repeated measures.
## The simple approach is to analyze each time point separately,
##  which avoids pseudo-replication, but also does not measure the effect of Time.

### Response Variable *****
Ycols <- c('grow01', 'grow12', 'grow23', 'grow13')
Y.col <- 'Prod01'     # Column to analyze as response variable           *****
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****

### container for results at each measurement time period
Y.results <- list()
Y.fits    <- list()
Y.CxP.mcp <- list()
Y.CxP.eff <- list()
Y.effects <- list()


### loop each variable / measurement time period
for (Y.col in Ycols)
{
  cat("\n\n\n======== Processing Variable:", Y.col, "========\n")

### Load default settings (based on response variable) *****
  source("./SECCanova/SECC - ANOVA settings.R", echo = FALSE) 

  ##================================================
  ## CUSTOM SETTINGS 
  ##================================================
  ## delete lines to use the defaults.

  ## Specify which treatment levels to include (by index is probably easiest)
  Time.use     <- levels(SECC$Time)               # samples collected from t3 AND t4 (t3 not included in SECC)
  Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include
  Position.use <- levels(SECC$Position)[c(1, 3)]   # Chamber treatments to include

  ## Define Labels
  Y.units <- bquote( .(Y.units) )     # store as quote(expression)  *****

  ## Output Results?
  Save.results  <- TRUE


### Load default Labels - dependent on above settings. *****
  source("./SECCanova/SECC - ANOVA labels.R", echo = FALSE) 

  ##================================================
  ## CUSTOM LABELS
  ##================================================


##################################################
### RUN nested model
##################################################
  if (T) {                               # mixed modelling for unbalanced data with nested factors
    source("./SECCanova/SECC - nested lm.R", echo = FALSE)
  } else {                               # Nested ANOVA won't work as is with unbalanced data 
    source("./SECCanova/SECC - nested ANOVA.R", echo = FALSE)
  }

  ## Save results for final graphs
  Y.fits[[Y.col]]    <- Yp.fit
  Y.CxP.mcp[[Y.col]] <- CxP.mcp
  Y.CxP.eff[[Y.col]] <- CxP.eff
  Y.effects[[Y.col]] <- Yp.effects

  ##   print( effect("Chamber:Position", Yp.fit) )
  print( effect("Chamber:Position", Y.fits[[Y.col]]) ) # works fine here

}



###===============================================
### Include Time as a factor in nested ANOVA
###===============================================
## Repeated measures of the same samples (patches) across time...


###===============================================
### Prepare results for graphing
###===============================================
## effect("Chamber:Position", Y.fits[["Prod01"]]) # Why the FUCK does this throw an error? Because Yp.fixed has changed since the object was saved (call).  GRRR :@
effect2df <- function(eff, column = "effect", 
                         response.var = "Response", fac1.label="fac1", fac2.label = "fac2")
{
  fac1.label <- names(eff$variables)[1]
  fac2.label <- names(eff$variables)[2]
  effw <- as.data.frame( summary(eff)[[column]] )
  fac2 <- colnames(effw)
  EffCols   <- paste(response.var, colnames(effw), sep=".")
  colnames(effw) <- EffCols
  effw[[fac1.label]] <- rownames(effw)
  effw$Var  <- response.var
  effl <- reshape(effw, varying = list(EffCols), direction = "long")
  effl[[fac2.label]] <- factor(effl$time, labels = fac2)
  effl$Response <- effl[, EffCols[1]]
  effl <- effl[, c("Var", fac1.label, fac2.label, "Response")]
  effl
}
effect2df(Y.CxP.eff[[1]], "effect", names(Y.CxP.eff)[1])

effect.to.df <- function(eff)
{
  vars <- names(eff$variables) 
  eff.df <- cbind(eff$x, effect = eff$fit, lower = eff$lower, upper = eff$upper)
  eff.df
}

CxP.data <- cbind(Var = names(Y.CxP.eff)[[1]], effect.to.df(Y.CxP.eff[[1]]) )
CxP.data <- CxP.data[NULL, ]           # basic structure, but empty
for(Var in names(Y.CxP.eff) )
{
  CxP.data <- rbind(CxP.data,
                    cbind( Var = Var, effect.to.df( Y.CxP.eff[[Var]] ) )
  )
}
rownames(CxP.data) <- NULL
CxP.labels <- attr(SECC, "labels")[which(names(attr(SECC, "labels")) %in% 
                                         unique(CxP.data$Var) )]
CxP.TimeP <- substring(CxP.labels, gregexpr("[0-9-]+ months", CxP.labels) )
Months <- c("August", "June", "August")
Months <- c(Months, rep( "", len = length(CxP.TimeP) - length(Months) ) )
CxP.data$Panel <- factor(CxP.data$Var, labels = paste(Months, CxP.TimeP, sep="\n") )
CxP.data$Chamber <- factor(CxP.data$Chamber, labels = c("Ambient", "Chamber"))



##################################################
### PUBLICATION GRAPHS
##################################################
Y.lim <- range(c(CxP.data$lower, CxP.data$upper))
Y.labels  <- gsub("\\s+[0-9-]+ months", "", CxP.labels) 
Y.plotlab <- bquote( .(Y.labels[1]) * " (" * .(Y.units) * ")" )
Plot.Title <- bquote("Patch means " %+-% "95% Confidence Intervals")
Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = c("Inner", "other", "Outer") )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
PositionPts  <- ggPts.SECC(Position.map, Position.label) 


CxP.plot <- ggplot(data = CxP.data, 
                   aes(x = Chamber, y = effect, 
                       group = Position, colour = Position, shape = Position)
                   ) + ylim(Y.lim) +
                    labs(x = attr(SECC, "labels")[["Chamber"]], y = Y.plotlab) +
                    opts(title = Plot.Title) +
                    geom_point(aes(group = Position), size = 3) +
                    geom_line(aes(group = Position), size = 0.8) +
                    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5)

CxP.plot <- CxP.plot + PositionPts +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label)
CxP.plot <- CxP.plot + facet_grid(.~ Panel) + jaw.ggplot()
print(CxP.plot)


if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  Save.final <- paste(SaveDir.plots(), "Figure - ", Y.labels[1], sep="")
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 3, scale = 1.5)
}
