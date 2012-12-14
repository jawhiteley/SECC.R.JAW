##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### basic analyses of experimental data
### Moss Growth / Productivity
### Jonathan Whiteley     R v2.12     2012-07-05
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

## Moss growth is a little different than other measured variables:
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
                                    "1='Inner'; 0='Outer'; 'c'='corridor'; 'c2'='corridor'; else='intermediate'", # **** non-standard!
                                    levels=c( "Inner", "intermediate", "Outer", "corridor"),
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
Y.col <- 'grow01'     # Column to analyze as response variable
Y.use <- 'Y'    # Which transformation is being used (for labels)? *****

### container for results at each measurement time period
Y.results <- list()
Y.fits    <- list()
Y.CxP.mcp <- list()
Y.CxP.eff <- list()
Y.effects <- list()
Time.eff  <- list()


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
  ## Including "intermediate" patches throws errors when fitting models (for grow01)?
  ## Include Partial Chambers throws errors when doing post-hoc multiple tests (grow01)?
  Time.use     <- levels(SECC$Time)               # samples collected from t3 AND t4 (t3 not normally included in SECC)
  Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include
  Position.use <- levels(SECC$Position)[c(1, 3)]  # Position treatments to include

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
  if (length(levels(SECCp$Time)) > 1)
  {
    Time.eff[[Y.col]] <- effect("Time:Chamber:Position", Yp.fit)
  }

  ##   print( effect("Chamber:Position", Yp.fit) )
  print( effect("Chamber:Position", Y.fits[[Y.col]]) ) # works fine here

}
## Check ANOVA outputs with backward-model selection ***
## Borderline Fragmentation effects in ts (summer 2)?


##================================================
## Prepare results for graphing
##================================================
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
Season <- c("Year 1", "Winter 2", "Summer 2", "Year 2")
Season <- c(Season, rep( "", len = length(CxP.TimeP) - length(Season) ) )    # fill panels with blank months?
CxP.data$Panel <- factor(CxP.data$Var, labels = paste(Season, CxP.TimeP, sep="\n") )
CxP.data$Chamber <- factor(CxP.data$Chamber, labels = c("Ambient", "Chamber"))


##================================================
## EXPLORATION: effect of "Time"?
##================================================
if (F)
{                                      # Hide prompts when source()ing
  ## Remember, Time is really just replication within blocks: there should be no effect.
  plot(Y.effects[["grow12"]])
  print(Time.eff[["grow12"]])
  plot(Time.eff[["grow12"]], multiline = TRUE)
  plot(Time.eff[["grow12"]], x.var = "Chamber", z.var = "Time")
  plot(Time.eff[["grow23"]], x.var = "Chamber", z.var = "Time")
  plot(Time.eff[["grow13"]], x.var = "Chamber", z.var = "Time")
}





##################################################
### PUBLICATION GRAPHS
##################################################
Y.lim <- range(c(CxP.data$lower, CxP.data$upper))
Y.labels  <- gsub("\\s+[0-9-]+ months", "", CxP.labels) 
Y.plotlab <- bquote( .(Y.labels[1]) * " (" * .(Y.units) * ")" )
Plot.Title <- bquote("Patch means " %+-% "95% Confidence Intervals")
Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC.full$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
PositionPts  <- ggPts.SECC(Position.map, Position.label) 


CxP.plot <- ggplot(data = CxP.data, 
                   aes(x = Chamber, y = effect, 
                       group = Position, colour = Position, fill = Position, size = Position, shape = Position)
                   ) + ylim(Y.lim) +
                    labs(x = attr(SECC, "labels")[["Chamber"]], y = Y.plotlab) +
                    opts(title = Plot.Title) +
                    geom_line(aes(group = Position, size = Position)) +
                    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
                    geom_point(aes(group = Position), size = 3)

CxP.plot <- CxP.plot + PositionPts +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label) +
scale_size_manual(name = Position.label,
                  values = Position.map$lwd, 
                  breaks = Position.map$label)
CxP.plot <- CxP.plot + facet_grid(.~ Panel) + jaw.ggplot()
print(CxP.plot)









##################################################
### Include Time as a factor in nested model (ANOVA)
##################################################
## Repeated measures of the same samples (patches) across years
cat("\n\n\n======== Nested analysis of Moss Productivity between Years ========\n")


if (F)
{ ## Specify which treatment levels to include (optional: hide to use same as previous analyses)
  Time.use     <- levels(SECC$Time)               # samples collected from t3 AND t4 (t3 not normally included in SECC)
  Chamber.use  <- levels(SECC$Chamber)[c(1, 3)]   # Chamber treatments to include (No Productivity measurements from Partial; no Moss_Biomass data for these treatments).
  Position.use <- levels(SECC$Position)[c(1:3)]  # Position treatments to include: intermediate patches normally excluded for consistency, but they are kind of interesting here ...
}


## Reshape columns of variables in a single column with a Year factor
SECC.grow <- reshape(SECC, idvar = "SampleID", varying = list(c("grow01", "grow13")), v.names = "grow",
                     timevar = "Year", times = c("grow01", "grow13"), direction = "long",
                     drop = setdiff(names(SECC), c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos', 'Position','grow01', 'grow13'))
                     )
SECC.prod <- reshape(SECC, idvar = "SampleID", varying = list(c("Prod01", "Prod13")),  v.names = "Prod",
                     timevar = "Year", times = c("Prod01", "Prod13"), direction = "long",
                     drop = setdiff(names(SECC), c('SampleID', 'Block', 'Time', 'Chamber', 'Frag', 'Pos', 'Position','Prod01', 'Prod13'))
                     )
if( all(SECC.grow$SampleID == SECC.prod$SampleID) ) SECC.prod <- cbind(SECC.prod, grow = SECC.grow[, "grow"])
SECC.prod$Year <- factor(SECC.prod$Year, levels = c("Prod01", "Prod13"), labels = c(1, 2))
rownames(SECC.prod) <- gsub("Prod01", 1, rownames(SECC.prod))
rownames(SECC.prod) <- gsub("Prod13", 2, rownames(SECC.prod))

## Filter data for analysis, according to settings of PREVIOUS ANALYSIS, above (see top of script).
SECC.prod <- SECC.prod[SECC.prod$Time     %in% Time.use     &
                       SECC.prod$Chamber  %in% Chamber.use  &
                       SECC.prod$Frag     %in% Frag.use     &
                       SECC.prod$Position %in% Position.use 
                       , ]
## drop unused factor levels (for plotting)
SECC.prod <- within( SECC.prod, 
                    {
                      Time     <- factor(Time,     exclude=NULL)
                      Chamber  <- factor(Chamber,  exclude=NULL)
                      Frag     <- factor(Frag,     exclude=NULL)
                      Position <- factor(Position, exclude=NULL)
                    })


##================================================
## check for Block effects: I can't include it in the full model, because of lack of replication
Prod.b   <- gls(Prod ~ Block*Chamber*Position,
                data = SECC.prod[SECC.prod$Year == 2, ],  # no productivity measurements from t3 in year 1 :(
                method="REML", na.action = na.omit
                )
anova(Prod.b)
Prod.b0 <- update(Prod.b , method = "ML")
Prod.b1 <- update(Prod.b0, . ~ . -Block:Chamber:Position)
anova(Prod.b0, Prod.b1)                # barely significant

Block.plot <- ggplot(SECC.prod, aes(y = Prod, x = Chamber, group = Position)) + 
                geom_point(aes(colour = Position, shape = Position)) +
                stat_summary(fun.y = mean, geom = "line", aes(colour = Position, size = Position)) +
                facet_wrap(~ Block)
print(Block.plot)                      # maybe something funny hapenning in block 4 or 8, but the pattern is pretty consistent.


##================================================
## FIT MODEL
##================================================
## Repeated measures in time (approximated via nesting structure)
lmd <- lmeControl()                    # save defaults
lmc <- lmeControl(niterEM = 500, msMaxIter = 100, opt="optim")
## Frag:Position:Year prevents solution in many cases; is the missing data really that unbalanced?
Yp.fixed  <- Prod ~ Chamber+Frag+Position+Year + 
                Chamber:Frag + Chamber:Position + Chamber:Year + Frag:Position + Frag:Year + Position:Year +
                Chamber:Frag:Position + Chamber:Frag:Year + Chamber:Position:Year # + Frag:Position:Year + Chamber:Frag:Position:Year,
Yp.random <- ~ 1 | Block/Chamber/Frag/Position/Year

## Main Model
## I *could* include all measurements from t3, which are only present in the second year, although I worry about what kind of comparison that creates (years aren't really comparable at that rate).
## t4 measurements alone appear to be sufficient.
## the results are basically the same, except the p-values for Frag drop a bit (still nowhere near 0.05).
Prod.lme <- lme(Yp.fixed, random = Yp.random, 
                data = SECC.prod[SECC.prod$Time == "t4", ],  # no productivity measurements from t3 in year 1 :(
                method="REML", control = lmc, na.action = na.omit)
Prod.lme2 <- lme(Yp.fixed, random = ~ 1 | Block/Time/Chamber/Frag/Position/Year, 
                data = SECC.prod,  # no productivity measurements from t3 in year 1 :(
                method="REML", control = lmc, na.action = na.omit)
## control for heterogeneity between blocks & years (see data exploration graphs)
Prod.hB  <- update(Prod.lme, weights = varIdent(form= ~1 | Block))
Prod.hY  <- update(Prod.lme, weights = varIdent(form= ~1 | Year))
Prod.hBY <- update(Prod.lme, weights = varIdent(form= ~1 | Block * Year))

Prod.gls <- gls(Yp.fixed,
                data = SECC.prod[SECC.prod$Time == "t4", ],  # no productivity measurements from t3 in year 1 :(
                method="REML", na.action = na.omit
                )


##================================================
## MODEL SELECTION?
##================================================
## Random Structure
anova(Prod.hBY, Prod.hB)
anova(Prod.hBY, Prod.hY)
anova(Prod.hBY, Prod.hY, Prod.lme, Prod.gls)              # random structure makes a BIG difference :)

## controlling for heterogeneity between years improves the model fit by AIC, but the residuals are no longer normally-distributed :(
## The results don't change much otherwise ...
## Prod.lme <- Prod.hY                   # replace with best model

if (FALSE)
{   ## Backwards-selection (using likelihood-ratio tests)
  ML <- update(Prod.lme, method = "ML")
  ML1 <- update(ML, . ~ . - Chamber:Frag:Position)
  ML2 <- update(ML, . ~ . - Chamber:Frag:Year)
  ML3 <- update(ML, . ~ . - Chamber:Position:Year)
  anova(ML, ML1)
  anova(ML, ML2)
  anova(ML, ML3)                         # ***
  ML4 <- update(ML1, . ~ . - Chamber:Frag:Year)
  ML5 <- update(ML1, . ~ . - Chamber:Position:Year)
  anova(ML1, ML4)
  anova(ML1, ML5)                        # ***
  ML6 <- update(ML4, . ~ . - Frag:Year)
  ML7 <- update(ML4, . ~ . - Frag:Position)
  ML8 <- update(ML4, . ~ . - Chamber:Frag)
  anova(ML4, ML6)
  anova(ML4, ML7)
  anova(ML4, ML8)
  ML9  <- update(ML6, . ~ . - Frag:Position)
  ML10 <- update(ML6, . ~ . - Chamber:Frag)
  anova(ML6, ML9)
  anova(ML6, ML10)
  ML11 <- update(ML9, . ~ . - Chamber:Frag)
  anova(ML9, ML11)
  ML12 <- update(ML11, . ~ . - Frag)
  anova(ML11, ML12)
  ## Basically, the anova output produces largely the same results (slightly different p-values, but same terms are significant, or not).
}


##================================================
## CHECK ASSUMPTIONS: analyse residuals, standard diagnostic plots
##================================================
## independence?
    ## experimental design: random position of Frag within Chambers.
    ## possibility of Block gradient (7&8 SW -> 1-6 E:N), which also corresponds roughly to order of samples.
Save.plots <- paste("./graphs/", "Results - Productivity - 12", ".pdf", sep = "")
if (Save.results == TRUE && is.null(Save.plots) == FALSE) pdf( file = Save.plots )

if (T)
{                                      #  PROMPT for plot too annoying
  ## trellis plots: any pattern across blocks, within frag & chambers?
  xyplot(Prod ~ Block | Frag + Chamber, data=SECC.prod, 
         pch=21, col="black", bg="grey", cex=0.8,
         main = Dataset.labels[1]
         )
  Yp.residuals <- resid(Prod.lme)
  par( mfrow=c(2,2), cex=0.8) # panel of figures: 2 rows & 2 columns
  ## homogeneity of variances?
  if (T) with( SECC.prod, plot(Prod ~ Block*Chamber*Frag*Position*Year) )    # fixed effects only, no nesting
  ##* PROMPT *##
  ## plot.new()  # put next plot in empty panel?

  ## normal distribution?
  ## controlling for heterogeneity might have a better fit AIC-wise, but the residuals are way less normal
  with(SECC.prod, qqnorm( Yp.residuals, main="Residuals" ) )   # are residuals normally distributed?
  qqline(Yp.residuals,  col="grey50")
  hist(Yp.residuals)  # plot residuals
  ## with(SECC.prod, shapiro.test( Yp.residuals ) )    # normality?
  par( mfrow=c(1,1) )
}


cat("\n\n")
if( isTRUE( paste("anova", class(Prod.lme), sep=".") %in% methods(anova) ) ) 
{
  Yp.summary <- anova(Prod.lme)
} else {
  Yp.summary <- summary(Prod.lme)
}
print( Yp.summary )        # summary output
Yp.mtab <- try( model.tables(Prod.lme, "means")   # effect sizes
               , silent = TRUE)        # wrapped in try() statement, because unbalanced designs throw errors :(

# Interaction Plots
par(mfrow=c(2,2))   # panel of figures: 2 rows & 2 columns
with( SECC.prod, interaction.plot(Frag, Chamber, Prod,
                              fun = function(x) mean(x, na.rm=TRUE),
                              ylab = paste("mean of", "Productivity")
                              )
     )
with( SECC.prod, interaction.plot(Position, Chamber, Prod,
                              fun = function(x) mean(x, na.rm=TRUE),
                              ylab=paste("mean of", "Productivity")
                              )
     )
with( SECC.prod, interaction.plot(Position, Frag, Prod,
                              fun = function(x) mean(x, na.rm=TRUE),
                              ylab=paste("mean of", "Productivity")
                              )
     )

if (Save.results == TRUE && is.null(Save.plots) == FALSE) dev.off()


##==============================================================
## unplanned Multiple Comparisons using multcomp package?
##   -> See `SECC - nested lm.R` script for experiments using multcomp.
if (T)
{                                      # NOT WORKING ?! :(
  ## estimates of DIFFERENCES and confidence intervals, adjusted for multiple comparisons.  This is close to what I want. (depends on model / contrast matrix K)
  ## fake factor for Chamber x Position x Year interaction
  SECC.prod$CPY <- with(SECC.prod, interaction(Chamber, Position, Year) )
  CPY.fit <- lme(Prod ~ CPY, random = Yp.random, data = SECC.prod[SECC.prod$Time == "t4", ], na.action = na.omit) #  with or without intercept?
  ## K <- mcp(CPY="Tukey")                  # differences between each combination of factors
  CPY.mcp <- glht(CPY.fit, linfct = mcp(CPY="Tukey"))
  summary(CPY.mcp)
  CPY.ci  <- confint(CPY.mcp)
  plot(CPY.ci)
}

##   -> comparison intervals for graphical display: using `effects` package for now.
  library(effects)    # attractive interaction plots.  Does not work with aovlist.

  Yp.effects <- allEffects(Prod.lme)
  plot(effect("Chamber:Position:Year", Prod.lme))
  print( plot(effect("Chamber:Position:Year", Prod.lme), 
              multiline = TRUE, x.var = "Chamber", z.var = "Position") )
  print( plot(effect("Chamber:Frag:Position", Prod.lme),
              multiline = TRUE, x.var = "Frag", z.var = "Position") )
  CPY.eff <- effect("Chamber:Position:Year", Prod.lme)
  CPYs    <- summary(CPY.eff)
  print( plot(CPY.eff) )




##================================================
## Prepare results for graphing
##================================================
CPY.data <- effect.to.df(CPY.eff)
CPY.data$Chamber <- factor(CPY.data$Chamber, labels = c("Ambient", "Chamber"))
CPY.data$Year    <- factor(CPY.data$Year,    labels = c("Year 1", "Year 2"))



##################################################
### PUBLICATION GRAPHS
##################################################
Y.lim <- range(c(CPY.data$lower, CPY.data$upper))
Y.plotlab <- bquote( .(attr(SECC, "labels")[["Prod"]]) * " (" * .(attr(SECC, "units")[["Prod"]]) * ")" )
Plot.Title <- bquote("Patch means " %+-% "95% Confidence Intervals")
Position.label <- "Patch\nPosition" # attr(SECC, "labels")[["Pos"]]
Position.map <- plotMap( factor = "Position", labels = levels(SECC.full$Position) )
Position.map <- Position.map[ levels(SECC$Position) %in% Position.use, ]
PositionPts  <- ggPts.SECC(Position.map, Position.label) 


CPY.plot <- ggplot(data = CPY.data, 
                   aes(x = Chamber, y = effect, 
                       group = Position, colour = Position, fill = Position, shape = Position)
                   ) + ylim(Y.lim) +
                    labs(x = attr(SECC, "labels")[["Chamber"]], y = Y.plotlab) +
                    opts(title = Plot.Title) +
                    geom_line(aes(group = Position, size = Position)) +
                    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
                    geom_point(aes(group = Position), size = 3)

CPY.plot <- CPY.plot + PositionPts +
scale_fill_manual(name = Position.label,
                  values = Position.map$bg, 
                  breaks = Position.map$label) +
scale_size_manual(name = Position.label,
                  values = Position.map$lwd, 
                  breaks = Position.map$label)
CPY.plot <- CPY.plot + facet_grid(.~ Year) + jaw.ggplot()
print(CPY.plot)


##################################################
### SAVE OUTPUTS
##################################################
if (Save.results == TRUE && is.null(Save.text) == FALSE) 
{
  Save.header  <- paste(  "Nested ANOVA Results for:", attr(SECC, "labels")[["Prod"]],
                        "\nTransformation used:     ", Y.use,
                        "\nExpt. Time:   ", paste(Time.use,     collapse = ", "),
                        "\nChamber:      ", paste(Chamber.use,  collapse = ", "),
                        "\nFragmentation:", paste(Frag.use,     collapse = ", "),
                        "\nPatches:      ", paste(Position.use, collapse = ", "),
                        paste("\n\n", date(), "\n\n", Save.divider, sep = "")
                        )
  capture.output(cat(Save.header, sep=""),
				 print(Yp.summary),                # results summary
				 cat("\n\n"),                      # for output
				 xtable(Yp.summary, digits = c(0, 0, 0, 1, 3)),
				 cat("\n\n"),                      # for output
                 summary(CPY.mcp),                 # multiple comparisons of Chamber x Position interaction
                 cat("\n\n"),                      # for output
				 print(Yp.effects),                # effect sizes
				 cat(Save.end),                    # END OUTPUT #
				 file = "./output/Results - Productivity - Years 1-2.txt"
				)
}

if (Save.results == TRUE && is.null(Save.final) == FALSE) {
  Save.final <- paste(SaveDir.plots(), "Figure - Moss Growth", sep="")
  ggsave(file = paste(Save.final, "- CxP.eps"), plot = CxP.plot, width = 6, height = 3, scale = 1.2)
  ggsave(file = paste(Save.final, "- CPY.eps"), plot = CPY.plot, width = 5, height = 4, scale = 1)
}
  
