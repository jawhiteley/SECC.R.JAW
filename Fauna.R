##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Analyses of Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2012-04-03
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
## SETTINGS
##################################################

## t4 samples only: I'm not really using the other data yet anyway.
Samples.fauna <- which(SECC.fauna$Time == 4)
Vars.fauna    <- setdiff(colnames(SECC.fauna), 
                         SECC.fauna.meta$ID[SECC.fauna.meta$Taxonomic.Group 
                                            %in% c("Other", "Prostigmata")]
)
Spp.fauna <- intersect(SECC.fauna.meta$ID, Vars.fauna)
Spp.cols  <- intersect(colnames(SECC.fauna), SECC.fauna.meta$ID)


##################################################
## CHECK & PROCESS DATA
##################################################
SECC.fauna.full <- SECC.fauna   # save a copy, just in case
## Filter data early on to maintain correspondence in row numbers 
##   between the various data frames.  
Fauna.sp <- SECC.fauna[Samples.fauna, Vars.fauna]
Fauna.meta <- SECC.fauna.meta[which(SECC.fauna.meta$ID %in% Spp.fauna), ]


##================================================
## Aggregate 'morphospecies' records into 'species'
##================================================
## Aggregate count data by morphospecies with "confidence" 
## (i.e. lumping things together that are probably the same)
SECC.sp <- Fauna.sp[, c('SampleID', Spp.fauna)]

if (F)
{                                      # now handled when loaded & data sets merged (JAW + ZL)
  SECC.sp <- data.frame(SampleID = Fauna.sp[['SampleID']])

  Taxa.groups <- rev( unique(Fauna.meta$sp_alias) )
  SECC.sp <- within(SECC.sp, 
                    {
                      for (taxa in Taxa.groups) {
                        taxa.sp <- Fauna.meta$ID[which(Fauna.meta$sp_alias == taxa)]
                        taxa.sp <- intersect(taxa.sp, colnames(Fauna.sp))
                        if (length(taxa.sp) > 0) {
                          assign( taxa, 
                                 if (length(taxa.sp) > 1) 
                                   apply(Fauna.sp[, taxa.sp], 1, sum) 
                                 else 
                                   Fauna.sp[, taxa.sp] 
                          )
                        }
                      }
                      rm(taxa, taxa.sp)
                    })
}

rownames(SECC.sp) <- SECC.sp$SampleID
SECC.sp <- SECC.sp[, -1]    # Data matrix only: remove SampleID column, leave IDs in rownames
  
## Data matrix only: more convenient name :P
Fauna <- SECC.sp  


##================================================
## Summary variables for univariate analyses
## suitable for merging into main SECC dataframe
##================================================
## Already calculated in`/lib/load_fauna.R`
SECC.sp.sum <- SECC.fauna.sum[Samples.fauna, ]

##================================================
## Species Richness & Diversity metrics
##================================================
## diversity() and specpool() in vegan
SECC.sp.sum <- within(SECC.sp.sum, {
    Richness <- apply(SECC.sp, 1, function(x) length(which(x>0)) )  # observed # spp.
    Evenness <- diversity(SECC.sp, index = "invsimpson")
})

attr(SECC.sp.sum, "labels") <- attr(SECC, "labels")        # starting point
attr(SECC.sp.sum, "units" ) <- attr(SECC, "units" )        # starting point
attr(SECC.sp.sum, "labels")[["Richness"]] <- "Morphospecies Richness"
attr(SECC.sp.sum, "units" )[["Richness"]] <- quote("(# observed taxa)")
attr(SECC.sp.sum, "labels")[["Evenness"]] <- "Inverse Simpson's Index"
attr(SECC.sp.sum, "units" )[["Evenness"]] <- quote("(# effective taxa)")
attr(SECC.sp.sum, "labels")[["Predators"]] <- "Predators (Mesostigmata + Prostigmata)"
attr(SECC.sp.sum, "units" )[["Predators"]] <- quote("#" %.% g^-1)
attr(SECC.sp.sum, "labels")[["Grazers"]]   <- "Grazers (fungivorous Acari + Collembola)"
attr(SECC.sp.sum, "units" )[["Grazers"]]   <- quote("#" %.% g^-1)



##################################################
## DATA EXPLORATION
##################################################
if (FALSE) 
{

  library(mgcv)                        # for GAMs
  Fauna.bip <- qplot(Grazers, Predators, data = SECC.fauna.sum,
                     ) + stat_smooth(aes(x = Grazers), method = "gam") + jaw.ggplot() # group = Chamber, shape = Pos, colour = Chamber
  print(Fauna.bip)
  Fauna.bip <- Fauna.bip + facet_grid(facets = Chamber~Pos)
  print(Fauna.bip)

  ## check variation & transformations
  boxplot( Fauna            , main = "raw data")
  boxplot( sqrt( Fauna )    , main = "sqrt")
  boxplot( Fauna ^0.25      , main = "4th rt")
  boxplot( log( Fauna +1 )  , main = "log_e(X +1)")
  boxplot( log( Fauna +1 , base = 2)  , main = "log_2(X +1)")
  boxplot( decostand( Fauna, method = 'log' )         , main = "log (decostand)") # ?? wtf?
  boxplot( decostand( Fauna, method = 'normalize' )   , main = "normalized") 

  with(SECC.sp.sum, plotMeans(Richness, Chamber, Pos, error.bars="conf.int", level=0.95) )
  with(SECC.sp.sum, plotMeans(Richness, Chamber, Frag, error.bars="conf.int", level=0.95) )

}


##################################################
### UNIVARIATE ANALYSES (ANOVAs)
##################################################

Y.col <- 'Richness' # Column to analyze as response variable           *****
Y.use <- 'Y'        # Which transformation is being used (for labels)? *****
Y.lim1 <- c(0, 20)  # consistent Y limits :/
source('Fauna-univariate.R')

Y.col <- 'Evenness'
Y.use <- 'Y'
Y.lim1 <- c(0, 10)  # consistent Y limits :/
source('Fauna-univariate.R')

Y.col <- 'Predators'
Y.use <- 'Y'
Y.lim1 <- c(0, 25)  # consistent Y limits :/
source('Fauna-univariate.R')

Y.col <- 'Grazers'
Y.use <- 'Y'
Y.lim1 <- c(0, 25)  # consistent Y limits :/
source('Fauna-univariate.R')




##################################################
### MULTIVARIATE ANALYSES
##################################################

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
## ANOSIM
##################################################
## Distance/Similarity matrix (included in MDS call)
Fauna.dist <- vegdist( Fauna.trans, method = "bray")

ANOSIM.results <- capture.output(      # to make saving easier @ end
## ANOSIM by single factor
  Fauna.Chamber.anosim <- anosim( Fauna.dist, Fauna.sp$Chamber )
, Fauna.Chamber.anosim # Look at results

,Fauna.Pos.anosim <- anosim( Fauna.dist, Fauna.sp$Pos )
,Fauna.Pos.anosim # Look at results

,Fauna.Frag.anosim <- anosim( Fauna.dist, Fauna.sp$Frag )
,Fauna.Frag.anosim # Look at results

## ANOSIM within main factors
,for (lvl in levels(Fauna.sp$Chamber)) {
  cat("Chamber:", lvl, "\n")
  Fauna.data <- Fauna.trans[which(Fauna.sp$Chamber == lvl), ]
  Fauna.d <- vegdist(Fauna.data, method = "bray")
  ## Chamber x Pos
  cat("Chamber:", lvl, "x Pos",  "\n")
  Fauna.CxP.anosim <- anosim( Fauna.d, Fauna.sp$Pos )
  print(Fauna.CxP.anosim) # Look at results
  ## Chamber x Frag
  cat("Chamber:", lvl, "x Frag",  "\n")
  Fauna.CxF.anosim <- anosim( Fauna.d, Fauna.sp$Frag )
  print(Fauna.CxF.anosim) # Look at results
}

,for (lvl in levels(Fauna.sp$Frag)) {
  cat("Frag:", lvl, "\n")
  Fauna.data <- Fauna.trans[which(Fauna.sp$Frag == lvl), ]
  Fauna.d <- vegdist(Fauna.data, method = "bray")
  ## Frag x Pos
  cat("Frag:", lvl, "x Pos",  "\n")
  Fauna.FxP.anosim <- anosim( Fauna.d, Fauna.sp$Pos )
  print(Fauna.FxP.anosim) # Look at results
  ## Frag x Chamber
  cat("Frag:", lvl, "x Chamber",  "\n")
  Fauna.FxC.anosim <- anosim( Fauna.d, Fauna.sp$Chamber )
  print(Fauna.FxC.anosim) # Look at results
}
)
cat( paste(ANOSIM.results, collapse = "\n") )


################################################################
### CONSTRAINED ORDINATION: FITTING ENVIRONMENTAL VARIABLES
################################################################
## Extract matching ordination rows from Environmental data, *in same order*
Fauna.sort <- match( rownames(Fauna.trans), SECC$SampleID )
Fauna.sort <- na.omit(Fauna.sort)
Fauna.env <- SECC[Fauna.sort, ]
## Process standardized environmental variables
Fauna.env <- within(Fauna.env,
                    {
                      Stigonema.g <- Stigonema / Cells.dwt
                      Nostoc.g <- Nostoc / Cells.dwt
                      Other.cells.g <- Other.cells / Cells.dwt
                      Trt <- paste(Chamber, Frag, Position, sep="\n")
                    })

## a posteriori fitting of environmental variables to nMDS (not really constrained ordination)
## Samples               : Same experimental subset, or ALL available?
## Response Variables    : Fauna species [SECC.fauna]
## Explanatory Variables : H2O, Temperature?, Patch.dwt ; Cells.g, Stigonema.g, Nostoc.g, ARA.g ; Decomposition, Productivity (all time periods), NH4, NO3, TAN
##    Experimental Factors (if subset of samples)?  <- partial these out?

Fauna.fit <- envfit(Fauna.mds ~ Block + Chamber + Frag + Position + 
                    H2O + Patch.dwt + Cells.g + Stigonema.g + Nostoc.g + ARA.g + 
                    Decomposition + Prod01 + Prod12 + Prod23 + 
                    NH4 + NO3 + TAN
                    , data = Fauna.env, perm = 1000, na.rm = TRUE)
Fauna.fit

## Add environmental variables to ordination
ordiplot( Fauna.mds, type="points", display = "sites")
ordispider(Fauna.mds, Fauna.env$Trt, label = TRUE, col = "red")
plot(Fauna.fit, p.max = 0.025, adj = 0.5, srt = 0, cex = 1)
ordiplot( Fauna.mds, type="text", display = "species", xlim = c(-1, 2), ylim = c(-1, 2) )
plot(Fauna.fit, p.max = 0.05, adj = 0.5, srt = 0, cex = 1)



##==============================================================
## CCA
##==============================================================
## Missing env values (productivity) 
Env.rows     <- apply(Fauna.env, 1, function (x) !any(is.na(x)) )
Fauna.envs   <- Fauna.env[which(Env.rows), ]
Fauna.subset <- Fauna.trans[which(Env.rows), ]

Fauna.cca <- cca(Fauna.subset ~ Block + Chamber + Frag + Position + 
                 H2O + Patch.dwt + Cells.g + Stigonema.g + Nostoc.g + ARA.g + 
                 Decomposition + Prod01 + Prod12 + Prod23 + 
                 NH4 + NO3 + TAN 
                 , data = Fauna.envs)
Fauna.cca
summary(Fauna.cca)
## First 2 axes don't seem to capture much
## the nMDS seems to be a better overall fit than the Gaussian model of CCA?
## Permutation-based tests: SLOW
CCA.anova.terms <- anova(Fauna.cca, by = "term",   step = 200) # Permutation-based tests
CCA.anova.mar   <- anova(Fauna.cca, by = "margin", step = 99 ) # *marginal* effects of each term (in a model with all other terms)
## Plots
ordiplot(Fauna.cca, type="text")
ordiplot(Fauna.cca, type="text", 
         display = c("species", "bp"), scaling = 2, 
         xlim = c(-2, 2), ylim = c(-2, 2)) # zoom in, but still messy.
refCol <- "#AAAAAA" 
refLty <- 3
rect(-0.5, -0.5, 0.5, 0.5, col = NA, border = refCol, lty = refLty)

## without experimental factors
Fauna.cca1 <- cca(Fauna.subset ~  
                 H2O + Patch.dwt + Cells.g + Stigonema.g + Nostoc.g + ARA.g + 
                 Decomposition + Prod01 + Prod12 + Prod23 + 
                 NH4 + NO3 + TAN 
                 , data = Fauna.envs)
Fauna.cca1
summary(Fauna.cca1)
ordiplot(Fauna.cca1, type="text")
ordiplot(Fauna.cca1, type="text", 
         display = c("species", "bp"), scaling = 2, 
         xlim = c(-2, 2), ylim = c(-2, 2)) # zoom in, but still messy.


##================================================
## Correlations
##================================================
library(plyr)

cols.numeric <- sapply(SECC.sp.sum, function (x) class(x) == "numeric")
Sum.cor <- cor(SECC.sp.sum[, cols.numeric], method = "pearson")
Spp.cor <- cor(SECC.sp)
summary(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE) # lower triangle of the correlation matrix (lower.tri) :D
Spp.meancor       <- mean(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE)
Spp.meancor.sd    <- sd(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE)
Spp.meancor.error <- sd(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE) / 
                    sqrt( length( which(!is.na(Spp.cor[lower.tri(Spp.cor)])) ) )
Spp.cortest       <- t.test( Spp.cor[lower.tri(Spp.cor)] )
## Most species are positively correlated :P
library(MASS)
truehist(Spp.cor[lower.tri(Spp.cor)], xlim = c(-1, 1), col = "#999999",
         xlab = "Species correlations", ylab = "frequency")
abline(v = 0, lty = "dashed", lwd = 3, col = "#333333")

## Correlations *within groups* (Predators, Grazers, etc.)?
## data.frame with:  column of correlation values, column of Group (including "Between Groups")
Spp.cort  <- Spp.cor
Spp.cort[upper.tri(Spp.cort, diag = TRUE)] <- NA
Spp.cordf <- reshape(data = as.data.frame(Spp.cort), 
                     varying = list(colnames(Spp.cor)),
                     times = colnames(Spp.cor),
                     timevar = "Species2",
                     idvar = "Species1",
                     ids = rownames(Spp.cor),
                     v.names = "cor",
                     direction = "long"
                     )
Spp.cordf <- subset(Spp.cordf, !is.na(cor) )
Spp.cordf <- within(Spp.cordf[, c("cor", "Species1", "Species2")], {
                    Group1 <- match(Species1, SECC.fauna.meta$ID)
                    Group1 <- SECC.fauna.meta$Taxonomic.Group[Group1]
                    Taxa <- match(Species1, SECC.fauna.meta$ID)
                    Taxa <- SECC.fauna.meta$Major.Taxa[Taxa]
                    Group1 <- ifelse(Taxa == "Uropodina", Taxa, Group1)
                    Group2 <- match(Species2, SECC.fauna.meta$ID)
                    Group2 <- SECC.fauna.meta$Taxonomic.Group[Group2]
                    Taxa <- match(Species2, SECC.fauna.meta$ID)
                    Taxa <- SECC.fauna.meta$Major.Taxa[Taxa]
                    Group2 <- ifelse(Taxa == "Uropodina", Taxa, Group2)
                    ## 1 column of between Group labels
                    Group <- ifelse(Group1 == Group2, Group1, "Between Groups")
                    Group <- ifelse(Group %in% c("Mesostigmata", "Prostigmata"), "Predators", Group)
                    Group <- ifelse(Group %in% c("Collembola", "Uropodina"), "Grazers", Group)
                    rm(Taxa)
                     })

## Predators are essentially uncorrelated (~0), Grazers & Between Groups are slightly +vely correlated (~0.05)
Graz.cortest <- t.test( subset(Spp.cordf, Group == "Grazers", select = "cor") )
Pred.cortest <- t.test( subset(Spp.cordf, Group == "Predators", select = "cor") )
Btwn.cortest <- t.test( subset(Spp.cordf, Group == "Between Groups", select = "cor") )
Spp.cortest
Graz.cortest
Pred.cortest
Btwn.cortest

## Stacked bar graph :D
Spp.corplot <- ggplot( Spp.cordf, aes(x = cor, fill = Group)) + 
    xlim(-1, 1) + xlab("Species correlation") + ylab("frequency") +
    geom_bar(binwidth = 0.1, colour = "#333333", lwd = 0.2) + 
    scale_fill_manual(values = c("#CCCCCC", "#999999", "#444444"), 
                      breaks = c("Predators", "Grazers", "Between Groups")) +
    geom_vline(xintercept = 0, colour = "#666666", lwd = 1, lty = "dashed") +
    jaw.ggplot()
print(Spp.corplot)


PredGraz.cor    <- cor(SECC.sp.sum$Predators, SECC.sp.sum$Grazers)
PredGraz.cortxt <- sprintf("r = %.3f", PredGraz.cor)

PredGraz.corplot <- ggplot(SECC.sp.sum, 
                           aes(x = Grazers, y = Predators)
                           ) + 
                 xlab(bquote("Grazers (" * .( attr(SECC.sp.sum, "units")[["Grazers"]]) * ")")) +
                 ylab(bquote("Predators (" * .( attr(SECC.sp.sum, "units")[["Predators"]]) * ")")) +
                 geom_point(aes(shape = Pos, colour = Chamber), size = 2.5) + 
                 scale_shape_manual(name = Position.label,
                                    values = Position.map$pch, 
                                    breaks = Position.map$label,
                                    labels = c("Inner (Wet)", "Outer (Dry)")
                                    ) +
                 scale_colour_manual(name = Chamber.label,
                                     values = Chamber.map$col, 
                                     breaks = Chamber.map$label,
                                     labels = c("Ambient", "Chamber")
                                     ) +
                 stat_quantile(quantiles = 0.5, colour = "#666666") +
                 geom_text(aes(x = max(Grazers), 
                               y = max(Predators),
                               label = PredGraz.cortxt
                               ), hjust = 1, vjust = 1) +
                 jaw.ggplot() + coord_equal()

print(PredGraz.corplot)


##================================================
## Composition
##================================================

SECC.sp.CxP <- checkSECCdata(SECC.sp.sum, 'SECC.sp.sum')
SECC.sp.CxP <- recodeSECC( SECC.sp.CxP )
SECC.sp.CxP <- reshape(data = SECC.sp.CxP, 
                       idvar = colnames(SECC.sp.sum)[!cols.numeric], 
                       varying = list(c("Predators", "Grazers")),
                       times = c("Predators", "Grazers"),
                       v.names = "Abundance", timevar = "Group",
                       direction = "long"
                       )
SECC.sp.CxP <- aggregate(Abundance ~ Chamber + Position + Group, 
                         data = SECC.sp.CxP, 
                         FUN = mean)

Compo.plot <- ggplot(SECC.sp.CxP, aes(x = Position, y = Abundance, fill = Group)) + 
    geom_bar(stat="identity", position = "fill", colour = NA) +
    scale_fill_manual(values = c("Predators" = "#444444", "Grazers" = "#AAAAAA"),
                      breaks = c("Predators", "Grazers")) +
    facet_wrap(~ Chamber) + 
    jaw.ggplot()

print(Compo.plot)




##################################################
### PUBLICATION GRAPHS
##################################################
## Prepare data to plot
stress.label <- paste( "Stress\n= ", round( Fauna.mds$stress, 1)/100 )	# Make text label of Stress value, rounded to 1 decimal place, and converted from % to a decimal number.
MDS.pts <- as.data.frame(Fauna.mds$points)
colnames(MDS.pts) <- c('x', 'y')
MDS.pts <- cbind(MDS.pts, Time = Fauna.sp$Time, Chamber = Fauna.sp$Chamber, 
                 Frag = Fauna.sp$Frag, Pos = Fauna.sp$Pos
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
                    xlab = NULL, ylab = NULL
                    )   # facets=Frag~Time
Fauna.plot <- Fauna.plot + scale_shape_manual(name = Position.label,
                                              values = Position.map$pch, 
                                              breaks = Position.map$label,
                                              labels = c("Inner (Wet)", "Outer (Dry)")
                                              )
Fauna.plot <- Fauna.plot + scale_colour_manual(name = Chamber.label,
                                               values = Chamber.map$col, 
                                               breaks = Chamber.map$label,
                                               labels = c("Ambient", "Chamber")
                                               )
Fauna.plot <- Fauna.plot + scale_x_continuous(breaks = -4:4, limits = c(-3.7, 1.3))
Fauna.plot <- Fauna.plot + scale_y_continuous(breaks = -4:4, limits = c(-2.7, 2.3))
## options: see theme_get() for available theme options.
Fauna.plot <- Fauna.plot + theme_bw() + coord_equal()
Fauna.plot <- Fauna.plot + opts(axis.ticks = theme_blank(), 
                                axis.text.x = theme_blank(), 
                                axis.text.y = theme_blank(), 
                                legend.key = theme_rect(colour = NA),
                                legend.position = "right",
                                legend.direction = "vertical"
                                )
Fauna.plot <- Fauna.plot + geom_text(aes(label = stress.label, x = max(x)*0.5, y = max(y) ), 
                                     size = I(4), colour = "black", hjust = 0, vjust = 1)
print(Fauna.plot)
grid.text(label = stress.label, x = 0.85, y = 0.85, just = "left", hjust = 0, vjust = 1)

## mtext( stress.label , side=3, line=0, adj=1 ) # kludge: trying to copy or save the graph now will make R crash.

Fauna.plot.frag <- Fauna.plot + facet_grid(facets = Frag~.)
print(Fauna.plot.frag)




##================================================
## Save Multivariate Results
##================================================
if (TRUE) {
  cat(  paste(ANOSIM.results, collapse = "\n")
      , "\n\n", "======== Environmental Variables a posteriori fit to nMDS ========", "\n"
      , paste(capture.output(print(Fauna.fit)), collape="\n")
      , "\n\n", "======== CCA ========", "\n"
      , paste(capture.output(summary(Fauna.cca)), collape="\n")
      , "\n\n", "======== Permutation tests of CCA terms ========", "\n"
      , paste(capture.output(print(CCA.anova.terms)), collape="\n")
      , file = "./output/Fauna - Multivariate.txt"
  )
  ggsave(filename="./graphs/Figure - Fauna-MDS.eps", plot = Fauna.plot, width = 6, height = 6)
}



##################################################
### CLEAN-UP / HOUSEKEEPING
##################################################
SECC.fauna <- SECC.fauna.full   # recover full data frame
