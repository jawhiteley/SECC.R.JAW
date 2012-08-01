##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Analyses of Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2012-07-31
###===============================================
### Species identified to morphospecies (usually family-level)
### Counts / sample converted to # / g dwt of moss substrate (using 'Patch.dwt' column)
##################################################
## INITIALISE
##################################################
if (FALSE) {  # do not run automatically
  ## Working Directory: see lib/init.R below
  setwd("./ SECC/")  # relative to my usual default wd in R GUI (Mac).
  setwd("..")                          # relative to this file (\rd in vim).
  getwd()  # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./lib/init.R')

## Load packages
library(vegan)
library(ggplot2)
library(xtable)


##################################################
## SETTINGS
##################################################
cat("STARTING Fauna Analysis ...\n")
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
Fauna.meta$Trophic <- ifelse(Fauna.meta$Taxonomic.Group == "Mesostigmata" | Fauna.meta$Taxonomic.Group == "Prostigmata", "Predator", "Grazer")
Fauna.meta$Trophic <- ifelse(Fauna.meta$Major.Taxa == "Uropodina",  "Grazer", Fauna.meta$Trophic)
Fauna.meta$Trophic <- ifelse(Fauna.meta$Taxonomic.Group == "Other", "Other",  Fauna.meta$Trophic)



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
## filter additional data frames to match
Fauna.xy <- SECC.xy[which(SECC.xy$SampleID %in% rownames(Fauna)), ]
SECC.env <- SECC[which(SECC$SampleID %in% rownames(Fauna)), ]
## This also keeps the same *order*, assuming the rownames all match
Fauna.xy <- SECC.xy[rownames(Fauna), ]
SECC.env <- SECC[rownames(Fauna), ]
SECC.env <- within(SECC.env,
                   {
                     Stigonema.g <- Stigonema / Cells.dwt
                     Nostoc.g <- Nostoc / Cells.dwt
                     Other.cells.g <- Other.cells / Cells.dwt
                     Prod13 <- Prod12 + Prod23
                     Trt <- paste(Chamber, Frag, Position, sep="\n")
                   })


##================================================
## Merge some tables ...
##================================================
## data frame with microarthropod and cyanobacteria sp - look for correlations?
if (all(rownames(SECC.sp) == rownames(SECC.env))) 
  Trophic.sp <- cbind(SECC.sp, SECC.env[, c("Stigonema.g", "Nostoc.g")]) # , "Cells.g"
colnames(Trophic.sp) <- gsub(".g", "", colnames(Trophic.sp), fixed = TRUE)



##================================================
## Species Richness & Diversity metrics
##================================================
cat("Calculating Diversity metrics\n")
## diversity() and specpool() in vegan
SECC.sp.sum <- within(SECC.sp.sum, {
    Richness <- apply(SECC.sp, 1, function(x) length(which(x>0)) )  # observed # spp.
    Evenness <- diversity(SECC.sp, index = "invsimpson")
})

attr(SECC.sp.sum, "labels") <- attr(SECC, "labels")        # starting point
attr(SECC.sp.sum, "units" ) <- attr(SECC, "units" )        # starting point
attr(SECC.sp.sum, "labels")[["Richness"]] <- "Morphospecies Richness"
attr(SECC.sp.sum, "units" )[["Richness"]] <- quote("# observed taxa")
attr(SECC.sp.sum, "labels")[["Evenness"]] <- "Inverse Simpson's Index"
attr(SECC.sp.sum, "units" )[["Evenness"]] <- quote("# effective taxa")
attr(SECC.sp.sum, "labels")[["Predators"]] <- "Predators (Mesostigmata + Prostigmata)"
attr(SECC.sp.sum, "units" )[["Predators"]] <- quote("#" %.% g^-1)
attr(SECC.sp.sum, "labels")[["Grazers"]]   <- "Grazers (fungivorous Acari + Collembola)"
attr(SECC.sp.sum, "units" )[["Grazers"]]   <- quote("#" %.% g^-1)

## Merge with SECC data for direct comparisons
SECC.sum <- merge(SECC.env, recodeSECC(SECC.sp.sum), all = TRUE)


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
  ##   boxplot( Fauna            , horizontal = TRUE, show.names = FALSE, main = "raw data")
  boxplot( Fauna            , main = "raw data")
  boxplot( sqrt( Fauna )    , main = "sqrt")
  boxplot( Fauna ^0.25      , main = "4th rt")
  boxplot( log( Fauna +1 )  , main = "log_e(X +1)")
  boxplot( log1p(Fauna)     , main = "log1p(X)")
  boxplot( log( Fauna +1 , base = 2)  , main = "log_2(X +1)")
  boxplot( decostand( Fauna, method = 'log' )         , main = "log (decostand)") # ?? wtf?
  boxplot( decostand( Fauna, method = 'max' )         , main = "max") 
  boxplot( decostand( Fauna, method = 'normalize' )   , main = "normalized") ## Chord transformation
  ##   boxplot( decostand( Fauna, method = 'standardize' )   , main = "standardize") # problems (warnings)
  boxplot( decostand( Fauna, method = 'hellinger' )   , main = "hellinger") 
  ## Chord transformation?

  with(SECC.sp.sum, plotMeans(Richness, Chamber, Pos, error.bars="conf.int", level=0.95) )
  with(SECC.sp.sum, plotMeans(Richness, Chamber, Frag, error.bars="conf.int", level=0.95) )

  pairplot(SECC.sp.sum[, which(sapply(SECC.sp.sum, class) == "numeric")])
  pairplot(SECC.sum[, c("Richness", "Evenness", "Predators", "Grazers", "fauna.jaw", "H2O")])

  pairplot(Trophic.sp)  # way too much going on: nothing seems strongly correlated with either cyanobacteria group
  
  ## Species Occurences
  sp.pr <- apply(Fauna > 0, 2, sum)
  sort(sp.pr)
  sp.refr <- 100*sp.pr/nrow(Fauna)
  par(mfrow = c(1, 2))
  hist(sp.pr, right = FALSE, main = "Species Occurences", xlab = "Number of occurences", ylab = "Number of species", col = "grey", border = NA)
  hist(sp.refr, right = FALSE, main = "Species Relative Frequencies", xlab = "Frequencey of occurences (%)", ylab = "Number of species", col = "grey", border = NA)
  par(mfrow = c(1,1))
  
  ## 1 single patch with no Predators!
  hist(SECC.sum$Predators, color = "grey")
  hist(SECC.sum$Grazers, color = "grey")
  sum(SECC.sum$Predators == 0)
  SECC.sum$SampleID[ which(SECC.sum$Predators == 0) ] # 14C-2.0
  sum(SECC.sum$Grazers == 0)

  ## Blocks 1 & 7 are kind of diversity hotspots
  qplot(xE, yN, data = Fauna.xy, size = SECC.sp.sum$Richness, shape = 20, alpha = 0.1 )
  with(SECC.sum, plot(Richness ~ Block) )
  with(SECC.sum, plot(Evenness ~ Block) )
  ## although it looks like drying may have been less severe in block 3 (remember, it got flooded during spring); fewer low-richness patches

  plot(specaccum(Fauna), ci.type = "polygon", ci.col = "lightgrey")

  plot(radfit(Fauna[1:12, ]))

  ## BD-EF relationships?
  library(mgcv)
  qplot(x = Richness, y = Decomposition, data = SECC.sum) + stat_smooth(method = "gam")
  qplot(x = H2O, y = Richness, data = SECC.sum) + stat_smooth(method = "gam")
  qplot(x = H2O, y = Decomposition, data = SECC.sum) + stat_smooth(method = "gam")
}


