##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Process Fauna data (microarthropod morphospecies counts) 
### for multivariate analyses with other variables
### Based on Fauna.R
### Jonathan Whiteley     R v2.12     2012-07-11
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
## source('./lib/init.R')

## Load packages
library(vegan)

##################################################
## SETTINGS
##################################################
cat("Processing Fauna data ...\n")
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

