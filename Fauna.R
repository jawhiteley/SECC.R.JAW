##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Analyses of Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2011-08-22
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


##################################################
## DATA EXPLORATION
##################################################

plot(SECC.fauna.sum[, which( sapply(SECC.fauna.sum, is.numeric) )])



##################################################
## CHECK & PROCESS DATA
##################################################
## Aggregate count data by morphospecies with "confidence" 
## (i.e. lumping things together that are probably the same)
SECC.fauna.sp <- data.frame(SampleID = SECC.fauna[['SampleID']])
Taxa.groups <- rev( unique(SECC.fauna.meta$sp_alias) )

SECC.fauna.sp <- within(SECC.fauna.sp, {
    for (taxa in Taxa.groups) {
      taxa.sp <- SECC.fauna.meta$ID[which(SECC.fauna.meta$sp_alias == taxa)]
      assign( taxa, 
             if (length(taxa.sp) > 1) 
               apply(SECC.fauna[, taxa.sp], 1, sum) 
             else 
               SECC.fauna[, taxa.sp] 
      )
    }
    rm(taxa, taxa.sp)
})
  
Fauna <- SECC.fauna.sp  # more convenient name :P



##################################################
## NMDS
##################################################








##################################################
### PUBLICATION GRAPHS
##################################################



##################################################
### CLEAN-UP / HOUSEKEEPING
##################################################
