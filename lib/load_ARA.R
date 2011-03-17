##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Load, clean & process ARA data loaded from "./data/"
# Jonathan Whiteley		R v2.12		2011-03-15
##################################################
## This script is run as part of `./lib/load.R`
require(car)		# load external package 'car', for recode()

##================================================
## CHECK DATA
##================================================
# str(SECC.ARA.t1)  # Should already be in memory.    str() produces output on source().

SECC.ARA.t1 <- SECCcolumnNames(SECC.ARA.t1)  # Standardize ID column names (but not values)

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
# Manually clean & prepare data for automatic checking.

SECC.ARA.t1 <- within( SECC.ARA.t1, {
  # Some Fragmentation entries are '-' and need to be recoded as 'NA'
# Frag <- recode(Frag,
#                "'1'='1' ; '2'='2' ; '3'='3' ; '4'='4' ; else=''",
#                as.factor.result = TRUE,
#                levels = c('1', '2', '3', '4')
#               )
})


##################################################
## CALCULATIONS
##################################################
## 

SECC.ARA.t1 <- within( SECC.ARA.t1, {
  
})

##################################################
## MERGE TIME POINTS
##################################################
## 

SECC.ARA <- SECC.ARA.t1

## Remove old objects from memory
rm(list=c('SECC.ARA.t1'))

##################################################
## CHECK DATA
##################################################

head(SECC.ARA)  # have a peek at the first 6 rows & columns: is this what you expected?
# str(SECC.ARA.t1)   # check structure: are the appropriate variables factors, numeric, etc.?

##################################################
## SAVE DATA
##################################################
# leave in memory
# only rows for samples, exclude controls
# SECC.ARA.t1 <- SECC.ARA.t1[SECC.ARA.t1$SampleControl=="Sample",]
