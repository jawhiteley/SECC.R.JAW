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
# str(SECC.cyanobacteria)  # Should already be in memory.  str() produces output on source().

SECC.cyanobacteria <- SECCcolumnNames(SECC.cyanobacteria)  # Standardize ID column names (but not values)

##================================================
## MANUALLY CLEAN & PROCESS DATA
##================================================
# Manually clean & prepare data for automatic checking.

SECC.cyanobacteria <- within( SECC.cyanobacteria, {
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

SECC.cyanobacteria <- within( SECC.cyanobacteria, {
  
})

##################################################
## CHECK DATA
##################################################

head(SECC.cyanobacteria)  # have a peek at the first 6 rows & columns: is this what you expected?
# str(SECC.cyanobacteria)   # check structure: are the appropriate variables factors, numeric, etc.?

##################################################
## SAVE DATA
##################################################
# leave in memory
