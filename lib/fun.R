##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Functions
# Jonathan Whiteley		R v2.12		2011-01-26
##################################################
source('./lib/SECC.functions.R')
source('./lib/jaw.graph_functions.R')
source('./lib/jaw.copied_functions.R')
source('./lib/jaw.misc_functions.R')

##################################################
# DATA PROCESSING FUNCTIONS
##################################################
# Produce a data frame with standard columns
# used for standardized analysis steps.

##================================================
## CHECK DATA
##================================================

checkSECCdata <- function (data=NULL) {
  ## Checks a data frame argument to make sure that 
  ## it conforms to the standards for the 
  ## Schefferville Experiment on Climate Change (SEC-C).
  # If possible, it also tries to fix or convert common deviations from the standard.
  
  # Check if first argument is actually a data frame.
  if ('data.frame' %in% class(data) == FALSE) stop(
      cat("ERROR: the first argument of the checkSECCdata function", 
          "must be an object of class \"data.frame\"."
         )
     )
  
}
