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

checkSECCdata <- function (data=NULL, DataName="data") {
  ## Checks a data frame argument to make sure that 
  ## it conforms to the standards for the 
  ## Schefferville Experiment on Climate Change (SEC-C).
  ## - If possible, it also tries to fix or convert common deviations from the standard.

  # Check if first argument is actually a data frame.
  if ('data.frame' %in% class(data) == FALSE) stop(
      paste("The first argument of the checkSECCdata() function", 
          "must be an object of class \"data.frame\"."
         )
     )
  
  ## REQUIRES:
  load('./save/SECC_factors.R')	# includes SECC.base data.frame, 
    # and other vectors of standard column names and levels.
  Pos_lvls <- Pos_all_lvls  # for convenience later.
  SampleID.synonyms <- c("SampleID", "PatchID", "ID")
  Block.synonyms   <- c("Block", "block")
  Time.synonyms    <- c("Time", "time.pt", "TimePt")
  Chamber.synonyms <- c("Chamber", "chamber", "Warming", "warming", "warm", "Chamber.trt", "Warming.trt")
  Frag.synonyms    <- c("Frag", "Fragmentation", "Frag.trt", "fragmentation", "frag", "frag.trt")
  Pos.synonyms     <- c("Pos", "Position", "position", "pos")
  
  ## Standard ID column names and types (based on template)
  ColNames_std <- colnames(SECC.base)	# colnames(SECC.base) (SampleID) or Trt_nest_order
  
  ## BEGIN checking data frame
  Data_Check <- TRUE	# default value: set to FALSE if a check fails.
  # data.new   <- data    # this is what will be returned.
  
  if ( min( colnames(DataObject) %in% ColNames_std ) == 0 )
  {
    # min(boolean_vector) == 0 means there was at least one FALSE result
    # this check is built into the checks for each column,
    # to put all operations in a single loop.
  }
    for (ColName in ColNames_std) 
    {
      ## Check for presence of standard ID columns.--------------------------------  
      ## Check for most likely non-standard names (synonyms)
      col.synonyms <- get(paste(ColName, ".synonyms", sep=""))
      ## Attempt to rename based on likely matches.
      # which synonym is present in the data frame?
      colnames.std <- colnames(data) %in% col.synonyms  
      # replace synonymous colnames in data frame with standard name (1st item)
      colnames(data)[colnames.std==TRUE] <- col.synonyms[1]
        
      # If column name still does not match, throw an Error message: this file needs Special Attention.
      if ( ColName %in% colnames(DataObject) == FALSE )
      {
        Data_Check <- FALSE
        stop(
          paste( "Column \'", ColName, "\' is missing from data frame:", DataName, "\n",
                 "      Automatic search & re-naming based on common synonyms failed.\n", 
                 "      This File needs Special Attention."
               )
        )
      }
      ##-------------------------------- 
      
      ## Standardize ID column types & values (levels) ---------------- 
       # Check that column types & values (levels) match template
      if ( typeof(data[[ColName]]) != typeof(SECC.base[[ColName]])
      {
        # Attempt to convert column data type if necessary
        data[[ColName]] <- as( data[[ColName]], typeof(SECC.base[[ColName]]) )
        # Check for unexpected Warnings: 
        # Unexpected Warnings means this file needs Special Attention.
      } 
      if ( class(data[[ColName]]) != class(SECC.base[[ColName]])
      {
        # Attempt to convert column class if necessary
        data[[ColName]] <- as( data[[ColName]], class(SECC.base[[ColName]]) )
        # Check for unexpected Warnings: 
        # Unexpected Warnings means this file needs Special Attention.
      } 

      if ( class(data[[ColName]]) == "factor" )
      {
        ## Standardize factor levels
        Col_lvls <- get( paste(ColName, "_lvls", sep="") )
        if ( min( levels(data[[ColName]]) %in% Col_lvls) == 0 )
        {
          # min(boolean_vector) == 0 means there was at least on FALSE result
          ## Catch common differences in factor levels
          if ( ColName == "Pos" && (levels(data[[ColName]]) %in% Pos_old_lvls) ) 
          {
            # rename factor levels from old to new; 'new'='old'
            # (omitted levels are dropped and replaced with empty strings)
            levels(data[[ColName]]) <- list(
              'I'='1', 'S'='S', 'W'='W', 'E'='E', 'N'='N', 'O'='0'
            )
          }
          # Anything else?  Automatically re-map other factor levels?
        }
        if ( min( levels(data[[ColName]]) %in% Col_lvls) == 0 )
        {
          stop( 
            paste( "Column \'", ColName, "\' in data frame:", DataName,
                   "Has non-standard levels. \n", 
                   "      Auto-correction failed.  This File needs Special Attention."
                 )
          )
        }
      }
      ##--------------------------------
    }
  
  ## Check that all IDs match EXACTLY with Template (or that all rows in data are represented in the template)
  # Check number of rows - data files may only be a subset
  if ( dim(data)[1] != dim(SECC.base)[1] )
  {
    # If fewer rows, are they a subset: can missing rows be replaced with NA's?
      if ( min( data[["SampleID"]] %in% SECC.base[["SampleID"]]) == 0 )
      {
        # If IDs still do not match, throw an Error message:
        # this file needs Special Attention.
        stop( 
          paste( "Non-standard IDs in data frame:", DataName,
                 "\nThis File needs Special Attention."
               )
        )
      }
  }
  
  # If all checks pass, continue.  
  # Otherwise, throw an Error message: this file needs Special Attention.
  if (Data_Check == FALSE) stop(
    paste("There is an unknown problem with data frame:", DataName)
  )

  ## RETURN cleaned version if there are no problems.
  return data
}

##
function recodeSECC (data=NULL)
{
  data.recoded <- within( data, {
    ## Rename columns / convert to standard informative names for the rest of this script
    Block   <- factor(Block)
    Time    <- factor(Time)
    Chamber <- factor(Chamber)
    Frag    <- factor(Frag)
    Pos     <- factor(Pos)
    ## rename and reorder factor levels the easy way - maintains empty values if empty factor specified, otherwise converts to 'NA'.  Requires package 'car'
    Chamber <- recode( Chamber, 
      "'A'='Ambient'; 'B'='Partial Chamber'; 'C'='Full Chamber'; else=''", 
      levels=c( "Ambient", "Partial Chamber", "Full Chamber" ),
      as.factor.result=TRUE
    )
    levels(Frag) <- list( 'Continuous'=1, 'Full Corridors'=2, 'Pseudo-Corridors'=3, 'Isolated'=4 )
  #  Position <- factor(Pos, levels=c('1', 'S', 'W', 'E', 'N', '0'))	# safely reorder factor levels
  #  levels(Position) <- list( 'I'='1', 'S'='S', 'W'='W', 'E'='E', 'N'='N', 'O'='0' )	# rename some factor levels (omitted levels are dropped and replaced with empty strings).
    # New factor with simplified recoded values for Patch Position
    position <- recode( Pos, 
      "'I'='Inner'; 'O'='Outer'; else='other'", 
      levels=c( "Inner", "other", "Outer" ),
      as.factor.result=TRUE
    )
  })
  return data.recoded
}
