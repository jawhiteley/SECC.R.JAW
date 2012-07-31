##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Functions
# Jonathan Whiteley		R v2.12		2011-03-23
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

SECCcolumnNames <- function(data=NULL, DataName="data") {
  ## Standardize column names
  # Check if first argument is actually a data frame.
  if (is.data.frame(data) == FALSE) stop(
      "The first argument of the SECCcolumns() function must be a \"data.frame\"."
     )

  ## REQUIRES:
  load('./save/SECC_factors.R')	# includes SECC.base data.frame, 
    # and other vectors of standard column names and levels.
  ## Standard ID column names and types (based on template)
   # colnames(SECC.base) includes SampleID ; Trt_nest_order does not
  ColNames_std <- colnames(SECC.base)
  SampleID.synonyms <- c("SampleID", "PatchID", "ID", "Sample", "Treatment")
  Block.synonyms   <- c("Block", "block")
  Time.synonyms    <- c("Time", "TimePt", "time", "time.point", "time.pt", "t.point")
  Chamber.synonyms <- c("Chamber", "chamber", "Warming", "warming", "warm", "Chamber.trt", "Warming.trt")
  Frag.synonyms    <- c("Frag", "Fragmentation", "Frag.trt", "fragmentation", "frag", "frag.trt")
  Pos.synonyms     <- c("Pos", "Position", "position", "pos")

  for (ColName in ColNames_std) {
    ## Check for presence of standard ID columns.--------------------------------  
    ## Check for most likely non-standard names (synonyms)
    col.synonyms <- get(paste(ColName, ".synonyms", sep=""))
    ## Attempt to rename based on likely matches.
    # which synonym is present in the data frame?
    colnames.std <- colnames(data) %in% col.synonyms  
    # replace synonymous colnames in data frame with standard name (1st item)
    colnames(data)[colnames.std==TRUE] <- col.synonyms[1]
      
    # If column name still does not match, throw an Error message: this file needs Special Attention.
    if ( ColName %in% colnames(data) == FALSE ) {
      Data_Check <- FALSE
      print('column synonyms:', quote = FALSE)
      print(col.synonyms, quote = FALSE)
      print('column names:', quote = FALSE)
      print(colnames(data), quote = FALSE)
      print('columns present:', quote = FALSE)
      print(colnames.std, quote = FALSE)
      stop(
        paste( "Column \'", ColName, "\' is missing from data frame: ", DataName, "\n",
               "Automatic search & re-naming based on common synonyms failed.\n", 
               "This File needs Special Attention.",
               sep=""
            )
      )
    }
    ##--------------------------------
  }
  return(data)
}


checkSECCdata <- function(data=NULL, DataName="-", CheckValues = TRUE, CheckDuplicates = TRUE) {
  ## Checks a data frame argument to make sure that 
  ## it conforms to the standards for the 
  ## Schefferville Experiment on Climate Change (SEC-C).
  ## - If possible, it also tries to fix or convert common deviations from the standard.

  ## CheckValues = TRUE argument means the function will also check factor levels,
  ## and standard IDs.  Set to FALSE to skip these steps.

  ## CheckDuplicates = TRUE checks for duplicate SampleIDs.  This is independent of CheckValues.

  # Check if first argument is actually a data frame.
  if ('data.frame' %in% class(data) == FALSE) stop(
      paste("The first argument of the checkSECCdata() function", 
          "must be an object of class \"data.frame\"."
         )
     )
  
  ## REQUIRES:
  load('./save/SECC_factors.R')	# includes SECC.base data.frame, 
    # and other vectors of standard column names and levels.
# Pos_lvls <- Pos_all_lvls  # for convenience later.
  ## Standard ID column names and types (based on template)
   # colnames(SECC.base) includes SampleID ; Trt_nest_order does not
  ColNames_std <- colnames(SECC.base)
  
  ## BEGIN checking data frame
  Data_Check <- TRUE	# default value: set to FALSE if a check fails.
  # data.new   <- data    # this is what will be returned.

  ## Strip blank (empty) rows with no data? values of NA or '' in all columns?
  data <- data[is.na(data[1])==FALSE,]
  data <- data[data[1]!='',]

  if ( min( ColNames_std %in% colnames(data) ) == 0 ) {
    # min(boolean_vector) == 0 means there was at least one FALSE result
    data <- SECCcolumnNames(data, DataName)  # standardize column names.
  }
  ## Strip rows with no valid time points?  Need standard column names first!
#  data <- data[(data$Time %in% levels(SECC.base$Time) ),]
  for (ColName in ColNames_std) {
    ## Strip rows with no valid levels?  Need standard column names first!
     # This does not remove factor levels that are no longer present:
     #   removing unused levels would require recoding or re-definition of factor column.
#    data <- data[(data[[ColName]] %in% levels(SECC.base[[ColName]])),]

    ## Standardize ID column types & values (levels) ---------------- 
     # Check that column types & values (levels) match template
    if ( class(data[[ColName]]) != class(SECC.base[[ColName]]) ) {
      # Attempt to convert column class if necessary
      if (is.factor(SECC.base[[ColName]]))
        data[[ColName]] <- factor(data[[ColName]])  # Convert to a factor.
      else {
        # as() will not remove extraneous factor levels.
        # as() will wipe out integer values when converting to a factor :(
        data[[ColName]] <- as( data[[ColName]], class(SECC.base[[ColName]]) )
      }
      # Check for unexpected Warnings: 
      # Unexpected Warnings means this file needs Special Attention.
    } 

    if ( CheckValues && class(data[[ColName]]) == "factor" ) {
      ## Standardize factor levels
       # Unused factor levels may still be present after filtering.
       # It might be more reliable to check for unique() values
       # and leave the factor levels to recoding after the merge.
      # Actual Factors should contain *ONLY* standard levels.
      # Otherwise, it is too difficult to automatically check for errors.
      data[[ColName]] <- factor(data[[ColName]])  # to remove unused factor levels
      Col_lvls <- levels( SECC.base[[ColName]] )
      if ( min( unique(data[[ColName]]) %in% Col_lvls ) == 0 ) {
        # min(boolean_vector) == 0 means there was at least one FALSE result
        ## Catch common differences in factor levels
        if ( ColName == "Pos" && min( c("0", "1") %in% levels(data[[ColName]])) == 1 ) {
          ## could also check for Pos_all_lvls,
          ## but the most important thing is to recode 1s & 0s
          # Rename factor levels from old to new: 'new'='old'
          # non-existent levels are ignored; unspecified levels are replaced with NA
          #  - levels with no 'new' value would be replaced with empty strings.
          levels(data[[ColName]]) <- list('I'='1',
                                          'S'='S',
                                          'W'='W',
                                          'E'='E',
                                          'N'='N',
                                          'O'='0'
                                          )
          ## Re-compute SampleID column based on new values?
          ## Original SampleID should be kept for reference & error-checking:
          ## I don't want to have to use it in the merge,
          ## but I may need it to resolve duplicates
          ## for things like controls in the source data frames.
          ## This should really be handled for each data file individually
          ## if this simple algorithm would create duplicates.
          data <- within( data,
                         SampleID <- paste(Block, Time, Chamber, "-", Frag, ".", Pos,
                                           sep=""
                                           )
                         )
        }
        else {
          # Automatically re-map other factor levels to standard values?

          # Automatically purge non-standard levels (replace with NA)?
#         non.std.values <- data[[ColName]] %in%  Col_lvls
          data[[ColName]] <- factor( data[[ColName]], levels = Col_lvls )
#         , exclude = non.std.values
          # Values not in levels are replaced with NA.  This is exactly what I want.
        }
        # Anything else?
      }

      ## If factor levels still do not match, this file needs Special Attention
       # Extra factor levels may not be a problem with the correct merge settings:
       # - un-matched rows can simply be omitted from the merge result.
       # Nevertheless, I need this to be able to check for mismatches in the combinations of all ID columns.
      if ( min( levels(data[[ColName]]) %in% Col_lvls) == 0 ) {
        print("Standard Levels", quote = FALSE)
        print(Col_lvls, quote = FALSE)
        print(paste(ColName, "levels"), quote = FALSE)
        print(levels(data[[ColName]]), quote = FALSE)
        stop( 
          paste( "Non-standard levels in column \'", ColName, "\' of data frame: ",
                 DataName,
                 "\nAuto-correction failed.  This File needs Special Attention.",
                 sep=""
               )
        )
      }
    }
    ## Check that the standard column is at least not empty (only NAs or '' )
    EmptyStrings <- ifelse( class(data[[ColName]]) %in% c("factor", "character"),
                           min( data[[ColName]][!is.na(data[[ColName]])] == '' ),
                           FALSE
                           )
    NAcolumn <- min( is.na(data[[ColName]]) )
    if ( NAcolumn == 1 || EmptyStrings == 1 ) {
      print("Standard Levels:", quote = FALSE)
      print(unique(SECC.base[[ColName]]), quote = FALSE)
      print(paste("Values in:", ColName), quote = FALSE)
      print(unique(data[[ColName]]), quote = FALSE)
      stop( 
        paste( "Column \'", ColName, "\' of data frame: ",
               DataName, " is EMPTY.",
               "\nAuto-correction failed.  This File needs Special Attention.",
               sep=""
             )
      )
    }
    ##--------------------------------
  }

  ## Check that all IDs match EXACTLY with Template (or that all rows in data are represented in the template)
  ## Un-matched rows will be omitted from merge *****
  ## But, I want to check for mismatched IDs (especially Position), that would lead to a row being mistakenly excluded.
  if (CheckValues) {
    # Are rows at least a subset of Template: missing rows can be replaced with NA's?
    DataIDs <- with( data, paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep="" ) )
    Std.IDs <- DataIDs %in% SECC.base[["SampleID"]]
    GasBlanks <- grep( "[1-8][1-4][ABC]-NA\\.NA", DataIDs, perl=TRUE )  # ARA Gas Controls will have BOTH Fragmentation & Position as 'NA'
    Std.IDs[GasBlanks] <- TRUE # ok
  # If I trust the Position data in SECC.base, I could use it to correct values in data automatically (then re-calculate SampleIDs & re-check).
    if ( min( Std.IDs ) == 0 ) {
      # If IDs still do not match, throw an Error message:
      # this file needs Special Attention.
      Pblm.msg <- paste("\nNon-Standard IDs:\n",
                        paste(DataIDs[!Std.IDs], collapse ="\n"),
                        sep=""
                        )
      stop( paste( "Non-standard IDs in data frame:", DataName,
                  "\nThis File needs Special Attention.",
                  Pblm.msg
                  )
           )
    }
  }

  if (CheckDuplicates) {
    SampleID.counts <- table(data$SampleID)  # Counts of each 'factor level'
    SampleID.wrong   <- SampleID.counts[SampleID.counts!=1] # IDs with counts other than 1
    if( length(SampleID.wrong)>0 ) {
      Duplicates <- paste(names(SampleID.counts[SampleID.counts>1]), collapse="\n")
      Missing    <- paste(names(SampleID.counts[SampleID.counts<1]), collapse="\n")
      text.empty <- "<none>"
      Pblm.msg <- paste("\nDuplicates:\n",
                        ifelse(Duplicates=="", text.empty, Duplicates),
                        "\nMissing:\n",
                        ifelse(Missing==""   , text.empty, Missing   ),
                        sep = ''
                        )
      stop( paste("There are duplicate or missing IDs",
                  Pblm.msg, sep='') )
    }
  }
  
  # If all checks pass, continue.  
  # Otherwise, throw an Error message: this file needs Special Attention.
  if (Data_Check == FALSE)
    stop( paste("There is an unknown problem with data frame:", DataName) )
  else ## RETURN cleaned version if there are no problems.
    cat(paste("All checks passed for Data Object:", DataName, "\n"))
    return(data)
}


SECCstr <- function (data) {
  ## see also ?replications
  ## replications(~ Block + Time + Chamber + Frag + Position, data = SECC)
  # Check if first argument is actually a data frame.
  if ('data.frame' %in% class(data) == FALSE)
    stop( "The first argument must be a \'data.frame\'." )
  ## REQUIRES:
  load('./save/SECC_factors.R')	# includes SECC.base data.frame, 
    # and other vectors of standard column names and levels.

  SECC.summary <- list()  # empty list to add to.
  
  ## Assemble a table of counts for all levels of each standard ID column
  ID_cols <- colnames(SECC.base)
  if("Position" %in% names(data))
     ID_cols <- c( ID_cols, "Position" )     # include re-coded Positions.
  ID_cols <- ID_cols[ID_cols!="SampleID"] # drop SampleID column!! (Loops crash with this insanely long vector).
  Data_cols <- colnames(data)
  Data_cols <- Data_cols[!(Data_cols %in% ID_cols)]

  ## Get counts of values of each level (excluding NAs)
  ## Essentially: length( na.omit( data[ data[[ID.col]] == ID.lvl, Data.col ] ) )
   # More efficient to use a version of apply() or by() instead?
  for (ID.col in ID_cols) {
    Col.lvls <- levels(data[[ID.col]])
    Col.totals  <- as.vector( table(data[[ID.col]]) )  # summary of this column
    Col.summary <- data.frame(Levels = Col.lvls, Total = Col.totals)
    for (Data.col in Data_cols) {
      Col.summary[[Data.col]] <- NA  # initialize column
      for (ID.lvl in Col.lvls) {
        Col.data <- data[data[[ID.col]]==ID.lvl, Data.col]
        Col.data <- na.omit(Col.data)
        Col.summary[Col.lvls == ID.lvl, Data.col] <- length( Col.data )
      }
    }
    SECC.summary[[ID.col]] <- Col.summary  # add ID column summary to summary object.
  }
  return(SECC.summary)
}


##==================================================
## Recode standardized factors to standard values prior to merging
recodeSECC <- function(data=NULL) {
  require(car)  # for recode()
  data.recoded <- within( data, {
    ## Rename columns / convert to standard informative names (already handled in SECCcolumns).
    ## Ensure columns that should be factors are, and drop unused levels.
    Block   <- factor(Block)
    Time    <- factor(Time)
    Chamber <- factor(Chamber)
    Frag    <- factor(Frag)
    Pos     <- factor(Pos)
    ## I can also get informative labels using the 'labels' argument of factor(),
    ## rather than use recode() in package 'car'.
    Chamber  <- factor(Chamber,
                      levels = c("A", "B", "C"),
                      labels = c("Ambient", "Partial Chamber", "Full Chamber")
                      )

    ## Rename Fragmentation codes
    Frag     <- factor(Frag,
                       levels = 1:4,
                       labels = c("Continuous", "Full Corridors",   # Continuous / Contiguous ?
                                  "Pseudo-Corridors", "Isolated")
                       )
    ## Recode Time to approximate dates: replace t3 with NA, since it has no date (yet)
    Time     <- factor(Time,
                       levels = c(1, 2, 4),  # Drop t3, replace with NA
                       labels = c( "12 months", "22 months", "24 months" )
                       )
    ## New factor with simplified recoded values for Patch Position
    ## Rename and reorder factor levels using recode()
     # - maintains empty values if empty factor specified,
     # otherwise converts to 'NA'.
     # Requires package 'car'
    Position <- recode(Pos, 
                       "'I'='Inner'; 'O'='Outer'; else='other'", 
                       levels=c( "Inner", "other", "Outer" ),
                       as.factor.result=TRUE
                       )
  })

  return(data.recoded)
}


