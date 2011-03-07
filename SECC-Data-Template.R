####################################################
# Template generator for SECC data entry
# With standard columns, all treatment level combinations,
# and ID codes: factor levels are coded using single characters,
# as used for sample labels.
# Jonathan Whiteley	2011-01-24
# R v2.12
####################################################
## INITIALISE
####################################################
# rm(list=ls())	# clear memory

## LOAD LIBRARIES
# library(car)	# load external package 'car', for recode()
library(reshape)	# sort_df (sort data frame) wrapper for order


####################################################
## GENERATE DATA FRAME COLUMNS
####################################################

Block_lvls    <- as.integer( seq(1, 8) )
Time_lvls     <- as.integer( seq(1, 4) )
Chamber_lvls  <- c("A", "B", "C")
Chamber_labels<- c("Ambient", "Partial Chamber", "Full Chamber")
Chamber.Heat  <- c("Ambient", "Warm", "Hot")
Frag_lvls     <- as.integer( seq(1, 4) )
Frag_labels   <- c("Contiguous", "Full Corridors", "Pseudo-Corridors", "Isolated")
Pos_old_lvls  <- c( 1 , "S", "W", "E", "N",  0 )  # I used to use 0 & 1 for 'Outer' & 'Inner', respectively.
Pos_all_lvls  <- c("I", "S", "W", "E", "N", "O")
Pos_all_sort  <- c("O", "I", "N", "S", "E", "W") # desired sort order.
Pos_lvls      <- c("I", "S/N", "W/E", "O")	# S/N is for alphabetical sorting.  Do not change.
Pos_sort      <- c("O", "I", "S/N", "W/E")	# S/N is for alphabetical sorting.  Do not change.
Pos_labels    <- c("Inner", "other", "Outer")
Pos.Precip    <- c("Wet", "mesic", "Dry")

Trt_nest_order <- c("Block", "Time","Chamber", "Frag", "Pos")
Trt_sort_order <- c("Time", "Block", "Chamber", "Frag") # not including "Pos" or "Position"
  # The nesting structure puts "Block" first (largest experimental unit).
  # I tend to sort by Time first, however,
  # To reflect the order in which samples were actually collected & processed.

# Make a list of objects to be saved at the end.
save.ls <- c(
  'Block_lvls',
  'Time_lvls',
  'Chamber_lvls', 'Chamber_labels', 'Chamber.Heat',
  'Frag_lvls', 'Frag_labels',
  'Pos_old_lvls', 'Pos_all_lvls', 'Pos_all_sort', 
  'Pos_labels', 'Pos.Precip', 
  'Trt_nest_order', 'Trt_sort_order'
)

# Generate all combinations of treatment levels (cycling through faster, then slower columns).
SECC.grid <- expand.grid( 
                          Block   = Block_lvls, 
                          Time    = Time_lvls, 
                          Chamber = Chamber_lvls, 
                          Frag    = Frag_lvls, 
                          Pos     = Pos_sort, 
                          KEEP.OUT.ATTRS = FALSE 
                        )

# SECC_sorted <- SECC.grid
SECC_sorted <- sort_df( SECC.grid, 
  vars=c(Trt_sort_order, "Pos") 
) 
  # useful for wholesale sorting, but I really want the Position column in a particular, non-alphabetical order (which is used for data entry).  
  # unless, I can use the factor codes for Pos, rather than the levels. ***
  # The columns will be re-ordered in the next step anyway.

str(SECC_sorted)
head(SECC_sorted)
  
SECC.base <- with(SECC_sorted, data.frame(
  SampleID = paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep=""),
  Block   = factor(Block),
  Time = factor(Time),
  Chamber = factor(Chamber),
  Frag    = factor(Frag),
  Pos     = factor(Pos, levels = Pos_lvls),
  stringsAsFactors = FALSE
)  )
# SECC.base <- within( SECC.base, SampleID <- as.character(SampleID) )	# convert factor to character

# at this point, the basic structure is ready,
# but it does not yet contain the actual values for Position (Pos):
# Inner, Outer, either North or South, either East or West
# These must be imported from the metadata for the experimental design.


##==================================================
## CHECK DATA
str(SECC.base)		# check structure: are the appropriate columns factors, numeric, etc.?
# invisible(edit(SECC.base))	# opens spreadsheet and returns changes invisibly.  
# Use fix() to make permanent changes, or assign the result of edit() to an object.
head(SECC.base)		# have a peek at the first 6 rows & columns: is this what you expected?



####################################################
## LOAD METADATA & PATCH IDs
####################################################
# Load MetaData containing actual (observed) patch positions.
PatchIDs <- read.csv("./data/SECC-PatchIDs.csv")
PatchIDs <- with( PatchIDs, data.frame(
  PatchID  = as.character(PatchID),
  Block    = factor(Block),
  Time     = factor(Time.point),
  Chamber  = factor(Chamber),
  Frag     = factor(Fragmentation),
  Position = factor(Position, levels=Pos_all_sort),  # levels in sort order
  stringsAsFactors = FALSE
) )
# Sort PatchIDs in the same order as the base template (or make sure it is).
# This is great in theory, except that I really want patch Position in the order specified in the base template, and I'm not sure how to do this automatically:
  # O(uter), I(nner), N/S, E/W
  # This is mostly to facilitate data processing & entry, which typically occurs in this order.
  # Sorting a data frame actually uses the numeric codes for factors, not the levels.  
  # So, I just have to re-order the levels to be in the order I want the values to be sorted in.
Patch_pos <- sort_df( PatchIDs, 
  vars=c(Trt_sort_order, "Position") 
)

str(Patch_pos)
head(Patch_pos)


####################################################
## UPDATE BASE WITH ACTUAL Position VALUES (from MetaData)
####################################################
#SECC.merged <- merge(
#  SECC.base, # [colnames(SECC.base)!=c("SampleID","Pos")], 
#  Patch_pos, 
#  sort=FALSE
#) # replicates each entry by 4??  Look at 'by' & 'all' args.

# all.equal(SECC.base, PatchIDs)
Rows_mismatched <- sum(
  SECC.base$Block   != Patch_pos$Block   |
  SECC.base$Time    != Patch_pos$Time |
  SECC.base$Chamber != Patch_pos$Chamber |
  SECC.base$Frag    != Patch_pos$Frag
)	# rows where columns are mismatched.  Should be 0

if (Rows_mismatched == 0) {
  SECC.merged <- within( SECC.base, {
    Pos <- factor(Patch_pos$Position, levels = Pos_all_lvls)  # levels in desired order.
    SampleID = paste(Block, Time, Chamber, "-", Frag, ".", Pos, sep="")
  })
} else stop("Rows do not all match.  Check row sorting before trying to merge.")

IDs_mismatched <- sum( SECC.merged$SampleID != Patch_pos$PatchID )
if (IDs_mismatched == 0) {
  cat("All IDs match.  Success!")
} else cat(IDs_mismatched, "IDs did not match!  Check inputs (especially sorting) and try again.")


## Check Data before export.
str(SECC.merged)
head(SECC.merged)
# invisible(edit(SECC.merged))



####################################################
## EXPORT TO FILE (END)
####################################################
# I usually just use a default working directory on each computer, but this can speed up the process, for projects with files in a different directory:
# setwd("../Data/Data Template/")	# Set Working Directory.
# setwd("./ SECC/")	# Project Directory.
getwd()	            # check that we're in the right place.

# Export objects to files
SECC.str  <- SECC.base    # for checking 
SECC.base <- SECC.merged  # for saving
write.csv( SECC.merged, 
  file="./save/SECC_base.csv", 
  row.names=FALSE
)
save( list=c( 'SECC.base', save.ls ), 
  file="./save/SECC_factors.R" 
)

rm(list=ls())	# clear memory
