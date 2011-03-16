##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# functions used in most analyses
# Jonathan Whiteley		R v2.10.1		2010-05-08
##################################################
## access these functions in another file by using: 
## 	 source("path/to/this/file.R")

##################################################
## ANOVA: standardized analysis of individual response variables
##################################################
# Standard univariate Analysis of Variance on NESTED treatments: 
# No fancy mixed effects or error distribution families (GLMM) yet.
# The same experimental design applies for each response variable, 
# therefore a standard analytical procedure applies.
# Major differences include details such as:
#   - which Transformation is appropriate
#   - Graphs of significant Main Effects

##================================================
## DATA EXPLORATION
##================================================
# See `/lib/fun.R` for data processing functions
# used to produce a data frame with standard columns
# for standardized analysis steps.



##================================================
## DEFINE MODEL FORMULA
##================================================




##================================================
## CHECK ASSUMPTIONS: residuals, standard diagnostic plots
##================================================




##================================================
## ANALYSIS: GET RESULTS
##================================================



##================================================
## PLOTS & GRAPHS
##================================================
chamberMap <- function ( labels=c("A", "B", "C") ) {
  Chamber.map <- data.frame( label=labels, 
    col = c("#000000","#000099","#990000"), 
    bg  = c("#FFFFFF","#FFFFFF","#FFFFFF"), 
    pch = c(21,23,18), lty = c(3,2,1) 
  )
    # A) Ambient = black, open circles with dotted line ; 
    # B) Partial = blue, open diamonds with dashed line ; 
    # C) Full    = red, solid diamond with solid line.
}


