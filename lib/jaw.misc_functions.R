##################################################
# Miscellaneous useful functions
# Jonathan Whiteley		R v2.12		2012-07-15
##################################################

strip_empty_dims  <- function( data = NULL, dim = c(1, 2), 
							  rows = NULL, cols = NULL, col.class = NULL ) {
  ## strip whole rows or columns if they are entirely empty (NA or NaN)
  ## `data` is the data frame to strip empty values from.
  ## `rows` is a vector of indices or logical values of 
  ##        which rows to check for empty values (default all)
  ## `cols` see `rows` above (but applied to columns)
  ## `dim`  the dimensions to strip if all empty.  1 = rows, 2 = cols
  ## `col.class` is a vector of classes: columns that match will be included
  ##             This applies to the `cols` argument.

  if (!is.data.frame(data) && !is.matrix(data)) stop(
	  "The first argument must be a \"data.frame\" or \"matrix\"."
     )

  cols <- as.vector(cols)
  if ( is.null(cols) ) {
	cols <-  1:ncol(data)
  }
  if ( !is.null(col.class) ) {
	check.cols <- which( lapply(data[, cols], class) %in% as.vector(col.class) )
	check.cols <- cols[check.cols]
  } else {
	check.cols <- cols
  }
  rows <- as.vector(rows)
  if ( is.null(rows) ) {
	rows <-  1:nrow(data)
  }  
  
  empty.rows <- which( apply( data[, check.cols], 1, function(x) all(is.na(x)) ) )
  empty.cols <- which( apply( data[rows,       ], 2, function(x) all(is.na(x)) ) )

  if (all(dim == 1)) {
	if (length(empty.rows) > 0) data <- data[-empty.rows, ]
  } else if (all(dim == 2))  {
	if (length(empty.cols) > 0) data <- data[, -empty.cols]
  } else if (all(dim %in% c(1, 2))) {
	if (length(empty.rows) > 0 && length(empty.cols) > 0) 
      data <- data[-empty.rows, -empty.cols]
  }
  
  return(data)
  
  if (FALSE) {  # Testing
	test.df <-  data.frame( )
	df.proc <- strip_empty_dims(test.df)
	df.proc
  }
}


###=============================================================
### glmulti wrapper functions for mixed effects modelling
## lmer.glmulti from ?glmulti examples
lmer.glmulti <- function (formula, data, random = "", ...) {
  lmer(paste(deparse(formula), random), data = data, REML=FALSE, ...)
}
# the fixed-effects are passed as formula, and the random effects are passed as "random"
# Here we could redefine the getfit function, but this is necessary only to use coef() or predict().
# 3.Last, we must provide the corresponding aicc method, since the default will not work with mer objects
lmer.aicc <- function()
{
  require(glmulti)
  setMethod('aicc', 'mer', function(object, ...) 
            {
              liliac<- logLik(object)
              k<-attr(liliac,"df")
              n= object@dims['n']
              return(-2*as.numeric(liliac[1]) + 2*k*n/max(n-k-1,0))
            })
}
## nlme equivalents - these don't have the same coef() methods, so they may not be useful, even if they work
  lme.glmulti <- function (formula, data, random, REML=FALSE, ...) {
    if (REML) method.gls <- "REML" else method.gls <- "ML"
    lme(formula, random = random, data = data, method = method.gls, ...)
  }
  gls.glmulti <- function (formula, data, REML=FALSE, ...) {
    if (REML) method.gls <- "REML" else method.gls <- "ML"
    gls(formula, data = data, method = method.gls, ...)
  }

### glmulti functions for extracting plot data
getCoef.glmulti <- function(glmObj, minImportance=0) {
  glm.coef <- as.data.frame(coef(glmObj))
  ## from coef: these are weights of model *coefficients*, NOT *model terms*...
  ## I think I want weights of model terms (variables, rather than levels of each factor)
  ## ggplot will order the bars by levels of the explanatory factor
  glm.coef$Term <- factor(row.names(glm.coef), levels=unique(row.names(glm.coef))) 
  glm.order <- order(glm.coef$Importance, decreasing=TRUE)
  ## glm.coef$Term <- factor(glm.coef$Term, levels=levels(glm.coef$Term)[glm.order])
  ## Drop terms below the threshold
  glm.minImp <- which(glm.coef$Importance >= minImportance)
  glm.coef <- glm.coef[glm.minImp, ]
  ## facilitate confidence intervals on Estimates
  glm.coef$Emax <- glm.coef$Estimate + glm.coef[, "+/- (alpha=0.05)"]
  glm.coef$Emin <- glm.coef$Estimate - glm.coef[, "+/- (alpha=0.05)"]
  glm.coef
  ## Process coef table afterwards (text replacement, etc.)
}
importance.glmulti <- function(x) {
  ## collect importances (see plot.glmulti( type="s" )
  ww = exp(-(x@crits - x@crits[1])/2)
  ww = ww/sum(ww)
  clartou = function(x) {
    pieces <- sort(strsplit(x, ":")[[1]])
    if (length(pieces) > 1) 
      paste(pieces[1], ":", pieces[2], sep = "")
    else x
  }
  tet = lapply(x@formulas, function(x) sapply(attr(delete.response(terms(x)), 
                                                   "term.labels"), clartou))
  allt <- unique(unlist(tet))
  imp <- sapply(allt, function(x) sum(ww[sapply(tet, function(t) x %in% t)]))
  order.imp <- order(imp)
  glm.imp <- data.frame(Term=allt[order.imp], Importance=imp[order.imp], stringsAsFactors=F)
  glm.imp$Term <- factor(glm.imp$Term, levels=glm.imp$Term) # sorted levels?
  glm.imp
}




###=============================================================
## Partial Regression formula processing
PartialFormula <- function (model = "", x.var = "", part = "both")
{
  RHS.part <- as.character(formula(get(model))[3]) 
  x.term <- gsub("([().])", "\\\\\\1", x.var)
  RHS.part <- gsub(sprintf("[^ ]*%s[^ ]*", x.term), "", RHS.part) # term (and any non-space characters around it)
  ##   RHS.part <- gsub("\\s\\:[^ ]+", "", RHS.part) # :term
  ##   RHS.part <- gsub("[^ ]+\\:\\s", "", RHS.part) # term:
  ##   RHS.part <- gsub("\\s[^ ]+\\:\\:[^ ]+", "", RHS.part) # term::term
  ##   RHS.part <- gsub("[^ ]+\\:$", "", RHS.part) # last term
  RHS.part1 <- RHS.part
  RHS.part  <- ""                      # to start the loop
  while (RHS.part1 !=RHS.part)
  {                                    # keep doing this until nothing changes!
    RHS.part  <- RHS.part1
    RHS.part1 <- gsub("[+*]\\s*([+*])", "\\1", RHS.part) # leftovers
  }
  RHS.part <- gsub("\\s*([+*])\\s*$", "", RHS.part) # leftovers
  Y.part <- sprintf("update(%s, .~ %s )",  model, RHS.part)
  X.part <- sprintf("update(%s, %s ~ %s )", model, x.var, RHS.part)
  ## If I return the raw text strings, it would be more "robust", and allow custom post-processing
  ## But, then it becomes harder to *use*: eval(parse(text = Y.part))
  ##   Y.part <- eval(parse(text=Y.part))
  ##   X.part <- eval(parse(text=X.part))
  Y.part <- parse(text=Y.part)
  X.part <- parse(text=X.part)
  out <- NULL
  if (part == "both") out <- list(y = Y.part, x = X.part)
  if (part %in% c("x", "X")) out <- X.part
  if (part %in% c("y", "Y")) out <- Y.part
  out
}


RegPlot.annote <- function(model, part = c("equation", "pvalue", "r.squared") )
{   # extract some common annotations to add to a regression plot - suitable for ggplot geom_text(..., parse = TRUE)
  mod.summary <- summary(model)
  mod.eq <- sprintf("y = %.1f %s paste(%.3f,x)", round(coef(model)[1], digits = 2), 
                      ifelse(coef(model)[2] > 0, "+", "-"), abs(coef(model)[2]) )
  mod.eq <- sub("-?0\\.0 ", "", mod.eq) # clean-up
  mod.eq <- sub("= \\+ ", "= ", mod.eq) # clean-up
  mod.eq <- gsub("([xy])", "italic(\\1)", mod.eq) # prep for expression
  mod.eq <- gsub(" = ", "~\"=\"~", mod.eq) # prep for expression
  Pvalue <- pf(mod.summary$fstatistic[1], mod.summary$fstatistic[2],
               mod.summary$fstatistic[3], lower.tail = FALSE) 
  mod.pv <- substitute(italic(F)[list(df1,df2)]~"="~Fval~", "~italic(p)~pchar~pval, 
                       list(df1 = mod.summary$fstatistic['numdf'],
                            df2 = mod.summary$fstatistic['dendf'],
                            Fval     = sprintf("%.2f", mod.summary$fstatistic['value']),
                            pchar    = ifelse(Pvalue < 0.001, "<", "="),
                            pval     = sprintf("%.3f", ifelse(Pvalue < 0.001, 0.001, Pvalue) )
                            )
  )
  mod.r2 <- substitute( italic(r)^2~"="~r2, list(r2 = sprintf("%.3f", summary(model)$r.squared)) )

  out <- c()
  if (length(grep("eq", part)) > 0) out <- c(out, mod.eq)
  if (length(grep("pv", part)) > 0) out <- c(out, as.character(as.expression(mod.pv)) )
  if (length(grep("r",  part)) > 0) out <- c(out, as.character(as.expression(mod.r2)))

  out
}


###=============================================================

gmean <- function (x) exp( mean( log(x) ) ) # Geometric mean == nth root of product of n values


###=============================================================
## modify summary.aov() output to include Effect Sizes for factors (SS_factor / SS_total)

summary_ES <- function (object)
{   # takes a `summary` object as the first argument, and returns a modified version with a column added for 'effect size' (ES)
  if (strsplit(class(object)[1], "\\.")[[1]][1] != "summary")
  {
    stop("The first argument to the summary.es() function must be the result of summary(). \nPlease run summary() on your fitted model object, and pass the result to summary.es().")
  }
  isaovlist <- if ("summary.aovlist" %in% class(object)) TRUE else FALSE
  if (isaovlist)
    object.flat <- unlist(object, recursive = FALSE)
  else
    object.flat <- object
  object.df <- rbind.fill(object.flat)  # make one big df: to calculate SS.total
  SS.total <- sum(object.df[["Sum Sq"]], na.rm=TRUE)    # First, I need the Total SS, which is the sum of ALL SS, over all levels.
  Insert.ES <- function (df)
  { # subroutine to use on anova data frames, wherever they happen to be
    df$ES <- df[["Sum Sq"]] / SS.total
    ESi <- which(colnames(df) == "ES")  # in case it's already there
    colI <- setdiff( 1:ncol(df), ESi )
    df <- df[, c(colI[1:3], ESi, colI[4:length(colI)]) ]   # re-arrange columns reproducibly
    df
  }
  for (i in 1:length(object))
  { # summary objects are lists, with an item for each level; loop through and add the "ES" column to each.
    if (isaovlist)
    {   ## Actually, `summary.aov` is a list of 1 data.frame; summary.aovlist is a list of `summary.aov` objects ... :/
      object[[i]][[1]] <- Insert.ES( object[[i]][[1]] )
    } else {
      object[[i]] <- Insert.ES( object[[i]] )
    }
  }
  if (FALSE)
  { # test code
    ## using examples from ?aov
    utils::data(npk, package="MASS")
    summary_ES(summary( aov(yield ~ block + N * P + K, npk) ))
    npk.aovE <- aov(yield ~  N*P*K + Error(block), npk)
    summary_ES( summary(npk.aovE) )
    ## lm : this will not work, as there is no "Sum Sq" column in the df (and the structure is totally different).
    ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
    ## Page 9: Plant Weight Data.
    ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
    trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
    group <- gl(2,10,20, labels=c("Ctl","Trt"))
    weight <- c(ctl, trt)
    lm.D9 <- lm(weight ~ group)
    summary_ES(summary(lm.D9))  # FAIL
    ## lme? ditto
    fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
    summary_ES(summary(fm1))    # FAIL
  }
  object
}

### re-define wrapper functions for summary() methods, to add ES automatically (saves having to rewrite old code! :D)
### I could just redefine summary(), but I haven't tested summary_ES() with all the different summary classes, so that would not be safe right now.
summary.aov <- function (object, ...)
{
  summary_ES( stats::summary.aov(object, ...) )
}

summary.aovlist <- function (object, ...)
{
  summary_ES( stats::summary.aovlist(object, ...) )
}
