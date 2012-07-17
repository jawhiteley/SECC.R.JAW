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
  RHS.part <- gsub("[+*]\\s*([+*])", "\\1", RHS.part) # leftovers
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
