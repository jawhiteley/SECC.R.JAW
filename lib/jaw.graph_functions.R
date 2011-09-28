##################################################
# functions for graphing with custom error bars
# Jonathan Whiteley     R v2.12     2011-08-21
##################################################
## access these functions in another file by using: 
## 	source("path/to/this/file.R")

##================================================
## MULTIPLE COMPARISONS
##================================================
MSD <- function( model, alpha=0.05, method=c("unplanned", "planned", "MSR", "LSD"), mode.df=c("pairwise", "MSE", "manual"), n=NULL ){
  ## Planned Multiple Comparisons using Least Significant Differences (LSD)
  ## -> comparison intervals for graphical display.
  ## Expects a model of class 'aov' and/or 'lm'? (aovlist)
  ## Design must be balanced (easier for now).

  method  <- match.arg(method)      # default "planned"; synonymous with "MSR"
  mode.df <- match.arg(mode.df)     # default "pairwise"
  conf.level <- 1 - (alpha/2)	# for a 2-tailed test

  ##   vars <- attr( attr(model$terms, "factors"), "dimnames")[[1]]  # dig in for the variable names ...
  ##   data <- data[, vars]  # keep only relevant columns
  ## re-factor columns to drop unused levels.
  ##   for( col in names(data) ) {
  ##     if (class(data[[col]]) == "factor") {
  ##       data[[col]] <- factor(data[[col]])
  ##     }
  ##   }

  if(mode.df == "manual"){
    n <- as.numeric(n)	# relevant group ?
  }
  if(mode.df == "pairwise"){
    ##  n <- replications( model$terms , data )	# sample sizes, according to model structure (formula).  
    ## I couldn't find an easy way to derive this directly from an aov object passed in, so this is the only reason that data is required as an argument (& formula?).  replications() returns a nasty list if the data are unbalanced :(
    tables <- model.tables(model, "means")
    means <- tables$tables
    n <- tables$n	# sample sizes, according to model structure (formula).  
    ### model.tables(model)$n  # does what I want!  Only with *full* model
    ## if the item name is a problem, use as.numeric() to convert to a pure number.  If no group specified, a (named) vector is produced with values for all treatment combinations.
  }

  model.ls <- model  # in case of aovlist

  msd.width <- c()  # container variable
  for (lvl in 1:length(n)) {  # loop through each treatment combination
    
    if ("aovlist" %in% class(model.ls)) {
      msd.lvl <- names(n)[lvl]
      model.lvl <- grep(msd.lvl, names(model.ls))  # find matching name
      ## use last level if no match
      if( length(model.lvl) == 0 ) model.lvl <- length(names(model.ls))
      model.lvl <- model.lvl[1]  # keep only first match
      model <- model.ls[[model.lvl]]
    }

    ## Calculate MSE & se for each level of nesting
    MSE <- sum(resid(model)^2)/model$df.residual	# MSE from model ( SS / df )
    msd.se <- sqrt(2*MSE/n)	# se of a difference (Crawley 2007, pg.465; Sokal & Rohlf 2003/1995, pg. 243)
    ## unbalanced: sqrt( (var[1]/n[1]) + (var[2]/n[2]) )
    
    ## A Minimum Significant Difference (MSD) = (critical value) * SE
    ## - Sokal & Rohlf, pg.243-244
    ## Student's t = a difference / se of the difference
    ## - Crawley, pg. 464
    if( method == "unplanned" || method == "MSR") {
      ngroups <-  length( means[[lvl +1]] )  # skip "Grand Mean"
      msd.df <- ngroups * (n -1)
      width <- qtukey(conf.level, ngroups, msd.df) * sqrt(MSE/n) 
      ## for unplanned comparisons, it should really be Minimum Significant Range
      ## MSR = Q[alpha,k, v] * sqrt(MSE/n)    (see Sokal & Rohlf, pg.244).
    } else {
      if(mode.df=="MSE") {
        msd.df <- model$df	# for MSE from fitted model?  This is probably completely incorrect.
      } else {
        msd.df <- (2*n)-2   # for differences between any 2 levels. 
      }

      width <- qt(conf.level, msd.df)*msd.se  # LSD based on error rate (alpha).
    }

    msd.width[lvl] <- width[lvl]
  }
  names(msd.width) <- names(n)
  return(msd.width)
}



##================================================
## PLOTS & GRAPHS
##================================================
##________________________________________________
## text_to_plot: mtext() wrapper.  Outputs text to a plotting device.
## useful for getting some text output into a pdf
## file with graphs.

text_to_plot <- function( txt = "", auto.size = TRUE,
                         cex = 0.6, family = "mono" , ... ) {
  ## can fit about 60 (56) letters @ cex = 1
  ## auto-word-wrap?  Insert "\n" every X characters, depending on cex?
  plot.new()
  mtext(txt,
        adj = 0, side = 3, line = 0,
        cex = cex, family = family, col = "black",
        ...
        )
}



##________________________________________________
## plotMeans: modified version from 'R commander' package (Rcmdr).  I just hate loading the entire package & GUI, just for this one convenient function.  
# I added an option for custom error bar widths: set 'error.bars' to "custom", and pass the custom width via the 'level' argument (sort of hack).  TO DO: allow a vector of upper & lower error bars for each point.
# I also added the ability to accept more graphing parameters: bg, cex, lwd.  Could access to these & other graphical parameters be provided more easily with the '...' (or '…') argument?
plotMeans <- function (response, factor1, factor2, 
  error.bars = c("se", "sd", "conf.int", "none", "custom"), 
  level = 0.95, xlab = deparse(substitute(factor1)), 
  ylab = paste("mean of", deparse(substitute(response))), 
  legend.lab = deparse(substitute(factor2)),
  main = "Plot of Means", pch = NULL, lty = NULL,
  col = palette(), bg=NULL, cex=2, lwd=NULL, ylim=NULL, ... ) 
{
    if (!is.numeric(response)) 
        stop(gettextRcmdr("Argument response must be numeric."))
    xlab
    ylab
    legend.lab
    error.bars <- match.arg(error.bars)
    if (missing(factor2)) {
        if (!is.factor(factor1)) 
            stop(gettextRcmdr("Argument factor1 must be a factor."))
        valid <- complete.cases(factor1, response)
        factor1 <- factor1[valid]
        response <- response[valid]
        means <- tapply(response, factor1, mean)
        sds <- tapply(response, factor1, sd)
        ns <- tapply(response, factor1, length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        if (error.bars == "custom") 
            sds <- level ## custom error bar widths.  Added by JAW.  
            # custom widths passed as an appropriate matrix, or a single value replicated: matrix( level, nrow = dim(sds)[1], ncol = dim(sds)[2] )
        sds[is.na(sds)] <- 0
        yrange <- if (error.bars != "none") 
            c(min(means - sds, na.rm = TRUE), max(means + sds, 
                na.rm = TRUE))
        else range(means, na.rm = TRUE)
        yrange <- if ( is.null(ylim)==FALSE ) ylim else yrange # added by JAW
        levs <- levels(factor1)
        n.levs <- length(levs)
        if (is.null(pch)) pch <- 1
        if (is.null(lty)) lty <- 1
        plot(c(1, n.levs), yrange, type = "n", xlab = xlab, ylab = ylab, 
            axes = FALSE, main = main, ... )	# `...` argument added by JAW
        points(1:n.levs, means, type = "b", pch = pch, cex = cex, bg = bg, lty=lty, lwd=lwd)	# originally: pch=16, cex=2.  Modified by JAW.
        box()
        axis(2)
        axis(1, at = 1:n.levs, labels = levs)
        if (error.bars != "none") 
            arrows(1:n.levs, means - sds, 1:n.levs, means + sds, 
                angle = 90, lty = 2, code = 3, length = 0.125)
    }
    else {
        if (!(is.factor(factor1) | is.factor(factor2))) 
            stop(gettextRcmdr("Arguments factor1 and factor2 must be factors."))
        valid <- complete.cases(factor1, factor2, response)
        factor1 <- factor1[valid]
        factor2 <- factor2[valid]
        response <- response[valid]
        means <- tapply(response, list(factor1, factor2), mean)
        sds <- tapply(response, list(factor1, factor2), sd)
        ns <- tapply(response, list(factor1, factor2), length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        if (error.bars == "custom") 
            sds <- level ## custom error bar widths.  Added by JAW.
        sds[is.na(sds)] <- 0
        yrange <- if (error.bars != "none") 
            c(min(means - sds, na.rm = TRUE), max(means + sds, 
                na.rm = TRUE))
        else range(means, na.rm = TRUE)
        yrange <- if ( is.null(ylim)==FALSE ) ylim else yrange # added by JAW
        levs.1 <- levels(factor1)
        levs.2 <- levels(factor2)
        n.levs.1 <- length(levs.1)
        n.levs.2 <- length(levs.2)
        if (is.null(pch)) pch <- 1:n.levs.2
        if (is.null(lty)) lty <- 1:n.levs.2
        if (length(pch) == 1) 
            pch <- rep(pch, n.levs.2)
        if (length(col) == 1) 
            col <- rep(col, n.levs.2)
        if (length(lty) == 1) 
            lty <- rep(lty, n.levs.2)
        ## More graphical parameters added by JAW.
        if (length(bg) == 1) 
            bg <- rep(bg, n.levs.2)
        if (length(lwd) == 1) 
            lwd <- rep(lwd, n.levs.2)
        if (length(cex) == 1) 
            cex <- rep(cex, n.levs.2)
        ## ##
        if (n.levs.2 > length(col)) 
            stop(sprintf(gettextRcmdr("Number of groups for factor2, %d, exceeds number of distinct colours, %d."), 
                n.levs.2, length(col)))
        plot(c(1, n.levs.1 * 1.4), yrange, type = "n", xlab = xlab, 
            ylab = ylab, axes = FALSE, main = main, ... )	# `...` argument added by JAW
        box()
        axis(2)
        axis(1, at = 1:n.levs.1, labels = levs.1)
        for (i in 1:n.levs.2) {
            points(1:n.levs.1, means[, i], type = "b", pch = pch[i], 
                cex = cex[i], col = col[i], lty = lty[i], lwd = lwd[i], bg = bg[i] )
            if (error.bars != "none") 
                arrows(1:n.levs.1, means[, i] - sds[, i], 1:n.levs.1, 
                  means[, i] + sds[, i], angle = 90, code = 3, 
                  col = col[i], lty = lty[i], length = 0.125)
        }
        x.posn <- n.levs.1 * 1.1
        y.posn <- sum(c(0.1, 0.9) * par("usr")[c(3, 4)])
        text(x.posn, y.posn, legend.lab, adj = c(0, -0.5))
        legend(x.posn, y.posn, levs.2, pch = pch, col = col, 
            lty = lty, pt.bg=bg, pt.lwd=lwd, lwd=lwd, pt.cex=cex*0.75 )
    }
    invisible(NULL)
}


##================================================
## PANEL FUNCTIONS
##================================================

## Coplots with linear fits (from Zuur et al. 2007 Chapter 22 R code)
panel.lines2 <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
    cex = 1, ...)
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)){
        tmp=lm(y[ok]~x[ok])
        abline(tmp)}

}


##================================================
## GGPLOT2 custom options
##================================================

jaw.ggplot <- function(x) {
  require(ggplot2)
  ## add common custom options to ggplot2 objects.
  ## I got tired of copying the same code over &over again.
  x <- x + theme_bw()   # change to white background colour scheme.
  x <- x + opts(legendkey = theme_rect(colour = NA))    # remove rectangle from legend keys.
  return(x)
}
