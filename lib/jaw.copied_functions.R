##################################################
# functions written by people other than me, which I copied
# (copied / acquired / adopted / etc.)
# Jonathan Whiteley		R v2.12		2011-01-26
##################################################

cleanVarName <- function(variable.name)
{
  # clean.var.name by John Myles White, package "ProjectTemplate"
  variable.name <- gsub('_', '.', variable.name, perl = TRUE)
  variable.name <- gsub('-', '.', variable.name, perl = TRUE)
  variable.name <- gsub('\\s+', '.', variable.name, perl = TRUE)
  return(variable.name)
}


################################################################
## Functions from "Numerical Ecology in R"  http://www.bio.umontreal.ca/numecolR/
################################################################
# evplot()
# Plot eigenvalues and percentages of variation
# Kaiser rule and broken stick model
#
# License: GPL-2 
# Author: Francois Gillet, Octber 2007

evplot = function(ev) {
	# Broken stick model (MacArthur 1957)
	n = length(ev)
	bsm = data.frame(j=seq(1:n), p=0)
	bsm$p[1] = 1/n
	for (i in 2:n) {
		bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
	}
	bsm$p = 100*bsm$p/n
	# Plot eigenvalues and % of variation for each axis
	op = par(mfrow=c(2,1))
	barplot(ev, main="Eigenvalues", col="bisque", las=2)
	abline(h=mean(ev), col="red")
	legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
	barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=T, 
		main="% variation", col=c("bisque",2), las=2)
	legend("topright", c("% eigenvalue", "Broken stick model"), 
		pch=15, col=c("bisque",2), bty="n")
	par(op)
}

'cleanplot.pca' <- function(res.pca, ax1=1, ax2=2, point=FALSE, ahead=0.07, cex=0.7) {

# A function to draw biplots from a PCA done with vegan.
# res.pca: an object of class "rda" (PCA or RDA result from vegan)
#
# License: GPL-2 
# Authors: Francois Gillet & Daniel Borcard, April 2010
# tweaked by Jonathan Whiteley (for personal preference)

  require("vegan")

  # Two PCA biplots: scaling 1 and scaling 2
  # ****************************************

  par(mfrow=c(1,2))
  p <- length(res.pca$CA$eig)

  # Scaling 1: "species" scores scaled to relative eigenvalues
  sit.sc1 <- scores(res.pca, display="wa", scaling=1, choices=c(1:p))
  spe.sc1 <- scores(res.pca, display="sp", scaling=1, choices=c(1:p))
  ## plot.cca() throws an error (can't find the function) - jaw
  plot(res.pca, choices=c(ax1,ax2), display=c("wa","sp"), type="n", 
	 main="PCA - scaling 1 (distance)", scaling=1)
  if (point) {
    points(sit.sc1[,ax1], sit.sc1[,ax2], pch=20)
    ##     text(res.pca, display="wa", choices=c(ax1,ax2), cex=cex, pos=3, scaling=1)
  }
  else {
    text(res.pca, display="wa", choices=c(ax1,ax2), cex=cex, scaling=1)
  }
  text(res.pca, display="sp", choices=c(ax1,ax2), cex=cex, pos=4, 
    col="red", scaling=1)
  arrows(0, 0, spe.sc1[,ax1], spe.sc1[,ax2], length=ahead, angle=20, col="red")
  pcacircle(res.pca)

  # Scaling 2: site scores scaled to relative eigenvalues
  sit.sc2 <- scores(res.pca, display="wa", choices=c(1:p))
  spe.sc2 <- scores(res.pca, display="sp", choices=c(1:p))
  plot(res.pca, choices=c(ax1,ax2), display=c("wa","sp"), type="n", 
  	main="PCA - scaling 2 (correlation)")
  if (point) {
    points(sit.sc2[,ax1], sit.sc2[,ax2], pch=20)
    ##     text(res.pca, display="wa", choices=c(ax1,ax2), cex=cex, pos=3)
  }
  else {
    text(res.pca, display="wa", choices=c(ax1,ax2), cex=cex)
  }
  text(res.pca, display="sp", choices=c(ax1,ax2), cex=cex, pos=4, col="red")
  arrows(0, 0, spe.sc2[,ax1], spe.sc2[,ax2], length=ahead, angle=20, col="red")
}



'pcacircle' <- function (pca) {

# Draws a circle of equilibrium contribution on a PCA plot 
# generated from a vegan analysis.
# vegan uses special constants for its outputs, hence 
# the 'const' value below.

    eigenv <- pca$CA$eig
    p <- length(eigenv)
    n <- nrow(pca$CA$u)
    tot <- sum(eigenv)
    const <- ((n - 1) * tot)^0.25
    radius <- (2/p)^0.5
    radius <- radius * const
    symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)    
}


`PCA` <- 
   function(Y, stand=FALSE)
# 
# Principal component analysis (PCA) with option for variable standardization
#
# stand = FALSE : center by columns only, do not divide by s.d.
# stand = TRUE  : center and standardize (divide by s.d.) by columns
#
# License: GPL-2 
# Author: Pierre Legendre, May 2006
{
   Y = as.matrix(Y)
   obj.names = rownames(Y)
   var.names = colnames(Y)
   size = dim(Y)
   Y.cent = apply(Y, 2, scale, center=TRUE, scale=stand)
   Y.cov = cov(Y.cent)
   Y.eig = eigen(Y.cov)
   k = length(which(Y.eig$values > 1e-10))
   U  = Y.eig$vectors[,1:k]
   F  = Y.cent %*% U
   U2 = U %*% diag(Y.eig$value[1:k]^(0.5))
   G  = F %*% diag(Y.eig$value[1:k]^(-0.5))
   rownames(F)  = obj.names
   rownames(U)  = var.names
   rownames(G)  = obj.names
   rownames(U2) = var.names
   axenames <- paste("Axis",1:k,sep=" ")
   colnames(F)  = axenames
   colnames(U)  = axenames
   colnames(G)  = axenames
   colnames(U2) = axenames

#
# Fractions of variance
   varY = sum(diag(Y.cov))
   eigval = Y.eig$values[1:k]
   relative = eigval/varY
   rel.cum = vector(length=k)
   rel.cum[1] = relative[1]
   for(kk in 2:k) { rel.cum[kk] = rel.cum[kk-1] + relative[kk] }
#
out <- list(total.var=varY, eigenvalues=eigval, rel.eigen=relative, 
       rel.cum.eigen=rel.cum, U=U, F=F, U2=U2, G=G, stand=stand, 
       obj.names=obj.names, var.names=var.names, call=match.call() )
class(out) <- "PCA"
out
}

`print.PCA` <-
    function(x, ...)
{
    cat("\nPrincipal Component Analysis\n")
    cat("\nCall:\n")
    cat(deparse(x$call),'\n')
    if(x$stand) cat("\nThe data have been centred and standardized by column",'\n')
    cat("\nTotal variance in matrix Y: ",x$total.var,'\n')
    cat("\nEigenvalues",'\n')
    cat(x$eigenvalues,'\n')
    cat("\nRelative eigenvalues",'\n')
    cat(x$rel.eigen,'\n')
    cat("\nCumulative relative eigenvalues",'\n')
    cat(x$rel.cum.eigen,'\n')
    invisible(x) 
}

`biplot.PCA` <-
    function(x, scaling=1, plot.axes=c(1,2), color.obj="black", color.var="red", ...)
# scaling = 1 : preserves Euclidean distances among the objects
# scaling = 2 : preserves correlations among the variables
{
    #### Internal function
	larger.frame <- function(mat, percent=0.07)
	# Produce an object plot 10% larger than strictly necessary
	{
	range.mat = apply(mat,2,range)
	z <- apply(range.mat, 2, function(x) x[2]-x[1])
	range.mat[1,]=range.mat[1,]-z*percent
	range.mat[2,]=range.mat[2,]+z*percent
	range.mat
	}
	####
	
	if(length(x$eigenvalues) < 2) stop("There is a single eigenvalue. No plot can be produced.")
	if(length(which(scaling == c(1,2))) == 0) stop("Scaling must be 1 or 2")

	par(mai = c(1.0, 0.75, 1.0, 0.5))

	if(scaling == 1) {

	# Distance biplot, scaling type = 1: plot F for objects, U for variables
	# This projection preserves the Euclidean distances among the objects
	
	lf.F = larger.frame(x$F[,plot.axes])
	biplot(x$F[,plot.axes],x$U[,plot.axes],col=c(color.obj,color.var), xlim=lf.F[,1], ylim=lf.F[,2], arrow.len=0.05, asp=1)
	title(main = c("PCA biplot","scaling type 1"), family="serif", line=3)

	} else {

	# Correlation biplot, scaling type = 2: plot G for objects, U2 for variables
	# This projection preserves the correlation among the variables
   
	lf.G = larger.frame(x$G[,plot.axes])
	biplot(x$G[,plot.axes],x$U2[,plot.axes],col=c(color.obj,color.var), xlim=lf.G[,1], ylim=lf.G[,2], arrow.len=0.05, asp=1)
	title(main = c("PCA biplot","scaling type 2"), family="serif", line=3)

	}
invisible()
}

# mite.hel = decostand(mite, "hel")
# mite.hel.D = dist(mite.hel)
# mite.correlog = mantel.correlog(mite.hel.D, XY=mite.xy, nperm=99, cutoff=FALSE)
# mite.correlog
# mite.correlog$mantel.res
# plot(mite.correlog)



