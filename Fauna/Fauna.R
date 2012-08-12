##################################################
### Schefferville Experiment on Climate Change (SEC-C)
### Analyses of Fauna data (microarthropod morphospecies counts)
### Jonathan Whiteley     R v2.12     2012-07-31
###===============================================
### Species identified to morphospecies (usually family-level)
### Counts / sample converted to # / g dwt of moss substrate (using 'Patch.dwt' column)
##################################################
## INITIALISE
##################################################
if (FALSE) {  # do not run automatically
  ## Working Directory: see lib/init.R below
  setwd("./ SECC/")                    # relative to my usual default wd in R GUI (Mac).
  setwd("..")                          # relative to this file (\rd in vim).
  getwd()                              # Check that we're in the right place
}

## Load data, functions, etc.  Includes rm(list=ls()) to clear memory
source('./Fauna/Fauna_setup.R')




##################################################
### UNIVARIATE ANALYSES (ANOVAs)
##################################################
cat("Univariate Analyses of Fauna data ...\n")
SECCxtract.aovlist <- function(x)
{
  if ("aovlist" %in% class(x) == FALSE) 
    stop("SECCxtract.aovlist() expects an object of class 'aovlist'")
  xs <- summary(x)
  for (i in names(xs))
  {
    lvl <- xs[[i]][[1]]
    lvl.out <- lvl[, c("Df", "F value", "Pr(>F)")] 
    res.row <- grep("Residuals", rownames(lvl.out))
    res.df <- lvl.out[res.row, "Df"]
    lvl.out$Df.resid <- res.df
    lvl.out <- lvl.out[, c("Df", "Df.resid", "F value", "Pr(>F)")]
    if (i == names(xs)[1])
      out <- lvl.out
    else
      out <- rbind(out, lvl.out)
  }
  ## cleanup
  empty.rows <- which( is.na(out[, "F value"]) )
  out <- out[-empty.rows, ]
  rownames(out) <- gsub("Frag", "Fragmentation", rownames(out))
  colnames(out) <- c("Df", "Df.resid", "F", "pF")
  return(out)
}
SECCcbind.aovlist <- function(FullANOVA = Fauna.ANOVA, newANOVA = SECCxtract.aovlist(Yp.aov) )
{
  FullANOVA
}

Y.col <- 'Richness' # Column to analyze as response variable           *****
Y.use <- 'Y'        # Which transformation is being used (for labels)? *****
Y.lim1 <- c(0, 20)  # consistent Y limits :/
source('Fauna/Fauna-univariate.R')
ANOVAtable <- SECCxtract.aovlist(Yp.aov)
Fauna.ANOVA <- ANOVAtable
Ycols <- grep("^p?F$", colnames(Fauna.ANOVA))
colnames(Fauna.ANOVA)[Ycols] <- paste(Y.col, colnames(Fauna.ANOVA)[Ycols], sep=".")

Y.col <- 'Evenness'
Y.use <- 'Y'
Y.lim1 <- c(0, 10)  # consistent Y limits :/
source('Fauna/Fauna-univariate.R')
ANOVAtable <- SECCxtract.aovlist(Yp.aov)
if (ANOVAtable$Df == Fauna.ANOVA$Df)
{
  Fauna.ANOVA <- cbind(Fauna.ANOVA, ANOVAtable[, c("F", "pF")])
  Ycols <- grep("^p?F$", colnames(Fauna.ANOVA))
  colnames(Fauna.ANOVA)[Ycols] <- paste(Y.col, colnames(Fauna.ANOVA)[Ycols], sep=".")
}

Y.col <- 'Predators'
Y.use <- 'Y'
Y.lim1 <- c(0, 25)  # consistent Y limits :/
source('Fauna/Fauna-univariate.R')
ANOVAtable <- SECCxtract.aovlist(Yp.aov)
if (ANOVAtable$Df == Fauna.ANOVA$Df)
{
  Fauna.ANOVA <- cbind(Fauna.ANOVA, ANOVAtable[, c("F", "pF")])
  Ycols <- grep("^p?F$", colnames(Fauna.ANOVA))
  colnames(Fauna.ANOVA)[Ycols] <- paste(Y.col, colnames(Fauna.ANOVA)[Ycols], sep=".")
}

## Chamber x Frag
Chamber.label <- "Chamber" # attr(SECC, "labels")[["Pos"]]
Chamber.map <- plotMap( factor = "Chamber", labels = levels(SECC$Chamber) )
Chamber.map <- Chamber.map[ levels(SECC$Chamber) %in% Chamber.use, ]
Chamber.map[2, "label"] <- "Chamber"
plot.means <- SECCplotDataANOVA(SECCp$Y.trans, 
                                list(Chamber=SECCp$Chamber, Frag=SECCp$Frag, 
                                     Time=SECCp$Time), 
                                error = msd["Chamber:Frag"]
                                )
levels(plot.means$Chamber)[2] <- "Chamber"

CxF.plot <- qplot(Frag, x, data = plot.means, group = Chamber, 
                  geom = "line", ylim = Y.lim, size = Chamber, 
                  colour = Chamber, fill = Chamber, 
                  shape = Chamber, # lty = Chamber,
                  main = Plot.Title, sub = Sub.msd,
                  xlab = attr(SECC, "labels")[["Chamber"]],
                  ylab = Y.plotlab,
                  legend = FALSE)
CxF.plot <- CxF.plot + geom_errorbar(aes(ymin = lower, ymax = upper), 
                                     width = 0.2, size = 0.5)
CxF.plot <- CxF.plot + geom_point(aes(group = Chamber), size = 3)
CxF.plot <- CxF.plot + scale_colour_manual(name = Chamber.label,
                                           values = Chamber.map$col, 
                                           breaks = Chamber.map$label)
CxF.plot <- CxF.plot + scale_fill_manual(name = Chamber.label,
                                         values = Chamber.map$bg, 
                                         breaks = Chamber.map$label)
CxF.plot <- CxF.plot + scale_shape_manual(name = Chamber.label,
                                          values = Chamber.map$pch, 
                                          breaks = Chamber.map$label)
CxF.plot <- CxF.plot + scale_size_manual(name = Chamber.label,
                                         values = Chamber.map$lwd*0.5, 
                                         breaks = Chamber.map$label)
## CxF.plot <- CxF.plot + scale_linetype_manual(name = Chamber.label,
##                                              values = Chamber.map$lty, 
##                                              breaks = Chamber.map$label)
CxF.plot <- CxF.plot + jaw.ggplot() 
CxF.plot <- CxF.plot + scale_x_discrete(labels = names(FragIconList), # c(1, 2, 4), 
                                        breaks = levels(plot.means$Frag)) +
     opts(axis.ticks.margin = unit(0.2, "lines"),
          axis.text.x = picture_axis(FragIconList, icon.size = unit(1.4, "lines")) 
     )
print(CxF.plot)

ggsave(file = paste(Save.final, "- CxF.eps"), plot = CxF.plot, width = 4, height = 4, scale = 1.2)


Y.col <- 'Grazers'
Y.use <- 'Y.sqrt'
Y.lim1 <- c(0, 25)  # consistent Y limits :/
source('Fauna/Fauna-univariate.R')
ANOVAtable <- SECCxtract.aovlist(Yp.aov)
if (ANOVAtable$Df == Fauna.ANOVA$Df)
{
  Fauna.ANOVA <- cbind(Fauna.ANOVA, ANOVAtable[, c("F", "pF")])
  Ycols <- grep("^p?F$", colnames(Fauna.ANOVA))
  colnames(Fauna.ANOVA)[Ycols] <- paste(Y.col, colnames(Fauna.ANOVA)[Ycols], sep=".")
}

if (TRUE)
{                                      # Totally NS: sample sizes too small to detect anything meaningful
  Y.col <- 'Richness' # Column to analyze as response variable           *****
  Y.use <- 'Y'        # Which transformation is being used (for labels)? *****
  Y.lim1 <- c(0, 26)  # consistent Y limits :/
  source('Fauna/Fauna-regional.R')
  Fauna.mcANOVA <- SECCxtract.aovlist(Ymc.aov)
}




##################################################
### MULTIVARIATE ANALYSES
##################################################
cat("Multivariate Analyses of Fauna data ...\n")
## Transformation
Fauna.log   <- log( Fauna +1 )         # or log1p
Fauna.max   <- decostand( Fauna, method = "max" ) 
Fauna.hell  <- decostand( Fauna, method = "hellinger" ) 
Fauna.chord <- decostand( Fauna, method = "normalize" ) ## Chord transformation
Fauna.pa    <- decostand( Fauna, method = "pa" ) 

Fauna.trans <- Fauna.log

## Distance/Similarity matrix (can be included in MDS call)
Fauna.bray  <- vegdist(Fauna.trans, method = "bray") # add `binary = TRUE` for Sorensen dissimilarity ;)
Fauna.jacc  <- vegdist(Fauna.trans, method = "jaccard", binary = TRUE) # presence/absence data
Fauna.hellD <- vegdist(Fauna.hell,  method = "euclidean")
Fauna.chorD <- vegdist(Fauna.chord, method = "euclidean") # Chord Distance
Fauna.euclD <- vegdist(Fauna, method = "euclidean") # euclidean Distance

Fauna.dist <- Fauna.bray               # default


##################################################
## NMDS
##################################################
## Distance/Similarity matrix (included in MDS call)
## Compute MDS
Fauna.mds  <- metaMDS( Fauna.trans, distance = "bray", autotransform = FALSE, k = 2)
Fauna.mds3 <- metaMDS( Fauna.trans, distance = "bray", autotransform = FALSE, k = 3)

## Check MDS quality (Stress)
## stress.label <- paste( "Stress = ", round( Fauna.mds$stress, 1), "%", sep="" )
## has vegan changed?  stress values now seem to be <1, rather than in % ?!?
stress.label <- paste( "Stress = ", round( Fauna.mds$stress, 2), sep="" )
stressplot(Fauna.mds, main=stress.label)

stress.label <- paste( "Stress = ", round( Fauna.mds3$stress, 2), sep="" )
stressplot(Fauna.mds3, main=stress.label)

## Look at the pretty graph
Fauna.plot <- ordiplot( Fauna.mds, type="text", main = "Bray-Curtis Distance" )

## Alternate transformations & distances
Jacc.mds   <- metaMDS( Fauna.pa,    distance = "jaccard",   autotransform = FALSE, k = 2)
Chord.mds  <- metaMDS( Fauna.chord, distance = "euclidean", autotransform = FALSE, k = 2)
HellD.mds  <- metaMDS( Fauna.hell,  distance = "euclidean", autotransform = FALSE, k = 2)
ordiplot( Jacc.mds, type="text" , main = "Jaccard Distance" )
ordiplot( Chord.mds, type="text", main = "Chord Distance" )
ordiplot( HellD.mds, type="text", main = "Hellinger Distance" )



##################################################
## ANOSIM
##################################################
## Chamber, Pos, Frag, 
## Frag:1x2, Frag:2x4, Frag:1x4, 
## Chamber:A x Pos, Chamber:A x Frag,
## Chamber:C x Pos, Chamber:C x Frag,
## Pos:I x Frag, Pos:O x Frag,
## X Chamber:A x Pos:I x Frag,
## Chamber:C x Pos:I x Frag,
## Chamber:C x Pos:O x Frag,
## split = ",\\s\\n?\\s*", fixed = FALSE)[[1]]

## Assemble a data frame of specific group comparisions x distance metrics
ANOSIM.results <- data.frame(Groups = c())
dist.methods <- c("bray", "bray", "euclidean", "euclidean", "jaccard")
dist.trans   <- c("",     "log",  "hell",      "chord",     "pa") # transformations
for (d in 1:length(dist.methods))
{
  dista <- dist.methods[d]
  trans <- dist.trans[d]
  col.label <- paste(dista, trans, sep="_")
  Fauna.t <- paste("Fauna", trans, sep=".")
  if (Fauna.t == "Fauna.") Fauna.t <- "Fauna"
  Fauna.t <- get( Fauna.t )
  cat("ANOSIM: using distance: \'", dista, "\', with transformation: \'", 
      ifelse(trans == "", "[none]", trans), "\'\n", 
      sep="")
  dist.g <- c()                        # group labels
  dist.R <- c()
  dist.p <- c()
  ## ANOSIM by single factor (basically useless, since Frag is NS as a main factor)
  cat("- Single Factor ANOSIMs\n")
  Fauna.d <- vegdist(Fauna.trans, method = dista)
  Fauna.Chamber.anosim <- anosim( Fauna.d, Fauna.sp$Chamber )
  Fauna.Pos.anosim     <- anosim( Fauna.d, Fauna.sp$Pos )
  Fauna.Frag.anosim    <- anosim( Fauna.d, Fauna.sp$Frag )
  dist.R <- c(dist.R, 
              Fauna.Chamber.anosim$statistic,
              Fauna.Pos.anosim$statistic,
              Fauna.Frag.anosim$statistic
              )
  dist.p <- c(dist.p, 
              Fauna.Chamber.anosim$signif,
              Fauna.Pos.anosim$signif,
              Fauna.Frag.anosim$signif
              )
  dist.g <- c(dist.g, "Chamber", "Pos", "Frag")
  ## Pairs of Frag levels
  cat("- among Frag levels\n")
  for (lvl in c("1x2", "2x4", "1x4")) 
  {
    g.rows <- which(Fauna.sp$Frag %in% strsplit(lvl, "x")[[1]])
    g.lvls <- Fauna.sp[g.rows, c("Frag")]
    Fauna.data <- Fauna.trans[g.rows, ]
    Fauna.d <- vegdist(Fauna.data, method = dista)
    Fauna.Fraglvl.anosim <- anosim( Fauna.d, g.lvls )
    dist.R <- c(dist.R, Fauna.Fraglvl.anosim$statistic)
    dist.p <- c(dist.p, Fauna.Fraglvl.anosim$signif)
    dist.g <- c(dist.g, paste("Frag:", lvl, sep="") )
  }
  ## ANOSIM within main factors (interactions)
  cat("- ANOSIM within levels of main factors (2-way interactions)\n")
  cat("  - ANOSIM within Chamber\n")
  for (lvl in levels(Fauna.sp$Chamber)) 
  {
    g.rows <- which(Fauna.sp$Chamber == lvl) 
    g.lvls <- Fauna.sp[g.rows, c("Chamber")] 
    Fauna.data <- Fauna.trans[g.rows, ]
    Fauna.d <- vegdist(Fauna.data, method = dista)
    ## Chamber x Pos
    Fauna.CxP.anosim <- anosim( Fauna.d, Fauna.sp$Pos[g.rows] )
    ## Chamber x Frag
    Fauna.CxF.anosim <- anosim( Fauna.d, Fauna.sp$Frag[g.rows] )
    dist.R <- c(dist.R, 
                Fauna.CxP.anosim$statistic,
                Fauna.CxF.anosim$statistic
                )
    dist.p <- c(dist.p, 
                Fauna.CxP.anosim$signif,
                Fauna.CxF.anosim$signif
                )
    dist.g <- c(dist.g, 
                paste("Chamber:", lvl, " x Pos",  sep = ""),
                paste("Chamber:", lvl, " x Frag", sep = "")
                )
  }
  cat("  - ANOSIM within Pos\n")
  for (lvl in levels(Fauna.sp$Pos)) 
  {
    g.rows <- which(Fauna.sp$Pos == lvl) 
    g.lvls <- Fauna.sp[g.rows, c("Pos")] 
    Fauna.data <- Fauna.trans[g.rows, ]
    Fauna.d <- vegdist(Fauna.data, method = dista)
    ## Pos x Frag
    Fauna.PxF.anosim <- anosim( Fauna.d, Fauna.sp$Frag[g.rows] )
    dist.R <- c(dist.R, 
                Fauna.PxF.anosim$statistic
                )
    dist.p <- c(dist.p, 
                Fauna.PxF.anosim$signif
                )
    dist.g <- c(dist.g, 
                paste("Pos:", lvl, " x Frag",  sep = "")
                )
  }
  cat("  - ANOSIM: Frag within Chamber x Pos\n")
  for (lvl in levels(Fauna.sp$Pos))
  {
    g.rows <- which(Fauna.sp$Pos == lvl & Fauna.sp$Chamber == "C") 
    g.lvls <- Fauna.sp[g.rows, c("Frag")] 
    Fauna.data <- Fauna.trans[g.rows, ]
    Fauna.d <- vegdist(Fauna.data, method = dista)
    ## Chamber x Pos x Frag
    Fauna.CxPxF.anosim <- anosim( Fauna.d, Fauna.sp$Frag[g.rows] )
    dist.R <- c(dist.R, 
                Fauna.CxPxF.anosim$statistic
                )
    dist.p <- c(dist.p, 
                Fauna.CxPxF.anosim$signif
                )
    dist.g <- c(dist.g, 
                paste("Chamber:C x Pos:", lvl, " x Frag",  sep = "")
                )
  }
  cat("Assembling ANOSIM results for method:", dista, "\n")
  if (nrow(ANOSIM.results) != length(dist.g)) 
    ANOSIM.results <- data.frame(Groups = dist.g)      # labels
  ##   ANOSIM.results$Groups <- dist.g      # labels
  ANOSIM.results <- cbind(ANOSIM.results, dist.R, dist.p)
  colnames(ANOSIM.results)[colnames(ANOSIM.results) == "dist.R"] <- 
    paste(col.label, ".R", sep="")
  colnames(ANOSIM.results)[colnames(ANOSIM.results) == "dist.p"] <- 
    paste(col.label, ".p", sep="")
}




##==============================================================
## SIMPER (slow)
## using my version of the simper algorithm in SECC.functions, which uses a few slow `for` loops
any( rownames(Fauna.trans) != rownames(Fauna.sp) ) # just checking: should be FALSE ;)
Chamber.simp <- simper(Fauna.trans, Fauna.sp$Chamber)
ChamberFrag <- paste(Fauna.sp$Chamber, Fauna.sp$Frag, sep=".")
ChamberPosn <- paste(Fauna.sp$Chamber, Fauna.sp$Pos , sep=".")
FragPosn    <- paste(Fauna.sp$Frag, Fauna.sp$Pos , sep=".")
ChamberFrag.simp <- simper(Fauna.trans, ChamberFrag)
ChamberPosn.simp <- simper(Fauna.trans, ChamberPosn)
Chamber.rows     <- which(Fauna.sp$Chamber == "C")
CxFragPos.simp   <- simper(Fauna.trans[Chamber.rows, ], FragPosn[Chamber.rows])



##==============================================================
## SPECIES OVERLAP / UNIQUENESS (BETA-diversity)
##==============================================================
## I just want to know how many species are unique to a group of sites, and how many are shared between two groups.  
## Surely, there is a package or routine that can do this automatically?
All.sp <- colnames(Fauna)
CO.rows <- which(Fauna.sp$Chamber == "C" & Fauna.sp$Pos == "O")
CO.sp <- apply(Fauna[CO.rows, ] > 0, 2, sum)
CO.sp <- names(CO.sp)[CO.sp > 0]
XCO.sp <- apply(Fauna[-CO.rows, ] > 0, 2, sum)
XCO.sp <- names(XCO.sp)[XCO.sp > 0]
COsp <- setdiff(CO.sp, XCO.sp)         # Isotomidae_sp.5 is found at a single site, in this group.
rownames(Fauna)[Fauna[[COsp]] > 0]

## betadisper test for differences in multivariate *dispersion* (variance)
CO.bdisp <- betadisper(betadiver(Fauna, "z"), ChamberPosn)
plot(CO.bdisp)
TukeyHSD(CO.bdisp)                     # Outer Chambers are not just different: they are more variable.

## What I want is a matrix of groups (as rows & columns), with cells that give:
## - # unique species in that group, or pairs of groups (triangular with diagonal)
## - # shared species in BOTH groups (triangular; diagonal is # species)
## - # species in a group NOT in the other (full matrix; diagonal is 0)
## * first two could be combined with one of the diagonals included in the third; 2 separate matrices / tables?





##==============================================================
## CA
##==============================================================
Fauna.ca <- cca(Fauna.trans)
ordiplot(Fauna.ca, scaling = 2, type = "text", xlim = c(-2, 12))
## Different species pop out vs. the SIMPER analyses; different distance metric
## Although the graphs are nice and show complex relationships easily,
## I'm losing confidence that CA / CCA is the most appropriate for this data,
## or, I need a different transformation for it.



##==============================================================
## PCoA
##==============================================================
## Yay, I get to use whatever distance measure I want :D
## Borcard et al. (2011) pg. 142
Fauna.pcoa <- cmdscale(Fauna.dist, k = 2, eig = TRUE)
Fauna.pcoa.spscores <- wascores(Fauna.pcoa$points[, 1:2], Fauna.trans) # species 'scores'

ordiplot(Fauna.pcoa, type = "text")
text(Fauna.pcoa.spscores, rownames(Fauna.pcoa.spscores), cex = 0.8, col = "red", srt = 20)




################################################################
## Correlations
################################################################
cat("Computing fauna species correllations\n")

cols.numeric <- sapply(SECC.sp.sum, function (x) class(x) == "numeric")
Sum.cor <- cor(SECC.sp.sum[, cols.numeric], method = "spearman")     # spearman rank correlation - non-parametric, doesn't have to be linear, just monotonic
Spp.cor <- cor(SECC.sp, method = "spearman")
summary(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE) # lower triangle of the correlation matrix (lower.tri) :D
Spp.meancor       <- mean(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE)
Spp.meancor.sd    <- sd(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE)
Spp.meancor.error <- sd(Spp.cor[lower.tri(Spp.cor)], na.rm = TRUE) / 
                    sqrt( length( which(!is.na(Spp.cor[lower.tri(Spp.cor)])) ) )
Spp.cortest       <- t.test( Spp.cor[lower.tri(Spp.cor)] )
## Most species are positively correlated :P
library(MASS)
truehist(Spp.cor[lower.tri(Spp.cor)], xlim = c(-1, 1), col = "#999999",
         xlab = "Species correlations", ylab = "frequency")
abline(v = 0, lty = "dashed", lwd = 3, col = "#333333")

##==============================================================
## Correlations *within groups* (Predators, Grazers, Cyanobacteria)?
## data.frame with:  column of correlation values, column of Group (including "Between Groups")
Spp.cort  <- cor(Trophic.sp, method = "spearman")
Spp.cort[upper.tri(Spp.cort, diag = TRUE)] <- NA
Spp.cordf <- reshape(data = as.data.frame(Spp.cort), 
                     varying = list(colnames(Spp.cort)),
                     times = colnames(Spp.cort),
                     timevar = "Species2",
                     idvar = "Species1",
                     ids = rownames(Spp.cort),
                     v.names = "cor",
                     direction = "long"
                     )
Spp.cordf <- subset(Spp.cordf, !is.na(cor) )
Spp.cordf <- within(Spp.cordf[, c("cor", "Species1", "Species2")], 
                    {
                      Group1 <- match(Species1, SECC.fauna.meta$ID)
                      Group1 <- SECC.fauna.meta$Trophic.Group[Group1]
                      Group1 <- ifelse(Species1 %in% c("Stigonema", "Nostoc"), "Cyanobacteria", Group1)
                      Group2 <- match(Species2, SECC.fauna.meta$ID)
                      Group2 <- SECC.fauna.meta$Trophic.Group[Group2]
                      Group2 <- ifelse(Species2 %in% c("Stigonema", "Nostoc"), "Cyanobacteria", Group2)
                      ## 1 column of between Group labels
                      Group <- ifelse(Group1 == Group2, Group1, 
                                      paste(substr(Group1, 1, 4), substr(Group2, 1, 4), sep="-") 
                      )
                      Group <- gsub("Grazer"    , "Grazers"   , Group , fixed = TRUE)
                      Group <- gsub("Predator"  , "Predators" , Group , fixed = TRUE)
                      Group <- gsub("Graz-Pred" , "Pred-Graz" , Group , fixed = TRUE)
                      Group <- gsub("Cyan-Graz" , "Graz-Cyan" , Group , fixed = TRUE)
                      Group <- gsub("Cyan-Pred" , "Pred-Cyan" , Group , fixed = TRUE)
                    })

## Predators are essentially uncorrelated (~0), Grazers & Between Groups are slightly +vely correlated (~0.05)
Graz.cortest <- t.test( subset(Spp.cordf, Group == "Grazers", select = "cor") )
Pred.cortest <- t.test( subset(Spp.cordf, Group == "Predators", select = "cor") )
## Cyan.cortest <- t.test( subset(Spp.cordf, Group == "Cyanobacteria", select = "cor") )  # only 1 correlation (>0)!
Btwn.cortest <- t.test( subset(Spp.cordf, ! Group %in% c("Grazers", "Predators", "Cyanobacteria"), select = "cor") )
PrGr.cortest <- t.test( subset(Spp.cordf, Group == "Pred-Graz", select = "cor") )
GrCy.cortest <- t.test( subset(Spp.cordf, Group == "Graz-Cyan", select = "cor") )
PrCy.cortest <- t.test( subset(Spp.cordf, Group == "Pred-Cyan", select = "cor") )
Spp.cortest
Graz.cortest
Pred.cortest
Btwn.cortest

## Stacked bar graph :D
Spp.corplot <- ggplot(subset(Spp.cordf, Group %in% c("Predators", "Grazers", "Pred-Graz")), 
                      aes(x = cor, fill = Group)) + 
    xlim(-1, 1) + xlab("Species correlation") + ylab("frequency") +
    geom_bar(binwidth = 0.1, colour = "#333333", lwd = 0.2) + 
    scale_fill_manual(values = c("#CCCCCC", "#999999", "#444444"), 
                      breaks = c("Predators", "Grazers", "Pred-Graz")) +
    geom_vline(xintercept = 0, colour = "#666666", lwd = 1, lty = "dashed") +
    jaw.ggplot()
print(Spp.corplot)

pairplot(SECC.sum[, c("Predators", "Grazers", "Cells.g", "H2O")])
## scatterplot of group abundances in Final Publication Graphs below

## Correlations between groups, within experimental treatment combinations
library(plyr)
cordf <- function(data, ...)
{
  ## aggregate over species groups: Predators, Grazers, cyanobacteria
  ## pass in SECC.sum: already aggregated.  Pull out relevant columns
  cols.numeric <- sapply(data, is.numeric)
  ## calculate pairwise correlations
  data.cor <- cor(data[, cols.numeric], ...)

  ## convert result to a data.frame for output
  colnames.cor <- paste(rep( dimnames(data.cor)[[1]], each  = length(dimnames(data.cor)[[2]]) ), 
                        rep( dimnames(data.cor)[[2]], times = length(dimnames(data.cor)[[1]]) ), sep=".")
  colnames.cor <- colnames.cor[lower.tri(data.cor)]
  data.cor <- data.cor[lower.tri(data.cor, diag = FALSE)]
  names(data.cor) <- colnames.cor

  data.cor
}

Trophic.cor <- ddply(.data = SECC.sum[, c("Chamber", "Frag", "Position", "Predators", "Grazers", "Cells.g")], 
                     .variables = c("Chamber", "Position"), .fun = cordf , method = "spearman")
Trophic.cor <- Trophic.cor[, c(1, 2, 3, 5, 4)] # rearrange columns: the number of columns MATTERS!
Trophic.cor$Chamber  <- factor(Trophic.cor$Chamber,  labels = c("Ambient", "Chamber"))
Trophic.cor$Position <- factor(Trophic.cor$Position)

## Correlations among **individual species**
corag <- function(data, ...)
{
  ## aggregate over species groups: Predators, Grazers, cyanobacteria
  ## pass in SECC.sum: already aggregated.  Pull out relevant columns
  cols.numeric <- sapply(data, is.numeric)
  ## calculate pairwise correlations
  data.cor <- cor(data[, cols.numeric], ...)

  ## convert result to a data.frame for output
  data.cor[upper.tri(data.cor, diag = TRUE)] <- NA
  data.cordf <- reshape(data = as.data.frame(data.cor), 
                       varying = list(colnames(data.cor)),
                       times = colnames(data.cor),
                       timevar = "Species2",
                       idvar = "Species1",
                       ids = rownames(data.cor),
                       v.names = "cor",
                       direction = "long"
                       )
  data.cordf <- subset(data.cordf, !is.na(cor) )
  data.cordf <- within(data.cordf[, c("cor", "Species1", "Species2")], 
                      {
                        Group1 <- match(Species1, SECC.fauna.meta$ID)
                        Group1 <- SECC.fauna.meta$Trophic.Group[Group1]
                        Group1 <- ifelse(Species1 %in% c("Stigonema", "Nostoc"), "Cyanobacteria", Group1)
                        Group2 <- match(Species2, SECC.fauna.meta$ID)
                        Group2 <- SECC.fauna.meta$Trophic.Group[Group2]
                        Group2 <- ifelse(Species2 %in% c("Stigonema", "Nostoc"), "Cyanobacteria", Group2)
                        ## 1 column of between Group labels
                        Group <- ifelse(Group1 == Group2, Group1, 
                                        paste(substr(Group1, 1, 4), substr(Group2, 1, 4), sep="-") 
                        )
                        Group <- gsub("Grazer"    , "Grazers"   , Group , fixed = TRUE)
                        Group <- gsub("Predator"  , "Predators" , Group , fixed = TRUE)
                        Group <- gsub("Graz-Pred" , "Pred-Graz" , Group , fixed = TRUE)
                        Group <- gsub("Cyan-Graz" , "Graz-Cyan" , Group , fixed = TRUE)
                        Group <- gsub("Cyan-Pred" , "Pred-Cyan" , Group , fixed = TRUE)
                      })

  cor.mean <- aggregate(data.cordf$cor, by = list(Group = data.cordf$Group), FUN = mean)
  cor.sd   <- aggregate(data.cordf$cor, by = list(Group = data.cordf$Group), FUN = sd  , simplify = TRUE)
  cordf <- cbind(cor.mean, cor.sd$x)
  colnames(cordf) <- c("Group", "cor.mean", "cor.sd")

  cordf
}

Trophicor <- ddply(.data = cbind(SECC.env[, c("Chamber", "Frag", "Position")], Trophic.sp), 
                     .variables = c("Chamber", "Position"), .fun = corag , method = "spearman")
Trophicor$Chamber  <- factor(Trophicor$Chamber,  labels = c("Ambient", "Chamber"))
Trophicor$Position <- factor(Trophicor$Position)
Trophicor$Group <- factor(Trophicor$Group, levels = c("Predators", "Grazers", "Cyanobacteria", "Pred-Graz", "Graz-Cyan", "Pred-Cyan"))
Trophicor$id     <- factor(paste(Trophicor$Chamber, Trophicor$Position, sep="_"))

## On average: 
## * Predators & Grazers are +vely correlated,
## * Grazers & Cells are weakly +vely correlated,
## * and Predators & Cells are weakly negatively / un- correlated.


##================================================
## Change in Abundance
##================================================
cat("Computing changes in fauna species densities\n")

Fauna.dab <- Fauna.sp[Fauna.sp$Frag == 4 & Fauna.sp$Time == 4, ]
Fauna.dab$Trt <- paste( Fauna.dab$Chamber, "_", Fauna.dab$Frag, ".", Fauna.dab$Pos, sep="" )
cols.sp <- sapply(Fauna.dab, is.numeric) & ( colnames(Fauna.dab) %in% colnames(Fauna) )
cols.sp <- names(cols.sp)[which(cols.sp)]
Fauna.dab <- aggregate(Fauna.dab[, cols.sp], 
                       by=list(Trt     = Fauna.dab$Trt), 
                       FUN=mean
                       )
row.names(Fauna.dab) <- Fauna.dab$Trt
Fauna.dab <- Fauna.dab[, -1]
Fauna.dab <- as.data.frame( t(Fauna.dab) )

axl <- c(0, 4)
plot(Fauna.dab$A_4.O, Fauna.dab$A_4.I, xlim = axl, ylim = axl )
abline(0, 1, lwd = 2, lty = 3)
plot(Fauna.dab$C_4.I, Fauna.dab$A_4.I, xlim = axl, ylim = axl )
abline(0, 1, lwd = 2, lty = 3)
plot(Fauna.dab$C_4.I, Fauna.dab$C_4.O, xlim = axl, ylim = axl )
abline(0, 1, lwd = 2, lty = 3)

with(SECC.sp.sum, plot(Richness, Predators/(Predators + Grazers), pch = 21, ylim = c(0,1) ))
with(SECC.sp.sum, plot(Richness, Mesostigmata/fauna, pch = 21, ylim = c(0,1) ))
with(SECC.sp.sum, plot(Richness, Collembola/fauna, pch = 21, ylim = c(0,1) ))

## Proportion (abundance) by Richness
library(mgcv)
PredProp.plot <- ggplot(SECC.sp.sum, aes(x = Richness, y = Predators/fauna)) + geom_point() + stat_smooth(method = "gam") + ylim(c(0,1))
print(PredProp.plot)
MesoProp.plot <- ggplot(SECC.sp.sum, aes(x = Richness, y = Mesostigmata/fauna)) + geom_point() + stat_smooth(method = "gam") + ylim(c(0,1))
print(MesoProp.plot)
CollProp.plot <- ggplot(SECC.sp.sum, aes(x = Richness, y = Collembola/fauna)) + geom_point() + stat_smooth(method = "gam") + ylim(c(0,1))
print(CollProp.plot)

## Proportion (# species) by Richness
SECC.sp.sum <- within(SECC.sp.sum, {
    PredRich <- apply(SECC.sp[, which(SECC.fauna.meta[match(colnames(SECC.sp), SECC.fauna.meta$ID), 
                                          "Trophic.Group"] == "Predator")], 
                          1, function(x) length(which(x>0)) 
                          )  # observed # spp.
    MesostigRich <- apply(SECC.sp[, which(SECC.fauna.meta[match(colnames(SECC.sp), SECC.fauna.meta$ID), 
                                          "Taxonomic.Group"] == "Mesostigmata")], 
                          1, function(x) length(which(x>0)) 
                          )  # observed # spp.
    CollemboRich <- apply(SECC.sp[, which(SECC.fauna.meta[match(colnames(SECC.sp), SECC.fauna.meta$ID), 
                                          "Taxonomic.Group"] == "Collembola")], 
                          1, function(x) length(which(x>0)) 
                          )  # observed # spp.
})

PredRich.plot <- ggplot(SECC.sp.sum, aes(x = Richness, y = PredRich/Richness)) + geom_point() + stat_smooth(method = "gam") + ylim(c(0,1))
print(PredRich.plot)
MesoRich.plot <- ggplot(SECC.sp.sum, aes(x = Richness, y = MesostigRich/Richness)) + geom_point() + stat_smooth(method = "gam") + ylim(c(0,1))
print(MesoRich.plot)
CollRich.plot <- ggplot(SECC.sp.sum, aes(x = Richness, y = CollemboRich/Richness)) + geom_point() + stat_smooth(method = "gam") + ylim(c(0,1))
print(CollRich.plot)




##================================================
## Treatment Presence/Absence ; Species Occurences
##================================================
cat("Computing fauna species occurences in treatments (presence/absence)\n")

Fauna.spr   <- recodeSECC(Fauna.sp)
Fauna.trtpa <- aggregate(Fauna.spr[, cols.sp], 
                         by=list(Time     = Fauna.spr$Time,
                                 Chamber  = Fauna.spr$Chamber,
                                 Frag     = Fauna.spr$Frag,
                                 Position = Fauna.spr$Position
                                 ),
                         FUN = sum
                         )
Fauna.trtpa[, cols.sp] <- decostand(Fauna.trtpa[, cols.sp], method = "pa")

xtable(t(Fauna.trtpa))



##================================================
## Composition
##================================================

SECC.sp.CxP <- checkSECCdata(SECC.sp.sum, 'SECC.sp.sum')
SECC.sp.CxP <- recodeSECC( SECC.sp.CxP )
SECC.sp.CxP <- reshape(data = SECC.sp.CxP, 
                       idvar = colnames(SECC.sp.sum)[!cols.numeric], 
                       varying = list(c("Predators", "Grazers")),
                       times = c("Predators", "Grazers"),
                       v.names = "Abundance", timevar = "Group",
                       direction = "long"
                       )
SECC.sp.CxP <- aggregate(Abundance ~ Chamber + Position + Group, 
                         data = SECC.sp.CxP, 
                         FUN = mean)

Compo.plot <- ggplot(SECC.sp.CxP, aes(x = Position, y = Abundance, fill = Group)) + 
    geom_bar(stat="identity", position = "fill", colour = NA) +
    scale_fill_manual(values = c("Predators" = "#444444", "Grazers" = "#AAAAAA"),
                      breaks = c("Predators", "Grazers")) +
    facet_wrap(~ Chamber) + 
    jaw.ggplot()

print(Compo.plot)




################################################################
### SPATIAL ANALYSIS
################################################################
cat("Spatial analysis of fauna data ...\n")
## Test for spatial trend
Fauna.rda <- rda(Fauna.hell, Fauna.xy[, c("xE", "yN")]) 
anova(Fauna.rda)
## Yes; not sure what this means, however
ordiplot(Fauna.rda, type = "t")

## Mantel test: Are similarities spatially correlated?
Fauna.xydist <- vegdist(Fauna.xy[, c("xE", "yN")], "euclidean", upper = TRUE)
Mantel.bray  <- mantel(Fauna.bray, Fauna.xydist, permutations = 1000)
Mantel.jacc  <- mantel(Fauna.jacc, Fauna.xydist, permutations = 1000)
Mantel.hell  <- mantel(Fauna.hellD, Fauna.xydist, permutations = 1000)
Mantel.chord <- mantel(Fauna.chorD, Fauna.xydist, permutations = 1000)
Mantel.euc   <- mantel(Fauna.euclD, Fauna.xydist, permutations = 1000)
## No? Yes, for Chord & Euclidean; No for everything else.

## correlgoram: At what scale are similarities correlated?
Fauna.det <- resid(lm(as.matrix(Fauna.chorD) ~ ., data = Fauna.xy[, c("xE", "yN")]))
Fauna.ddist <- vegdist(Fauna.det, method = "euclidean")
( Fauna.correlog <- mantel.correlog(Fauna.ddist, XY = Fauna.xy[, c("xE", "yN")], nperm = 99) )
plot(Fauna.correlog, main = "Mantel correlogram (Chord distance matrix)")
## strongly negative correlation around distance class 2 (index 13).  Otherwise, nothing.
## This corresponds to 9--18 m (Fauna.correlog$break.pts)

## Partial Mantel Test: differences among treatment groups, after removing effect of space?
## Not much point, if there seems to be very little effect of space.
Mantelp.design <- sapply( SECC.env[, c("Chamber", "Frag", "Pos")], as.numeric )
Mantelp.design <- sapply( as.data.frame(Mantelp.design), factor )
Mantelp.design <- sapply( as.data.frame(Mantelp.design), as.numeric ) # use smallest numbers, rather than leftover factor levels from reduced set.
Mantelp.dist   <- vegdist(Mantelp.design, method = "manhattan")

## Partial out TREATMENTS: is there an effect of space?
Mantelp.bray  <- mantel.partial(Fauna.bray, Fauna.xydist, Mantelp.dist, permutations = 1000)
Mantelp.chord <- mantel.partial(Fauna.chorD, Fauna.xydist, Mantelp.dist, permutations = 1000)
## Partial out SPACE: are there treatment differences? This isn't really the best way to answer this question (rather than a mixed-effects model that includes spatial autocorrelation)
Mantels.bray  <- mantel.partial(Fauna.bray, Mantelp.dist, Fauna.xydist, permutations = 1000)
Mantels.chord <- mantel.partial(Fauna.chorD, Mantelp.dist, Fauna.xydist, permutations = 1000)




################################################################
### CONSTRAINED ORDINATION: FITTING ENVIRONMENTAL VARIABLES
################################################################
cat("Constrained ordinations of fauna data ...\n")
## Extract matching ordination rows from Environmental data, *in same order*
Fauna.sort <- match( rownames(Fauna.trans), SECC$SampleID )
Fauna.sort <- na.omit(Fauna.sort)
Fauna.env  <- SECC[Fauna.sort, ]
## Process standardized environmental variables - should be able to substitue with SECC.env, which already has this
Fauna.env <- within(Fauna.env,
                    {
                      Stigonema.g <- Stigonema / Cells.dwt
                      Nostoc.g <- Nostoc / Cells.dwt
                      Other.cells.g <- Other.cells / Cells.dwt
                      Prod13 <- Prod12 + Prod23
                      Trt <- paste(Chamber, Frag, Position, sep="\n")
                    })



## Post Hoc fit of environmental variables to nMDS (not really constrained ordination)
## Samples               : Same experimental subset, or ALL available?
## Response Variables    : Fauna species [SECC.fauna]
## Explanatory Variables : H2O, Temperature?, Patch.dwt ; Cells.g, Stigonema.g, Nostoc.g, ARA.g ; Decomposition, Productivity (all time periods), NH4, NO3, TAN
##    Experimental Factors (if subset of samples)?  <- partial these out?

Fauna.fit <- envfit(Fauna.mds ~ Block + Chamber + Frag + Position + 
                    H2O + Patch.dwt + Cells.g + Stigonema.g + Nostoc.g + ARA.g + 
                    Decomposition + Prod01 + Prod12 + Prod23 + 
                    NH4 + NO3 + TAN
                    , data = Fauna.env, perm = 1000, na.rm = TRUE)
Fauna.fit

## Add environmental variables to ordination
ordiplot( Fauna.mds, type="points", display = "sites")
ordispider(Fauna.mds, Fauna.env$Trt, label = TRUE, col = "red")
plot(Fauna.fit, p.max = 0.025, adj = 0.5, srt = 0, cex = 1)
ordiplot( Fauna.mds, type="text", display = "species", xlim = c(-1, 2), ylim = c(-1, 2) )
plot(Fauna.fit, p.max = 0.05, adj = 0.5, srt = 0, cex = 1)



##==============================================================
## CCA
##==============================================================
## Missing env values (productivity) 
Env.rows     <- apply(Fauna.env, 1, function (x) !any(is.na(x)) )
Fauna.envs   <- Fauna.env[which(Env.rows), ]
Fauna.subset <- Fauna.trans[which(Env.rows), ]

Fauna.cca <- cca(Fauna.subset ~ Block + Chamber + Frag + Position + 
                 H2O + Patch.dwt + ARA.g + Cells.g + # Stigonema.g + Nostoc.g + 
                 Decomposition + Prod01 + Prod13 + # Prod12 + Prod23 + 
                 TAN # + NH4 + NO3
                 , data = Fauna.envs)
Fauna.cca
summary(Fauna.cca)
## First 2 axes don't seem to capture much
## the nMDS seems to be a better overall fit than the Gaussian model of CCA?
## Permutation-based tests: SLOW
CCA.anova.terms <- anova(Fauna.cca, by = "term",   step = 200) # Permutation-based tests
CCA.anova.mar   <- anova(Fauna.cca, by = "margin", step = 99 ) # *marginal* effects of each term (in a model with all other terms)
## Plots
ordiplot(Fauna.cca, type="text")
ordiplot(Fauna.cca, type="text", 
         display = c("species", "bp"), scaling = 2, 
         xlim = c(-2, 2), ylim = c(-2, 2)) # zoom in, but still messy.
refCol <- "#AAAAAA" 
refLty <- 3
rect(-0.5, -0.5, 0.5, 0.5, col = NA, border = refCol, lty = refLty)

## without experimental factors
Fauna.cca1 <- cca(Fauna.subset ~  
                 H2O + Patch.dwt + ARA.g + Cells.g + # Stigonema.g + Nostoc.g + 
                 Decomposition + Prod01 + Prod13 + # Prod12 + Prod23 + 
                 TAN # + NH4 + NO3
                 , data = Fauna.envs)
Fauna.cca1
summary(Fauna.cca1)
ordiplot(Fauna.cca1, type="text")
ordiplot(Fauna.cca1, type="text", 
         display = c("species", "bp"), scaling = 2, 
         xlim = c(-2, 2), ylim = c(-2, 2)) # zoom in, but still messy.





##################################################
### PUBLICATION GRAPHS
##################################################
cat("Producing Figures for Publication\n")

## plot symbol map
## Inner = full
## Outer = open
## Chamber = Red
## Ambient = Black
## Pos: pch (shape); Chamber: colour; Frag: shape?
Position.map <- plotMap( factor = "Position" )
Position.map <- Position.map[ Position.map$label %in% levels(SECC.fauna$Pos), ]
Chamber.map  <- plotMap( factor = "Chamber"  )
Chamber.map  <- Chamber.map[ Chamber.map$label %in% levels(SECC.fauna$Chamber), ]

Position.label <- attr(SECC, "labels")[["Pos"]] # "Patch\nPosition"
Chamber.label  <- "Chamber"         # attr(SECC, "labels")[["Chamber"]]

## Prepare data to plot
## has vegan changed?  stress values now seem to be <1, rather than in % ?!?
## stress.label <- paste( "Stress\n= ", round( Fauna.mds$stress, 1)/100 )	# Make text label of Stress value, rounded to 1 decimal place, and converted from % to a decimal number.
stress.label <- paste( "Stress\n= ", round( Fauna.mds$stress, 2) )	# Make text label of Stress value.
MDS.pts <- as.data.frame(Fauna.mds$points)
colnames(MDS.pts) <- c('x', 'y')
MDS.pts <- cbind(MDS.pts, Time = Fauna.sp$Time, Chamber = Fauna.sp$Chamber, 
                 Frag = Fauna.sp$Frag, Pos = Fauna.sp$Pos
)

## Data for spider plot
MDS.spider <- MDS.pts[which(MDS.pts$Chamber == "C" & MDS.pts$Pos == "O"), ]
MDS.centroids <- MDS.spider[0, ]
for(f in levels(MDS.spider$Frag))
{
  g.pts <- MDS.spider[MDS.spider$Frag == f,]
  centrx <- mean(g.pts$x)
  centry <- mean(g.pts$y)
  MDS.centroids <- rbind(MDS.centroids, 
                         data.frame(x = centrx, y = centry,
                                    Time = 4, Chamber = "C",
                                    Frag = f, Pos = "O")
  )
}
c.rows <- match(MDS.spider$Frag, MDS.centroids$Frag)
MDS.spider$xend <- MDS.centroids$x[c.rows]
MDS.spider$yend <- MDS.centroids$y[c.rows]

## Build plot in ggplot
## including grouping variables (shape, colour, etc.) in the default aes() (or qplot) causes geom_text() to look REALLY UGLY
Fauna.plot <- ggplot(data = MDS.pts, aes(x = x, y = y) ) +  
                 ## Add spider plots of outer chamber patches (under the points)
                 geom_segment(aes(x = x, y = y, xend = xend, yend = yend, 
                                  colour = Chamber), data = MDS.spider) +  # , alpha = 0.6
                 ## I want to overlay on the spider plot (behind the points):
                 ## - a hexagon with a red outline ('outer chamber')
                 ## - an icon of the frag treatment (a la grimport)
                 ## complicated examples:
                 ## http://stackoverflow.com/questions/2181902/how-to-use-an-image-as-a-point-in-ggplot
                 ## http://stackoverflow.com/questions/8905101/how-can-i-use-a-graphic-imported-with-grimport-as-axis-tick-labels-in-ggplot2-u
                 ## I don't have time to figure this out now, so I might just add those manually in Illustrator
                 ## Add colour-coded points on top of spider plot :D
                 geom_point(aes(group = Chamber, shape = Pos, colour = Chamber), 
                           size = 3, lwd = 1.5) +
                 xlab(NULL) + ylab(NULL)
Fauna.plot <- Fauna.plot + scale_colour_manual(name = Chamber.label,
                                               values = Chamber.map$col, 
                                               breaks = Chamber.map$label,
                                               labels = c("Ambient", "Chamber")
                                               )
Fauna.plot <- Fauna.plot + scale_shape_manual(name = Position.label,
                                              values = Position.map$pch, 
                                              breaks = Position.map$label,
                                              labels = c("Inner (Wet)", "Outer (Dry)")
                                              )
Fauna.plot <- Fauna.plot + 
                scale_x_continuous(breaks = -4:4, limits = c(-3.7, 1.3)) +
                scale_y_continuous(breaks = -4:4, limits = c(-2.7, 2.3))
## options: see theme_get() for available theme options.
Fauna.plot <- Fauna.plot + jaw.ggplot() + coord_equal()
Fauna.plot <- Fauna.plot + opts(axis.ticks = theme_blank(), 
                                axis.text.x = theme_blank(), 
                                axis.text.y = theme_blank(), 
                                legend.key = theme_rect(colour = NA),
                                legend.position = "right",
                                legend.direction = "vertical"
                                )
Fauna.plot <- Fauna.plot + geom_text(aes(label = stress.label, x = -3.5, y = max(y) ), 
                                     size = 4, colour = "black", hjust = 0, vjust = 1)
print(Fauna.plot)

## this is where I really want the stress label, but ggplot doesn't seem to be able to do this easily :(
## grid.text(label = stress.label, x = 0.85, y = 0.85, just = "left", hjust = 0, vjust = 1)
## mtext( stress.label , side=3, line=0, adj=1 ) # kludge: trying to copy or save the graph now will make R crash.

Fauna.plot.frag <- Fauna.plot + facet_grid(facets = Frag~.)
print(Fauna.plot.frag)



##==============================================================
## Species correlations
## Scatter plot; fitted line should be a Model II regression (Major Axis)
library(lmodel2)
PredGraz.lm2 <- lmodel2(Predators ~ Grazers, data = SECC.sp.sum, nperm = 999)
if (FALSE)
{
  ## Model 2 regression: axes / order of variables does not matter!
  GrazPred.lm2 <- lmodel2(Grazers ~ Predators, data = SECC.sp.sum, nperm = 999)
  plot(PredGraz.lm2)
  plot(GrazPred.lm2)
}
PredGraz.abline <- PredGraz.lm2$regression.results
lm2.row <- which( PredGraz.abline == "MA" ) 
PredGraz.abline <- PredGraz.abline[lm2.row, ]
PredGraz.abline <- within(PredGraz.abline,
                          {
                            b.lower <- PredGraz.lm2$confidence.intervals[lm2.row, 4]
                            a.lower <- mean(PredGraz.lm2$y) - b.lower * mean(PredGraz.lm2$x)
                            b.upper <- PredGraz.lm2$confidence.intervals[lm2.row, 5]
                            a.upper <- mean(PredGraz.lm2$y) - b.upper * mean(PredGraz.lm2$x)
                          })
PredGraz.ribbon <- data.frame(x = c(0, mean(PredGraz.lm2$x), 40),
                              y = c(PredGraz.abline[1, 2],
                                    PredGraz.abline[1, 2] +
                                    mean(PredGraz.lm2$x) * PredGraz.abline[1, 3],
                                    PredGraz.abline[1, 2] +
                                    40 * PredGraz.abline[1, 3]
                                    ),
                              lower = c(mean(PredGraz.lm2$y) - 
                                        PredGraz.abline[1, "b.upper"] * mean(PredGraz.lm2$x),
                                        mean(PredGraz.lm2$y),
                                        mean(PredGraz.lm2$y) + 
                                        PredGraz.abline[1, "b.lower"] * 
                                        (40 - mean(PredGraz.lm2$x))
                                        ),
                              upper = c(mean(PredGraz.lm2$y) - 
                                        PredGraz.abline[1, "b.lower"] * mean(PredGraz.lm2$x),
                                        mean(PredGraz.lm2$y),
                                        mean(PredGraz.lm2$y) + 
                                        PredGraz.abline[1, "b.upper"] * 
                                        (40 - mean(PredGraz.lm2$x))
                                        )
                              )

PredGraz.cor    <- cor(SECC.sp.sum$Predators, SECC.sp.sum$Grazers)
PredGraz.cortxt <- sprintf("r = %.2f", PredGraz.cor)
Chamber.label  <- "Chamber"         # attr(SECC, "labels")[["Chamber"]]

PredGraz.corplot <- ggplot(SECC.sp.sum, 
                           ## aes(x = Grazers, y = Predators)
                           ) + 
                geom_ribbon(aes(x = x, ymin = lower, ymax = upper), data = PredGraz.ribbon,
                            fill = "#999999", alpha = 0.5
                ) +
                geom_line(aes(x = x, y = y),
                          data = PredGraz.ribbon,
                          colour = "#666666", size = 0.4
                          ) +
                geom_point(aes(x = Grazers, y = Predators,
                               shape = Pos, colour = Chamber), 
                            size = 2.5) + 
                scale_colour_manual(name = Chamber.label,
                                    values = Chamber.map$col, 
                                    breaks = Chamber.map$label,
                                    labels = c("Ambient", "Chamber")
                                    ) +
                scale_shape_manual(name = Position.label,
                                   values = Position.map$pch, 
                                   breaks = Position.map$label,
                                   labels = c("Inner (Wet)", "Outer (Dry)")
                                   ) +
                xlab(bquote("Grazers (" * .( attr(SECC.sp.sum, "units")[["Grazers"]]) * ")")) +
                ylab(bquote("Predators (" * .( attr(SECC.sp.sum, "units")[["Predators"]]) * ")")) +
                geom_text(aes(x = max(Grazers), 
                              y = max(Predators),
                              label = PredGraz.cortxt
                              ), data = SECC.sp.sum, hjust = 1, vjust = 1) +
                jaw.ggplot() + coord_equal()
## WHY isn't the legend drawing properly??!? might be due to a different version of Chamber.map or Position.map previously used from the ANOVA scripts.  Seems to work fine here...
print(PredGraz.corplot)

if (FALSE)
{
                  geom_abline(intercept = PredGraz.abline[1, 2],
                              slope = PredGraz.abline[1, 3],
                              colour = "#666666", size = 0.4
                              ) +
                                 geom_abline(intercept = PredGraz.abline[1, 6],
                                             slope = PredGraz.abline[1, 7],
                                             colour = "#999999", size = 0.4, lty = "dashed"
                                             ) +
                                 geom_abline(intercept = PredGraz.abline[1, 8],
                                             slope = PredGraz.abline[1, 9],
                                             colour = "#999999", size = 0.4, lty = "dashed"
                                             ) +
                xlim(c(0, 40)) 
}



##==============================================================
## Correlations between hypothesized "trophic groups"
## Uses Chamber.map & Position.map from above
if (FALSE)
{   ## reorganize data for plotting
Trophicor <- reshape(data = Trophic.cor, 
                     varying = list(colnames(Trophic.cor)[3:5]),
                     times = colnames(Trophic.cor)[3:5],
                     timevar = "Groups",
                     ##                      idvar = "Trt",
                     ##                      ids = paste(Trophic.cor$Chamber, Trophic.cor$Position, sep="_"),
                     v.names = "cor",
                     direction = "long"
                     )
}
## Modify labels (now done in ggplot code directly)
## Trophicor$Groups <- gsub("Cells\\.g", "cyanobacteria", Trophicor$Groups)
## Trophicor$Groups <- gsub("\\.", " %prop% \n", Trophicor$Groups)

Trophic.corplot  <- ggplot(data = droplevels( subset(Trophicor, Group %in% c('Pred-Graz', 'Graz-Cyan', 'Pred-Cyan')) ), 
                           aes(x = Group, y = cor.mean, group = id)) +
                    geom_line(aes(colour = Chamber, size = Position, lty = Chamber)) +
                    ##                     geom_errorbar(aes(ymax = cor.mean + cor.sd, ymin = cor.mean - cor.sd, 
                    ##                                       colour = Chamber, size = Position, lty = Chamber), alpha = 0.3) +
                    geom_point(aes(shape = Position, colour = Chamber, fill = Chamber), size = 3) +
                    geom_hline(aes(x = 0)) + ylim(c(-1, 1)) +
                    xlab("Trophic level comparisons") + ylab("mean spearman correlation coefficient") +
                    scale_x_discrete(labels = c('Pred-Graz' = "Predators &\nGrazers",
                                                'Graz-Cyan' = "Grazers &\nCyanobacteria",
                                                'Pred-Cyan' = "Predators &\nCyanobacteria"
                                                )
                    ) + 
                    scale_shape_manual(name = Position.label,
                                       values = Position.map$pch) +
                    scale_colour_manual(name = Chamber.label,
                                        values = Chamber.map$col) +
                    scale_fill_manual(name = Chamber.label,
                                      values = Chamber.map$bg) +
                    scale_linetype_manual(name = Chamber.label,
                                          values = c(3, 1)) +
                    scale_size_manual(name = Position.label,
                                      values = Position.map$lwd*0.5) +
                    jaw.ggplot()

print(Trophic.corplot)

Within.corplot  <- ggplot(data = droplevels( subset(Trophicor, Group %in% c('Predators', 'Grazers', 'Cyanobacteria')) ), 
                           aes(x = Group, y = cor.mean, group = id)) +
                    geom_line(aes(colour = Chamber, size = Position, lty = Chamber)) +
                    ##                     geom_errorbar(aes(ymax = cor.mean + cor.sd, ymin = cor.mean - cor.sd, 
                    ##                                       colour = Chamber, size = Position, lty = Chamber), alpha = 0.3) +
                    geom_point(aes(shape = Position, colour = Chamber, fill = Chamber), size = 3) +
                    geom_hline(aes(x = 0)) + ylim(c(-1, 1)) +
                    xlab("Trophic level comparisons") + ylab("mean spearman correlation coefficient") +
                    scale_shape_manual(name = Position.label,
                                       values = Position.map$pch) +
                    scale_colour_manual(name = Chamber.label,
                                        values = Chamber.map$col) +
                    scale_fill_manual(name = Chamber.label,
                                      values = Chamber.map$bg) +
                    scale_linetype_manual(name = Chamber.label,
                                          values = c(3, 1)) +
                    scale_size_manual(name = Position.label,
                                      values = Position.map$lwd*0.5) +
                    jaw.ggplot()

print(Within.corplot)


##==============================================================
## Change in Abundance between Treatments (Gonzalez-style)
## prepare data for plotting
Fauna.dab <- within(Fauna.dab, 
                    {
                      Trophic <- Fauna.meta$Trophic[ match(rownames(Fauna.dab), rownames(Fauna.meta)) ]
                    })

axl <- c(0, 4)
Faunachange.plot  <- ggplot(data = Fauna.dab, aes(x = Fauna.dab$C_4.I, y = Fauna.dab$C_4.O)) +
                    geom_point(aes(shape = Trophic), colour = "black", size = 3) +
                    geom_abline(aes(intercept = 0, slope = 1), size = 0.5, lty = "dashed", colour = "#666666") +
                    xlab(expression(atop( "Species density in " , italic("Isolated ") * bolditalic("Inner ") * italic("Chamber") * " patches (#" %.% g^-1 * ")" ))) + 
                    ylab(expression(atop( "Species density in " , italic("Isolated ") * bolditalic("Outer ") * italic("Chamber") * " patches (#" %.% g^-1 * ")" ))) +
                    scale_shape_manual(name = "Trophic group",
                                       values = c(19, 21)) +
                    xlim(axl) + ylim(axl) +
                    coord_equal() + jaw.ggplot() + opts( plot.margin = unit(c(0, 0, 0, 1), "lines") )

print(Faunachange.plot)





##================================================
## Save (Multivariate) Results
##================================================
cat("Saving Results...\n")
if (TRUE) 
{
  library(xtable)
  ANOVA.table <- xtable(Fauna.ANOVA, digits = c(0, 0, 0, 2, 3, 2, 3, 2, 3, 2, 3))
  ANOSIM.table <- xtable(ANOSIM.results[, c("Groups", "bray_.p", 
                                                   "bray_log.p", "jaccard_pa.p")], 
                                digits = 3) # bray, bray_log, jaccard **** 
  capture.output(  cat("FAUNA results\n", date(), "\n")
      , cat("\n\n", "======== ANOVA results summary ========", "\n", sep="")
      , print(Fauna.ANOVA), cat("\n\n")
      , print(ANOVA.table), cat("\n\n")
      , cat("\n\n", "======== Species Densities by Trt (for plotting changes in densities) ========", "\n", sep="")
      , print(Fauna.dab),   cat("\n\n")
      , cat("\n\n", "======== Species Occurence by Trt (presence/absence) ========", "\n", sep="")
      , print(Fauna.trtpa), cat("\n\n")
      , print(xtable(t(Fauna.trtpa))), cat("\n\n")
      , cat("\n\n", "======== Trophic Group Correlations by Trt ========", "\n", sep="")
      , print(Trophic.cor), cat("\n\n")
      , cat("\n\n", "======== ANOSIM results ========", "\n", sep="")
      , print(ANOSIM.results), cat("\n\n")
      , print(ANOSIM.table), cat("\n\n")
      , cat("\n\n", "======== Environmental Variables Post Hoc fit to nMDS ========", "\n", 
            sep="")
      , print(Fauna.fit)
      , cat("\n\n", "======== CCA ========", "\n", sep="")
      , summary(Fauna.cca)
      , cat("\n\n", "======== Permutation tests of CCA terms ========", "\n", sep="")
      , print(CCA.anova.terms)
      , file = "./output/Fauna - Main Results.txt"
  )
  ggsave(filename="./graphs/Figure - Fauna-MDS.eps", plot = Fauna.plot, width = 6, height = 6)
  ggsave(filename="./graphs/Figure - Fauna-Correlation.pdf", plot = PredGraz.corplot, width = 6, height = 2.8) # short height to push x-axis label a little closer to the axis (coord_equal bug?)
  ## semi-transparency (`alpha = `) is not supported by eps.  But it is by pdf ;)
  ggsave(filename="./graphs/Figure - Trophic-Correlations.eps", plot = Trophic.corplot, width = 6, height = 4, scale = 1)
  ggsave(filename="./graphs/Figure - Density Changes.eps", plot = Faunachange.plot, width = 5, height = 4, scale = 1)
}



##################################################
### CLEAN-UP / HOUSEKEEPING
##################################################
cat("======== Analysis of Fauna data FINISHED ========\n")
SECC.fauna <- SECC.fauna.full   # recover full data frame
