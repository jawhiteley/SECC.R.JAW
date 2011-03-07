##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# Response Variable(s) @ time #s
# R v2.10.1
##################################################
## INITIALISE
##################################################
rm(list=ls())	# clear memory
setwd("/Users/jonathan/Documents/ My Documents/PhD/Analysis")	# Set Working Directory: replace with a path in quotes "".
library(car)	# for recode
library(lattice)	# mostly for xyplot
library(effects)	# for plotting effect sizes?
# library(nlme)	# for mixed-effects models

#########################
## LOAD & Process DATA
#########################
# N-fixation
NfixFull <- read.csv('ARAdata.csv')
head(NfixFull)
str(NfixFull)
NfixFull <- within( NfixFull, {
	Block 	<- as.factor(Block)
	Warming <- as.factor(Warming)
	Frag <- as.factor(Frag)
	pos <- as.factor(pos)
})
NfixFull$ARA2_md <- as.numeric(NfixFull$ARA2_md)
## new factor column with simplified values
NfixFull$Chamber <- as.factor(NfixFull$Warming)
	levels(NfixFull$Chamber) <- c("", "Ambient", "Partial Chamber", "Full Chamber")	 ## empty values = empty factor level
levels(NfixFull$Chamber) <- list( "Ambient"='A', "Partial Chamber"='B', "Full Chamber"='C' )	## empty values converted to 'NA'
## the easy way - maintains empty values as "NA"
NfixFull$Chamber <- recode(NfixFull$Warming, "'A'='Ambient'; 'B'='Partial Chamber'; 'C'='Full Chamber'", 
					as.factor.result=TRUE, levels=c("Ambient", "Partial Chamber", "Full Chamber"))
NfixFull$Frag <- recode(NfixFull$Frag, "1='Continuous'; 2='Full Corridors'; 3='Pseudo-Corridors'; 4='Isolated'", 
					as.factor.result=TRUE, levels=c("Continuous", "Full Corridors", "Pseudo-Corridors", "Isolated"))
NfixFull$pos <- factor(NfixFull$pos, levels=c("1", "S", "W", "E", "N", "0"))	# safely reorder factor levels
## OR
	levels(NfixFull$pos) <- list( 'I'='1', 'S'='S', 'W'='W', 'E'='E', 'N'='N', 'O'='0' )
## new factor column with combined factor levels
NfixFull$pos1 <- as.factor(NfixFull$pos)
levels(NfixFull$pos1)[which( levels(NfixFull$pos1)=="I" )] <- "Inner"
levels(NfixFull$pos1)[1] <- "Inner"
levels(NfixFull$pos1)[c(2, 3, 4, 5)] <- "Other"
levels(NfixFull$pos1)[3] <- "Outer"
## the easy way
NfixFull$pos1 <- recode(NfixFull$pos, "'I'='Inner'; 'O'='Outer'; else='Other'",
					as.factor.result=TRUE, levels=c("Inner", "Other", "Outer"))
NfixFull$WxP <- paste(NfixFull$Warming, NfixFull$pos1)	# Warming by position
# filter samples only
Nfix <- NfixFull[NfixFull$SampleControl=="Sample",]	# if a filter is needed
# less confusing column names
Nfix$ARA <- Nfix$ARA2_md
with( Nfix, {
	ARA <- replace( ARA, which( ARA<0 ), 0 )
})
Nfix$ARA <- with( Nfix, replace( ARA, which( ARA<0 ), 0 ) ) ## replace values less than 0 with 0 - necessary for transformations.
# controls
ARAc1 <- NfixFull[NfixFull$SampleControl=="control1",]	# filter
ARAc2 <- NfixFull[NfixFull$SampleControl=="Control2",]	# filter

#########################
## EXPLORE - PLOTS
#########################
# plot.design() might be useful, but it does not accept nested designs (nested errors)
## trellis plots
xyplot( ARA ~ Block | Frag + Chamber, data=Nfix, pch=21, col="black", bg="grey", cex=0.8)
## check normality & transform data if necessary
Nfix$ARA.ln <- log( Nfix$ARA +1 )	# defaults to base e=exp(1).
Nfix$ARA.log <- log10( Nfix$ARA +1 )	# base 10.
Nfix$ARA.trans <- ( Nfix$ARA ) ^ (1/4)	# fourth-root - NaNs produced by -ve values (-14)?.
Nfix$ARA.std <- ( Nfix$ARA - mean(Nfix$ARA) ) / var(Nfix$ARA)	# standardized - doesn't change shape: useless
par(mfrow=c(2,2), cex=0.6)	# panel of figures: 2 rows & 2 columns
hist( Nfix$ARA )
hist( Nfix$ARA.log )
hist( Nfix$ARA.ln )
hist( Nfix$ARA.trans )
qqnorm( Nfix$ARA, main="untransformed" )
qqline( Nfix$ARA )
qqnorm( Nfix$ARA.log, main="log-10 transformed" )
qqline( Nfix$ARA.log )
qqnorm( Nfix$ARA.ln, main="ln transformed" )
qqline( Nfix$ARA.ln )
qqnorm( Nfix$ARA.trans, main="4th-root" )
qqline( Nfix$ARA.trans )

#########################
## MODELS
#########################
# with position as a 'fixed effect' (N,S,E,W recoded as 'other')
ARA.model <- with( Nfix, ARA ~ Warming*Frag*pos1+Error(Block/Warming/Frag) )
# ignoring effect of position: 'regional' effects only
ARA.model2 <- with( Nfix, ARA ~ Warming*Frag+Error(Block/Warming/Frag) )
ARA.ln.model <- with( Nfix, ARA.ln ~ Warming*Frag*pos1+Error(Block/Warming/Frag) )

#########################
## ANALYSIS
#########################
## with position as a 'fixed effect' (N,S,E,W recoded as 'other')
ARA.anova <- aov(ARA.ln.model)
summary(ARA.anova)
par(mfrow=c(2,2))	# panel of figures: 2 rows & 2 columns
with( Nfix, interaction.plot( Frag, Chamber, ARA ) )
with( Nfix, interaction.plot( pos1, Chamber, ARA ) )
with( Nfix, interaction.plot( pos1, Frag, ARA ) )

## lm does not accept nested error terms - see nlme package 
# ARA.lm <- with( ARA, lm( ARA2_md ~ Warming*Frag*pos1+Error(Block/Warming/Frag) ) )
# anova(ARA.lm)
## ignoring effect of position: 'regional' effects only
ARA.anova2 <- with( Nfix, aov( ARA2_md ~ Warming*Frag+Error(Block/Warming/Frag) ) )
summary(ARA.anova2)
## random effects using lme
ARA.lme <- with( Nfix, lme( fixed = ARA2_md ~ Warming*Frag*pos1, random = ~1|Block/Warming/Frag ) )
summary( ARA.lme )

##################################################
## CHECK ASSUMPTIONS (analyse residuals, standard diagnostic plots)
##################################################
par(mfrow=c(3,2))	# panel of figures: 3 rows & 2 columns
#aov
names(ARA.anova)
plot(ARA.anova)	# are residuals well-behaved?
hist(resid(ARA.anova))	# plot residuals?
with(Nfix, plot(ARA.anova, ARA~fitted(.) ) )	# is response a reasonably linear combination of fitted values?
with(Nfix, qqnorm(ARA.lme, ~resid(.)|Block ) )	# are residuals normally distributed within blocks?

#lme
plot(ARA.lme)	# are residuals well-behaved?
hist(resid(ARA.lme))	# plot residuals
with(Nfix, plot(ARA.lme, ARA2_md~fitted(.) ) )	# is response a reasonably linear combination of fitted values?
with(Nfix, qqnorm(ARA.lme, ~resid(.)|Block ) )	# are residuals normally distributed within blocks?

#########################
## GET RESULTS OF ANALYSIS
#########################
summary(ARA.lme)
## factor-level comparisons (multiple comparisons?)
with( Nfix, pairwise.t.test(ARA, Block) )
with( Nfix, pairwise.t.test(ARA, Frag) )
with( Nfix, pairwise.t.test(ARA, Warming) )
with( Nfix, pairwise.t.test(ARA, pos1) )
with( Nfix, pairwise.t.test(ARA, WxP) )	#Warming by position
library(multcomp)	# multiple comparisons package
ARA.mc <- glht(ARA.lme)	# General Linear Hypotheses?

#########################
## PLOT RESULTS
#########################
par(mfrow=c(2,2))	# panel of figures: 2 rows & 2 columns
with( Nfix, interaction.plot(Frag, Chamber, ARA2_md) )
with( Nfix, interaction.plot(pos1, Chamber, ARA2_md) )
ARA.effects <- allEffects(ARA.anova)

Chamber.map <-	data.frame( label=levels(SECC$Chamber), col=c("#000000","#000099","#990000"), bg=c("#FFFFFF","#FFFFFF","#FFFFFF"), pch=c(21,23,18), lty=c(3,2,1) )	
## matplot is convenient if tapply() is used to generate summary of means.
means <- with( SECC, 
	tapply( cbind(Y.log), 
		list(pos=pos1, Chamber=Chamber), 
		mean 
	) 
)
matplot( means, type='l', 
	xlab=levels(SECC$pos1),
	col=as.character(Chamber.map$col), 
	lty=Chamber.map$lty, lwd=2 
)	# draw lines
	# col=Chamber.map$col[match(plot.means$Chamber, label)])
matpoints( means, pch=Chamber.map$pch, col=as.character(Chamber.map$col), bg=as.character(Chamber.map$bg), cex=2, lty=Chamber.map$lty, lwd=2 )	# add points **

#########################
## 95% CI
tcrit <- abs( qt(p=0.025, df=336) )	# 2-tailed?
WarmingA <- Nfix$ARA2_md[Nfix$Warming %in% "A"]
WarmingA.mean <- mean(WarmingA)
WarmingA.df <- length(WarmingA)-1
WarmingA.se <- sqrt( var(WarmingA)/length(WarmingA) )	#se()
WarmingA.mean + ( abs( qt(p=0.025, df= WarmingA.df) ) * WarmingA.se )