#---------------------------------------------------------------------------------------------
# Comparing result of  flowDensity with manual gating
# Author: Mehrnoush MAlek
#Created: April 2014
#Last modified: August 2014
#This script is only made to reproduce the Figure 1 in the flowDensity Manuscript
#---------------------------------------------------------------------------------------------
library(flowCore)
library(flowDensity)
f <- read.FCS("~/Fig1-a_Sample.fcs")
message("tube: ", f@description$'TUBE NAME')


message("all")
channels <- c("FSC-A", "DAPI-A")
live <- flowDensity(f, channels=channels, position=c(T,T), upper=c(NA,F))


##############################################################################################
#                                     Gating Granulocytes
##############################################################################################
channels <- c("Pacific Blue-A","PE-Texas Red-A")
lf <- live
tmp.lf <- flowDensity(lf@flow.frame, channels=channels, position=c(NA,T))
hgate <- tmp.lf@gates[2]
granulo <- flowDensity(tmp.lf@flow.frame, channels=channels, position=c(T,NA), ellip.gate=T)
tmp.lf <- flowDensity(lf@flow.frame, channels=channels, position=c(NA,F))
hgate <- tmp.lf@gates[2]
CD5 <- flowDensity(tmp.lf@flow.frame, channels=channels, position=c(T,F), upper=c(NA,T),
                   ellip.gate=T, scale=0.9,)
CD5@proportion <- CD5@cell.count/lf@cell.count*100
CD5@gates[2] <- hgate
tmp.lf <- notSubFrame(lf@flow.frame, channels=channels, position=c(T,T),gates=granulo@gates)
not.granCD5 <- notSubFrame(tmp.lf@flow.frame, channels=channels, gates=CD5@gates, 
                           position=c(T,F))

##############################################################################################
#                                   Gating Monocytes
##############################################################################################
message("NOT(Granulo or CD5)")
channels <- c("BV570-A","PE-Texas Red-A")
lf <- not.granCD5
tmp.lf <- flowDensity(lf@flow.frame, channels=channels, position=c(T,T))
mono <- flowDensity(tmp.lf@flow.frame, channels=channels, position=c(T,T), ellip.gate=T)
mono@proportion <- mono@cell.count/lf@cell.count*100
not.mono <- notSubFrame(lf@flow.frame, channels=channels, gates=mono@gates, position=c(T,T))

##############################################################################################
#                                 Gating Macro
##############################################################################################
message("NOT(Mono)")
channels <- c("PE-A", "PE-Texas Red-A")
lf <- not.mono
tmp.lf <- flowDensity(lf@flow.frame, channels=channels, position=c(T,NA))
macro <- flowDensity(tmp.lf@flow.frame, channels=channels, position=c(T,NA), ellip.gate=T)
macro@proportion <- macro@cell.count/lf@cell.count*100
not.macro <- notSubFrame(lf@flow.frame, channels=channels, gates=macro@gates, position=c(T,NA))

##############################################################################################
#                               Gating NK cells
##############################################################################################
message("NOT(Macro)")
channels <- c("FITC-A", "APC-A")
lf <- not.macro
nk <- flowDensity(lf@flow.frame, channels=channels, position=c(NA,T), scale=0.9, ellip.gate=T)

##############################################################################################
#                           Plotting the populations
##############################################################################################
x11()
par(mfrow=c(4,4),mar=c(4,5,4,3))

plotDens(live@flow.frame, c("Pacific Blue-A","PE-Texas Red-A"), main="Singlets")
lines(granulo@filter,type="l")
points(granulo@flow.frame@exprs[, c("Pacific Blue-A", "PE-Texas Red-A")], pch=".", 
       col=colors()[53])
text(x=3, y=3, labels="Granulo")

## NOT(Granulo OR CD5)
plotDens(not.granCD5@flow.frame, c("BV570-A","PE-Texas Red-A"), main="NOT(Granulo or CD5)")
lines(mono@filter,type="l")
points(mono@flow.frame@exprs[, c("BV570-A","PE-Texas Red-A")], pch=".", col=colors()[41])
text(x=3.5, y=3, labels="Mono")

## NOT(Mono)
channels <- c("PE-A", "PE-Texas Red-A")
plotDens(not.mono@flow.frame, c("PE-A", "PE-Texas Red-A"), main="NOT(Mono)")
points(macro@flow.frame@exprs[, c("PE-A", "PE-Texas Red-A")], pch=".", col=colors()[76])
lines(macro@filter,type="l")
text(x=3.5, y=1, labels="Macro")

## NOT(Macro)
plotDens(not.macro@flow.frame, c("FITC-A", "APC-A"), main="NOT(Macro)")
lines(nk@filter,type="l")
points(nk@flow.frame@exprs[, c("FITC-A", "APC-A")], pch=".", col=colors()[116])
text(x=3, y=3.5, labels="NK cells")

#---------------------------------------------------------------------------------------------
# Comparing result of SamSPECTRAL, flowMeans, XCYTE with flowDensity 
#---------------------------------------------------------------------------------------------
#FCS file downloaded from flowRepository mentioned in the mani manuscript
frame <- read.FCS("~/Fig1-d-TransformedSample.fcs")

##############################################################################################
#Running flowDensity
##############################################################################################
temp <- flowDensity(frame,c("FSC-A","SSC-A"),position=c(T,NA),percentile=c(.05,NA),
                    use.percentile=c(T,F))
fs.gate <-  deGate(temp@flow.frame,"FSC-A",tinypeak.removal=1/50)
lymph.1 <- flowDensity(temp,c("FSC-A","SSC-A"),position=c(F,F),gates=c(fs.gate,NA),
                       ellip.gate=T, scale=.98)
lymph.1@proportion <- (lymph.1@cell.count/nrow(frame))*100
Tcell <- flowDensity(lymph.1@flow.frame,c(11,7),position=c(T,F),percentile=c(NA,.97))
CD4 <- flowDensity(Tcell@flow.frame,c(9,14),position=c(T,F))
CD8 <- flowDensity(Tcell@flow.frame,c(9,14),position=c(F,T))
##############################################################################################
#Running flowMeans
##############################################################################################
library(flowMeans)
#Only considering channels used in flowDensity for pre-gating
fmean.6<- flowMeans(x=frame,colnames(frame)[c(1,4,7,9,11,14)],MaxN=14)
##############################################################################################
#Running flowMeans
##############################################################################################
#We ran SamSPECTRAL several times to find the optimal value of normal.sigma
#Colors might not be identical to the Figure on the manuscript, as 
#SamSPECTRAL performs faithful sampling and assign random colors to clusters.
library(SamSPECTRAL)
s.result<-SamSPECTRAL(data.points=exprs(frame),
                      dimensions=c(1:ncol(frame)),
                      normal.sigma=700,
                      separation.factor=1.4,
                      talk=TRUE,
                      return_only.labels=T)
##############################################################################################
#Running XCyt 
##############################################################################################
#As explained in supplementary material of the X-Cyt manuscript,
#We also tried downsampling of the data before running X-Cyt
####To use X-Cyt algorithm:
#1)Download the package from http://www.broadinstitute.org/mpg/xcyt/
#2)Install EMMIX-skew Package and source Emmixskew.R from X-Cyt package
seed = 1111
init <- NULL
data <- exprs(frame)
obj<-EmSkew(data,g=10,distr="mvn",ncov=3,seed=seed,init=init)
clust <- obj$clust

##############################################################################################
#Plotting the comparison
#Creating labels for flowDensity populations
ind.1<-which(!is.na(exprs(Tcell@flow.frame)[,1]))
ind.2<-which(!is.na(exprs(CD4@flow.frame)[,1]))
ind.3<-which(!is.na(exprs(CD8@flow.frame)[,1]))
lbls<- rep(1,nrow(frame))

lbls[ind.1]<-4
lbls[ind.2]<-2
lbls[ind.3]<-3

#flowDensity
plot(exprs(frame)[,c(9,14)],pch=".",col=lbls, xlab="CD4", ylab="CD8",
     main="All Cells",cex=1.7,cex.lab=1.6)
legend("topright","flowDensity" ,cex=2,box.lty=2)
points(CD4@filter,type="l",col=2,lwd=4)
points(CD8@filter,type="l",col=2,lwd=4)
#Fmeans
plot(exprs(frame)[,c(9,14)],pch=".",col=fmean.6@Label,xlab="CD4",
     ylab="CD8", main="All Cells",cex=1.7,cex.lab=1.6)
legend("topright", "flowMeans",cex=1.52,box.lty=2)
#SamSPEC
plot(exprs(frame)[,c(9,14)],pch=".",col=s.result,xlab="CD4", ylab="CD8",
     main="All Cells",cex=1.7,cex.lab=1.6)
legend("topright", "SamSPECTRAL",cex=1.5,box.lty=2)
#XCYTE
plot(exprs(frame)[,c(9,14)],pch=".",col=clust,xlab="CD4", ylab="CD8",
     main="All Cells",cex=1.7,cex.lab=1.6)
legend("topright", "XCYTE",cex=1.5,box.lty=2)
