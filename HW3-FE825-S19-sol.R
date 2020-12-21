
############################
# PROBLEM 1
#1) In sample Optimization, Standard Use of the Frontier
#
setwd("/Users/ivy 1/Desktop/825hw4/")
indus<-read.csv("47indus-day.csv")

nind <-length(indus[1,])-3;nind 
rets <-cbind(indus[,1],indus[,2:(nind+1)]-indus[,nind+3], indus[,nind+2])

ybeg<-2009; yend<-2013
rets<-rets[rets[,1]>ybeg*10000&rets[,1]<(yend+1)*10000,]
dpy<-length(rets[,1])/(yend-ybeg+1);dpy

# Statistics of the 47 industries, Market,  and EW

muE <-mean(apply(rets[,2:(nind+1)],1,mean))*dpy
sigE<-sd(apply(rets[,2:(nind+1)],1,mean))*sqrt(dpy)
muM <-mean(rets[,nind+2])*dpy
sigM<-sd(rets[,nind+2])*sqrt(dpy)
ones  <- rep(1,nind)
muvec <- apply(rets[,2:(nind+1)],2,mean)*dpy
sigvec<- apply(rets[,2:(nind+1)],2,sd)*sqrt(dpy)
covmat<- cov(rets[,2:(nind+1)])*dpy
covinv<- solve(covmat)

# Quantities of the Frontier

AA<- sum(covinv)
BB<- t(muvec)%*%covinv%*%ones
CC<- t(muvec)%*%covinv%*%muvec
DEL<-AA*CC-BB^2
AA
BB
CC
DEL

# Global MVP

wG  <- covinv%*%ones/AA
muG <- BB/AA
sigG<- sqrt(1/AA)

# Frontier equation mean/sig starting from Global MVP

CC/DEL
-2*BB/DEL
AA/DEL

mufront  <-seq(muG,max(muvec),length=200)
sigfront <-sqrt( CC - 2*BB*mufront + AA*mufront^2)/sqrt(DEL) 

#
# M1*:	Min Var ST mu = muM
#

lams<-solve(matrix(c(AA,BB,BB,CC),ncol=2))%*%c(1,muM); lams
w1  <- covinv %*% (lams[1]*ones+lams[2]*muvec)
mu1 <- t(w1)%*%muvec
sig1<- sqrt(t(w1)%*%covmat%*%w1)

#
# M2*:  Min Var ST beta = 1
# Standard Min Var with w*mu=mup constraint replaced by
# w*bet=betp
# B, C, must be computed with betvec not muvec!!

betvec<-lsfit(rets[,nind+2],rets[,2:(nind+1)])$coefficients[2,]
mean(betvec)
BBbet<-t(betvec)%*%covinv%*%ones
CCbet<-t(betvec)%*%covinv%*%betvec
lams<-solve(matrix(c(AA,BBbet,BBbet,CCbet),ncol=2))%*%c(1,1); lams
w2  <- covinv %*% (lams[1]*ones+lams[2]*betvec)
mu2 <- t(w2)%*%muvec
sig2<- sqrt(t(w2)%*%covmat%*%w2)

betvec%*%w2		# It worked !

# By the way, what is the beta of the Global MVP?
betvec%*%wG
# 

#
# P:  MAX CE St
#
gam <-4
wP  <- wG + covinv%*%(muvec-muG*ones)/gam
muP <- t(wP)%*%muvec
sigP<-sqrt(t(wP)%*%covmat%*%wP)

#
# Tangency portfolio
#
wT  <-covinv%*%muvec;wT<-wT/sum(wT)
muT <- t(wT)%*%muvec
sigT<-sqrt(t(wT)%*%covmat%*%wT)

#
# The Table
#
zmuze<-c(muG,muM,muE,mu1,mu2,muP,muT)
zsigs<-c(sigG,sigM,sigE,sig1,sig2,sigP,sigT)
zSharpes<-zmuze/zsigs

write.csv(round(cbind(zmuze,zsigs,zSharpes),2),"tab1.csv")

# The Figure

par(mfrow=c(1,1),mgp=c(1.5,0.5,0),mar=c(2.6,2.6,1.6,0.1))
plot(sigfront,mufront,ylim=range(c(0,muvec)),xlim=range(0,max(sigvec))
,xlab=expression(sigma),ylab=expression(paste(mu,"-Rf")),type="l")
points(sigvec,muvec)
points(c(sigE,sigM,sigG,sig1,sig2,sigP),c(muE,muM,muG,mu1,mu2,muP),
pch="*",col="blue")
segments(0,0,sigT,muT)
text(c(sigE,sigM,sigG,sig1,sig2,sigP),c(muE,muM,muG,mu1,mu2,muP),
labels=c("E","M","G","M1","M2","P*"),col="red",pos=c(1,2,1,4,4,2))
title(paste("Figure 1: ",bquote(.(nind)),"Industries, Market, and Frontier"), line=0.2)


# Of course we could do this for every portfolio: sum(wG[wG<0])
# This may be faster:

zw  <-cbind(wG,w1,w2,wP,wT)
zneg<- zw<0		# A logical map with T/F
apply(zw*zneg,2,sum)

# Discussion
# Compare G , M2, M Sharpe ratios
# The Sharpe ratio allows us to compute the benefits of leverage

abline(0,muG/sigG,lty=2)

# Large amounts of shorts driven by one or two positions?
# A nice visual check:

plot(scale(sigvec),wG)
plot(scale(muvec),wT)

# Other helpful visual: Are two portfolio different in weights

plot(wP,wG);abline(0,1)

# But remember, most of the action is in the N(N-1) covariances

#
# Problem 4. Investing with these weights - tracking performance
#

# Let's get our 2014-16 returns
rets <-cbind(indus[,1],indus[,2:(nind+1)]-indus[,nind+3], indus[,nind+2])

ybeg<-2014; yend<-2014
rets<-rets[rets[,1]>ybeg*10000&rets[,1]<(yend+1)*10000,]
dpy<-length(rets[,1])/(yend-ybeg+1);dpy

muM4 <- mean(rets[,nind+2])*dpy
sigM4<- sd(rets[,nind+2])*sqrt(dpy)
muvec4 <- apply(rets[,2:(nind+1)],2,mean)*dpy
sigvec4<- apply(rets[,2:(nind+1)],2,sd)*sqrt(dpy)

indr <-rets[,2:(nind+1)]	 	# Just the industries
wmat <-cbind(wG,rep(1/47,47),w1,w2,wP,wT) 	# All 6 Portfolio weights
pret <-as.matrix(indr)%*%wmat	# All 6 portfolio returns

pmeans4<-apply(pret,2,mean)*dpy
psigs4 <-apply(pret,2,sd)*sqrt(dpy)
Sharp4 <-pmeans4/psigs4
write.csv(round(cbind(pmeans4,psigs4,Sharp4),2),"tab.csv")

# The realized performance Figure

par(mfrow=c(1,1),mgp=c(1.5,0.5,0),mar=c(2.6,2.6,1.6,0.1))
plot(sigvec4,muvec4,ylim=range(c(0,muvec)),xlim=range(0,max(sigvec))
,xlab=expression(sigma),ylab=expression(paste(mu,"-Rf")))
points(c(psigs4,sigM4),c(pmeans4,muM4),pch="*",col="blue")
text(c(psigs4,sigM4),c(pmeans4,muM4),
labels=c("G","E","M1","M2","P*","T","M"),col="red",pos=c(3,2,1,4,2,2,3))
title(paste("Figure 2: ",bquote(.(nind)),"Industries, 2014 Performance"), line=0.2)

#
# Doing it again for 2014-16 as a whole
#
rets <-cbind(indus[,1],indus[,2:(nind+1)]-indus[,nind+3], indus[,nind+2])

ybeg<-2014; yend<-2016
rets<-rets[rets[,1]>ybeg*10000&rets[,1]<(yend+1)*10000,]
dpy<-length(rets[,1])/(yend-ybeg+1);dpy

muM6 <- mean(rets[,nind+2])*dpy
sigM6<- sd(rets[,nind+2])*sqrt(dpy)
muvec6 <- apply(rets[,2:(nind+1)],2,mean)*dpy
sigvec6<- apply(rets[,2:(nind+1)],2,sd)*sqrt(dpy)

indr <-rets[,2:(nind+1)]	 	# Just the industries
wmat <-cbind(wG,rep(1/47,47),w1,w2,wP,wT) 	# All 6 Portfolio weights
pret <-as.matrix(indr)%*%wmat	# All 6 portfolio returns

pmeans6<-apply(pret,2,mean)*dpy
psigs6 <-apply(pret,2,sd)*sqrt(dpy)
Sharpe6<-pmeans6/psigs6
write.csv(round(cbind(pmeans6,psigs6,Sharpe6),2),"tab.csv")

# The realized performance Figure

par(mfrow=c(1,1),mgp=c(1.5,0.5,0),mar=c(2.6,2.6,1.6,0.1))
plot(sigvec6,muvec6,ylim=range(c(0,muvec)),xlim=range(0,max(sigvec))
,xlab=expression(sigma),ylab=expression(paste(mu,"-Rf")))
points(c(psigs6,sigM6),c(pmeans6,muM6),pch="*",col="blue")
text(c(psigs6,sigM6),c(pmeans6,muM6),
labels=c("G","E","M1","M2","P*","T","M"),col="red",pos=c(1,2,3,3,2,2,3))
title(paste("Figure 3: ",bquote(.(nind)),"Industries, 2014-16 Performance"), line=0.2)

# Let's try even later. 2015-2018

rets <-cbind(indus[,1],indus[,2:(nind+1)]-indus[,nind+3], indus[,nind+2])

ybeg<-2015; yend<-2018
rets<-rets[rets[,1]>ybeg*10000&rets[,1]<(yend+1)*10000,]
dpy<-length(rets[,1])/(yend-ybeg+1);dpy

muM8 <- mean(rets[,nind+2])*dpy
sigM8<- sd(rets[,nind+2])*sqrt(dpy)
muvec8 <- apply(rets[,2:(nind+1)],2,mean)*dpy
sigvec8<- apply(rets[,2:(nind+1)],2,sd)*sqrt(dpy)

indr <-rets[,2:(nind+1)]	 	# Just the industries
wmat <-cbind(wG,rep(1/47,47),w1,w2,wP,wT) 	# All 6 Portfolio weights
pret <-as.matrix(indr)%*%wmat	# All 6 portfolio returns

pmeans8<-apply(pret,2,mean)*dpy
psigs8 <-apply(pret,2,sd)*sqrt(dpy)
Sharpe8<-pmeans8/psigs8
write.csv(round(cbind(pmeans8,psigs8,Sharpe8),2),"tab.csv")

# The realized performance Figure

par(mfrow=c(1,1),mgp=c(1.5,0.5,0),mar=c(2.6,2.6,1.6,0.1))
plot(sigvec8,muvec8,ylim=range(c(0,muvec)),xlim=range(0,max(sigvec))
,xlab=expression(sigma),ylab=expression(paste(mu,"-Rf")))
points(c(psigs8,sigM8),c(pmeans8,muM8),pch="*",col="blue")
text(c(psigs8,sigM8),c(pmeans8,muM8),
labels=c("G","E","M1","M2","P*","T","M"),col="red",pos=c(1,2,3,3,2,2,3))
title(paste("Figure 4: ",bquote(.(nind)),"Industries, 2015-18 Performance"), line=0.2)


#
# Make some functions
#
glob<-function(muvec,covmat){
	ones<-rep(1,length(muvec))
	wmvp<-solve(covmat)%*%ones;wmvp<-wmvp/sum(wmvp)
	mumvp<-t(wmvp)%*%muvec
	sdmvp<-sqrt(t(wmvp)%*%covmat%*%wmvp)
	globout<-list(wmvp,mumvp,sdmvp)
	return(globout)
}


