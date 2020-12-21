library(quadprog)
library(tseries)
setwd("/Users/ivy 1/Desktop/825hw4/")
df <- read.csv("47indus-day.csv")


rets <- df
rets[2:48] <- rets[2:48]-rets$rf
rets <- rets[rets$X>'20090000'&rets$X<='20131231',1:49]
rets[,2:49] <- rets[,2:49]/100

# Problem1
muall <- apply(rets[,2:48],2,mean)*252
sdall <- apply(rets[,2:48],2,sd)*sqrt(252)
covall <- cov(rets[,2:48])*252 # 47*47
plot(sdall,muall,xlim = c(0.1,max(sdall)),xlab="Standard Deviation",ylab="Mu",main="Figure 1 Mean vs Standard Deviation (2009-2013) (47 industries)",pch=20)

# Market
mumkt <- mean(rets$xrm)*252
sdmkt <- sd(rets$xrm)*sqrt(252)
points(sdmkt,mumkt,col='purple')
text(sdmkt,mumkt,labels = "M",col="purple",pos=3)

# Gn
# use solve.QP
# solve.QP(Dmat,dvec,Amat,bvec,meq,factorized)
Amat <- cbind(rep(1,47),diag(1,47))
bvec <- c(1,rep(0,47))
wgn <- solve.QP(covall,rep(0,47),Amat,bvec,meq=1,factorized=FALSE)$solution
muGn <- as.numeric(t(wgn)%*% muall)
sdGn <- sqrt(as.numeric((t(wgn)%*%covall%*%wgn)))
shGn <- muGn/sdGn
points(sdGn,muGn,pch=2,col='red')
text(sdGn,muGn,labels = "Gn",col="red",pos=2)
a <- sapply(wgn, function(x) {all.equal(x,0)})


# Frontier with no short sale
mufront <- seq(muGn ,max(muall),length = 100)
sdfront <- rep(0,length(mufront))
# formulate standard QP
Dmat <- covall
dvec <- rep(0, 47)
Amat <- cbind(muall,rep(1,47),diag(1,47))
count <- 0
for(mu in mufront){
  count <- count + 1
  bvec <- c(mu,1,rep(0,47))
  soln <- solve.QP(Dmat, dvec, Amat, bvec, meq=2, factorized=FALSE)
  weight <- soln$solution
  sdfront[count] <- sqrt(t(weight) %*% covall %*% weight)
}
lines(sdfront,mufront,type = 'l')

# Pn
# Gamma = 1
# solve.QP(Dmat,dvec,Amat,bvec,meq,factorized)
gamma <- 100
Amatp <- cbind(rep(1,47),diag(1,47))
wp <- solve.QP(covall*gamma,muall,Amatp,c(1,rep(0,47)),meq=1)$solution
meanp <- as.numeric(t(wp)%*%muall)
sigp <- sqrt(as.numeric(t(wp)%*%covall%*%wp))
sh_p <- meanp/sigp
points(sigp,meanp,pch=4,col='skyblue')
text(sdGn,muGn,labels = "Pn",col="skyblue",pos=3)
a <- sapply(wp, function(x) {all.equal(x,0)})


# Tangency portfolio
Amatt <- cbind(muall,diag(1,47))
wt <- solve.QP(covall,rep(0,47),Amatt,c(0.15,rep(0,47)),meq=1)$solution
wt <- wt/sum(wt)
muTn <- as.numeric(t(wt)%*%muall)
sdTn <- as.numeric(sqrt(t(wt)%*%covall%*%wt))
sh_t <- muTn/sdTn
points(sdTn,muTn,col="pink",pch=4)
text(sdTn,muTn,labels = "Tn",col="pink",pos=3)
abline(a=0,b=muTn/sdTn)
a <- sapply(wt, function(x) {all.equal(x,0)})


# mu=mumkt
mumkt/sh_t 
# sd=sdmkt
sdmkt *sh_t



# Problem 2
df2 <- df
df2[2:48] <- df2[2:48]-df2$rf
df2 <- df2[df2$X>"20140000"&df2$X<"20150000",1:49]
df2[2:49] <- df2[2:49]/100
muall <- apply(df2[,2:48],2,mean)*252
sdall <- apply(df2[,2:48],2,sd)*sqrt(252)
covall <- cov(df2[,2:48])*252 # 47*47
plot(sdall,muall,xlim = c(0.1,max(sdall)),xlab="Standard Deviation",ylab="Mu",main="Figure 2 Mean vs Standard Deviation (2014) (47 industries)",pch=20)

mumkt <- mean(df2$xrm)*252
sdmkt <- sd(df2$xrm)*sqrt(252)
points(sdmkt,mumkt,col='purple')
text(sdmkt,mumkt,labels = "M",col="purple",pos=3)


muGn <- as.numeric(t(wgn)%*% muall)
sdGn <- sqrt(as.numeric((t(wgn)%*%covall%*%wgn)))
shGn <- muGn/sdGn
points(sdGn,muGn,pch=2,col='red')
text(sdGn,muGn,labels = "Gn",col="red",pos=2)

meanp <- as.numeric(t(wp)%*%muall)
sigp <- sqrt(as.numeric(t(wp)%*%covall%*%wp))
sh_p <- meanp/sigp
points(sigp,meanp,pch=4,col='skyblue')
text(sigp,meanp,labels = "Pn",col="skyblue",pos=3)

muTn <- as.numeric(t(wt)%*%mum)
sdTn <- as.numeric(sqrt(t(wt)%*%covm%*%wt))
sh_t <- muTn/sdTn
points(sdTn,muTn,col="pink",pch=4)
text(sdTn,muTn,labels = "Tn",col="pink",pos=3)
