#Problem 1
#a
p <- c(0.3,0.6,0.1)
wealth <- c(125000,100000,40000)
EUW <- sum(log(wealth)*p)
RP <- sum(wealth*p)-exp(EUW)
#b
mW <- sum(p*wealth)
varW <- sum(p*(wealth-mW)^2)
sdW <- sqrt(varW)
PA_RP <- 1/2*varW*(1/mW)
#c
CV <- sdW/mW
#It is not a good estimate, RP is 3938, PA_RP is 2685. The risk is around 25% of the wealth
#It is not a small fraction of wealth

#Problem 2
#a
rf <- 0.02
rM <- 0.13
sigM <- 0.22
ri <- 0.18
sigi <- 0.34
Sharpe_M <- (rM-rf)/sigM
Sharpe_i <- (ri-rf)/sigi
#fund M's sharpe ratio is higher
#By shortselling asset i we should be able to increase the Sharpe ratio of the whole portfolio
#b
rho <- 0.7
wM <- (sigi^2-sigM*sigi*rho)/(sigi^2+sigM^2-2*rho*sigM*sigi)
wi <- 1-wM
#c
minsd <- sqrt(wM^2*sigM^2+wi^2*sigi^2+2*wM*wi*rho*sigM*sigi)
plot(c(sigM,sigi,minsd),c(rM,ri,wM*rM+wi*ri),xlab='standard deviation',pch=18,ylab='mean',main='standard deviation vs mean',col=c('red','green','blue'),ylim=c(0,0.2),xlim=c(0,0.35))
legend("bottomright", legend=c('M','i','MVP'), col=c('red','green','blue'),pch=18)
#d
#We can see the Sharpe ratio does increase after shorselling some i
covM <- matrix(c(sigM^2,sigi*sigM*rho,sigi*sigM*rho,sigi^2),ncol=2)
wMS <- solve(covM)%*%c(rM-rf,ri-rf)/sum(solve(covM)%*%c(rM-rf,ri-rf))
points(sqrt(t(wMS)%*%covM%*%wMS),t(wMS)%*%c(rM,ri),pch=18,col='yellow')
legend("bottomright", legend=c('M','i','MVP','MaxSharpe'), col=c('red','green','blue','yellow'),pch=18)
#e
#beta
betai <- rho*sigi/sigM
#alpha
ri-rf-betai*(rM-rf)
#epsilon
sigi^2-betai^2*sigM^2
#Rsquared
rho^2

#Problem 3
setwd("F:/大学/研究生/Semester 2/FE825/Assign 3")
dt <- read.csv("47indus-day.csv")
#1
nrfdt <- dt
nrfdt[2:48] <- nrfdt[2:48]-nrfdt$rf
nrfdt <- nrfdt[nrfdt$X>'20090000'&nrfdt$X<='20131231',1:49]
nrfdt[,2:49] <- nrfdt[,2:49]/100

covm <- cov(nrfdt[,2:48])*252
mum <- apply(nrfdt[,2:48],2,mean)*252
plot(apply(nrfdt[,2:48],2,sd)*sqrt(252),mum,pch=18,xlab='standard deviation',ylab='mu',xlim = c(0,sqrt(max(diag(covm)))),ylim=c(min(mum),0.9))
points(sd(nrfdt$xrm)*sqrt(252),mean(nrfdt$xrm)*252,pch=1)
#Equal weights
we <- rep(1/47,47)
points(sqrt(t(we)%*%covm%*%we),mean(apply(nrfdt[,2:48],2,mean))*252,pch=2)
#G
wg <- (solve(covm)%*%rep(1,47))/as.numeric((t(rep(1,47))%*%solve(covm)%*%rep(1,47)))
points(sqrt(t(wg)%*%covm%*%wg),as.numeric(t(wg)%*%mum),pch=3)
A=as.numeric(t(rep(1,47))%*%solve(covm)%*%rep(1,47))
B=as.numeric(t(mum)%*%solve(covm)%*%rep(1,47))
C=as.numeric(t(mum)%*%solve(covm)%*%mum)
delta=A*C-B^2
mup <- seq(as.numeric(t(wg)%*%mum),max(mum),0.001)
sigp <- sqrt((C-2*B*mup+A*mup^2)/delta)
lines(sigp,mup,lty=1)
#M1
lambda1 <- (C-B*mean(nrfdt$xrm)*252)/delta
lambda2 <- (A*mean(nrfdt$xrm)*252-B)/delta
wm1 <- lambda1*solve(covm)%*%rep(1,47)+lambda2*solve(covm)%*%mum
points(sqrt(t(wm1)%*%covm%*%wm1),as.numeric(t(wm1)%*%mum),pch=4)
#M2, basically replace mup with beta
betam <- t(cov(nrfdt$xrm,nrfdt[,2:48])/var(nrfdt$xrm))
Bb <- as.numeric(t(betam)%*%solve(covm)%*%rep(1,47))
Cb <- as.numeric(t(betam)%*%solve(covm)%*%betam)
deltab <- A*Cb-Bb^2
lambda1b <- (Cb-Bb)/deltab
lambda2b <- (A-Bb)/deltab
wm2 <- lambda1b*solve(covm)%*%rep(1,47)+lambda2b*solve(covm)%*%betam
points(sqrt(t(wm2)%*%covm%*%wm2),as.numeric(t(wm2)%*%mum),pch=5)
#3 ???STRANGE RESULT?
gamma <- 4
wp <- 1/gamma*solve(covm)%*%(mum-as.numeric(t(wg)%*%mum)*rep(1,47))+wg
points(sqrt(t(wp)%*%covm%*%wp),as.numeric(t(wp)%*%mum),pch=6)
#4 Also a little bit strange?
wt <- solve(covm)%*%mum/as.numeric(sum(solve(covm)%*%mum))
points(sqrt(t(wt)%*%covm%*%wt),as.numeric(t(wt)%*%mum),pch=7)
x <- seq(0,0.5,0.01)
y <- as.numeric(as.numeric(t(wt)%*%mum)/sqrt(t(wt)%*%covm%*%wt))*x
lines(x,y,lty=1)
#The market portfolio is far from efficient
legend("topleft",legend=c("Market","EW","MVP","M*1","M*2","P*","T"),pch = c(1,2,3,4,5,6,7),cex=0.5)

#Problem 4
dt2 <- dt
dt2[2:48] <- dt2[2:48]-dt2$rf
dt2 <- dt2[dt2$X>'20140000'&dt2$X<='20141231',1:49]
dt2[,2:49] <- dt2[,2:49]/100
dt3 <- dt
dt3[2:48] <- dt3[2:48]-dt3$rf
dt3 <- dt3[dt3$X>'20140000'&dt3$X<='20161231',1:49]
dt3[,2:49] <- dt3[,2:49]/100

#Table and Figure 2
covm <- cov(dt2[,2:48])*252
mum <- apply(dt2[,2:48],2,mean)*252
plot(apply(dt2[,2:48],2,sd)*sqrt(252),mum,pch=18,xlab='standard deviation',ylab='mu',xlim = c(0,0.7),ylim=c(min(mum),0.4))
points(sd(dt2$xrm)*sqrt(252),mean(dt2$xrm)*252,pch=1)
points(sqrt(t(we)%*%covm%*%we),mean(apply(dt2[,2:48],2,mean))*252,pch=2)
points(sqrt(t(wg)%*%covm%*%wg),as.numeric(t(wg)%*%mum),pch=3)
points(sqrt(t(wm1)%*%covm%*%wm1),as.numeric(t(wm1)%*%mum),pch=4)
points(sqrt(t(wm2)%*%covm%*%wm2),as.numeric(t(wm2)%*%mum),pch=5)
points(sqrt(t(wp)%*%covm%*%wp),as.numeric(t(wp)%*%mum),pch=6)
points(sqrt(t(wt)%*%covm%*%wt),as.numeric(t(wt)%*%mum),pch=7)
x <- seq(0,0.5,0.01)
y <- as.numeric(as.numeric(t(wt)%*%mum)/sqrt(t(wt)%*%covm%*%wt))*x
lines(x,y,lty=1)
A=as.numeric(t(rep(1,47))%*%solve(covm)%*%rep(1,47))
B=as.numeric(t(mum)%*%solve(covm)%*%rep(1,47))
C=as.numeric(t(mum)%*%solve(covm)%*%mum)
delta=A*C-B^2
mup <- seq(as.numeric(t(wg)%*%mum),max(mum),0.001)
sigp <- sqrt((C-2*B*mup+A*mup^2)/delta)
lines(sigp,mup,lty=1)
legend("topright",legend=c("Market","EW","MVP","M*1","M*2","P*","T"),pch = c(1,2,3,4,5,6,7),cex=0.5)
rbind(c("G","M","EW","M1","M2","P*","T"),
      c(as.numeric(t(wg)%*%mum),mean(dt2$xrm)*252,mean(apply(dt2[,2:48],2,mean))*252,
        as.numeric(t(wm1)%*%mum),as.numeric(t(wm2)%*%mum),as.numeric(t(wp)%*%mum),as.numeric(t(wt)%*%mum)),
      c(sqrt(t(wg)%*%covm%*%wg),sd(dt2$xrm)*sqrt(252),sqrt(t(we)%*%covm%*%we),sqrt(t(wm1)%*%covm%*%wm1),
        sqrt(t(wm2)%*%covm%*%wm2),sqrt(t(wp)%*%covm%*%wp),sqrt(t(wt)%*%covm%*%wt)),
      c(as.numeric(t(wg)%*%mum)/sqrt(t(wg)%*%covm%*%wg),(mean(dt2$xrm)*252)/(sd(dt2$xrm)*sqrt(252)),
        as.numeric(t(we)%*%mum)/sqrt(t(we)%*%covm%*%we),as.numeric(t(wm1)%*%mum)/sqrt(t(wm1)%*%covm%*%wm1),
        as.numeric(t(wm2)%*%mum)/sqrt(t(wm2)%*%covm%*%wm2),as.numeric(t(wp)%*%mum)/sqrt(t(wp)%*%covm%*%wp),
        as.numeric(t(wt)%*%mum)/sqrt(t(wt)%*%covm%*%wt)),
      c(sum(wg[wg>0]),1,1,sum(wm1[wm1>0]),sum(wm2[wm2>0]),sum(wp[wp>0]),sum(wt[wt>0])),
      c(sum(wg[wg<0]),0,0,sum(wm1[wm1<0]),sum(wm2[wm2<0]),sum(wp[wp<0]),sum(wt[wt<0])))
#Table and Figure 3
covm <- cov(dt3[,2:48])*252
mum <- apply(dt3[,2:48],2,mean)*252
plot(apply(dt3[,2:48],2,sd)*sqrt(252),mum,pch=18,xlab='standard deviation',ylab='mu',xlim = c(0,0.74),ylim=c(min(mum),0.3))
points(sd(dt3$xrm)*sqrt(252),mean(dt3$xrm)*252,pch=1)
points(sqrt(t(we)%*%covm%*%we),mean(apply(dt3[,2:48],2,mean))*252,pch=2)
points(sqrt(t(wg)%*%covm%*%wg),as.numeric(t(wg)%*%mum),pch=3)
points(sqrt(t(wm1)%*%covm%*%wm1),as.numeric(t(wm1)%*%mum),pch=4)
points(sqrt(t(wm2)%*%covm%*%wm2),as.numeric(t(wm2)%*%mum),pch=5)
points(sqrt(t(wp)%*%covm%*%wp),as.numeric(t(wp)%*%mum),pch=6)
points(sqrt(t(wt)%*%covm%*%wt),as.numeric(t(wt)%*%mum),pch=7)
x <- seq(0,0.5,0.01)
y <- as.numeric(as.numeric(t(wt)%*%mum)/sqrt(t(wt)%*%covm%*%wt))*x
lines(x,y,lty=1)
A=as.numeric(t(rep(1,47))%*%solve(covm)%*%rep(1,47))
B=as.numeric(t(mum)%*%solve(covm)%*%rep(1,47))
C=as.numeric(t(mum)%*%solve(covm)%*%mum)
delta=A*C-B^2
mup <- seq(as.numeric(t(wg)%*%mum),max(mum),0.001)
sigp <- sqrt((C-2*B*mup+A*mup^2)/delta)
lines(sigp,mup,lty=1)
legend("bottomright",legend=c("Market","EW","MVP","M*1","M*2","P*","T"),pch = c(1,2,3,4,5,6,7),cex=0.5)
rbind(c("G","M","EW","M1","M2","P*","T"),
      c(as.numeric(t(wg)%*%mum),mean(dt3$xrm)*252,mean(apply(dt3[,2:48],2,mean))*252,
        as.numeric(t(wm1)%*%mum),as.numeric(t(wm2)%*%mum),as.numeric(t(wp)%*%mum),as.numeric(t(wt)%*%mum)),
      c(sqrt(t(wg)%*%covm%*%wg),sd(dt3$xrm)*sqrt(252),sqrt(t(we)%*%covm%*%we),sqrt(t(wm1)%*%covm%*%wm1),
        sqrt(t(wm2)%*%covm%*%wm2),sqrt(t(wp)%*%covm%*%wp),sqrt(t(wt)%*%covm%*%wt)),
      c(as.numeric(t(wg)%*%mum)/sqrt(t(wg)%*%covm%*%wg),(mean(dt3$xrm)*252)/(sd(dt3$xrm)*sqrt(252)),
        as.numeric(t(we)%*%mum)/sqrt(t(we)%*%covm%*%we),as.numeric(t(wm1)%*%mum)/sqrt(t(wm1)%*%covm%*%wm1),
        as.numeric(t(wm2)%*%mum)/sqrt(t(wm2)%*%covm%*%wm2),as.numeric(t(wp)%*%mum)/sqrt(t(wp)%*%covm%*%wp),
        as.numeric(t(wt)%*%mum)/sqrt(t(wt)%*%covm%*%wt)),
      c(sum(wg[wg>0]),1,1,sum(wm1[wm1>0]),sum(wm2[wm2>0]),sum(wp[wp>0]),sum(wt[wt>0])),
      c(sum(wg[wg<0]),0,0,sum(wm1[wm1<0]),sum(wm2[wm2<0]),sum(wp[wp<0]),sum(wt[wt<0])))
