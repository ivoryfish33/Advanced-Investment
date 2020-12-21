#Problem 3
setwd("/Users/ivy 1/Desktop/825hw5/")
indus <- read.csv("47indus-day.csv")
indus[,2:50] <- indus[,2:50]/100
nind <-47
rets <-cbind(indus[,1],indus[,2:(nind+1)]-indus[,nind+3], indus[,nind+2])
year  <-floor(rets[,1]/10000)	
ivec  <-rep(1,nind)			# Vector of ones
#weights matrix
wG1<- matrix(0,ncol=nind,nrow = 2018-1992)
wG2<- wG1; wGL<-wG1; wT <- wG1; wM  <- wT 	
wB <- array(0,dim=c(nind,2018-1992,5)); wBs<-wB
targbeta<-c(0,0.5,1,1.5,2)
dim(rets)

##########################################################
for (i in 1992:2017){		# Begin Year Loop for portfolio design
  Zrets <-rets[year==i,2:(nind+1)]		# The returns
  Zfac  <-rets[year==i,nind+2]	# The market
  Zmeans<-apply(Zrets,2,mean)	# Used for momentum strategy
  Zvars <-apply(Zrets,2,var) 
  Sigmat<-cov(Zrets)				# Covariance matrix
  vinv  <-solve(Sigmat)			# inverse COV
  
  # Global MVP 
  wG1[i-1991,] <- vinv%*%ivec			
  wG1[i-1991,] <- wG1[i-1991,]/sum(wG1[i-1991,])		
  
  # 1 Factor based Covariance Matrix and betas
  # Get betas, and residual s^2, shrunk betas
  # BB'varM + diag(residual vars)
  
  varm   <- var(Zfac)				# Index variance
  Zregs  <- lsfit(Zfac,Zrets) 			# 47 regressions
  Zbetas <- Zregs$coefficients[2,]		# 47 Betas, don't need the intercepts 
  Zbetcvar<- var(Zbetas) 						# Cross-sectional variance of 47 betas 
  Zbetsds <- ls.diag(Zregs)$std.err[2,]^2		# OLS standard deviations of OLS Betas
  Zresids<- Zregs$residuals					# 47 residuals 
  Dmat   <- diag(diag(var(Zresids)))  		# Diagonal part of residual cov. mat.
  Sigmatf<- Dmat + varm * Zbetas %*% t(Zbetas) 	# Constrained Covmat 
  vinvf  <- solve(Sigmatf)					# Inverse
  Zbetass <- (Zbetcvar*Zbetas+Zbetsds)/(Zbetcvar+Zbetsds)
  
  G1beta <- Zbetas%*%wG1[i-1991,]
  wGL[i-1991,] <- wG1[i-1991,]*1.5/sum(G1beta)
  wG2[i-1991,]<- vinvf %*%ivec				# Global MVP 1-Factor Cov Matrix
  wG2[i-1991,]<- wG2[i-1991,] / sum(wG2[i-1991,])
  
  # Now code for MinVar S.T beta = 0, 0.5, 1, 1.5, 2
  # Do it first for OLS betas
  
  AA  <- sum(vinv)  # i' vinv i
  BB  <- t(Zbetas) %*% vinv %*% ivec
  CC  <- t(Zbetas) %*% vinv %*% Zbetas
  DEL <- AA*CC - BB^2
  for ( ii in 1:5){		# 5 target betas
    zbet<-targbeta[ii]
    lam1<-drop(CC- BB * zbet)
    lam2<-drop(AA*zbet-BB)
    wB[,i-1991,ii] <- (lam1*vinv %*% ivec + lam2*vinv %*% Zbetas)/as.numeric(DEL)
  }
 
  # Now do the same for shrinkage betas
  AA  <- sum(vinv) 
  BB  <- t(Zbetass) %*% vinv %*% ivec
  CC  <- t(Zbetass) %*% vinv %*% Zbetass
  DEL <- AA*CC - BB^2
  for ( ii in 1:5){		# 5 target betas
    zbet<-targbeta[ii]
    lam1<-drop(CC- BB * zbet)
    lam2<-drop(AA*zbet-BB)
    wBs[,i-1991,ii] <- (lam1*vinv %*% ivec + lam2*vinv %*% Zbetass)/as.numeric(DEL)
  }
  
  # Now Max Sharpe Ratio (tangency)
  # using expected returns proportional to market beta
  
  wT[i-1991,] <- vinv %*% Zbetas
  wT[i-1991,] <- wT[i-1991,] / sum(wT[i-1991,])
  # Momentum
  # Make a (non-optimized) momentum portfolio with weights summing to 1
  
  # In www, put a vector with 0.5 for the top 4 means, and -0.25 for the
  # bottom 4 means. 
  # Use rank(Zmeans), you want rank>43, rank<5, zero weight on the others
  # Hint, use logic:  rank(Zmeans)>43  etc...
  wM[i-1991,] <-(rank(Zmeans)>43)*0.5+(rank(Zmeans)<5)*-0.25  
    
}	# END YEAR LOOP FOR PORTFOLIO DESIGN

# Now have the weights 
# Compute realized returns -  Put them ALL in one matrix portret

portret<- matrix(0,nrow=length(year),ncol=15) #15 portfolios

for (i in 1993:2018){
  portret[year==i,1]<- ts(rets[year==i,2:48])%*% wG1[i-1992,]
  portret[year==i,2]<- ts(rets[year==i,2:48])%*% wG2[i-1992,]
  portret[year==i,3]<- ts(rets[year==i,2:48])%*% wGL[i-1992,]-(sum(wGL[i-1992,])-1)*ts(indus[year==i,50])
  portret[year==i,4]<- ts(rets[year==i,2:48])%*% wB[,i-1992,1]
  portret[year==i,5]<- ts(rets[year==i,2:48])%*% wB[,i-1992,2]
  portret[year==i,6]<- ts(rets[year==i,2:48])%*% wB[,i-1992,3]
  portret[year==i,7]<- ts(rets[year==i,2:48])%*% wB[,i-1992,4]
  portret[year==i,8]<- ts(rets[year==i,2:48])%*% wB[,i-1992,5]
  portret[year==i,9]<- ts(rets[year==i,2:48])%*% wBs[,i-1992,1]
  portret[year==i,10]<- ts(rets[year==i,2:48])%*% wBs[,i-1992,2]
  portret[year==i,11]<- ts(rets[year==i,2:48])%*% wBs[,i-1992,3]
  portret[year==i,12]<- ts(rets[year==i,2:48])%*% wBs[,i-1992,4]
  portret[year==i,13]<- ts(rets[year==i,2:48])%*% wBs[,i-1992,5]
  portret[year==i,14]<- ts(rets[year==i,2:48])%*% wT[i-1992,]
  portret[year==i,15]<- ts(rets[year==i,2:48])%*% wM[i-1992,]
}


# Compute statistics
portret <- portret[portret[,1]!=0,]
mupor <- apply(portret,2,mean)*252
sdpor <- apply(portret,2,sd)*sqrt(252)

# Market
muM<- mean(indus[indus$X>19930000,49])*252
sdM<- sd(indus[indus$X>19930000,49])*sqrt(252)

# Figure 1, All industries 1993-2018
muvec<-apply(indus[indus$X>19930000,2:48],2,mean)*252
sdvec<-apply(indus[indus$X>19930000,2:48],2,sd)*sqrt(252)
  
plot(sdvec,muvec,xlab="Standard Deviation",ylab="Mean",xlim=c(0.1,0.9),ylim=c(0.03,0.35),main="Figure 1 Mean vs. Standard Deviation (1993-2018)")
points(sdM,muM,pch=2)
abline(a=0,b=muM/sdM)
points(sdpor,mupor,pch=3:17)
legend("bottomright",legend=c("Market","G1","G2","GL","B0","B1","B2","B3","B4",
                          "S0","S1","S2","S3","S4","T","M"),pch = 2:17,cex=0.6)

#Table 1
preg <- lsfit(indus[indus$X>19930000,49],portret)
cbind(mupor,sdpor,preg$coefficients[2,],mupor/sdpor,preg$coefficients[1,],ls.diag(preg)$std.err[1,])
ttest <- preg$coefficients[1,]/ls.diag(preg)$std.err[1,]
qnorm(0.025/15)

#Figure 2
plot(preg$coefficients[2,],mupor,xlab='beta',ylab='mu(excess)',ylim=c(0,0.35),main="Figure 2 Mu(excess) vs. Beta (1993-2018)")
abline(a=0,b=muM)
abline(lsfit(preg$coefficients[2,],mupor))
      
      