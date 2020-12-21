setwd("/Users/ivy 1/Desktop/825hw5/")
indus<-read.csv("47indus-day.csv")
nind <-47
rets <-cbind(indus[,1],(indus[,2:(nind+1)]-indus[,nind+3]), indus[,nind+2])
rets[,2:49] <- rets[,2:49]/100
rets <- rets[rets[,1]>'19920000'&rets[,1]<='20181231',]
rets[,1] <- floor(rets[,1]/10000)	
nday <- length(rets[,1])/(2018-1992)

# Create a year index to roll by year
year  <-floor(rets[,1]/10000)	
ivec  <-rep(1,nind)			# Vector of ones

# store annual portfolio weights in  #years x48 Matrices
# Or you could have just one array  #years x 48 x #portfolios 

wG<- matrix(0,ncol=nind,nrow = 2018-1992)	 # 1st weight used for 1993
wG2<- wG; wGL<-wG; wT <- wG;wM  <- wT 	
wB <- array(0,dim=c(nind,2018-1992,5)); wBs<-wB # weight arrays for beta ports
wM <- wG1
targbeta<-c(0,0.5,1,1.5,2)	# Target Betas



for(i in year){
  mu_1y <- apply(rets[rets[,1]>=i&rets[,1]<i+1,2:(nind+1)],2,mean)*nday
  sd_1y <- apply(rets[rets[,1]>=i&rets[,1]<i+1,2:(nind+1)],2,sd)*sqrt(nday)
  cov_1y <- cov(rets[rets[,1]>=i&rets[,1]<i+1,2:(nind+1)])*nday
  D_1y <- diag(cov_1y)
  mumkt <- mean(rets[rets[,1]>=i&rets[,1]<i+1,nind+2])*nday
  sdmkt <- sd(rets[rets[,1]>=i&rets[,1]<i+1,nind+2])*sqrt(nday)
  wG[i,] <- (solve(cov_1y)%*%ivec)/as.numeric((t(ivec)%*%solve(cov_1y)%*%ivec))
  cov_G2 <- 
}
