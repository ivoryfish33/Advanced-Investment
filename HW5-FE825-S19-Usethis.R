#
#	you CAN use this as a skeleton if you want
#   you are professionals, so YOU are responsible for your
#   code, and for checking that everything you use works
#   as it should !!!
#	Question marks in the matrices are obvious things you need
#   to fill in yourself
#	Inputs: Matrix of Returns, or Mu, VCOV
#
indus<-read.csv("47indus-day.csv")
nind <-47
rets <-cbind(indus[,1],indus[,2:(nind+1)]-indus[,nind+3], indus[,nind+2])


# Create a year index to roll by year

year  <-floor(rets[,1]/10000)	
ivec  <-rep(1,nind)			# Vector of ones

# store annual portfolio weights in  #years x48 Matrices
# Or you could have just one array  #years x 48 x #portfolios 


wG1<- matrix(0,ncol=nind,nrow = 2018-199)	 # 1st weight used for 1993
wG2<- wG; wGL<-wG; wT  <- wG;wM  <- wT 	
wB <- array(0,dim=c(nind,2018-1992,5)); wBs<-wB # weight arrays for beta ports
wM <- wG1
targbeta<-c(0,0.5,1,1.5,2)	# Target Betas

##########################################################
# First estimation with 1992 goes in first row, invest in 1993
# Start with 1992, then loop around i

for (i in 1992:201?){		# Begin Year Loop for portfolio design

Zrets <-rets[year==i,]		# The returns
Zfac  <-rets[year==?,nind+2]	# The market
Zmeans<-apply(Zrets,2,mean)	# Used for momentum strategy
Zvars <-apply(Zrets,2,var) 
Sigmat<-cov(?)				# Covariance matrix
vinv  <-solve(?)			# inverse COV

# Global MVP 

wG1[i-199?,] <- vinv%*%ivec			
wG1[i-199?,] <- wG1[i-19?,]/sum(?)		

# 1 Factor based Covariance Matrix and betas
# Get betas, and residual s^2, shrunk betas
# BB'varM + diag(residual vars)

varm   <- var(Zfac)				# Index variance
Zregs  <- lsfit(Zfac,?) 			# 47 regressions
Zbetas <- Zregs$coefficients[?,]		# 47 Betas, don't need the intercepts 
Zbetcvar<- 						# Cross-sectional variance of 47 betas 
Zbetsds <- 						# OLS standard deviations of OLS Betas
Zresids<- Zregs$?					# 47 residuals 
Dmat   <- diag(diag(var(?)))  		# Diagonal part of residual cov. mat.
Sigmatf<- Dmat + varm * ?? %*% t(?) 	# Constrained Covmat 
vinvf  <- solve(?)					# Inverse

wBs<-	?

wG2[?,]<- ? %*%ivec				# Global MVP 1-Factor Cov Matrix
wG2[?,]<- wG2[?,] / sum(wG2[ ? ,])

# Now code for MinVar S.T beta = 0, 0.5, 1, 1.5, 2
# Do it first for OLS betas

AA  <- sum( ? ) 
BB  <- t(Zbetas) %*% vinv %*% ? 
CC  <-  ?  %*%vinv%*% ?
DEL <- AA*CC - BB^2

for ( ii in 1:5){		# 5 target betas
	zbet<-targbeta[ii]
	lam1<-drop(CC- BB * ? )
	lam2<-drop( ? )
	wB[?,?] <- (lam1*vinv%*%ivec + ? *vinv%*% ? ) / ?
}

# Now do the same for shrinkage betas

?


# Now Max Sharpe Ratio (tangency)
# using expected returns proportional to market beta

wT[?,] <- ? %*% Zbetas
wT[?,] <- wT[?,] / sum (?)


# Momentum
# Make a (non-optimized) momentum portfolio with weights summing to 1

# In www, put a vector with 0.5 for the top 4 means, and -0.25 for the
# bottom 4 means. 
# Use rank(Zmeans), you want rank>43, rank<5, zero weight on the others
# Hint, use logic:  rank(Zmeans)>43  etc...

wMom[?,] <- 

}	# END YEAR LOOP FOR PORTFOLIO DESIGN

# Now have the weights 
# Compute realized returns -  Put them ALL in one matrix portret

portret<- matrix(0,nrow=length(year),ncol=?) # You made ? portfolios

for (i in 1993:2018){
portret[year==i,1]<- ts(allrets[year==i,])%*% wG1[i-199?,]
portret[year==i,2]<- ts(allrets[year==i,])%*% wG2[?,]
# and so on
}

# Compute statistics
mupor <- apply(portret[?,],?,mean)*252
sdpor <- 

# Market
muM<-
sdM<-


# All industries 1993-2018
muvec<-apply(
sdvec<-apply(

plot(sdvec,muvec,xlab="Std.Dev.",ylab="Mean",ylim=c(),xlim=c())
# And the rest


# Realized Betas for beta targets

bees<-cov(portret[dayvec,],factor[?]) / var( ? )
bees

# "Realized" security market line plot (return vs beta)



