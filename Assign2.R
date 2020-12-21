#define a dirac delta function
dirac<-function(n) {
  y <- rep(0, length(n))
  y[n==0] = 1
  return(y)
}


#define a function that performs fft on Heston process
Heston_fft<-function(alpha, n, B, K, params) {
  r<-params[7]
  T<-params[8]
  S0<-params[6]
  N<-2^n
  Eta<-B/N
  Lambda_Eta<-2*pi/N
  Lambda<-Lambda_Eta/Eta
  J = 1:N
  vj = (J-1)*Eta
  m = 1:N
  Beta = log(S0) - Lambda*N / 2
  km = Beta + (m-1)*Lambda
  i<-complex(real=0, imaginary=1)
  #calculate values of characteristic function
  Psi_vj = rep(0, length(J))
  for (zz in 1:N) {
    u <- (vj[[zz]] - (alpha+1.0)*i)
    numer <- Heston_cf(u, params) 
    denom <- ((alpha+i*vj[[zz]])*(alpha+1.0+i*vj[[zz]]))
    Psi_vj[[zz]] <- (numer / denom)
  }
  #compute fft
  XX = (Eta/2)*Psi_vj*exp(-i*Beta*vj)*(2-dirac(J-1))
  ZZ = fft(XX)
  #calculate option prices
  Multiplier = exp(-alpha*km)/pi
  ZZ2 <- Multiplier*ZZ
  Km <- exp(km)
  #discard strikes that are 0 or infinity to avoid errors in interpolation
  inds <- (!is.infinite(Km) & !is.nan(Km) & (Km > 1e-16) & (Km < 1e16))
  px_interp <- approx(Km[inds], Re(ZZ2[inds]), method = "linear", xout=K)
  fft_price = Re(exp(-r*T)*px_interp$y)
  return(fft_price)
}


Heston_cf<-function(u, params) {
  sigma<-params[1]
  v0<-params[2]
  k<-params[3]
  rho<-params[4]
  theta<-params[5]
  S0<-params[6]
  r<-params[7]
  T<-params[8]
  q<-0
  i<-complex(real=0, imaginary=1)
  lemda<-sqrt(sigma^2*(u^2+i*u)+(k-i*rho*sigma*u)^2)
  w<-exp(i*u*log(S0)+i*u*(r-q)*T+k*theta*T*(k-i*rho*sigma*u)/sigma^2)/(cosh(lemda*T/2)+(k-i*rho*sigma*u)/lemda*sinh(lemda*T/2))^(2*k*theta/sigma^2)
  phi<-w*exp(-(u^2+i*u)*v0/(lemda/tanh(lemda*T/2)+k-i*rho*sigma*u))
  return(phi)
}



sigma = 0.2
v0 = 0.08
k = 0.7
rho = -0.4
theta = 0.1
S0 = 250
r = 0.02
T = 6/12
params <- c(sigma,v0,k,rho,theta,S0,r,T)
K = seq(230, 270, by=2.5)
t1 = Sys.time()
alpha <- 1
n <- 14
B <- 300.0
Heston_fft(alpha, n, B, K, params)[9]
t2 = Sys.time()
print(call_px)
put_px = vg_fft(-alpha, n, B, K, params )
print(put_px)


K = seq(230, 270, by=2.5)
alpha <- 1
temp = c()
for (n in c(10,11,12,13,14)){
  for (B in c(150,200,250,300)){
    t1 = Sys.time()
    temp2 = (Heston_fft(alpha, n, B, K, params )[9]-21.26887)^2
    t2 = Sys.time()
    temp = c(temp,1/(temp2*as.numeric(t2-t1)))
  }
}
plot(temp,type='l',ylab = "Efficiency")


alpha <- 1
n <- 14
B <- 200.0
Heston_fft(alpha, n, B, 260, params )

alpha <- 1
temp = c()
for (n in c(9,10,11,12,13)){
  for (B in c(100,150,200,250)){
    t1 = Sys.time()
    temp2 = (Heston_fft(alpha, n, B, 260, params )-16.77443)^2
    t2 = Sys.time()
    temp = c(temp,1/(temp2*as.numeric(t2-t1)))
  }
}
plot(temp,type = "l",ylab = "Efficiency")