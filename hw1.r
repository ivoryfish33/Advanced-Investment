# FE825 HW1

library(forecast)
library(sandwich)
library(tseries)
library(moments)
setwd("/Users/ivy/Documents/02学校/02BU/2019Spring/FE825_Advanced_Investment/Assignment/825hw1/")

# Problem 1
month_data <- read.csv("monthly.csv")
annual_data <- read.csv("annually.CSV")

month_data$Mkt <- month_data$Mkt.RF + month_data$RF
annual_data$Mkt <- annual_data$Mkt.RF + annual_data$RF


# a
m_data <- month_data[month_data[,1] >= 197012 & month_data[,1] <= 201712,]
a_data <- annual_data[annual_data[,1] >= 1970 & annual_data[,1] <= 2017,]

mrow <- dim(m_data)[1]
arow <- dim(a_data)[1]

m_data$mkt_dollar <- m_data$Mkt
m_data$mkt_dollar[1] <- 1
m_data$rf_dollar <- m_data$RF
m_data$rf_dollar[1] <- 1

a_data$mkt_dollar <- a_data$Mkt
a_data$mkt_dollar[1] <- 1
a_data$rf_dollar <- a_data$RF
a_data$rf_dollar[1] <- 1

for(i in 2:mrow){
  m_data$mkt_dollar[i] = m_data$mkt_dollar[i-1] * (1 + m_data$Mkt[i]/100)
  m_data$rf_dollar[i] = m_data$rf_dollar[i-1] * (1 + m_data$RF[i]/100)
}

for(i in 2:arow){
  a_data$mkt_dollar[i] = a_data$mkt_dollar[i-1] * (1 + a_data$Mkt[i]/100)
  a_data$rf_dollar[i] = a_data$rf_dollar[i-1] * (1 + a_data$RF[i]/100)
}


# Figure 1
ts.plot(cbind(m_data$mkt_dollar,m_data$rf_dollar), col=c("blue","red"), xlab = "Month passed", ylab="Total Wealth",main="Figure 1 Accumulated Wealth Beginning with $1 (100% $)")
ts.plot(cbind(a_data$mkt_dollar,a_data$rf_dollar), col=c("blue","red"), xlab = "Year passed", ylab="Total Wealth", main = "Figure 1 Accumulated Wealth Beginning with $1 (100% $)")


head(m_data)
head(a_data)
tail(m_data)
tail(a_data)

# annualized return
(m_data$mkt_dollar[mrow])^(1/47) -1
(m_data$rf_dollar[mrow])^(1/47) -1
((a_data$mkt_dollar[arow])^(1/47) - 1) 
((a_data$rf_dollar[arow])^(1/47) - 1 ) 



# b 

# table 1
mm_data <- month_data[month_data[,1] <= 201712,]
aa_data <- annual_data[annual_data[,1] <= 2017,]

table1 <- matrix(nrow=2, ncol=6, dimnames=list(c("Tbill","Market"),c("LogRet,Ann.miu","LogRet,Ann.std_mon","LogRet,Ann.std_ann","LogRet,ACF(1)_mon","AnnDis,miu","AnnDis,std")))

month_data_log <- log(1 + mm_data/100)
month_data_log[,1] <- mm_data[,1]

annual_data_log <- log(1 + aa_data/100)
annual_data_log[,1] <- aa_data[,1]


head(month_data_log)
mean(month_data_log$Mkt)

# log data
# mu
table1[1,1] <- mean(month_data_log$RF) * 12
table1[2,1] <- mean(month_data_log$Mkt) * 12
# std from m
table1[1,2] <- sd(month_data_log$RF) * (12^0.5)
table1[2,2] <- sd(month_data_log$Mkt) * (12^0.5)
# std from a
table1[1,3] <- sd(annual_data_log$RF)
table1[2,3] <- sd(annual_data_log$Mkt) 
# ACF
Acf(month_data_log$RF)[1]
Acf(month_data_log$Mkt)[1]
table1[1,4] <- 0.976
table1[2,4] <- 0.101
# discrete data
# mu
table1[1,5] <- mean(annual_data$RF /100)
table1[2,5] <- mean(annual_data$Mkt /100)
# std
table1[1,6] <- sd(annual_data$RF /100)
table1[2,6] <- sd(annual_data$Mkt /100)

table1
paste("1b")
# correlation 
cor(mm_data$Mkt, mm_data$RF)



# 1c

m_data$flag <- m_data$Mkt > m_data$RF
a_data$flag <- a_data$Mkt > a_data$RF
m_data$best <- m_data$Mkt
a_data$best <- a_data$Mkt

for(i in 1:mrow){ 
  if(m_data$flag[i] == FALSE){ 
    m_data$best[i] <- m_data$RF[i]
    }
}


for(i in 1:arow){
  if(a_data$flag[i] == FALSE){
    a_data$best[i] <- a_data$RF[i]
  }
}

m_data
a_data


m_data$best_dollar <- m_data$best
m_data$best_dollar[1] <- 1

a_data$best_dollar <- a_data$best
a_data$best_dollar[1] <- 1

for(i in 2:mrow){
  m_data$best_dollar[i] <- m_data$best_dollar[i-1] * (1 + m_data$best[i] / 100)
}

for(i in 2:arow){
  a_data$best_dollar[i] <- a_data$best_dollar[i-1] * (1 + a_data$best[i] / 100)
}

m_final_wealth <- m_data$best_dollar[mrow]
a_final_wealth <- a_data$best_dollar[arow]
tail(m_data)

# d
m_data0918 <- month_data[month_data[,1] >= 200901 & month_data[,1] <= 201811,]

cor(m_data0918$Mkt,m_data0918$SMB) * sd(m_data0918$SMB) / sd(m_data0918$Mkt)
cor(m_data0918$Mkt,m_data0918$HML) * sd(m_data0918$HML) / sd(m_data0918$Mkt)



# Problem 2
library(tseries)

com <- data.frame('Date' = c(2019:2059), "0fee" = 0, "20bp" = 0,"1%" = 0,"2%"= 0)
com <- data.frame(matrix(rep(10000,40*4),40,4))


comf <- data.frame(matrix(rep(10000,40*4),40,4))
for (i in 2:40){
  comf[i,] <- 10000 + comf[i-1,]*c(1.05,1.048,1.04,1.03)
}
plot(2019:2058,comf$X1,type='l',col='red',xlab='years',ylab='wealth',main='Power of Compounding',ylim=)
lines(2019:2058,comf$X2,col='blue')
lines(2019:2058,comf$X3,col='green')
lines(2019:2058,comf$X4,col='yellow')


nrow <- dim(com)[1]

for(i in 2:nrow){
  com[i,2] = 10000 + com[i-1,2]*(1+0.05)
  com[i,3] = 10000 + com[i-1,3]*(1+0.05-0.20/100)
  com[i,4] = 10000 + com[i-1,4]*(1+0.05-0.01)
  com[i,5] = 10000 + com[i-1,5]*(1+0.05-0.02)
}


#plot(2019:2058,com$X0fee,type='l',col='red',xlab='years',ylab='wealth',main='Power of Compounding',ylim=)
#lines(2019:2058,com$X20bp,col='blue')
#lines(2019:2058,com$X1.,col='green')
#lines(2019:2058,com$X2.,col='yellow')
ts.plot(cbind(com$X0fee,com$X20bp,com$X1.,com$X2.),col=c("black","red","blue","orange"),xlab="Period Passed",ylab='Wealth in $',main="Figure Plot")
# ts.plot(ts(sub_vix[,2],frequency=12,start=c(2007,1)),ylab="VIX",xlab="Day")

finalwealth <- com[nrow,]
finalwealth





