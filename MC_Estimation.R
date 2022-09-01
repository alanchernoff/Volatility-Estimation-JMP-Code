
#Model Parameters

t=23400 #number of seconds in a 6.5 hour trading day
m = 0.05 
phi = 5  
v = 0.04 
eta = 0.5 
rho = -0.5 
del = (1/252)/t 
sdel = sqrt(del)
srho = sqrt(1-rho^2)

#Simulation set up 

n = 10000 #number of simulations

blanks=matrix(0,nrow=n,ncol=(t-1))
y = cbind(matrix(m,nrow=n,ncol=1),blanks)
S2 = cbind(matrix(v,nrow=n,ncol=1),blanks)
W1 = matrix(rnorm(n*t),nrow=n,ncol=t)
n1 = matrix(rnorm(n*t),nrow=n,ncol=t)
W2 = (n1*srho +rho*W1)

#To test for the correlation of -0.5 between brownian motions W1 and W2
#use cor(W1[i,],W2[i,]) for any i < number of simulations


#Heston Volatility Equations

for (i in 1:(t-1)){
  S2[,i+1] <- S2[,i] +phi*(v - S2[,i])*del + eta*sqrt(S2[,i])*W2[,i]*sdel
  y[,i+1] <- y[,i] + (m - (S2[,i])/2 )*del + sqrt(S2[,i])*W1[,i]*sdel 
}

#Generate data for simulated 10, and 5 minute intervals

y_10 = y[,seq(1,ncol(y), 600)]
y_5 = y[,seq(1,ncol(y), 300)]
y_2.5 = y[,seq(1,ncol(y), 150)]
y_1 = y[,seq(1,ncol(y), 60)]


#Volatility Calculation

RV = (rowSums(S2))*del
RV_new <- matrix(RV,nrow=100,byrow=TRUE)

write.csv(y_5,"C:\\Users\\chernofa\\Downloads\\MC Data no jumps 5min.csv")
write.csv(y_2.5,"C:\\Users\\chernofa\\Downloads\\MC Data no jumps 2.5min.csv")
write.csv(y_1,"C:\\Users\\chernofa\\Downloads\\MC Data no jumps 1min.csv")
write.csv(RV,"C:\\Users\\chernofa\\Downloads\\MC Data no jumps 5min.csv")
