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

n = 100 #number of simulations

#Large Jumps No Noise

blanks=matrix(0,nrow=n,ncol=(t-1))
y = cbind(matrix(m,nrow=n,ncol=1),blanks)
S2 = cbind(matrix(v,nrow=n,ncol=1),blanks)
S2_jumps = matrix(0,nrow=n,ncol=t)
W1 = matrix(rnorm(n*t),nrow=n,ncol=t)
n1 = matrix(rnorm(n*t),nrow=n,ncol=t)
W2 = (n1*srho +rho*W1)

#To test for the correlation of -0.5 between brownian motions W1 and W2
#use cor(W1[i,],W2[i,]) for any i < number of simulations

#Jump generation

jmin = 0 ;
jmax = 1/252 ;
jp = jmin + runif(n)*(jmax - jmin);
jt = round(t*jp*252);
ser = 1:t ;
jump=matrix(0,nrow=t,ncol=n)
for (i in 1:t){
  jump[i,] = ser[i];
}

for (j in 1:n){
  for (i in 1:t){
    if (jump[i,j] == jt[j]){
      jump[i,j] = 1;
    }else{
      jump[i,j] = 0;
    }
  }
}

j_price = t(matrix(rnorm(n*t,0,0.02^2),nrow=t,ncol=n)*jump) #rnorm(n*t,0,0.02^2)

j_vol = t(matrix(exp(rnorm(n*t,-5,1)),nrow=t,ncol=n)*jump)

#Heston Volatility Equation pre-jumps (calculated separately for ease of IV calculation)

for (i in 1:(t-1)){
  S2[,i+1] <- S2[,i] +phi*(v - S2[,i])*del + eta*sqrt(S2[,i])*W2[,i]*sdel
}

#Heston Volatility Equation including jumps

for (i in 0:(t-1)){
  S2_jumps[,i] <- S2[,i] + j_vol[,i]
}

#Heston Volatility Price Equation

for (i in 1:(t-1)){
  y[,i+1] <- y[,i] + (m - (S2[,i])/2 )*del + sqrt(S2[,i])*W1[,i]*sdel + j_price[,i]
}

#Generate data for simulated 10, and 5 minute intervals

y_10 = y[,seq(1,ncol(y), 600)]
y_5 = y[,seq(1,ncol(y), 300)]
y_2.5 = y[,seq(1,ncol(y), 150)] 
y_1 = y[,seq(1,ncol(y), 60)] 

#Volatility Calculation

RV = (rowSums(S2))*del

write.csv(y_5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps 5min.csv", row.names = TRUE)
write.csv(y_2.5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps 2.5min.csv", row.names = TRUE)
write.csv(y_1,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps 1min.csv", row.names = TRUE)
write.csv(RV,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps RV.csv", row.names = TRUE)

#Small Jumps

blanks=matrix(0,nrow=n,ncol=(t-1))
y = cbind(matrix(m,nrow=n,ncol=1),blanks)
S2 = cbind(matrix(v,nrow=n,ncol=1),blanks)
S2_jumps = matrix(0,nrow=n,ncol=t)
W1 = matrix(rnorm(n*t),nrow=n,ncol=t)
n1 = matrix(rnorm(n*t),nrow=n,ncol=t)
W2 = (n1*srho +rho*W1)

#To test for the correlation of -0.5 between brownian motions W1 and W2
#use cor(W1[i,],W2[i,]) for any i < number of simulations

#Jump generation

jmin = 0 ;
jmax = 1/252 ;
jp = jmin + runif(n)*(jmax - jmin);
jt = round(t*jp*252);
ser = 1:t ;
jump=matrix(0,nrow=t,ncol=n)
for (i in 1:t){
  jump[i,] = ser[i];
}

for (j in 1:n){
  for (i in 1:t){
    if (jump[i,j] == jt[j]){
      jump[i,j] = 1;
    }else{
      jump[i,j] = 0;
    }
  }
}

j_price = t(matrix(rnorm(n*t,0,0.005^2),nrow=t,ncol=n)*jump) #rnorm(n*t,0,0.02^2)

j_vol = t(matrix(exp(rnorm(n*t,-5,1)),nrow=t,ncol=n)*jump)

#Heston Volatility Equation pre-jumps (calculated separately for ease of IV calculation)

for (i in 1:(t-1)){
  S2[,i+1] <- S2[,i] +phi*(v - S2[,i])*del + eta*sqrt(S2[,i])*W2[,i]*sdel
}

#Heston Volatility Equation including jumps

for (i in 0:(t-1)){
  S2_jumps[,i] <- S2[,i] + j_vol[,i]
}

#Heston Volatility Price Equation

for (i in 1:(t-1)){
  y[,i+1] <- y[,i] + (m - (S2[,i])/2 )*del + sqrt(S2[,i])*W1[,i]*sdel + j_price[,i]
}

#Generate data for simulated 10, and 5 minute intervals

y_10 = y[,seq(1,ncol(y), 600)]
y_5 = y[,seq(1,ncol(y), 300)]
y_2.5 = y[,seq(1,ncol(y), 150)] 
y_1 = y[,seq(1,ncol(y), 60)] 

#Volatility Calculation

RV = (rowSums(S2))*del

write.csv(y_5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps 5min.csv", row.names = TRUE)
write.csv(y_2.5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps 2.5min.csv", row.names = TRUE)
write.csv(y_1,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps 1min.csv", row.names = TRUE)
write.csv(RV,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps RV.csv", row.names = TRUE)

#Large Jumps Large Noise

blanks=matrix(0,nrow=n,ncol=(t-1))
y = cbind(matrix(m,nrow=n,ncol=1),blanks)
S2 = cbind(matrix(v,nrow=n,ncol=1),blanks)
S2_jumps = matrix(0,nrow=n,ncol=t)
W1 = matrix(rnorm(n*t),nrow=n,ncol=t)
n1 = matrix(rnorm(n*t),nrow=n,ncol=t)
W2 = (n1*srho +rho*W1)

#To test for the correlation of -0.5 between brownian motions W1 and W2
#use cor(W1[i,],W2[i,]) for any i < number of simulations

#Noise
noise = t(matrix(rnorm(n*t,0,0.001^2),nrow=t,ncol=n))
#Jump generation

jmin = 0 ;
jmax = 1/252 ;
jp = jmin + runif(n)*(jmax - jmin);
jt = round(t*jp*252);
ser = 1:t ;
jump=matrix(0,nrow=t,ncol=n)
for (i in 1:t){
  jump[i,] = ser[i];
}

for (j in 1:n){
  for (i in 1:t){
    if (jump[i,j] == jt[j]){
      jump[i,j] = 1;
    }else{
      jump[i,j] = 0;
    }
  }
}

j_price = t(matrix(rnorm(n*t,0,0.02^2),nrow=t,ncol=n)*jump)

j_vol = t(matrix(exp(rnorm(n*t,-5,1)),nrow=t,ncol=n)*jump)

#Heston Volatility Equation pre-jumps (calculated separately for ease of IV calculation)

for (i in 1:(t-1)){
  S2[,i+1] <- S2[,i] +phi*(v - S2[,i])*del + eta*sqrt(S2[,i])*W2[,i]*sdel
}

#Heston Volatility Equation including jumps

for (i in 0:(t-1)){
  S2_jumps[,i] <- S2[,i] + j_vol[,i]
}

#Heston Volatility Price Equation

for (i in 1:(t-1)){
  y[,i+1] <- y[,i] + (m - (S2[,i])/2 )*del + sqrt(S2[,i])*W1[,i]*sdel + j_price[,i]  + noise[,i]
}

#Generate data for simulated 10, and 5 minute intervals

y_10 = y[,seq(1,ncol(y), 600)]
y_5 = y[,seq(1,ncol(y), 300)]
y_2.5 = y[,seq(1,ncol(y), 150)] 
y_1 = y[,seq(1,ncol(y), 60)]

#Volatility Calculation

RV = (rowSums(S2))*del

write.csv(y_5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps large noise 5min.csv", row.names = TRUE)
write.csv(y_2.5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps large noise 2.5min.csv", row.names = TRUE)
write.csv(y_1,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps large noise 1min.csv", row.names = TRUE)
write.csv(RV,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data large jumps large noise RV.csv", row.names = TRUE)

#Small Jumps Small Noise

blanks=matrix(0,nrow=n,ncol=(t-1))
y = cbind(matrix(m,nrow=n,ncol=1),blanks)
S2 = cbind(matrix(v,nrow=n,ncol=1),blanks)
W1 = matrix(rnorm(n*t),nrow=n,ncol=t)
n1 = matrix(rnorm(n*t),nrow=n,ncol=t)
W2 = (n1*srho +rho*W1)

#To test for the correlation of -0.5 between brownian motions W1 and W2
#use cor(W1[i,],W2[i,]) for any i < number of simulations

#Noise
noise = t(matrix(rnorm(n*t,0,0.0005^2),nrow=t,ncol=n))

#Heston Volatility Equation pre-jumps (calculated separately for ease of IV calculation)

for (i in 1:(t-1)){
  S2[,i+1] <- S2[,i] +phi*(v - S2[,i])*del + eta*sqrt(S2[,i])*W2[,i]*sdel
}

#Heston Volatility Equation including jumps

for (i in 0:(t-1)){
  S2_jumps[,i] <- S2[,i] + j_vol[,i]
}

#Heston Volatility Price Equation

for (i in 1:(t-1)){
  y[,i+1] <- y[,i] + (m - (S2[,i])/2 )*del + sqrt(S2[,i])*W1[,i]*sdel + j_price[,i]  + noise[,i]
}

#Generate data for simulated 10, and 5 minute intervals

y_10 = y[,seq(1,ncol(y), 600)]
y_5 = y[,seq(1,ncol(y), 300)]
y_2.5 = y[,seq(1,ncol(y), 150)] 
y_1 = y[,seq(1,ncol(y), 60)]

#Volatility Calculation

RV = (rowSums(S2))*del
RV_new <- matrix(RV, nrow = 100, byrow = TRUE)

write.csv(y_5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps small noise 5min.csv", row.names = TRUE)
write.csv(y_2.5,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps small noise 2.5min.csv", row.names = TRUE)
write.csv(y_1,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps small noise 1min.csv", row.names = TRUE)
write.csv(RV,"C:\\Users\\19084\\My Backup Files\\Data\\MC Reg Data\\MC Data small jumps small noise RV.csv", row.names = TRUE)
