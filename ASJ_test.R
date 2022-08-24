setwd("/Users/19084/My Backup Files/Data/Data")

jump_calc <- function (df,interval,yr){
dfraw=read.csv(file=paste0(df,interval,yr,"e.csv"), header=FALSE, sep=",")
dfmat=data.matrix(dfraw, rownames.force = NA)
dfmat=t(dfmat)
dfmat=log(dfmat)
dfdiff=diff(dfmat)
delta=1/nrow(dfmat)
alpha=matrix(0, (nrow(dfmat)-1), ncol(dfmat))
for (j in 1: ncol(dfmat)){for (i in 1: (nrow(dfmat)-1)){if (abs(dfdiff[i,j]) <= sqrt(delta)){alpha[i,j]=abs((dfdiff[i,j]))^2} else {alpha[i,j]=0}}} 
alph_fin=5*sqrt(colSums(alpha)) 
omega=0.47 
BPD = colSums((abs(dfdiff))^4) 
kfreq = (nrow(dfmat)+1)/2 
data_10 = matrix(0, kfreq, ncol(dfmat)) 
for (j in 1: ncol(dfmat)){for (i in 1:kfreq){data_10[i,j] = dfmat[(i-1)*2+1,j]}} 
BPK = colSums((abs(diff(data_10)))^4) 
SPK = BPK/BPD 
trun_4 = matrix(0,nrow(dfmat)-1,ncol(dfmat)) 
for (j in 1: ncol(dfmat)){for (i in 1:(nrow(dfmat)-1)){ 
  if (abs(dfdiff[i,j]) <= alph_fin[j]*delta^omega){trun_4[i,j] = abs((dfdiff[i,j]))^4}
  else {trun_4[i,j] = 0}}} 
mp = pi^(-0.5)*4*gamma(5/2) 
AP = (delta^(-1)/mp)*colSums(trun_4) 
trun_8 = matrix(0, (nrow(dfmat)-1), ncol(dfmat)) 
for (j in 1: ncol(dfmat)){for (i in 1: (nrow(dfmat)-1)){ 
  if (abs(dfdiff[i,j]) <= alph_fin[j]*delta^omega){trun_8[i,j] = abs((dfdiff[i,j]))^8} 
  else {trun_8[i,j] = 0}}} 
mp_8 = pi^(-0.5)*16*gamma(9/2) 
AP_8 = (delta^(-3)/mp_8)*colSums(trun_8) 
Var = (delta* AP_8*160)/(3*AP^2) 
ASJ = (2 - SPK)/sqrt(Var) 
normASJ=pnorm(ASJ)
jumpday = ifelse(normASJ > 0.95, 1, 0)
sum(jumpday/length(jumpday))
}

jump_calc(comp_list[1],interval_list[2],yr[2])
  

comp_list = c("HD","IBM","aapl","msft") 
interval_list = c("1min","2.5min","5min")
yr= c("2012-2016","2017-2019","2012-2019")

df <- data.frame(matrix(ncol = length(yr), nrow = (length(interval_list)*length(comp_list))))
colnames(df) <- yr
comp1 <- lapply(comp_list, function(x) paste(x,interval_list[1], sep="_"))
comp2 <- lapply(comp_list, function(x) paste(x,interval_list[2], sep="_"))
comp3 <- lapply(comp_list, function(x) paste(x,interval_list[3], sep="_"))
comp = do.call(c, list(comp1,comp2,comp3))
comp = as.character(comp)
rownames(df) <- comp


for (i in 1:ncol(df)){
  for (j in 1:length(interval_list)){
    for (k in 1:length(comp_list))
    df[(j-1)*length(comp_list)+k,i] <- jump_calc(comp_list[k],interval_list[j],yr[i])
  }
}
