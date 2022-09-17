setwd("/Users/19084/My Backup Files/Data/Data") #double data folder
library(dplyr)
library(TTR)
library(tidyverse)

Vol_calc <- function (df,interval,yr){
raw=read.csv(file=paste0(df,interval,yr,"e.csv"), header=FALSE, sep=",") 
mat=data.matrix(raw, rownames.force = NA)
mat=t(mat)
mat=log(mat)
#calculate RV
dif=diff(mat)
RV=colSums((dif)^2)
#calculate BPV
dif1=dif[1:(nrow(mat)-2),1:ncol(mat)] 
dif2=dif[2:(nrow(mat)-1),1:ncol(mat)] 
BPV=(sqrt(2)/sqrt(pi))^(-2)*colSums(abs(dif1)*abs(dif2))
#calculate TPV
dif11=dif[1:(nrow(mat)-3), 1:ncol(mat)] 
dif22=dif[2:(nrow(mat)-2), 1:ncol(mat)] 
dif33=dif[3:(nrow(mat)-1), 1:ncol(mat)] 
cons=(((2^(1/3))*gamma(5/6))/gamma(1/2))^(-3) 
TPV=cons*colSums((abs(dif11)^(2/3))*(abs(dif22)^(2/3))*(abs(dif33)^(2/3)))
#calculate MinRV
minvec=matrix(0, (nrow(mat)-2), ncol(mat)) 
for (j in 1:ncol(mat)){for (i in 1:(nrow(mat)-2)){minvec[i,j]=(min( abs(dif1[i,j]), abs(dif2[i,j]) ))^2}} 
MinRV=(pi/(pi-2))*(nrow(mat)/(nrow(mat)-1))*colSums(minvec) 
#calculate MedRV
medvec=matrix(0, (nrow(mat)-3), ncol(mat)) 
for (j in 1:ncol(mat)){for (i in 1:(nrow(mat)-3)){medvec[i,j]=(median(c(abs(dif11[i,j]),abs(dif22[i,j]),abs(dif33[i,j]))))^2}} 
MedRV=(pi/(6-4*sqrt(3)+pi))*(nrow(mat)/(nrow(mat)-2))*colSums(medvec)
#calculate TRV
delta=1/nrow(mat) 
omega=0.47 
alpha=matrix(0, (nrow(mat)-1), ncol(mat)) 
for (j in 1: ncol(mat)){for (i in 1: (nrow(mat)-1)){if (abs(dif[i,j]) <= sqrt(delta)){alpha[i,j]=abs((dif[i,j]))^2} else {alpha[i,j]=0}}} 
alph_fin=5*sqrt(colSums(alpha)) 
trun=matrix(0, (nrow(mat)-1), ncol(mat))
for (j in 1: ncol(mat)){for (i in 1: (nrow(mat)-1)){if (abs(dif[i,j]) <= alph_fin[j]*delta^omega){trun[i,j]=abs((dif[i,j]))^2} else {trun[i,j]=0}}} 
TRV=colSums(trun) 
vol_all = bind_cols(RV,BPV,TPV,MinRV,MedRV,TRV)
mean_vol=rowMeans(vol_all)
#vol_3 = bind_cols(RV,BPV,TPV,MinRV,MedRV,TRV)
#mean_3vol=rowMeans(vol_3)
vol_j = bind_cols(RV,BPV,TPV,MinRV,MedRV,TRV)
mean_jvol=rowMeans(vol_j)
result = data.frame(RV = RV,BPV=BPV,TPV=TPV,MinRV=MinRV,MedRV=MedRV,TRV=TRV,mean_vol=mean_vol,mean_jvol=mean_jvol)
result
}

lagpad <- function(x) { return (c(x[(1+1) : length(x)], rep(NA, 1)))}

sharp_Calc = function(vols,rf3mo,close,sharps,j){
  for (i in 1:ncol(vols)){
    movavg5=SMA(vols[,i],5)
    movavg20=SMA(vols[,i],20)
    signal = ifelse(movavg5 < movavg20, 1, 0)
    returns=lagpad(ROC(close))*signal
    returns=na.omit(returns-rf3mo)
    #sharps[nrow(sharps)+1,]<- c((mean(returns))/sd(returns))
    sharps[j,i+1] <- c((mean(returns))/sd(returns))
    #names(sharps)[length(names(sharps))]<-names(vols)[i] 
  }
  sharps
}

convert_table <-function(sharps,vol_names,interval_list){
  vol_names1 <- lapply(vol_names, function(x) paste(x,interval_list[1], sep="_"))
  vol_names2 <- lapply(vol_names, function(x) paste(x,interval_list[2], sep="_"))
  vol_names3 <- lapply(vol_names, function(x) paste(x,interval_list[3], sep="_"))
  vol_names = do.call(c, list(vol_names1,vol_names2,vol_names3))
  vol_names = as.character(vol_names)
  names(vol_names) <- paste0("V",2:(length(vol_names)+1))
  sharps <- rename_with(.data = sharps, .cols = starts_with("V"), .fn = function(x){vol_names[x]})
  sharps
}

sharp_table<-function(comp_list,vol_names,interval_list,yr){
  sharps = data.frame(comp_list)
  for(i in 1:length(comp_list)){
    df = comp_list[1] #i
    #yr1 = "2012-2016"
    #yr2 = "2017-2019"
    
    rf3mo = read.csv(file=paste("3mo",yr,".csv",sep=""),header=FALSE)
    rf3mo=rf3mo[,1]
    rf3mo=rf3mo/360
    close=read.csv(file=paste(df,"close",yr,"e.csv",sep=""),header=FALSE) #e
    close = close[,1]
    
    vols1=Vol_calc(df,interval_list[1],yr)
    vols2=Vol_calc(df,interval_list[2],yr)
    vols3=Vol_calc(df,interval_list[3],yr)
    
    vols_t <- merge(vols1,vols2, by = 'row.names',all = TRUE)
    vols_t$Row.names = as.numeric(vols_t$Row.names)
    vols_t = vols_t[order(vols_t$Row.names),]
    rownames(vols_t) <- NULL
    vols_t <- subset(vols_t, select = -c(Row.names))
    
    vols <- merge(vols_t,vols3, by = 'row.names', all = TRUE)
    vols$Row.names = as.numeric(vols$Row.names)
    vols = vols[order(vols$Row.names),]
    rownames(vols) <- NULL
    vols <- subset(vols, select = -c(Row.names))
    
    sharps=sharp_Calc(vols,rf3mo,close,sharps,i)
  }
  
  sharps = convert_table(sharps,vol_names,interval_list)
  sharps
  
} 
  
comp_list = c("HD","IBM","aapl","msft") 
vol_names =c("RV","BPV","TPV","MinRV","MedRV","TRV","mean_vol","mean_jvol")
interval_list = c("1min","2.5min","5min")
yr="2017-2019"

sharps = sharp_table(comp_list,vol_names,interval_list,yr) #<-run this

View(sharps)