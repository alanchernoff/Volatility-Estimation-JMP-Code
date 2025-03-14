setwd("C:/Users/acher/JMP/data/")
library(dplyr)
library(TTR)
library(tidyverse)
library(matrixStats)
library(Cairo)
library(ggplot2)

#Data Parameters
comp_list = c("HD","IBM","aapl","msft") 
yr="2017-2019"
dates = paste0("dates", yr,".csv")
interval_list = c("1min","2.5min","5min")
datefile=read.csv(dates, header=FALSE, sep=",")
datefile = datefile[,1]

#Volatility Calculation Functions
rv_calc<- function (raw){
  mat=data.matrix(raw, rownames.force = NA)  
  mat=t(mat)
  #mat=log(mat)
  #calculate RV
  dif=diff(mat)
  RV=colSums((dif)^2)
  RV
}
Vol_calc<- function (raw){
  #Vol_calc <- function (df,interval,yr){
  #raw=read.csv(file=paste0(df,interval,yr,"e.csv"), header=FALSE, sep=",")
  mat=data.matrix(raw, rownames.force = NA)  
  mat=t(mat)
  #mat=log(mat)
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
  #calculate RK
  kfunTH <- function(x) (sin((((1 - x)^2) * pi) / 2))^2
  kfunB <- function(x) (1-x)
  kfun2O <- function(x) (1-2*x+x^2)
  kfunE <- function(x) (1-x^2)
  kfunC <- function(x) (1-3*x^2+2*x^2)
  H=3
  gam_1 <- colSums(dif[1:(nrow(mat)-2),] * dif[2:(nrow(mat)-1),])
  gam_2 <- colSums(dif[1:(nrow(mat)-3),] * dif[3:(nrow(mat)-1),])
  gam_3 <- colSums(dif[1:(nrow(mat)-4),] * dif[4:(nrow(mat)-1),])
  gam <- cbind(gam_1, gam_2, gam_3)
  j_values <- 1:3 
  RK_TK = RV + colSums(t(matrix(kfunTH((j_values -1)/H), nrow = ncol(mat), ncol = 3, byrow = TRUE) * (2 * gam[,j_values])))
  RK_B = RV + colSums(t(matrix(kfunB((j_values -1)/H), nrow = ncol(mat), ncol = 3, byrow = TRUE) * (2 * gam[,j_values])))
  RK_2O = RV + colSums(t(matrix(kfun2O((j_values -1)/H), nrow = ncol(mat), ncol = 3, byrow = TRUE) * (2 * gam[,j_values])))
  RK_E = RV + colSums(t(matrix(kfunE((j_values -1)/H), nrow = ncol(mat), ncol = 3, byrow = TRUE) * (2 * gam[,j_values])))
  RK_C = RV + colSums(t(matrix(kfunC((j_values -1)/H), nrow = ncol(mat), ncol = 3, byrow = TRUE) * (2 * gam[,j_values])))
  #TSRV
  n <- 5
  K <- (nrow(mat)-1)/(n+1)
  sub_samp<- array(, dim = c(n, ncol(mat), K))
  for (j in 1:K) {for (i in 1:n) {
    sub_samp[i,,j] <- (mat[ i * K + j,] - mat[ (i - 1) * K + j,])^2
  }  }
  sub_vols <- matrix(,K, ncol(mat))
  
  for (j in 1:ncol(mat)) {for (i in 1:K) {
    sub_vols[i, j] <- sum(sub_samp[1:n, j, i])
  }  }
  avg_vols <- 1/K * colSums(sub_vols)
  TSRV = (1-(n-K+1)/(K*n))^(-1)*(avg_vols - (n-K+1)/(K*n) * RV)
  #combine columns
  vol_all=data.frame(RV,BPV,TPV,TRV,RK_TK,RK_B,RK_2O,RK_E,RK_C,TSRV,MinRV,MedRV)
  mean_vol=rowMeans(vol_all)
  mean_vol
}

#Choose Colors for Charts
c2 <- "#000000"
c1 <- "#5D8AA8" #69b3a2 #rgb(0.2, 0.6, 0.9, 1)
time1 <- 2
coeff_list = c(0.000006,0.000015,0.00005,0.00005,
               0.0000045,0.000015,0.00005,0.00003,
               0.000006,0.000020,0.000005,0.000011)

#Generate Charts
for(j in 1:length(interval_list)){
  for(i in 1:length(comp_list)){
    comp = comp_list[i]
    coeff = coeff_list[((j-1)*4+i)]
    raw1=read.csv(file=paste0(comp,interval_list[j],yr,"e.csv"), header=FALSE, sep=",")
    rv=rv_calc(log(raw1))
    meanvol=Vol_calc(log(raw1))
    close=read.csv(file=paste0(comp,"close",yr,"e.csv"),header=FALSE)
    close = close[,1]
    
    merged1 = data.frame(datefile,close,rv)
    merged1$dates = as.Date(merged1$datefile, "%m/%d/%Y")
    
    merged2 = data.frame(datefile,close,meanvol)
    merged2$dates = as.Date(merged2$datefile, "%m/%d/%Y")
    
    fig1 = ggplot(merged1, aes(x=dates)) +
      geom_bar( aes(y=rv/coeff), stat="identity", size=.1, fill=c1, color=c1, alpha=.1) + 
      geom_line( aes(y=close),size=.1, color=c2 ) +
      scale_y_continuous(
        # Features of the first axis
        name = "Closing Price",
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~.*coeff, name="Volatility"),
        expand = c(0.005,0)
      ) + 
      scale_x_date(date_breaks="1 year",date_labels="%Y", expand = c(0.005,0)) +
      theme_classic()+
      theme(
        aspect.ratio = 1,
        axis.title.y = element_text(color = c2, size=13),
        axis.title.y.right = element_text(color = c1, size=13),
        axis.title.x=element_blank(),
        axis.line.y=element_blank(),
        axis.line.y.right=element_blank(),
        axis.line.x=element_blank(),
        panel.border=element_rect(color="black",fill=NA,size=0.4)
      )
    CairoPNG(file=paste0(comp,interval_list[j],"rv.png"),width = 3.5,height = 2.8,units="in",dpi=96) 
    print(fig1)
    Sys.sleep(time1)
    print(fig1)
    dev.off()
    
    fig2 = ggplot(merged2, aes(x=dates)) +
      geom_bar( aes(y=meanvol/coeff), stat="identity", size=.1, fill=c1, color=c1, alpha=.1) + 
      geom_line( aes(y=close),size=.1, color=c2 ) +
      scale_y_continuous(
        # Features of the first axis
        name = "Closing Price",
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~.*coeff, name="Volatility"),
        expand = c(0.005,0)
      ) + 
      scale_x_date(date_breaks="1 year",date_labels="%Y", expand = c(0.005,0)) +
      theme_classic()+
      theme(
        aspect.ratio = 1,
        axis.title.y = element_text(color = c2, size=13),
        axis.title.y.right = element_text(color = c1, size=13),
        axis.title.x=element_blank(),
        axis.line.y=element_blank(),
        axis.line.y.right=element_blank(),
        axis.line.x=element_blank(),
        panel.border=element_rect(color="black",fill=NA,size=0.4)
      )
    CairoPNG(file=paste0(comp,interval_list[j],"meanvol.png"),width = 3.5,height = 2.8,units="in",dpi=96) 
    print(fig2)
    Sys.sleep(time1)
    print(fig2)
    dev.off()
    
  }
}

#Generating a single chart for testing purposes
j=3
i=1

comp = comp_list[i]
coeff = coeff_list[((j-1)*4+i)]
raw1=read.csv(file=paste0(comp,interval_list[j],yr,"e.csv"), header=FALSE, sep=",")
rv=rv_calc(log(raw1))
#meanvol=Vol_calc(log(raw1))
close=read.csv(file=paste0(comp,"close",yr,"e.csv"),header=FALSE)
close = close[,1]

merged1 = data.frame(datefile,close,rv)
merged1$dates = as.Date(merged1$datefile, "%m/%d/%Y")

#merged2 = data.frame(datefile,close,meanvol)
#merged2$dates = as.Date(merged2$datefile, "%m/%d/%Y")

fig1 = ggplot(merged1, aes(x=dates)) +
  geom_bar( aes(y=rv/coeff), stat="identity", size=.1, fill=c1, color=c1, alpha=.1) + 
  geom_line( aes(y=close),size=.1, color=c2 ) +
  scale_y_continuous(
    # Features of the first axis
    name = "Closing Price",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*coeff, name="Volatility"),
    expand = c(0.005,0)
  ) + 
  scale_x_date(date_breaks="1 year",date_labels="%Y", expand = c(0.005,0)) +
  theme_classic()+
  theme(
    title = element_blank(),
    aspect.ratio =1,
    axis.title.y = element_text(color = c2, size=13),
    axis.title.y.right = element_text(color = c1, size=13),
    axis.title.x=element_blank(),
    axis.line.y=element_blank(),
    axis.line.y.right=element_blank(),
    axis.line.x=element_blank(),
    panel.border=element_rect(color="black",fill=NA,size=0.4)
  )
CairoPNG(file=paste0(comp,interval_list[j],"rv.png"),width = 3.5,height = 2.8,units="in",dpi=96) 
print(fig1)
fig1
#Sys.sleep(time1)
#print(fig1)
dev.off()
