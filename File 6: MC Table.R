library(dplyr)
library(matrixStats)
library(glmnet)

jump_ind =c("large jumps","small jumps","large jumps noise","small jumps noise")
interval_list = c("1min","2.5min","5min")
vol_names =c("RV","BPV","TPV","TRV","RK_TK","RK_B","RK_2O","RK_E","RK_C","TSRV","MinRV","MedRV")

Vol_calc <- function (raw){
  mat=data.matrix(raw, rownames.force = NA) 
  mat=t(mat)
  #calculate RV
  dif=diff(mat)
  RV=colSums((dif)^2)
  #calculate BPV
  dif1=dif[1:(nrow(mat)-2),] 
  dif2=dif[2:(nrow(mat)-1),] 
  BPV=(sqrt(2)/sqrt(pi))^(-2)*colSums(abs(dif1)*abs(dif2))
  #calculate TPV
  dif11=dif[1:(nrow(mat)-3),] 
  dif22=dif[2:(nrow(mat)-2),] 
  dif33=dif[3:(nrow(mat)-1),] 
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
  vol_df=data.frame(RV,BPV,TPV,TRV,RK_TK,RK_B,RK_2O,RK_E,RK_C,TSRV,MinRV,MedRV)
  vol_df
} #used in MC_table & MC_table_clean

MC_reg_estimation<-function(jump_ind,vol_names,interval_list,i){
  
  setwd("/Users/19084/My Backup Files/Data/MC Data")

  #read in pseudo true volatiltiy from larger data set
  PT_RV=read.csv(file=paste0("MC Data ",jump_ind[i]," RV.csv"), header=TRUE, sep=",", row.names=1)
  
  #read in sub-sampled data
  jumpind = jump_ind[i] #i
  raw1=read.csv(file=paste0("MC Data ",jump_ind[i]," ",interval_list[1],".csv"), header=TRUE, sep=",", row.names=1)
  raw2=read.csv(file=paste0("MC Data ",jump_ind[i]," ",interval_list[2],".csv"), header=TRUE, sep=",", row.names=1)
  raw3=read.csv(file=paste0("MC Data ",jump_ind[i]," ",interval_list[3],".csv"), header=TRUE, sep=",", row.names=1)
  
  vols1=Vol_calc(raw1)
  vols2=Vol_calc(raw2)
  vols3=Vol_calc(raw3)
  
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
  
  vol_names1 <- lapply(vol_names, function(x) paste(x,interval_list[1], sep="_"))
  vol_names2 <- lapply(vol_names, function(x) paste(x,interval_list[2], sep="_"))
  vol_names3 <- lapply(vol_names, function(x) paste(x,interval_list[3], sep="_"))
  vol_names_e = do.call(c, list(vol_names1,vol_names2,vol_names3))
  colnames(vols)<-vol_names_e
  
  vols['RV'] = PT_RV
  
  
  
  vols
  
}

model_reg_full <- function(vol_names,interval_list,vols){
  vol_names1 <- lapply(vol_names, function(x) paste(x,interval_list[1], sep="_"))
  vol_names2 <- lapply(vol_names, function(x) paste(x,interval_list[2], sep="_"))
  vol_names3 <- lapply(vol_names, function(x) paste(x,interval_list[3], sep="_"))
  vol_names_e = do.call(c, list(vol_names1,vol_names2,vol_names3))
  form = "RV~0"
  len = length(vol_names_e)
  for (i in 1:len){
    form =paste0(form,"+",vol_names_e[i])
  }
  model1 <-lm(eval(parse(text=form)),data=vols)
  model1
  
}

model_reg <- function(vol_names,interval_list,i,vols){
  vol_names <- lapply(vol_names, function(x) paste(x,interval_list[i], sep="_"))
  vol_names_e = do.call(c, list(vol_names))
  form = "RV~0"
  len = length(vol_names_e)
  for (i in 1:len){
    form =paste0(form,"+",vol_names_e[i])
  }
  model1 <-lm(eval(parse(text=form)),data=vols)
  model1
  
}

model_lasso <-function(vol_names,interval_list,vols){
  vol_names1 <- lapply(vol_names, function(x) paste(x,interval_list[1], sep="_"))
  vol_names2 <- lapply(vol_names, function(x) paste(x,interval_list[2], sep="_"))
  vol_names3 <- lapply(vol_names, function(x) paste(x,interval_list[3], sep="_"))
  vol_names_e = do.call(c, list(vol_names1,vol_names2,vol_names3))
  form = "c('"
  len = length(vol_names_e)
  for (i in 2:len-1){
    form =paste0(form,vol_names_e[i],"','")
  }
  form=paste0(form,vol_names_e[len],"')")
  
  y <- vols$RV
  x <- data.matrix(vols[,eval(parse(text=form))])
  
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x, y, alpha = 1)
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  best_lambda
  
  #produce plot of test MSE by lambda value
  #plot(cv_model) 
  
  best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda,intercept=FALSE)
  best_model  
}

MC_table_clean<-function(jump_ind,vol_names,interval_list,i){
  
  setwd("/Users/19084/My Backup Files/Data/MC Reg Data")
  
  #generate volatilities for OLS and LASSO estimation
  df = MC_reg_estimation(jump_ind,vol_names,interval_list,i)
  
  #estimate OLS and LASSO models 
  ols_full <-model_reg_full(vol_names,interval_list,df)
  ols1 <-model_reg(vol_names,interval_list,1,df)
  ols2.5 <-model_reg(vol_names,interval_list,2,df)
  ols5 <-model_reg(vol_names,interval_list,3,df)
  las <-model_lasso(vol_names,interval_list,df)
  
  setwd("/Users/19084/My Backup Files/Data/MC Data")
  
  #read in pseudo true volatiltiy from larger data set
  PT_RV=read.csv(file=paste0("MC Data ",jump_ind[i]," RV.csv"), header=TRUE, sep=",", row.names=1)
  
  #read in sub-sampled data
  jumpind = jump_ind[i] #i
  raw1=read.csv(file=paste0("MC Data ",jump_ind[i]," ",interval_list[1],".csv"), header=TRUE, sep=",", row.names=1)
  raw2=read.csv(file=paste0("MC Data ",jump_ind[i]," ",interval_list[2],".csv"), header=TRUE, sep=",", row.names=1)
  raw3=read.csv(file=paste0("MC Data ",jump_ind[i]," ",interval_list[3],".csv"), header=TRUE, sep=",", row.names=1)
  
  #calculate sub-sampled volatilities
  vols1=Vol_calc(raw1)
  vols2=Vol_calc(raw2)
  vols3=Vol_calc(raw3)
  
  #merge sub-sampled volatilities
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
  
  #rename sub-sampled volatilities
  vol_names1 <- lapply(vol_names, function(x) paste(x,interval_list[1], sep="_"))
  vol_names2 <- lapply(vol_names, function(x) paste(x,interval_list[2], sep="_"))
  vol_names3 <- lapply(vol_names, function(x) paste(x,interval_list[3], sep="_"))
  vol_names_e = do.call(c, list(vol_names1,vol_names2,vol_names3))
  colnames(vols)<-vol_names_e
  
  #calculaute mean volatilities
  vols$MeanVol_1min <-rowMeans(as.matrix(vols[,1:6]))
  vols$MeanVol_2.5min <-rowMeans(as.matrix(vols[,7:12]))
  vols$MeanVol_5min <-rowMeans(as.matrix(vols[,13:18]))
  vols$MeanVol_all <-rowMeans(as.matrix(vols[,1:18]))
  vols$MeanJVol_1min <-rowMeans(as.matrix(vols[,c(2,3,6)]))
  vols$MeanJVol_2.5min <-rowMeans(as.matrix(vols[,c(8,9,12)]))
  vols$MeanJVol_5min <-rowMeans(as.matrix(vols[,c(14,15,18)]))
  vols$MeanJVol_all <-rowMeans(as.matrix(vols[,c(2,3,6,8,9,12,14,15,18)]))
  
  #calculate OLS and LASSO volatilities
  vols$OLS_Vol_full <- predict(ols_full,vols)
  vols$OLS_Vol_1 <- predict(ols1,vols)
  vols$OLS_Vol_2.5 <- predict(ols2.5,vols)
  vols$OLS_Vol_5 <- predict(ols5,vols)
  vols$LASSO_Vol <- predict(las,as.matrix(vols[,1:18]))
  
  vols['RV'] = PT_RV
  
  vols
  
} 

for (k in 1:4){

  #k = 1

  df=MC_table_clean(jump_ind,vol_names,interval_list,k) # edit 

  RMSE <- data.frame(matrix(ncol = ncol(df)-1, nrow = nrow(df)))
  colnames(RMSE) <- colnames(df)[1:ncol(df)-1]


  for(i in 1:ncol(RMSE)){
  RMSE[,i]=(df[,i]-df[,ncol(df)])^2
}

  n=100
  RMSE_sub <- data.frame(matrix(ncol = ncol(df)-1, nrow = nrow(df)/n))
  colnames(RMSE_sub) <- colnames(df)[1:ncol(df)-1]

  for (j in 1:n){
  for (i in 1:ncol(RMSE_sub)){
    RMSE_sub[j,i] = sqrt(mean(RMSE[(1+(j-1)*n):(j*n),i]))
  }
  }
  
  RMSE_sub = RMSE_sub*10^6

  RMSE_table <- data.frame(matrix(ncol = ncol(df)-1, nrow = 2))
  colnames(RMSE_table) <- colnames(df)[1:ncol(df)-1]

  for (i in 1:ncol(RMSE_table)){
    RMSE_table[1,i] <- format(round(mean(RMSE_sub[,i]), 4), nsmall = 4)
    RMSE_table[2,i] <- toString(format(round(quantile(RMSE_sub[,i], c(.05,.95)), 4), nsmall = 4))
  }
  
  RMSE_table = t(RMSE_table)
  colnames(RMSE_table) <- c("Mean RMSE","5th and 95th percentiles")

  
  write.csv(RMSE_table,paste0("C:\\Users\\19084\\My Backup Files\\Data\\RMSE_table_",jump_ind[k],".csv"), row.names = TRUE)
} 
