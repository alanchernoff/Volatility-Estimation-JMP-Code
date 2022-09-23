library(dplyr)
library(matrixStats)
library(glmnet)

setwd("/Users/19084/My Backup Files/Data/MC Reg Data")

jump_ind =c("0","1","2","3","4","5")
interval_list = c("1min","2.5min","5min")
vol_names =c("RV","BPV","TPV","MinRV","MedRV","TRV")

Vol_calc <- function (raw){
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
  vol_all = bind_cols(RV,BPV,TPV,MinRV,MedRV,TRV)
  #mean_vol=rowMeans(vol_all)
  #vol_j = bind_cols(RV,BPV,TPV,MinRV,MedRV,TRV)
  #mean_jvol=rowMeans(vol_j)
  #result = data.frame(RV = RV,BPV=BPV,TPV=TPV,MinRV=MinRV,MedRV=MedRV,TRV=TRV,mean_vol=mean_vol,mean_jvol=mean_jvol)
  #result = data.frame(RV = RV,BPV=BPV,TPV=TPV,MinRV=MinRV,MedRV=MedRV,TRV=TRV)
  #result
  #vol_all
  vol_df=data.frame(RV,BPV,TPV,MinRV,MedRV,TRV)
  #vol_df=as.data.frame(vol_df)
  
  vol_df
} #used in MC_table & MC_table_clean

MC_table<-function(jump_ind,vol_names,interval_list,i){
  PT_RV=read.csv(file=paste0("T",jump_ind[i]," RV.csv"), header=TRUE, sep=",", row.names=1)
  
  jumpind = jump_ind[i] #i
  raw1=read.csv(file=paste0("T",jumpind," ",interval_list[1],".csv"), header=TRUE, sep=",", row.names=1)
  raw2=read.csv(file=paste0("T",jumpind," ",interval_list[2],".csv"), header=TRUE, sep=",", row.names=1)
  raw3=read.csv(file=paste0("T",jumpind," ",interval_list[3],".csv"), header=TRUE, sep=",", row.names=1)
  
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

model_reg <- function(vol_names,interval_list,vols){
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

df1 = MC_table(jump_ind,vol_names,interval_list,1)
ols1 <-model_reg(vol_names,interval_list,df1)
las1 <-model_lasso(vol_names,interval_list,df1)

df2 = MC_table(jump_ind,vol_names,interval_list,2)
ols2 <-model_reg(vol_names,interval_list,df2)
las2 <-model_lasso(vol_names,interval_list,df2)

df3 = MC_table(jump_ind,vol_names,interval_list,3)
ols3 <-model_reg(vol_names,interval_list,df3)
las3 <-model_lasso(vol_names,interval_list,df3)

df4 = MC_table(jump_ind,vol_names,interval_list,4)
ols4 <-model_reg(vol_names,interval_list,df4)
las4 <-model_lasso(vol_names,interval_list,df4)

df5 = MC_table(jump_ind,vol_names,interval_list,5)
ols5 <-model_reg(vol_names,interval_list,df5)
las5 <-model_lasso(vol_names,interval_list,df5)

df6 = MC_table(jump_ind,vol_names,interval_list,6)
ols6 <-model_reg(vol_names,interval_list,df6)
las6 <-model_lasso(vol_names,interval_list,df6)

setwd("/Users/19084/My Backup Files/Data/MC")
jump_ind_table =c("no jumps","T1 jumps","t2 jumps","no jumps noise","T1 jumps noise","T2 jumps noise")

MC_table_clean<-function(jump_ind,vol_names,interval_list,i){
  PT_RV=read.csv(file=paste0("MC Data ",jump_ind[i]," RV.csv"), header=TRUE, sep=",", row.names=1)
  
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
  
  vols$MeanVol_1min <-rowMeans(as.matrix(vols[,1:6]))
  vols$MeanVol_2.5min <-rowMeans(as.matrix(vols[,7:12]))
  vols$MeanVol_5min <-rowMeans(as.matrix(vols[,13:18]))
  vols$MeanVol_all <-rowMeans(as.matrix(vols[,1:18]))
  vols$MeanJVol_1min <-rowMeans(as.matrix(vols[,c(2,3,6)]))
  vols$MeanJVol_2.5min <-rowMeans(as.matrix(vols[,c(8,9,12)]))
  vols$MeanJVol_5min <-rowMeans(as.matrix(vols[,c(14,15,18)]))
  vols$MeanJVol_all <-rowMeans(as.matrix(vols[,c(2,3,6,8,9,12,14,15,18)]))
  
  
  vols$OLS_Vol_1 <- predict(ols1,vols)
  vols$LAS_Vol_1 <- predict(las1,as.matrix(vols[,1:18]))
  vols$OLS_Vol_2 <- predict(ols2,vols)
  vols$LAS_Vol_2 <- predict(las2,as.matrix(vols[,1:18]))
  vols$OLS_Vol_3 <- predict(ols3,vols)
  vols$LAS_Vol_3 <- predict(las3,as.matrix(vols[,1:18]))
  vols$OLS_Vol_4 <- predict(ols4,vols)
  vols$LAS_Vol_4 <- predict(las4,as.matrix(vols[,1:18]))
  vols$OLS_Vol_5 <- predict(ols5,vols)
  vols$LAS_Vol_5 <- predict(las5,as.matrix(vols[,1:18]))
  vols$OLS_Vol_6 <- predict(ols6,vols)
  vols$LAS_Vol_6 <- predict(las6,as.matrix(vols[,1:18]))
  
  vols['RV'] = PT_RV
  
  vols
  
} 

for (k in 1:6){

  #k = 1

  df=MC_table_clean(jump_ind_table,vol_names,interval_list,k) # edit 

  RMSE <- data.frame(matrix(ncol = ncol(df)-1, nrow = nrow(df)))
  colnames(RMSE) <- colnames(df)[1:ncol(df)-1]

  MAE <- data.frame(matrix(ncol = ncol(df)-1, nrow = nrow(df)))
  colnames(MAE) <- colnames(df)[1:ncol(df)-1]

  for(i in 1:ncol(RMSE)){
  RMSE[,i]=(df[,i]-df[,ncol(df)])^2
  MAE[,i]=abs(df[,i]-df[,ncol(df)])
}

  n=100
  RMSE_sub <- data.frame(matrix(ncol = ncol(df)-1, nrow = nrow(df)/n))
  MAE_sub <- data.frame(matrix(ncol = ncol(df)-1, nrow = nrow(df)/n))
  colnames(RMSE_sub) <- colnames(df)[1:ncol(df)-1]
  colnames(MAE_sub) <- colnames(df)[1:ncol(df)-1]

  for (j in 1:n){
  for (i in 1:ncol(RMSE_sub)){
    RMSE_sub[j,i] = sqrt(mean(RMSE[(1+(j-1)*n):(j*n),i]))
    MAE_sub[j,i] = mean(MAE[(1+(j-1)*n):(j*n),i])
  }
}

  RMSE_table <- data.frame(matrix(ncol = ncol(df)-1, nrow = 3))
  MAE_table <- data.frame(matrix(ncol = ncol(df)-1, nrow = 3))
  colnames(RMSE_table) <- colnames(df)[1:ncol(df)-1]
  colnames(MAE_table) <- colnames(df)[1:ncol(df)-1]

  for (i in 1:ncol(RMSE_table)){
    RMSE_table[1,i] <- format(round(mean(RMSE_sub[,i]), 4), nsmall = 4)
    RMSE_table[2,i] <- toString(format(round(quantile(RMSE_sub[,i], c(.10,.90)), 4), nsmall = 4))
    RMSE_table[3,i] <- toString(format(round(quantile(RMSE_sub[,i], c(.05,.95)), 4), nsmall = 4))
    MAE_table[1,i] <- format(round(mean(MAE_sub[,i]), 4), nsmall = 4)
    MAE_table[2,i] <- toString(format(round(quantile(MAE_sub[,i], c(.10,.90)), 4), nsmall = 4))
    MAE_table[3,i] <- toString(format(round(quantile(MAE_sub[,i], c(.05,.95)), 4), nsmall = 4))
  }
  

  write.csv(RMSE_table,paste0("C:\\Users\\19084\\My Backup Files\\Data\\RMSE_table_",toString(k),".csv"), row.names = TRUE)
  write.csv(MAE_table,paste0("C:\\Users\\19084\\My Backup Files\\Data\\MAE_table_",toString(k),".csv"), row.names = TRUE)
}
