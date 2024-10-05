setwd("/Users/acher/JMP/")
comp_list = c("AAPL","AXP","BA","CAT","CSCO", "CVX","DIS","HD","IBM",
              "INTC","JNJ","JPM","KO","MCD","MMM","MRK","MSFT","NKE",
              "PFE","UNH","UTX","VZ","WMT","XOM") 

yr="2017-2019"
yr_pre="2012-2016"
Vol_calc <- function (raw){
  mat=data.matrix(raw, rownames.force = NA) 
  mat=t(mat)
  #mat=log(mat)
  #calculate RV
  dif=diff(mat)
  RV=colSums((dif)^2)

  vol_df=data.frame(RV)
  #vol_df=as.data.frame(vol_df)

  vol_df
} #used in MC_table & MC_table_clean
for(i in 1:length(comp_list)){
  comp=comp_list[i]
  raw1_pre=read.csv(file=paste0(comp,"1sec",yr_pre,"e.csv"), row.names = 1, header=FALSE, sep=",")

  vols1=Vol_calc(log(raw1_pre))
  write.csv(vols1,file=paste0("C:\\Users\\19084\\My Backup Files\\Data\\Data\\pre",comp,"RV.csv"), row.names = TRUE)
