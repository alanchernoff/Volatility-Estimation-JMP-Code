setwd("/Users/acher/JMP/")
library(dplyr)
comp_list = c("AAPL","AXP","BA","CAT","CSCO", "CVX","DIS","HD","IBM",
              "INTC","JNJ","JPM","KO","MCD","MMM","MRK","MSFT","NKE",
              "PFE","UNH","UTX","VZ","WMT","XOM") 
start_date ="2012-01-03" #"2017-01-03" 
end_date="2016-12-30" #"2019-12-31" 
years = "2012-2016" #"2017-2019" 

for (i in 1:length(comp_list)){
  ticker = comp_list[i]
  df=read.csv(paste("/Users/acher/JMP/",ticker,'_1sec_Clean.csv',sep=""),header=TRUE)
  start=which(df$DATE==start_date)
  stop=which(df$DATE==end_date)
  df=df[start:stop,]
  datedf=df[,1]
  
  #export 1sec data
  write.table(df,paste("C:\\Users\\acher\\JMP\\Data\\"
                       ,ticker,"1sec",years,"e.csv",sep=""),sep=",", row.names = FALSE,col.names = FALSE,)
  #export close
  close=df[,ncol(df)]
  write.table(close,paste("C:\\Users\\acher\\JMP\\Data\\"
                          ,ticker,"close",years,"e.csv",sep=""),sep=",", row.names = FALSE,col.names = FALSE,)
  
  
  #1min
  minutedatadf=df[,grepl("00$",names(df))]
  write.table(minutedatadf,paste("C:\\Users\\acher\\JMP\\Data\\"
                                 ,ticker,"1min",years,"e.csv",sep=""),sep=",", row.names = FALSE,col.names = FALSE,)
  
  #5min
  df5min=minutedatadf[,!grepl("1.00$",names(minutedatadf))]
  df5min=df5min[,!grepl("2.00$",names(df5min))]
  df5min=df5min[,!grepl("3.00$",names(df5min))]
  df5min=df5min[,!grepl("4.00$",names(df5min))]
  df5min=df5min[,!grepl("6.00$",names(df5min))]
  df5min=df5min[,!grepl("7.00$",names(df5min))]
  df5min=df5min[,!grepl("8.00$",names(df5min))]
  df5min=df5min[,!grepl("9.00$",names(df5min))]
  #newdf=bind_cols(datedf,df5min)
  write.table(df5min,paste("C:\\Users\\acher\\JMP\\Data\\"
                           ,ticker,"5min",years,"e.csv",sep=""),sep=",", row.names = FALSE,col.names = FALSE,)
  
  #2.5min
  tensecdatadf=df[,grepl("0$",names(df))]
  df2halfmin=tensecdatadf[,!grepl("10$",names(tensecdatadf))]
  df2halfmin=df2halfmin[,!grepl("20$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("40$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("50$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("1.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("2.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("3.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("4.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("6.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("7.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("8.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("9.00$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("0.30$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("1.30$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("3.30$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("4.30$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("5.30$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("6.30$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("8.30$",names(df2halfmin))]
  df2halfmin=df2halfmin[,!grepl("9.30$",names(df2halfmin))]
  write.table(df2halfmin,paste("C:\\Users\\acher\\JMP\\Data\\"
                               ,ticker,"2.5min",years,"e.csv",sep=""),sep=",", row.names = FALSE,col.names = FALSE,)
  
}
