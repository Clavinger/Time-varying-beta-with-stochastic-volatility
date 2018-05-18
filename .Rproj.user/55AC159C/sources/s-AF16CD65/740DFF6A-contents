#skypt tworzy baze danych do badania

#ladowanie wymaganych pakietow
library(quantmod)
library(PerformanceAnalytics)
library(RCurl)


#rm(list = ls())
setwd("C:/Users/user/Dropbox/phd/Skrypty do R/Leverage-effect/Dane")

#funkcja do siagnia danych z stooq.pl
getStooqData <- function(asset_code,rodzaj) {
  w <- getURL(paste("https://stooq.pl/q/d/l/?s=",asset_code,"&i=",rodzaj,sep=""), 
              ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
  write(w,file="data_temp.csv")
  stooq_data <- read.csv("data_temp.csv")
  
  stooq_data
}

#wig
x=getStooqData("wig","d") 
ceny = x[,5]
daty = as.Date(x[,1])
x.xts =xts(ceny, daty)
zwroty=CalculateReturns(x.xts, method="log")
zwroty= zwroty[-1,]*100
write.zoo(zwroty,file="wig_zwroty.csv" ,sep=',')


#wig20
x=getStooqData("wig20","d") 
ceny = x[,5]
daty = as.Date(x[,1])
x.xts =xts(ceny, daty)
zwroty=CalculateReturns(x.xts, method="log")
zwroty= zwroty[-1,]*100
write.zoo(zwroty,file="wig20_zwroty.csv" ,sep=',')