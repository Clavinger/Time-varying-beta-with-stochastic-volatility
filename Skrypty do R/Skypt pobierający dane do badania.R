#skrypt pobiera nazwy spolek z WIG30 ze strony wikipedii, następnie ze stooq.pl pobiera dane 
#i zapisuje na dysku twardym

#pakiety
library(quantmod)
library(PerformanceAnalytics)
library(RCurl)
library("rvest")
library("tidyr")
library("dplyr")


#funkcja do siagnia danych z stooq.pl
getStooqData <- function(asset_code,rodzaj) {
  w <- getURL(paste("https://stooq.pl/q/d/l/?s=",asset_code,"&i=",rodzaj,sep=""), 
              ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
  write(w,file="data_temp.csv")
  stooq_data <- read.csv("data_temp.csv")
  
  stooq_data
}

 
#sciagnie nazw spolek z wikipedii
wikiPL <- "https://pl.wikipedia.org/wiki/WIG30"
webpage <- read_html(wikiPL,encoding = "CP-1250")
table_links <- html_nodes(webpage, '.wikitable')
tab <- html_table(table_links, fill=TRUE)
head(tab)
colnames(df)[1]
df <- tab[[1]]
tikery<-as.vector(as.matrix(select(df,  colnames(df)[2])))
nazwy.spolek<-as.vector(as.matrix(select(df,  colnames(df)[1] )))


#zapisanie tikerow do pliku
write.csv2(tikery,file="tikery.csv", quote = TRUE)
#zapisanie nazw spolek do pliku
write.csv2(nazwy.spolek,file="nazwy_spolek.csv", quote = TRUE)

####################################################################################
#sciaganie danych na dysk twardy
setwd("C:/Users/user/Documents/github/Time-varying-beta-with-stochastic-volatility/Dane")
for(i in 1:length(tikery)){
  x=getStooqData(tikery[i],"d") 
  ceny = x[,5]
  daty = as.Date(x[,1])
  x.xts =xts(ceny, daty)
  zwroty=CalculateReturns(x.xts, method="log")
  zwroty= zwroty[-1,]*100
  write.zoo(zwroty,file=paste(tikery[i],"_zwroty.csv",sep=""),sep=',')
  assign(paste(tikery[i],".zwroty",sep=""),value=zwroty)
  
}

x=getStooqData("wig","d") 
ceny = x[,5]
daty = as.Date(x[,1])
x.xts =xts(ceny, daty)
zwroty=CalculateReturns(x.xts, method="log")
zwroty= zwroty[-1,]*100
write.zoo(zwroty,file="wig_zwroty.csv" ,sep=',')


setwd("C:/Users/user/Documents/github/Time-varying-beta-with-stochastic-volatility")
####################################################################################




#czyszczenie pamieci podrecznej
#rm(list = ls())


#przyklad jak odczytywać dane
#y=read.csv.zoo("ALR_zwroty.csv",sep=',')
#y=as.xts(y)
#plot.xts(y)
