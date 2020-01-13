###clean every thing###
rm(list=ls())  
ls() 

###import libraries###
library(raster)
library(rgdal)
library(maptools)
library(RColorBrewer)
library(plotrix)
library(sp)
library(lattice)
library(rgeos)
require(RgoogleMaps)
library(shapefiles)
library(xts)
require(plyr)
require(tseries)
require(ggplot2)
require(lme4)
require(scales)
require(rpart)
require(stargazer)
require(randomForest)
require(dplyr)
require(xtable)
require(glmnet)
require(e1071)
require(chron)
require(mgcv)
require(Metrics)
require(scam)
require(data.table)

setwd("./")

###mobility data####
#source('mobility_mat.R')
source('function_code_new.R')


direct<-"fever_onset_files" #fever_onset_files #reporting_date_files

### code ###
town.1<-read.csv(paste("data/",direct,"/ISB1.csv",sep=""))
town.2<-read.csv(paste("data/",direct,"/ISB2.csv",sep=""))
town.3<-read.csv(paste("data/",direct,"/ISB3.csv",sep=""))
town.4<-read.csv(paste("data/",direct,"/ISB4.csv",sep=""))
town.5<-read.csv(paste("data/",direct,"/RK1.csv",sep=""))
town.6<-read.csv(paste("data/",direct,"/RWP1.csv",sep=""))
town.7<-read.csv(paste("data/",direct,"/RWP2.csv",sep=""))
town.8<-read.csv(paste("data/",direct,"/RWP3.csv",sep=""))
town.9<-read.csv(paste("data/",direct,"/RWP4.csv",sep=""))
town.10<-read.csv(paste("data/",direct,"/RWP5.csv",sep=""))
town.11<-read.csv(paste("data/",direct,"/RWP6.csv",sep=""))
town.12<-read.csv(paste("data/",direct,"/RWP7.csv",sep=""))
town.13<-read.csv(paste("data/",direct,"/RWP8.csv",sep=""))
town.14<-read.csv(paste("data/",direct,"/W1.csv",sep=""))

town.pop<-read.csv("data/rwp_population_area_towns.csv")
weather.rwp<-read.csv("data/weather_rwp.csv")


#town.pop<-adjust.pop(town.pop)

weather<-data.frame(date=weather.rwp$date, max_temperature=weather.rwp$max_temperature, min_temperature=weather.rwp$min_temperature, min_rain_fall=weather.rwp$min_rain_fall)
weather$max_temp<-weather$max_temperature
weather$max_temperature<-(weather$max_temp+weather$min_temperature)/2
weather$min_rain_fall[weather$min_rain_fall>=1]<-14
weather$min_rain_fall[weather$min_rain_fall<1]<-0

weather<-aggregate.bi.weekly(weather,mean)

#models_repl<-vector("list",100)
orig_df_repl<-list()
sel_df_repl<-list()


repl<-1

while (repl <=100){
  #lag=0
  report.rate<-100/100
  
  t1<-calc.feature(town.1,1,weather,town.pop,report.rate)
  t2<-calc.feature(town.2,2,weather,town.pop,report.rate)
  t3<-calc.feature(town.3,3,weather,town.pop,report.rate)
  t4<-calc.feature(town.4,4,weather,town.pop,report.rate)
  t5<-calc.feature(town.5,5,weather,town.pop,report.rate)
  t6<-calc.feature(town.6,6,weather,town.pop,report.rate)
  t7<-calc.feature(town.7,7,weather,town.pop,report.rate)
  t8<-calc.feature(town.8,8,weather,town.pop,report.rate)
  t9<-calc.feature(town.9,9,weather,town.pop,report.rate)
  t10<-calc.feature(town.10,10,weather,town.pop,report.rate)
  t11<-calc.feature(town.11,11,weather,town.pop,report.rate)
  t12<-calc.feature(town.12,12,weather,town.pop,report.rate)
  t13<-calc.feature(town.13,13,weather,town.pop,report.rate)
  t14<-calc.feature(town.14,14,weather,town.pop,report.rate)
  
  
  fin.data<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14)
  #fin.data<-rbind(t6,t7,t8,t11)
  
  fin.data.o<-fin.data
  fin.data<- remove.bad(fin.data)
  
  fin.data<-seperate.xi(fin.data)
  
  #fin.data<-fin.data[fin.data$lag1.irs==0,]
  
  #fin.data<-fin.data[(years(fin.data$date)==2013 | years(fin.data$date)==2016),]
  #fin.data<-subset(fin.data,as.Date(date,format="%Y-%m-%d")<as.Date('2017-06-30'))
  
  print(repl)
  
  orig_df_repl[[repl]]<-fin.data.o
  sel_df_repl[[repl]]<-fin.data
  
  repl<-repl+1
  
}


fin.data<-rbindlist(sel_df_repl)

#save(fin.data,file='rwp_single_replication_models_dfs_all_years_norm_area.Rdata')

lg=0

if (lg==1){
  pred_df=data.frame(y=fin.data$y,x=fin.data$x, y2=fin.data$y2, y3=fin.data$y3, tot.IoN=fin.data$tot.IoN, log_S=fin.data$log_S, lag.Fogging=fin.data$lag1.Fogging,lag.irs=fin.data$lag1.irs,lag.Larviciding=fin.data$lag2.Larviciding,
                     lag.FishSeed=fin.data$lag2.FishSeed,lag.Taps=fin.data$lag2.Taps,lag.Tyres=fin.data$lag2.Tyres,lag.Dewatering=fin.data$lag2.Dewatering,
                     lag.max_temperature=fin.data$lag2.max_temperature,lag.min_rain_fall=fin.data$lag2.min_rain_fall,pop.density=fin.data$pop.density,
                     x1=fin.data$x1,x2=fin.data$x2,x3=fin.data$x3,x4=fin.data$x4,x5=fin.data$x5,x6=fin.data$x6,x7=fin.data$x7,x8=fin.data$x8,x9=fin.data$x9,x10=fin.data$x10,
                     x11=fin.data$x11,x12=fin.data$x12,x13=fin.data$x13,x14=fin.data$x14)
  
  print (paste0('lg',lg))
  
} else {
  pred_df=data.frame(y=fin.data$y,x=fin.data$x, y2=fin.data$y2, y3=fin.data$y3, tot.IoN=fin.data$tot.IoN, log_S=fin.data$log_S, lag.Fogging=fin.data$lag0.Fogging,lag.irs=fin.data$lag0.irs,lag.Larviciding=fin.data$lag1.Larviciding,
                     lag.FishSeed=fin.data$lag1.FishSeed,lag.Taps=fin.data$lag1.Taps,lag.Tyres=fin.data$lag1.Tyres,lag.Dewatering=fin.data$lag1.Dewatering,
                     lag.max_temperature=fin.data$lag1.max_temperature,lag.min_rain_fall=fin.data$lag1.min_rain_fall,pop.density=fin.data$pop.density,
                     x1=fin.data$x1,x2=fin.data$x2,x3=fin.data$x3,x4=fin.data$x4,x5=fin.data$x5,x6=fin.data$x6,x7=fin.data$x7,x8=fin.data$x8,x9=fin.data$x9,x10=fin.data$x10,
                     x11=fin.data$x11,x12=fin.data$x12,x13=fin.data$x13,x14=fin.data$x14)
  
  print (paste0('lg',lg))
  
}

#model1<-gam(y~x+s(lag.Fogging)+s(lag.irs)+s(lag.Larviciding)+ s(lag.FishSeed)+s(lag.Taps)+s(lag.Tyres)+s(lag.Dewatering)+s(lag.max_temperature)+s(lag.min_rain_fall)+s(pop.density),data=pred_df,method="ML")



#model1 <- scam(y3~x+tot.IoN+s(lag.irs,bs="mdcx")+s(lag.Larviciding,bs="mdcx")+s(lag.Fogging,bs="mdcx")+s(lag.Tyres,bs="mdcx")+s(lag.max_temperature,bs="micv")+s(lag.min_rain_fall,bs="micv")+s(lag.Dewatering,bs="mdcx")+s(lag.FishSeed,bs="mdcx")+s(lag.Taps,bs="mdcx") + s(pop.density,bs="mpi") -1,data=pred_df)



#models_repl  <-try({ scam(y~x+s(lag.irs,bs="mdcx")+s(lag.Larviciding,bs="mdcx")+s(lag.Fogging,bs="mdcx")+s(lag.Tyres,bs="mdcx")+s(lag.max_temperature,bs="micv")+s(lag.min_rain_fall,bs="micv")+s(lag.Dewatering,bs="mdcx")+s(lag.FishSeed,bs="mdcx")+s(lag.Taps,bs="mdcx") + s(pop.density,bs="mpi") -1,data=pred_df)},silent = TRUE)

#model_single_x  <- scam(y~x+s(lag.irs,bs="mdcx")+s(lag.Larviciding,bs="mdcx")+s(lag.Fogging,bs="mdcx")+s(lag.Tyres,bs="mdcx")+s(lag.max_temperature,bs="micv")+s(lag.min_rain_fall,bs="micv")+s(lag.Dewatering,bs="mdcx")+s(lag.FishSeed,bs="mdcx")+s(lag.Taps,bs="mdcx") + s(pop.density,bs="mpi") -1,data=pred_df)
#model_single_x_no_contain<-scam(y~x+s(lag.max_temperature,bs="micv")+s(lag.min_rain_fall,bs="micv")+ s(pop.density,bs="mpi") -1,data=pred_df)

model_seperate_x<-scam(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+s(lag.irs,bs="mdcx")+s(lag.Larviciding,bs="mdcx")+s(lag.Fogging,bs="mdcx")+s(lag.Tyres,bs="mdcx")+s(lag.max_temperature,bs="micv")+s(lag.min_rain_fall,bs="micv")+s(lag.Dewatering,bs="mdcx")+s(lag.FishSeed,bs="mdcx")+s(lag.Taps,bs="mdcx") + s(pop.density,bs="mpi") -1,data=pred_df)

model_seperate_x_no_contain<-scam(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+s(lag.max_temperature,bs="micv")+s(lag.min_rain_fall,bs="micv")+ s(pop.density,bs="mpi") -1,data=pred_df)

#model1 <- scam(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+s(lag.irs,bs="mpd")+s(lag.Larviciding,bs="mpd")+s(lag.Fogging,bs="mpd")+s(lag.Tyres,bs="mpd")+s(lag.max_temperature,bs="mpi")+s(lag.min_rain_fall,bs="mpi")+s(lag.Dewatering,bs="mpd")+s(lag.FishSeed,bs="mpd")+s(lag.Taps,bs="mpd") + s(pop.density,bs="mpi") -1,data=pred_df)
#model1 <- scam(y2~x+log_S+s(lag.irs,bs="mdcx")+s(lag.Larviciding,bs="mdcx")+s(lag.Fogging,bs="mdcx")+s(lag.Tyres,bs="mdcx")+s(lag.max_temperature,bs="micv")+s(lag.min_rain_fall,bs="micv")+s(lag.Dewatering,bs="mdcx")+s(lag.FishSeed,bs="mdcx")+s(lag.Taps,bs="mdcx") + s(pop.density,bs="mpi") -1,data=pred_df)


#summary(model1)
save(model_seperate_x,model_seperate_x_no_contain,fin.data,orig_df_repl,sel_df_repl,file=paste0('lg',lg,'_single_replication_models_dfs_all_years.Rdata'))
#source('plot_data_repl_new.R')