###source functions###

##weekely##
aggregate.weekly <- function(data, method=sum) {
  data$date<-as.Date(data$date,"%m/%d/%Y")
  zz<-xts(data,order.by = data$date)
  aa<-apply.weekly(x=zz, FUN=function(v) {method(as.double(v[,3]))})
  zzs<-data.frame(date=index(aa))
  
  for (j in names(data)[-c(1,1)]) {
    zzs[,j]<-apply.weekly(x=zz, FUN=function(v) {method(as.double(v[,j]))})
  }
  
  zzs<-data.frame(date=zzs$date,data.matrix(zzs[,-1]))
  
  return(zzs)
}

aggregate.bi.weekly <- function(data,method=sum){
  data$date<-as.Date(data$date,"%m/%d/%Y")
  TimeStep<-"2 weeks"
  BINS<-seq(data$date[1],data$date[length(data$date)],by=TimeStep)
  
  zzs<-data.frame(date=BINS)
  vals<-matrix(0,nrow=length(BINS),ncol=(ncol(data)-1))
  
  for (i in 1:(length(BINS)-1)){
    tp<-data$date>=BINS[i] & data$date<BINS[i+1]
    
    for (j in 1:ncol(vals)){
      
      vals[i,j]<-method(data[tp,(j+1)])
      
    }
    
    
  }
  
  dff<-cbind(zzs,vals)
  names(dff)<-names(data)
  
  return(dff)
}



### RMSE ###
rmse <- function(x, y) { sqrt(mean((x-y)^2)) }



## spatial effect##
add.spatialeffect.old<-function(data,mobility){
  data2<-data
  data2$old_x<-data2$x
  data2$x<-0
  
  for (i in 1:length(data$date)){
    
    select<-which(data$date==data$date[i])
    temp.data<-data[select,]
    
    town<-data$town[i]
    
    tot.x<-as.vector(mobility[,town])*as.vector(temp.data$x)
    data2$x[i]<-sum(tot.x)
    
    
  }
  
  return(data2)
  
  
}

spatial.effect.cols<-function(data, mobility){
  
  data.o<-data
  
  len<-length(data$irs)
  #m.data<-data.frame(lag.irs.m=rep(0,len),lag.Larviciding.m=rep(0,len),lag.Dewatering.m=rep(0,len),lag.FishSeed.m=rep(0,len),lag.Fogging.m=rep(0,len),lag.Taps.m=rep(0,len),lag.Tyres.m=rep(0,len),lag.max_temperature.m=rep(0,len),lag.min_rain_fall.m=rep(0,len),pop.density.m=rep(0,len))
  #s.data<-data.frame(lag.irs.s=rep(0,len),lag.Larviciding.s=rep(0,len),lag.Dewatering.s=rep(0,len),lag.FishSeed.s=rep(0,len),lag.Fogging.s=rep(0,len),lag.Taps.s=rep(0,len),lag.Tyres.s=rep(0,len),lag.max_temperature.s=rep(0,len),lag.min_rain_fall.s=rep(0,len),pop.density.s=rep(0,len))
  #data<-data.frame(data.o,m.data,s.data)
  
  activ<- c("irs","Fogging","Larviciding","FishSeed","Taps","Tyres", "Dewatering","max_temperature","min_rain_fall")
  
  for (ac in activ ){
    
    data[paste("lag0.",ac,".m",sep="")]<-0 
    data[paste("lag1.",ac,".m",sep="")]<-0 
    data[paste("lag2.",ac,".m",sep="")]<-0 
    data[paste("lag3.",ac,".m",sep="")]<-0 
    
    data[paste("lag0.",ac,".s",sep="")]<-0 
    data[paste("lag1.",ac,".s",sep="")]<-0 
    data[paste("lag2.",ac,".s",sep="")]<-0 
    data[paste("lag3.",ac,".s",sep="")]<-0 
    
  }
  
  data["pop.density.m"]<-0
  data["pop.density.s"]<-0
  
  
  
  for (i in 1:length(data$date)){
    
    select<-which(data$date==data$date[i])
    temp.data<-data[select,]
    
    temp.data$I.m<-0
    
    for (k in 1:length(temp.data$I)){
      
      temp.data$I.m[k]<- sum(as.vector(temp.data$I)*as.vector(mobility[,k])) 
    }
    
    
    town<-data$town[i]
    
    pop.out<-data$S[i]*as.vector(mobility[town,])
    SoN<-as.vector(pop.out)/as.vector(temp.data$m.pop)
    SIoN<-as.vector(SoN)*temp.data$I.m
    
    
    for (j in 1:length(SIoN)){
      
      
      activ<- c("irs","Fogging","Larviciding","FishSeed","Taps","Tyres", "Dewatering","max_temperature","min_rain_fall")
      
      for (ac in activ ){
        
        data[i,paste("lag0.",ac,".m",sep="")]<-data[i,paste("lag0.",ac,".m",sep="")] + SIoN[j]*temp.data[j,paste("lag0.",ac,sep="")]
        data[i,paste("lag1.",ac,".m",sep="")]<-data[i,paste("lag1.",ac,".m",sep="")] + SIoN[j]*temp.data[j,paste("lag1.",ac,sep="")] 
        data[i,paste("lag2.",ac,".m",sep="")]<-data[i,paste("lag2.",ac,".m",sep="")] + SIoN[j]*temp.data[j,paste("lag2.",ac,sep="")] 
        data[i,paste("lag3.",ac,".m",sep="")]<-data[i,paste("lag3.",ac,".m",sep="")] + SIoN[j]*temp.data[j,paste("lag3.",ac,sep="")]
        
      }
      
      data$pop.density.m[i]<- data$pop.density.m[i] + SIoN[j]*temp.data$m.pop.density[j]
      
    }
    
    
  }
  
  fn<-0.1
  activ<- c("irs","Fogging","Larviciding","FishSeed","Taps","Tyres", "Dewatering","max_temperature","min_rain_fall")
  for (ac in activ ){
    
    data[paste("lag0.",ac,".s",sep="")]<-(1-fn)*as.vector(data[paste("lag0.",ac,".m",sep="")]) + fn*as.vector(data$x)*as.vector(data[paste("lag0.",ac,sep="")])
    data[paste("lag1.",ac,".s",sep="")]<-(1-fn)*as.vector(data[paste("lag1.",ac,".m",sep="")]) + fn*as.vector(data$x)*as.vector(data[paste("lag1.",ac,sep="")])
    data[paste("lag2.",ac,".s",sep="")]<-(1-fn)*as.vector(data[paste("lag2.",ac,".m",sep="")]) + fn*as.vector(data$x)*as.vector(data[paste("lag2.",ac,sep="")])
    data[paste("lag3.",ac,".s",sep="")]<-(1-fn)*as.vector(data[paste("lag3.",ac,".m",sep="")]) + fn*as.vector(data$x)*as.vector(data[paste("lag3.",ac,sep="")])
    
  }
  
  data["pop.density.s"]<-(1-fn)*as.numeric(data$pop.density.m) + fn*as.vector(data$x)*as.vector(data["pop.density"])
  
  
  return(data)
  
  
}


## residual effect##
residual.effect.inner<-function(activ,wks,decay){
  temp<-rep(0,length(activ))
  
  for (i in 1:length(activ)){
    
    val<-activ[i]
    
    for (j in 1:(wks-1)){
      
      if((i+j)<=length(activ)){
        
        temp[i+j]<- temp[i+j] + val*decay
        val<-val*decay
      }
    }
    
  }
  
  temp<-temp + activ
  
  return(temp)
  
  
}





residual.effect<-function(data){
  
  
  data$irs.o<- data$irs
  data$irs<-residual.effect.inner(data$irs.o,6,1)   #12 weeks divide by 2 , bi-weekly data
  
  data$Larviciding.o<- data$Larviciding
  data$Larviciding<-residual.effect.inner(data$Larviciding.o,3,1)
  
  data$FishSeed.o<- data$FishSeed
  data$FishSeed<-residual.effect.inner(data$FishSeed.o,200,1)   #guppyfish 9 years life.. growth rate not considered as adult each children
  
  data$Dewatering.o<- data$Dewatering
  data$Dewatering<-residual.effect.inner(data$Dewatering.o,2,0.5)
  
  data$Taps.o<- data$Taps
  data$Taps<-residual.effect.inner(data$Taps.o,2,1)
  
  data$Tyres.o<- data$Tyres
  data$Tyres<-residual.effect.inner(data$Tyres.o,2,1)
  
  data$min_rain_fall.o<- data$min_rain_fall
  #data$min_rain_fall<-residual.effect.inner(data$min_rain_fall.o,4,1)
  
  
  
  
  return(data)
}

norm.area<-function(data,area){
  data[,3:15]=data[,3:15]/area
  
  return(data)
  
}

## reporting rate ##

reporting_rate_construct_birth_rate_method<- function(data,pop,report.rate){
  orig_data<-data
  birth.r.y<-(22.5)/1000   # rate per year
  birth.r.w<- ((1+birth.r.y)^(1/(52/2))-1) # rate per week # updated with bi-week
  data$pop[1]<-pop
  data$tot_rep_cases[1]<-0
  data$birth<-0
  data$tot_birth<-0
  
  
  for (i in 2:nrow(data)){
    data$birth[i]<-birth.r.w*data$pop[i-1]
    data$pop[i]<-data$pop[i-1]+data$birth[i]
    data$tot_birth[i]<-sum(data$birth[1:i])
    data$tot_rep_cases[i]<-sum(data$confirmed[1:i])
    
    
  }
  
  mdl<-glm(tot_birth~tot_rep_cases,data = data)
  
  TSs<-data$confirmed
  RR<-as.numeric(mdl$coefficients)[2]
  #RR<-report.rate
  
  mySSTSM = smooth.spline(1:length(TSs), TSs)$y
  mySSTSMSum = mySSTSM * RR
  mySSTSMSum[mySSTSMSum<0] = 0
  PoisConstruct = rpois(length(mySSTSM),mySSTSMSum)
  
  orig_data$constructed_confirmed<-PoisConstruct
  return(orig_data)
  
  
}

reporting_rate_construct_actual<- function(data,pop,report.rate){
  orig_data<-data
  
  TSs<-data$confirmed
  
  mySSTSM = smooth.spline(1:length(TSs), TSs)$y
  mySSTSMSum = mySSTSM*(1/0.27)*2
  mySSTSMSum[mySSTSMSum<0] = 0
  PoisConstruct = rpois(length(mySSTSM),mySSTSMSum)
  
  #PoisConstruct =PoisConstruct*(1/0.27)
  #PoisConstruct = rpois(length(PoisConstruct),PoisConstruct)
  
  #PoisConstruct =PoisConstruct*(2)
  #PoisConstruct = rpois(length(PoisConstruct),PoisConstruct)
  
  orig_data$constructed_confirmed<-PoisConstruct
  return(orig_data)
  
  
}

## SIR ##
SIR.cols<- function(data, pop,area,town.num,cases_2011) {
  
  #reporting_rate<-report.rate
  #immunity<-0.1
  immune_pop<-cases_2011[town.num]*(1/0.27)*2
  #print(immune_pop)
  
  birth.r.y<-(22.5)/1000   # rate per year
  birth.r.w<- ((1+birth.r.y)^(1/(52/2))-1) # rate per week # updated with bi-week
  
  death.r.y<-(-6.49)/1000
  death.r.w<-((1+death.r.y)^(1/(52/2))-1)  # rate per week # updated with bi-week
  
  
  
  data$S<-0
  data$I<-data$constructed_confirmed
  #data$I<-(1/reporting_rate)*data$I
  data$R<-0
  data$pop<-0
  #data$m.pop<-0
  tot.pop<-pop
  data$tot.I<-0
  
  #data$m.pop[1]<-sum(as.vector(tot.pop)*as.vector(mobility[,town.num]))
  
  data$pop[1]<-pop[town.num]
  
  #data$R[1]<-immunity*pop[town.num]
  data$R[1]<-immune_pop
  
  data$S[1]<-data$pop[1]-data$R[1]-as.numeric(data$I[1])
  for (i in 2:length(data$I)){
    data$S[i]<-as.numeric(data$S[i-1])-as.numeric(data$I[i]) + birth.r.w* as.numeric(data$pop[i-1]) + death.r.w*as.numeric(data$S[i-1])
    data$R[i]<-data$R[i-1]+data$I[i-1] + death.r.w*as.numeric(data$R[i-1])
    data$pop[i]<- data$S[i]+ data$I[i] + data$R[i]
    
    tot.pop<- tot.pop +  birth.r.w*tot.pop +death.r.w*tot.pop 
    #data$m.pop[i]<-sum(as.vector(tot.pop)*as.vector(mobility[,town.num]))
    data$tot.I<-sum(data$I[1:i])
  }
  
  data$pop.density<-data$pop/area[town.num]
  #data$m.pop.density<-data$m.pop/area[town.num]
  
  data2<-data[-length(data$irs),]
  #data2<-data2[-length(data2$irs),]   #2 
  
  data2$St<-data2$S
  data2$SoNt<-data2$S/data2$pop
  data2$It<-data2$I
  
  data<-data[-1,]
  #data<-data[-1,]
  data2$Itp1<-data$I
  
  #data3<-remove.bad(data2)
  data3<-data2
  
  #data3$y<-data3$Itp1
  #data3$x<-data3$It*data3$SoNt
  
  data3$y<-log(data3$Itp1) - log(data3$SoNt) 
  data3$x<-log(data3$It)
  
  data3$y2<-log(data3$Itp1) + log(data3$pop)
  data3$log_S<-log(data3$S)
  
  data3$y3<-log(data3$Itp1)
  data3$tot.IoN<--1*data$tot.I/data3$pop
  
  
  return(data3) 
  
}


## adjust population
adjust.pop<- function (data){
  yr <- 2012
  growth<- 2.65 # %age per year
  #growth.week<-1.0189
  
  for (i in 2010:yr){
    data$population<-data$population*(1+growth/100)
  }
  
  return(data)
}

## remove zeros ##
remove.bad<- function (data){
  
  Bad<-data$x==(-Inf) | data$y==(-Inf)
  #Bad<-which(as.vector(data$lag0.max_temperature.s)==0 & as.vector(data$y)==0)
  #Bad<-which(as.vector(data$lag0.max_temperature.s)==0)
  
  data2<-data[!Bad,]
  
  return(data2)
  
}

## lag ##
add.lag<-function(data){
  
  #lag<-0
  data<-data[-length(data$irs),]
  
  temp1.1<-data
  temp1.2<-data
  #temp1.3<-data
  temp2<-data
  
  
  for (i in 1:2){
    temp2<-temp2[-1,]   
  }
  
  
  for (i in 1:1){
    temp1.1<-temp1.1[-1,]   #lag + 1 activity
  }
  
  for (i in 1:1){
    temp1.1<-temp1.1[-length(temp1.1$irs),]   #lag + 1 activity
  }
  
  #for (i in 1:1){
  #  temp1.2<-temp1.2[-1,]   #lag + 1 activity
 # }
  
  for (i in 1:2){
    temp1.2<-temp1.2[-length(temp1.2$irs),]   #lag + 1 activity
  }
  
  
  #for (i in 1:3){
  #  temp1.3<-temp1.3[-length(temp1.3$irs),]   #lag + 1 activity
 # }
  
  
  
  activ<- c("irs","Fogging","Larviciding","FishSeed","Taps","Tyres", "Dewatering","max_temperature","min_rain_fall")
  
  for (ac in activ ){
    
    temp2[paste("lag0.",ac,sep="")]<-temp2[ac] 
    temp2[paste("lag1.",ac,sep="")]<-temp1.1[ac] 
    temp2[paste("lag2.",ac,sep="")]<-temp1.2[ac] 
    #temp2[paste("lag3.",ac,sep="")]<-temp1.3[ac] 
    
  }
  
  #z=data.frame(t1$irs,t1$lag0.irs,t1$lag1.irs,t1$lag2.irs,t1$lag3.irs)
  
  return(temp2)
  
  
  
}

add.lag.normalize<-function(data){
  
  #lag<-0
  data<-data[-length(data$irs),]
  
  temp1.1<-data
  temp1.2<-data
  #temp1.3<-data
  temp2<-data
  
  
  for (i in 1:2){
    temp2<-temp2[-1,]   
  }
  
  
  for (i in 1:1){
    temp1.1<-temp1.1[-1,]   #lag + 1 activity
  }
  
  for (i in 1:1){
    temp1.1<-temp1.1[-length(temp1.1$irs),]   #lag + 1 activity
  }
  
  #for (i in 1:1){
  #  temp1.2<-temp1.2[-1,]   #lag + 1 activity
  # }
  
  for (i in 1:2){
    temp1.2<-temp1.2[-length(temp1.2$irs),]   #lag + 1 activity
  }
  
  
  #for (i in 1:3){
  #  temp1.3<-temp1.3[-length(temp1.3$irs),]   #lag + 1 activity
  # }
  
  
  
  activ<- c("irs","Fogging","Larviciding","FishSeed","Taps","Tyres", "Dewatering")
  
  for (ac in activ ){
    
    temp2[paste("lag0.",ac,sep="")]<-temp2[ac]/temp2$constructed_confirmed 
    temp2[paste("lag1.",ac,sep="")]<-temp1.1[ac]/temp2$constructed_confirmed 
    temp2[paste("lag2.",ac,sep="")]<-temp1.2[ac]/temp2$constructed_confirmed 
    #temp2[paste("lag3.",ac,sep="")]<-temp1.3[ac] 
    
  }
  
  clim<- c("max_temperature","min_rain_fall")
  
  for (ac in clim ){
    
    temp2[paste("lag0.",ac,sep="")]<-temp2[ac] 
    temp2[paste("lag1.",ac,sep="")]<-temp1.1[ac] 
    temp2[paste("lag2.",ac,sep="")]<-temp1.2[ac] 
    #temp2[paste("lag3.",ac,sep="")]<-temp1.3[ac] 
    
  }
  
  
  #z=data.frame(t1$irs,t1$lag0.irs,t1$lag1.irs,t1$lag2.irs,t1$lag3.irs)
  
  return(temp2)
  
  
  
}


calc.feature<-function(town.d,town.num,weather,town.pop,report.rate){
  
  t<-aggregate.bi.weekly(town.d)
  t<-reporting_rate_construct_actual(t,town.pop$population[town.num],report.rate)
  #print(report.rate)
  t$town<-town.num
  t<-cbind(t,weather[,-1])
  t<-residual.effect(t)
  
  #t<-norm.area(t,as.numeric(town.pop$area.sq.kms.[town.num]))
  t<-add.lag(t)
  #t<-add.lag.normalize(t)
  
  t<-SIR.cols(t,town.pop$population,town.pop$area.sq.kms.,town.num,town.pop$X2011.cases)
  
  return(t)
}

### mobility matrix ###
#zqq<- runif(100,min=0,max=0.1)
#mobility= matrix (  runif(100,min=0,max=0.5),nrow=10,ncol=10,byrow = TRUE)

#for (i in 1:10) {mobility[i,i]=1}


seperate.xi<-function(data){
  
  for (i in 1:14){
  
    data[paste("x",i,sep="")]<-0 
    data[data$town==i,paste("x",i,sep="")]<-1
    data[paste("x",i,sep="")]<-data[paste("x",i,sep="")]*data$x
  }
  
  return(data)
  
}

