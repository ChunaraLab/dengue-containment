rm(list = ls())


library(PBSmapping)
library(IDSpatialStats)
library(surveillance)
library(ggplot2)
library(matrixStats)

setwd('./')

app_data<-read.csv("../data/activities_2012_now.csv",header = TRUE,stringsAsFactors = FALSE)
portal_data<-read.csv("../data/dengue_patients_2012_now.csv",header = TRUE)

xx<-cbind(X=app_data$lng,Y=app_data$lat)
attr(xx,"projection")<-"LL"
coords.UTM<-convUL(xx,km=FALSE)
app_data$X<-coords.UTM[,1]
app_data$Y<-coords.UTM[,2]

app_data$created_at<-as.Date(app_data$created_at,format = "%Y-%m-%d")
app_data$created.days<-app_data$created_at - as.Date("2012-01-01",format = "%Y-%m-%d")

#containment_data<-subset(app_data,tag!="irs_patient")
containment_data<-app_data
containment_data<-subset(containment_data,tag!="Awareness" & tag!="Adult_Mosquitoes" & tag!="Helpline" & tag!="mowing" & tag!="Ovi_Trap" & tag!="Irrigation")
patient_tag_data<-subset(app_data,tag=="irs_patient")

portal_data$Reporting.Date<-as.Date(portal_data$Reporting.Date,format = "%Y-%m-%d")
portal_data$Onset.Date.Fever<-as.Date(portal_data$Onset.Date.Fever,format = "%Y-%m-%d")

portal_data$Reporting.days<-portal_data$Reporting.Date - as.Date("2012-01-01",format = "%Y-%m-%d")
portal_data$Fever.days<-portal_data$Onset.Date.Fever - as.Date("2012-01-01",format = "%Y-%m-%d")

merged_data<-merge(portal_data,patient_tag_data,by.x="Patient.ID",by.y="patient_id",all=TRUE)

#merged_data.city<-merged_data[merged_data$district=='Lahore' | merged_data$district=='LAHORE'| merged_data$District=='Lahore' | merged_data$District=='LAHORE',]
#merged_data.city<-merged_data[merged_data$District=='Lahore' | merged_data$District=='LAHORE',]
merged_data.city<-subset(merged_data,district=='Rawalpindi' | district=='RAWALPINDI'| district=='ISLAMABAD' | District=='Rawalpindi' | District=='RAWALPINDI' | District=='ISLAMABAD' | district=='Lahore' | district=='LAHORE'| District=='Lahore' | District=='LAHORE')
#merged_data.city<-subset(merged_data,district=='Rawalpindi' | district=='RAWALPINDI'| district=='ISLAMABAD' | District=='Rawalpindi' | District=='RAWALPINDI' | District=='ISLAMABAD' | district=='Lahore' | district=='LAHORE'| District=='Lahore' | District=='LAHORE')
#merged_data.city<-merged_data
merged_data.city<-subset(merged_data.city,place!='workplace')
merged_data.city<-subset(merged_data.city,!is.na(lng) & !is.na(Fever.days))
merged_data.city<-subset(merged_data.city,!duplicated(Patient.ID))

merged_data.city<-subset(merged_data.city,Onset.Date.Fever>as.Date("2014-01-01",format = "%Y-%m-%d"))

#old.merged_data.city<-merged_data.city
#merged_data.city<-subset(merged_data.city,Y>3712500 & Y<3727500 & X>312500 & X<329000)


### simple TAU ###

# hag.dat<-cbind(as.numeric(merged_data.city$X),as.numeric(merged_data.city$Y),merged_data.city$Fever.days)
# 
# colnames(hag.dat)<-c("x","y","t")
# 
# hagg.func<-function(a,b,tlimit=14){
#   if(abs(a[3]-b[3])<=tlimit){rc=1}
#   else{rc=2}
#   return(rc)
# }
# 
# wind.h=500
# r.max.h<-seq(100,10000,100)
# r.min.h<-r.max.h-wind.h
# r.min.h[which(r.min.h<0)]<-0
# r.mid.h<-(r.max.h+r.min.h)/2
# 
# tau.hagg<-get.tau(hag.dat,hagg.func,r=r.max.h,r.low=r.min.h,comparison.type = "independent")
# tau.ci<-get.tau.ci(hag.dat,hagg.func,r=r.max.h,r.low=r.min.h,boot.iter=100,comparison.type = "independent")
# 
# xpoly<-c(r.mid.h,rev(r.mid.h))
# ypoly<-c(tau.ci[1,],rev(tau.ci[2,]))
# 
# pdf("rawalpindi_2014-2016_14days_500m_10km_new.pdf", width=5, height=5)
# par(mar=c(3,3,3,1))
# plot(r.mid.h, tau.hagg ,ylim=c(0.7,5), type="l",xlim=c(0,9750),yaxt="n",axes=F)
# axis(2,las=1,cex.axis=0.75,at=c(0.8,1.5,2,2.5,3,4,5))
# axis(1,las=1,cex.axis=0.75,at=c(0,500,1000,5000,10000))
# #axis(1,las=1,cex.axis=0.75,at=c(0,1000,3450,5000,10000))
# mtext(text='distance(m)',side=1,line=2)
# mtext(text='tau',side=2,line=2)
# abline(h=1,lty=2,col=1)
# abline(v=500,lty=2,col=1)
# #abline(v=3450,lty=2,col=1)
# polygon(xpoly,ypoly,col=rgb(t(col2rgb("dark green")/255),alpha=0.6),border=NA)
# lines(r.mid.h, tau.hagg,lty=1,col="dark green")
# dev.off()

### end ###

pst<-30  # activity days from past in window
fut<-0   # activity days from future in window
rad_m<-20 # radius for activity presence in meters 
tim_break<-30

patient_rad_m<-500
patient_pst<-30
patient_fut<-0

patient_data<- merged_data.city

figure_df<-matrix(0,ncol=length(unique(containment_data$tag)),nrow=nrow(patient_data))
activ_name<-sort(unique(containment_data$tag))
#activ_name<-c(activ_name,'irs_patient')
patient_epi_context_count<-rep(0,nrow(patient_data))
patient_epi_context_count2<-rep(0,nrow(patient_data))

containment_context_count<-rep(0,nrow(patient_data))

tp_df<-rep(0,nrow(patient_data))
tp_df_tag<-rep(0,nrow(patient_data))

for (i in 1:nrow(patient_data)){
  
  selected_cont_data<-subset(containment_data,created.days>=(patient_data$Fever.days[i]-pst) & created.days<=(patient_data$Fever.days[i]+fut))
  
  restricted_selected_data<-subset(selected_cont_data,   sqrt( (patient_data$X[i]-X)^2 +(patient_data$Y[i]-Y)^2 )<rad_m )
  
  
  tp_df[i]<-nrow(restricted_selected_data)
  tp_df_tag[i]<-length(unique(restricted_selected_data$tag))
  
  unq<-unique(restricted_selected_data$tag)
  
  if (length(unq>0)){
    
    for (j in 1:8){
      
      for (u in 1:length(unq)){
        
        if (unq[u]==activ_name[j]){
          
          figure_df[i,j]<-1
          
        }    
        
        
      }
    }
  }
  
  
  
  selected_patient_data<-subset(patient_data,Fever.days>=(patient_data$Fever.days[i]-patient_pst) & Fever.days<=(patient_data$Fever.days[i]+patient_fut))
  
  restricted_selected_patient_data<-subset(selected_patient_data,   sqrt( (patient_data$X[i]-X)^2 +(patient_data$Y[i]-Y)^2 )<patient_rad_m )
  
  patient_epi_context_count[i]<-nrow(restricted_selected_patient_data)-1
  patient_epi_context_count2[i]<-(nrow(restricted_selected_patient_data)-1)/nrow(selected_patient_data)
  
  sel_cont_context_data<-subset(containment_data,   sqrt( (patient_data$X[i]-X)^2 +(patient_data$Y[i]-Y)^2 )<patient_rad_m)
  
  containment_context_count[i]<-nrow(sel_cont_context_data)
  
    print(i)
  
}

# pdf("lahore_2013_contain_activ_p90+f14.pdf", width=5, height=5)
# hist(tp_df_tag,main="containment activ. in past 90 days + future 14 days in 500m radius",ylab="Total cases",xlab="Total Distinct Type of Activities",cex.main=0.75)
# dev.off()


figure_df_new<-data.frame(figure_df,stringsAsFactors = FALSE)
names(figure_df_new)<-activ_name

count_df<-rep(0,length(activ_name))
activ_type_patient<-rep(0,nrow(figure_df_new))

for (j in 1:8){
  
  for (i in 1:nrow(figure_df_new)){
    
    if(figure_df_new[i,j]==1 & sum(figure_df_new[i,])==1){
      
      count_df[j]<-count_df[j]+1
      activ_type_patient[i]<-j
    }
    
    if (sum(figure_df_new[i,])>1){
      activ_type_patient[i]<--1
      
    }
    
  }
  
  
}


patient_data$Fever.month<-as.numeric(strftime(patient_data$Onset.Date.Fever, "%m"))

#save(patient_data,activ_type_patient,pst,fut,rad_m,tim_break,file=paste("df_p",pst,"_f",fut,"_rad",rad_m,".Rdata",sep = ''))
#load(paste("df_p",pst,"_f",fut,"_rad",rad_m,".Rdata",sep = ''))
#load('df_pst_15.Rdata')
#tim_break<-14


labels<-list()



source('TauRatios-with-ci-future-matching.R')
#source('TauRatios-with-ci-future-matching_effect30_60_days.R')


# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# labels[[1]]<-as.numeric(patient_data$labels)
# 
# 
# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# slx<- patient_data$Y<3712500 | patient_data$Y>3727500 | patient_data$X<312500 | patient_data$X>329000
# patient_data$labels[slx]<--2
# labels[[2]]<-as.numeric(patient_data$labels)
# 
# 
# 
# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# slx<- patient_data$Y<3720000 | patient_data$Y>3724000 | patient_data$X<317500 | patient_data$X>322500
# patient_data$labels[slx]<--2
# labels[[3]]<-as.numeric(patient_data$labels)
# 
# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# slx<- patient_data$Y<3712500 | patient_data$Y>3727500 | patient_data$X<312500 | patient_data$X>329000
# patient_data$labels[!slx]<--2
# labels[[4]]<-as.numeric(patient_data$labels)

#activ_type_patient_orig<-activ_type_patient
#activ_type_patient<-sample(-1:8, length(activ_type_patient), replace=T)

#quant<-quantile(patient_epi_context_count)



# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# patient_data$labels[patient_data$Fever.month>9 & patient_data$labels!=0]<--2
# labels[[1]]<-as.numeric(patient_data$labels)
# 
# 
# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# patient_data$labels[patient_data$Fever.month<=9 & patient_data$labels!=0]<--2
# labels[[2]]<-as.numeric(patient_data$labels)
# 
# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# labels[[3]]<-as.numeric(patient_data$labels)


 patient_data$labels<-activ_type_patient
 patient_data$labels[patient_data$labels==-1]<-9
 labels[[1]]<-as.numeric(patient_data$labels)




# 
# 
# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# patient_data$labels[patient_epi_context_count<=quant[3] |patient_epi_context_count>quant[4]]<--2
# labels[[3]]<-as.numeric(patient_data$labels)
# 
# 
# patient_data$labels<-activ_type_patient
# patient_data$labels[patient_data$labels==-1]<-9
# patient_data$labels[patient_epi_context_count<=quant[4]]<--2
# labels[[4]]<-as.numeric(patient_data$labels)


matching_dist<-1000
fil_nam<-paste(matching_dist,"-season-split-matched-future_",tim_break,"_RWP-ISB-LHR-combined_tauA-tauB_ratio_p",pst,"_f",fut,"_rad",rad_m,".pdf",sep = '')
plt_titl<-paste("Season Split Matched (",matching_dist," m) Rawalpindi and Islamabad data",sep='')

source('quant-plot_code_single_new.R')



plot_df<-plot_df[-1,]
clrs<-c("red", "grey", "green","orange","purple","pink", "brown", "blue")

pdf(paste("combined2_tauA-tauB_ratio_p",pst,"_f",fut,".pdf",sep = ''), width=7, height=5)
ggplot()+geom_line(data=plot_df,aes(x=dist,y=mean,color=activ)) +geom_ribbon(data=plot_df,aes(ymax=upper_ci,ymin=lower_ci,x=dist,fill=activ),alpha=0.3) + theme_bw()+ scale_colour_manual(values = clrs) + scale_fill_manual(values = clrs) +scale_y_log10(limits = c(0.2,2),breaks =c(0.2,0.5,1,2)) + ylab("ratio Tau(A)/Tau(B)") + xlab("distance(meters)") +ggtitle(paste("past ",pst," days and ",fut," future days)",sep = '')) +theme(plot.title = element_text(hjust = 0.5))
dev.off()
# pdf(paste("test-rawalpindi_containment_ratio_punjab-2014-2016-p",pst,"_f",fut,"_",rad_m,"m.pdf",sep = ''), width=10, height=7)
# ggplot()+ geom_line(data=plot_df,aes(x=dist, y=tau_val,colour=activ),size=0.7) + theme_bw()  + ylab("ratio P(A)/P(B)") + xlab("distance(meters)") +ggtitle(paste("Rawalpindi 2014-2016 (",rad_m,"m and past ",pst," days and ",fut," future days)",sep = ''))+theme(plot.title = element_text(size=15,hjust = 0.5))
# dev.off()


#jpeg(paste(patient_pst," days_test-rawalpindi_tauA-tauB_ratio_punjab-2014-2016-p",pst,"_f",fut,"_",rad_m,"m.jpeg",sep = ''), width = 10, height = 7, units = 'in',  res = 300)
#print(
 # ggplot()+ geom_line(data=plot_df,aes(x=dist, y=tau_val,colour=activ),size=0.7) + theme_bw()  + ylab("ratio P(A)/P(B)") + xlab("distance(meters)") +ggtitle(paste("Rawalpindi 2014-16 (",rad_m,"m and past ",pst," days and ",fut," future days) EpiContext ",j,", based on patients in past ",patient_pst," days",sep = ''))+theme(plot.title = element_text(size=10,hjust = 0.5)) + scale_y_log10(limits = c(0.05,5),breaks =c(0.2,1,5))
#)
#dev.off()




old.merged_data.city2<-merged_data.city


pdf("selected_area.pdf",width=7, height=7)
ggplot()+ geom_point(data=old.merged_data.city2,aes(x=X, y=Y),size=0.5) + geom_rect(aes(xmin=312500, xmax=329000, ymin=3712500, ymax=3727500), fill='blue', alpha=0.4)+ theme_bw()  + ylab("Y") + xlab("X")+theme(plot.title = element_text(size=15,hjust = 0.5))
dev.off()

patient_data2<-patient_data[patient_data$labels!=-2,]
patient_data2$type<-''
patient_data2$type[patient_data2$labels>0]<-'containment'
patient_data2$type[patient_data2$labels==0]<-'control'
patient_data2$type[patient_data2$labels<0]<-'multiple'
patient_data2$type[patient_data2$labels==5]<-'multiple'

pdf("colored_selected_area.pdf",width=9, height=7)
ggplot()+ geom_point(data=patient_data2,aes(x=X, y=Y,color=type),size=0.5)+ geom_rect(aes(xmin=317500, xmax=322500, ymin=3720000, ymax=3724000), fill='grey', alpha=0.4) + theme_bw()  + ylab("Y") + xlab("X")+theme(plot.title = element_text(size=15,hjust = 0.5)) +guides(color = guide_legend(override.aes = list(size=5)))
dev.off()


require(rgdal)
require(maptools)
require(ggplot2)


fn='Rawalpindi.kml'
kml <- readOGR(fn)
slot(kml, "polygons") <- lapply(slot(kml, "polygons"), checkPolygonsHoles)
kml.points = fortify(kml, region="Description")

xx<-cbind(X=kml.points$long,Y=kml.points$lat)
attr(xx,"projection")<-"LL"
coords.UTM<-convUL(xx,km=FALSE)
kml.points$X<-coords.UTM[,1]
kml.points$Y<-coords.UTM[,2]

pdf("tstt.pdf",width=19, height=17)
ggplot()+ geom_polygon(data=kml.points, aes(x=X,y=Y,group=id),color='black')+ geom_point(data=old.merged_data.city2,aes(x=lng, y=lat),size=0.5,color='red')

