require(foreach)
require(doParallel)
require(matrixStats)

locs<-as.matrix(cbind(patient_data$X,patient_data$Y))
dates<-as.numeric(patient_data$Fever.days)


rearrange_activ<-c(3,4,6,1,7,8)
activ_name_fixed<-c("Dewatering","Fish Seeding","Fogging","Indoor Residual Spray", "irs_patient", "Larviciding", "Tap Fixing", "Tire Shredding","multiple")
#for (i in c(3,4,5)){
plot_df<-matrix('',nrow=(length(dbreaks.mid)-1),ncol=length(rearrange_activ))
plot_df<-plot_df_orig<-as.data.frame(plot_df)
names(plot_df)<-names(plot_df_orig)<-activ_name_fixed[rearrange_activ]
  #data.frame(mean=0, dist=0, upper_ci=0, lower_ci=0, activ=" ", context=0,stringsAsFactors = FALSE)

sav_list<-list()

pdf(fil_nam, width=6, height=8)
#par(mfrow=c(4,2))
m<-matrix(c(1,2,3,4,5,6),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m)
par(oma=c(2,1,0,0))
par(mar=c(5,4,3,2))

for (i in rearrange_activ){
    
  mean_val<-list()
  xpoly<-list()
  ypoly<-list()
  
  print(paste0(activ_name_fixed[i],':',sum(labels[[1]]==i)))
  
  for (k in 1:1){
  
    temp_list<-list(c())
    
    # cores=detectCores()
    # cl <- makeCluster(3) #not to overload your computer
    # registerDoParallel(cl)
    
    
    for (j in 1:100){
      #tau<-tauRatioFunc(dates, locs, dbreaks.max,dbreaks.min, tbreak=tim_break,labels, i, 0)
      tau<-tauRatioFunc(dates, locs, dbreaks.max,dbreaks.min, tbreak=tim_break,labels[[k]], i, 0)

      temp_list[[j]]<-tau$tauRat
      #print(j)
    }
    
    
    # temp_df <- foreach(i=1:100, .combine=cbind) %dopar% {
    #   tempMatrix = tauRatioFunc(dates, locs, dbreaks.max,dbreaks.min, tbreak=tim_break,labels[[k]], i, 0)
    #   
    #   tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    # }
    # 
    # stopCluster(cl)
    
    temp_df<-as.matrix(as.data.frame(temp_list))
    mean_val[[k]]<-rowMeans(temp_df)
    mxx_val<-apply(temp_df,1, max, na.rm = TRUE)
    mnn_val<-apply(temp_df,1, min, na.rm = TRUE)
    
    temp_quant<-rowQuantiles(temp_df,probs = c(0.025,0.975),na.rm = TRUE)
    
    
    temp_ci<-matrix(0,nrow = 2,ncol = length(mean_val[[k]]))
    
    temp_ci[1,]<-temp_quant[,2]
    temp_ci[2,]<-temp_quant[,1]
    
    r.mid.h<-dbreaks.mid[-length(dbreaks.mid)]
    xpoly[[k]]<-c(r.mid.h,rev(r.mid.h))
    ypoly[[k]]<-c(temp_ci[1,],rev(temp_ci[2,]))
    
    plot_df[,activ_name_fixed[i]]<-paste0(round(mean_val[[k]],digits = 2),' (',round(temp_ci[2,],digits = 2),',',round(temp_ci[1,],digits = 2),')')
    plot_df_orig[,activ_name_fixed[i]]<-paste0(mean_val[[k]],' (',temp_ci[2,],',',temp_ci[1,],') (',mnn_val,',',mxx_val,')')
    #print(mean_val[[k]])
    sav_list[[i]]<-temp_df
    
  }
  #pdf("random_tau.pdf", width=5, height=5)
  #jpeg(paste(activ_name[i],"_tauA-tauB_ratio_p",pst,"_f",fut,".jpeg",sep = ''), width = 7, height = 7, units = 'in',  res = 300)
  #print(
  #par(mar=c(3,3,3,1))
  plot(r.mid.h, mean_val[[1]] ,ylim=c(0.8,1.1), type="l",log="y", xlim=c(0,max(r.mid.h)),yaxt="n",axes=F,ann=FALSE)
  axis(2,las=1,cex.axis=1.3,at=c(0.8,1,1.1))
  axis(1,las=1,cex.axis=1.3)
  #axis(1,las=1,cex.axis=0.75,at=c(0,1000,3450,5000,10000))
  mtext(text='Distance from index case (m)',side=1,line=2.2,cex=0.9)
  mtext(text=expression(xi["a"]),side=2,line=2.4,cex=1.2)
  abline(h=1,lty=2,col=1)
  #abline(v=500,lty=2,col=1)
  #abline(v=3450,lty=2,col=1)
  polygon(xpoly[[1]],ypoly[[1]],col=rgb(t(col2rgb("dark green")/255),alpha=0.4),border=NA)
  lines(r.mid.h, mean_val[[1]],lty=1,col="dark green")
  #polygon(xpoly[[2]],ypoly[[2]],col=rgb(t(col2rgb("red")/255),alpha=0.6),border=NA)
  #lines(r.mid.h, mean_val[[2]],lty=1,col="red")
  #polygon(xpoly[[3]],ypoly[[3]],col=rgb(t(col2rgb("blue")/255),alpha=0.6),border=NA)
  #lines(r.mid.h, mean_val[[3]],lty=1,col="blue")
  #polygon(xpoly[[4]],ypoly[[4]],col=rgb(t(col2rgb("brown")/255),alpha=0.6),border=NA)
  #lines(r.mid.h, mean_val[[4]],lty=1,col="brown")
  #title(activ_name[i])
  mtext(activ_name_fixed[i],side=1,line=3.8,cex=1.2)
  #)
  #dev.off()
  
  
  
  #tmp_tau<-data.frame(tau_val=tau$tauRat,dist=dbreaks.mid[-1],activ=activ_name[i])
  
  
}
#mtext(plt_titl, side = 3, outer = TRUE,cex=1.25)
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#legend(x = "top",inset = 0, ncol=3,legend = c("Cases before 1st Oct.","Cases after 1st Oct.", "Total"),fill = c("dark green","red","blue"),xpd = TRUE, horiz = TRUE, bty = "n",title="Context",cex = 1.5)

dev.off()

plot_df<-cbind(r.mid.h,plot_df)
write.csv(plot_df,paste0(matching_dist,"--",tim_break,'-temp.csv'))
write.csv(plot_df_orig,paste0(matching_dist,"--",tim_break,'-temp_orig.csv'))
save(sav_list,rearrange_activ,activ_name_fixed,file =paste0(matching_dist,"--",tim_break,'-data_frames.Rdata') )
