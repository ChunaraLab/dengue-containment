### Tau(A)/Tau(B) with confidence interval sampling, modified to monitor affect only in future. abs() removed from original and t() added###

dbreaks.max=seq(100,1500,100) # Maximum distances
window=500 # Window size for distances
#window=dbreaks.max[length(dbreaks.max)] # Window size for distances
dbreaks.min=dbreaks.max-window # minimum distances for windows
dbreaks.min[which(dbreaks.min<0)]<-0
dbreaks.mid<-(dbreaks.min+dbreaks.max)/2 # Middle distance for plotting

tauRatioFunc<-function(dates, locs, dbreaks.max,dbreaks.min, tbreak=30,labels, labelA, labelB){
  
  an<-c("Dewatering","Fish Seeding","Fogging","Indoor Residual Spray", "irs_patient", "Larviciding", "Tap Fixing", "Tire Shredding","multiple")
  
	dates2<-dates+runif(length(dates),-0.001,0.001)
	locs2<-locs
	locs2[,1]<-locs2[,1]+runif(nrow(locs2),-0.001,0.001)
	
	samp<-sample(nrow(locs),nrow(locs),replace=T)
	locs2<-locs2[samp,]
	dates2<-dates2[samp]
	labels<-labels[samp]
	
	match_dmat<-as.matrix(dist(locs2))<=matching_dist2 & as.matrix(dist(locs2))>matching_dist1
#	match_tmat<-t(outer(dates2,dates2,"-")<=tbreak & outer(dates2,dates2,"-")>=0)
	match_tmat<-t(abs(outer(dates2,dates2,"-"))<=(30/2) )
	
	match_mat<-match_dmat*match_tmat
	
	activ_sel_mat<-match_mat *(labels==labelA)
	
	activ_cont_sel_mat<-t( t(activ_sel_mat)*(labels==labelB)   )
	
	diag(activ_cont_sel_mat)<-0
	
	indA_new<-which(rowSums(activ_cont_sel_mat,na.rm = TRUE)>0)
	#prn<-paste(an[labelA],":",length(indA_new),sep='')
	#print(prn)
	
	indB_new<-c()
	
	for (x in indA_new){
	  
	  mtch_cont_list<-which(activ_cont_sel_mat[x,]>0)
	  
	  if (length(mtch_cont_list)==1){
	    indB_new<-c(indB_new,as.numeric(mtch_cont_list))
	    
	  } else {
	    sel_cont<-sample(mtch_cont_list,1)
	    indB_new<-c(indB_new,sel_cont)
	  }
	  
	}
	
	labels[indA_new]<-101
	labels[indB_new]<-100
	
	#samp<-sample(nrow(locs),nrow(locs),replace=T)
	
	tmat<-t(outer(dates2,dates2,"-")<=tbreak & outer(dates2,dates2,"-")>=0)
	dmat<-as.matrix(dist(locs2))
	diag(tmat)<-diag(dmat)<-NA
	tmat[which(tmat==0)]<-NA
  
	#dmat<-dmat[samp,samp]
	#tmat<-tmat[samp,samp]
	
	indA<-which(labels==101)
	indB<-which(labels==100)

	num.matA<-(tmat*dmat)[indA,]
	num.matB<-(tmat*dmat)[indB,]

	A.numnum1<-cumsum(hist(num.matA,breaks=c(0,dbreaks.max,1e10),plot=F)$counts)
	A.numnum2<-cumsum(hist(num.matA,breaks=c(0,dbreaks.min,1e10),plot=F)$counts)
	A.dennum<-sum(tmat[indA,],na.rm=T)
	A.numden1<-cumsum(hist(dmat[indA,],breaks=c(0,dbreaks.max,1e10),plot=F)$counts)
	A.numden2<-cumsum(hist(dmat[indA,],breaks=c(0,dbreaks.min,1e10),plot=F)$counts)
	A.denden<-sum(dmat[indA,]>=0,na.rm=T)
	tauA<-((A.numnum1-A.numnum2)/A.dennum)/((A.numden1-A.numden2)/A.denden)

	B.numnum1<-cumsum(hist(num.matB,breaks=c(0,dbreaks.max,1e10),plot=F)$counts)
	B.numnum2<-cumsum(hist(num.matB,breaks=c(0,dbreaks.min,1e10),plot=F)$counts)
	B.dennum<-sum(tmat[indB,],na.rm=T)
	B.numden1<-cumsum(hist(dmat[indB,],breaks=c(0,dbreaks.max,1e10),plot=F)$counts)
	B.numden2<-cumsum(hist(dmat[indB,],breaks=c(0,dbreaks.min,1e10),plot=F)$counts)
	B.denden<-sum(dmat[indB,]>=0,na.rm=T)
	tauB<-((B.numnum1-B.numnum2)/B.dennum)/((B.numden1-B.numden2)/B.denden)

	tauRat<-tauA/tauB

	outDat<-list(tauA=tauA[-c(1,length(tauRat))],tauB=tauB[-c(1,length(tauRat))],tauRat=tauRat[-c(1,length(tauRat))])
	return(outDat)
}
