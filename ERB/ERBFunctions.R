

#calculates index score
ERB=function(D1,D2,MIC,DIA,MICBrkpt1,MICBrkpt2,VM1,M1,m1,VM2,M2,m2,withinOneVector){

  tempMIC=MIC
  tempDIA=DIA
  
  temp=max(VM1,M1,m1,VM2,M2,m2)
  VM1=temp/VM1; M1=temp/M1; m1=temp/m1
  VM2=temp/VM2; M2=temp/M2; m2=temp/m2
  
  ###outside one intermediate range
  II=0; SS=0; VM=0; RR=0; M=0; m=0
  MIC=tempMIC[withinOneVector==0]
  DIA=tempDIA[withinOneVector==0]
  n=length(MIC)
  for (i in 1:n){
    #SS
    if(MIC[i]<=MICBrkpt1 & DIA[i]>=D2)SS=SS+1
    #intermediate
    else if ((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & (DIA[i]>D1 & DIA[i]<D2)) II=II+1
    #resistant
    else if(MIC[i]>=MICBrkpt2 & DIA[i]<=D1)RR=RR+1
    #I MIC S DIA
    else if((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & DIA[i]>=D2)m=m+1
    #R MIC S DIA
    else if(MIC[i]>=MICBrkpt2 & DIA[i]>=D2)VM=VM+1
    #S MIC I DIA
    else if(MIC[i]<=MICBrkpt1 & (DIA[i]>D1 & DIA[i]<D2))m=m+1
    #R MIC I DIA
    else if(MIC[i]>=MICBrkpt2 & (DIA[i]>D1 & DIA[i]<D2))m=m+1
    #S MIC R DIA
    else if(MIC[i]<=MICBrkpt1 & DIA[i]<=D1)M=M+1
    #I MIC R DIA
    else if((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & DIA[i]<=D1)m=m+1
  }
        
  index1=VM2*VM/n+M2*M/n+m2*m/n
  CorrectPerc1=sum(SS,RR,II)/n*100
  VMPerc1=VM/n*100
  MPerc1=M/n*100
  mPerc1=m/n*100
  SS1=SS; RR1=RR; II1=II
  VMCount1=VM; MCount1=M; mCount1=m

  ###inside one intermediate range
  II=0; SS=0; VM=0; RR=0; M=0; m=0
  if(sum(withinOneVector)>0){
    MIC=tempMIC[withinOneVector==1]
    DIA=tempDIA[withinOneVector==1]
    n=length(MIC)
    for (i in 1:n){
      #SS
      if(MIC[i]<=MICBrkpt1 & DIA[i]>=D2)SS=SS+1
      #intermediate
      else if ((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & (DIA[i]>D1 & DIA[i]<D2)) II=II+1
      #resistant
      else if(MIC[i]>=MICBrkpt2 & DIA[i]<=D1)RR=RR+1
      #I MIC S DIA
      else if((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & DIA[i]>=D2)m=m+1
      #R MIC S DIA
      else if(MIC[i]>=MICBrkpt2 & DIA[i]>=D2)VM=VM+1
      #S MIC I DIA
      else if(MIC[i]<=MICBrkpt1 & (DIA[i]>D1 & DIA[i]<D2))m=m+1
      #R MIC I DIA
      else if(MIC[i]>=MICBrkpt2 & (DIA[i]>D1 & DIA[i]<D2))m=m+1
      #S MIC R DIA
      else if(MIC[i]<=MICBrkpt1 & DIA[i]<=D1)M=M+1
      #I MIC R DIA
      else if((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & DIA[i]<=D1)m=m+1
    }
    index2=VM1*VM/n+M1*M/n+m1*m/n
    CorrectPerc2=sum(SS,RR,II)/n*100
    VMPerc2=VM/n*100
    MPerc2=M/n*100
    mPerc2=m/n*100
    SS2=SS; RR2=RR; II2=II
    VMCount2=VM; MCount2=M; mCount2=m
  }else{
    index2=0
    CorrectPerc2=100
    VMPerc2=0
    MPerc2=0
    mPerc2=0
    SS2=0; RR2=0; II2=0
    VMCount2=0; MCount2=0; mCount2=0
  }


  return(list(idx=index1+index2,CorPerc1=CorrectPerc1,VMPerc1=VMPerc1,MPerc1=MPerc1,mPerc1=mPerc1,CorrectCount1=sum(SS1,RR1,II1),
           VMCount1=VMCount1,MCount1=MCount1,mCount1=mCount1,
           CorPerc2=CorrectPerc2,VMPerc2=VMPerc2,MPerc2=MPerc2,mPerc2=mPerc2,
           CorrectCount2=sum(SS2,RR2,II2),VMCount2=VMCount2,MCount2=MCount2,mCount2=mCount2))
}


#finds optimum DIA given breakpoints M1 and M2 error rate bounded method
findBrkptsERB=function(MIC,DIA,VM1=10,M1=10,m1=40,VM2=2,M2=2,m2=5,
  MICBrkptL,MICBrkptU,minWidth=4,maxWidth=20){
  
  MICBrkpt1=MICBrkptL
  MICBrkpt2=MICBrkptU
  
  
  #find optimal
  parms=findBrkptsERBC(MIC,DIA,VM1,M1,m1,VM2,M2,m2,MICBrkpt1,MICBrkpt2,minWidth,maxWidth)
  D1=parms$D1; D2=parms$D2
  
  
  #find information for plotting and display information
  N=length(MIC)
  withinOne=rep(0,length(MIC))
  withinOne[which(MIC>=MICBrkpt1 & MIC<=MICBrkpt2)]=1
  numwithinOne=sum(withinOne)

  cat('Optimal DIA Breakpoints for ERB:',D1,D2,'\n')
  cat('Number of Isolates: ',N,'\n')
  cat('Number Observed Outside One of Intermediate Range: ',N-numwithinOne,'\n')
  cat('Number Observed Within One of Intermediate Range: ',numwithinOne,'\n')
  temp=matrix(nrow=2,ncol=5)
  parms=ERB(D1,D2,MIC,DIA,MICBrkpt1,MICBrkpt2,VM1,M1,m1,VM2,M2,m2,withinOne)
  cat('Index Score = ',parms$idx,'\n \n')
  cat('Count (%) \n')
  temp[1,2:5]=c(paste(parms$CorrectCount2,' (',round(parms$CorPerc2,digits=2),')',sep=''),
                paste(parms$VMCount2,' (',round(parms$VMPerc2,digits=2),')',sep=''),
                paste(parms$MCount2,' (',round(parms$MPerc2,digits=2),')',sep=''),
                paste(parms$mCount2,' (',round(parms$mPerc2,digits=2),')',sep=''))
  temp[2,2:5]=c(paste(parms$CorrectCount1,' (',round(parms$CorPerc1,digits=2),')',sep=''),
                paste(parms$VMCount1,' (',round(parms$VMPerc1,digits=2),')',sep=''),
                paste(parms$MCount1,' (',round(parms$MPerc1,digits=2),')',sep=''),
                paste(parms$mCount1,' (',round(parms$mPerc1,digits=2),')',sep=''))
  temp[1:2,1]=c('Within 1','Outside 1')
  temp=as.data.frame(temp)
  names(temp)=c('Range','Agree','Very Major','Major','Minor')
  name.width <- max(sapply(names(temp), nchar))
  names(temp) <- format(names(temp), width = name.width, justify = "centre")
  print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
  if(parms$idx==0) cat('\n \n NOTE: When the index score equals 0, there can be a range of optimal DIA breakpoints.  We present the smallest breakpoint in this range.')
  return(list(D1=D1,D2=D2))
}



ERBGivenDIA=function(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,DIABrkptL,DIABrkptU,
  VM1=10,M1=10,m1=40,VM2=2,M2=2,m2=5,flipGraph=0){
  
  MICBrkpt1=MICBrkptL
  MICBrkpt2=MICBrkptU
  DIABrkpt1=DIABrkptL
  DIABrkpt2=DIABrkptU
    
  D1=DIABrkpt1; D2=DIABrkpt2
  
  if(D2<=D1){
    stop('Lower DIA Breakpoint must be less than Upper DIA Breakpoint.')
  }
  
  N=length(MIC)
  withinOne=rep(0,length(MIC))
  
  #observations within one of intermediate range
  withinOne=rep(0,length(MIC))
  withinOne[which(MIC>=MICBrkpt1 & MIC<=MICBrkpt2)]=1
  numwithinOne=sum(withinOne)
  
  cat('Classification for DIA Breakpoints for ERB:',D1,D2,'\n')
  cat('Number of Isolates: ',N,'\n')
  cat('Number Observed Outside One of Intermediate Range: ',N-numwithinOne,'\n')
  cat('Number Observed Within One of Intermediate Range: ',numwithinOne,'\n')
  temp=matrix(nrow=2,ncol=5)
  parms=ERB(D1,D2,MIC,DIA,MICBrkpt1,MICBrkpt2,VM1,M1,m1,VM2,M2,m2,withinOne)
  cat('Index Score = ',parms$idx,'\n \n')
  cat('Count (%) \n')
  temp[1,2:5]=c(paste(parms$CorrectCount2,' (',round(parms$CorPerc2,digits=2),')',sep=''),
                paste(parms$VMCount2,' (',round(parms$VMPerc2,digits=2),')',sep=''),
                paste(parms$MCount2,' (',round(parms$MPerc2,digits=2),')',sep=''),
                paste(parms$mCount2,' (',round(parms$mPerc2,digits=2),')',sep=''))
  temp[2,2:5]=c(paste(parms$CorrectCount1,' (',round(parms$CorPerc1,digits=2),')',sep=''),
                paste(parms$VMCount1,' (',round(parms$VMPerc1,digits=2),')',sep=''),
                paste(parms$MCount1,' (',round(parms$MPerc1,digits=2),')',sep=''),
                paste(parms$mCount1,' (',round(parms$mPerc1,digits=2),')',sep=''))
  temp[1:2,1]=c('Within 1','Outside 1')
  temp=as.data.frame(temp)
  names(temp)=c('Range','Agree','Very Major','Major','Minor')
  name.width <- max(sapply(names(temp), nchar))
  names(temp) <- format(names(temp), width = name.width, justify = "centre")
  print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
  invisible()
}


PlotBrkptsERB2=function(MIC,DIA,xcens,ycens,VM1=10,M1=10,m1=40,VM2=2,M2=2,m2=5,
                        MICBrkptL,MICBrkptU,minWidth=4,maxWidth=20,flipGraph,logConvert){
  
  MICBrkpt1=MICBrkptL
  MICBrkpt2=MICBrkptU
  
  
  #find optimal
  parms=findBrkptsERBC(MIC,DIA,VM1,M1,m1,VM2,M2,m2,MICBrkpt1,MICBrkpt2,minWidth,maxWidth)
  D1=parms$D1; D2=parms$D2

  fit=plotBrkPtsERB(MIC,DIA,xcens,ycens,MICBrkpt1,MICBrkpt2,D1,D2,flipGraph,logConvert)
  return(fit)
}

PlotBrkptsERBGiven=function(MIC,DIA,xcens,ycens,VM1=10,M1=10,m1=40,VM2=2,M2=2,m2=5,
                        MICBrkptL,MICBrkptU,D1,D2,minWidth=4,maxWidth=20,flipGraph,logConvert){
  
  MICBrkpt1=MICBrkptL
  MICBrkpt2=MICBrkptU
  
  fit=plotBrkPtsERB(MIC,DIA,xcens,ycens,MICBrkpt1,MICBrkpt2,D1,D2,flipGraph,logConvert)
  return(fit)
}


#plot single scatterplot 
plotBrkPtsERB=function(MIC,DIA,xcens,ycens,MICBrkpt1,MICBrkpt2,DIABrkpt1,DIABrkpt2,flipGraph,logConvert){

  MIC1=MIC
  DIA1=DIA
  MIC[xcens==1 & MIC==max(MIC)]=max(MIC)+1
  MIC[xcens==-1 & MIC==min(MIC)]=min(MIC)-1
  DIA[ycens==1 & DIA==max(DIA)]=max(DIA)+1
  DIA[ycens==-1 & DIA==min(DIA)]=min(DIA)-1
  M1=MICBrkpt1; M2=MICBrkpt2
  D1=DIABrkpt1; D2=DIABrkpt2
  a1=data.frame(MIC,DIA)
  a1=count(a1,c('MIC','DIA'))
  Freq=a1$freq; MIC=a1$MIC; DIA=a1$DIA
  
  MICBrkpt1=MICBrkpt1+.5
  MICBrkpt2=MICBrkpt2-.5
  DIABrkpt1=DIABrkpt1+.5
  DIABrkpt2=DIABrkpt2-.5
  
  n=length(MIC); classification=rep(NA,n)
  for (i in 1:n){
    #SS
    if(MIC[i]<=M1 & DIA[i]>=D2)classification[i]='Correct'
    #intermediate
    else if ((MIC[i]>M1 & MIC[i]<M2) & (DIA[i]>D1 & DIA[i]<D2)) classification[i]='Correct'
    #resistant
    else if(MIC[i]>=M2 & DIA[i]<=D1)classification[i]='Correct'
    #I MIC S DIA
    else if((MIC[i]>M1 & MIC[i]<M2) & DIA[i]>=D2)classification[i]='Minor'
    #R MIC S DIA
    else if(MIC[i]>=M2 & DIA[i]>=D2)classification[i]='Very Major'
    #S MIC I DIA
    else if(MIC[i]<=M1 & (DIA[i]>D1 & DIA[i]<D2))classification[i]='Minor'
    #R MIC I DIA
    else if(MIC[i]>=M2 & (DIA[i]>D1 & DIA[i]<D2))classification[i]='Minor'
    #S MIC R DIA
    else if(MIC[i]<=M1 & DIA[i]<=D1)classification[i]='Major'
    #I MIC R DIA
    else if((MIC[i]>M1 & MIC[i]<M2) & DIA[i]<=D1)classification[i]='Minor'
  }
  a1=data.frame(MIC,DIA,Freq,classification,stringsAsFactors=FALSE)
#   a1=rbind(a1,c(999,999,0,'Minor'))
#   a1=rbind(a1,c(999,999,0,'Major'))
#   a1=rbind(a1,c(999,999,0,'Very Major'))
  a1[,1:3] = apply(a1[,1:3], 2, function(x) as.numeric(x));
  a1$classification=factor(a1$classification,levels=c("Correct", "Minor", "Major","Very Major"))
  
  if(flipGraph=='Yes' && logConvert==TRUE){
    fit=ggplot(a1,aes(MIC,DIA))+geom_text(aes(label=Freq,color=factor(classification)),size=4,show_guide=FALSE)+
      geom_point(aes(group=factor(classification),color=factor(classification)),size=0)+
      geom_vline(xintercept=MICBrkpt1,lty=2,alpha=.4)+
      geom_vline(xintercept=MICBrkpt2,lty=2,alpha=.4)+
      geom_hline(yintercept=DIABrkpt1,lty=2,alpha=.4)+
      geom_hline(yintercept=DIABrkpt2,lty=2,alpha=.4)+
      labs(x='MIC (Dilution Test in log(ug/mL))',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(breaks = seq(min(MIC1)-1,max(MIC1)+1,by=1),
          labels = c(paste("<",min(MIC1),sep=''),seq(min(MIC1),max(MIC1),by=1), paste(">",max(MIC1),sep='')),
          limits = c(min(MIC1)-1,max(MIC1)+1))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_color_manual(values=c('Correct'='Black','Minor'='blue','Major'='#CC9900','Very Major'='red'))+
      theme(
          legend.position='top',
          legend.title=element_blank(),
          legend.key=element_rect(fill="white",colour="white"),
          legend.text = element_text(size = 14))+
          guides(colour = guide_legend(override.aes = list(size=3,alpha = 1)))
  }
  if(flipGraph=='Yes' && logConvert==FALSE){
    MICBrkpt1=2^MICBrkpt1
    MICBrkpt2=2^MICBrkpt2
    MIC2=MIC1
    a1$MIC=2^a1$MIC
    MICTemp=c(min(MIC1)-1,min(MIC1):max(MIC1),max(MIC1)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC1):max(MIC1))
    
    fit=ggplot(a1,aes(MIC,DIA))+geom_text(aes(label=Freq,color=factor(classification)),size=4,show_guide=FALSE)+
      geom_point(aes(group=factor(classification),color=factor(classification)),size=0)+
      geom_vline(xintercept=MICBrkpt1,lty=2,alpha=.4)+
      geom_vline(xintercept=MICBrkpt2,lty=2,alpha=.4)+
      geom_hline(yintercept=DIABrkpt1,lty=2,alpha=.4)+
      geom_hline(yintercept=DIABrkpt2,lty=2,alpha=.4)+
      labs(x='MIC (Dilution Test in ug/mL)',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp,
                         labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_color_manual(values=c('Correct'='Black','Minor'='blue','Major'='#CC9900','Very Major'='red'))+
      theme(
        legend.position='top',
        legend.title=element_blank(),
        legend.key=element_rect(fill="white",colour="white"),
        legend.text = element_text(size = 14))+
      guides(colour = guide_legend(override.aes = list(size=3,alpha = 1)))
  }
  if(flipGraph=='No' && logConvert==TRUE){
    fit=ggplot(a1,aes(DIA,MIC))+geom_text(aes(label=Freq,color=factor(classification)),size=4,show_guide=FALSE)+
      geom_point(aes(group=factor(classification),color=factor(classification)),size=0)+
      geom_hline(yintercept=MICBrkpt1,lty=2,alpha=.4)+
      geom_hline(yintercept=MICBrkpt2,lty=2,alpha=.4)+
      geom_vline(xintercept=DIABrkpt1,lty=2,alpha=.4)+
      geom_vline(xintercept=DIABrkpt2,lty=2,alpha=.4)+
      labs(y='MIC (Dilution Test in log(ug/mL))',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(breaks = seq(min(MIC1)-1,max(MIC1)+1,by=1),
          labels = c(paste("<",min(MIC1),sep=''),seq(min(MIC1),max(MIC1),by=1), paste(">",max(MIC1),sep='')),
          limits = c(min(MIC1)-1,max(MIC1)+1))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_color_manual(values=c('Correct'='Black','Minor'='blue','Major'='#CC9900','Very Major'='red'))+
      theme(
          legend.position='top',
          legend.title=element_blank(),
          legend.key=element_rect(fill="white",colour="white"),
          legend.text = element_text(size = 14))+
          guides(colour = guide_legend(override.aes = list(size=3,alpha = 1)))
  }
  if(flipGraph=='No' && logConvert==FALSE){
    MICBrkpt1=2^MICBrkpt1
    MICBrkpt2=2^MICBrkpt2
    MIC2=MIC1
    a1$MIC=2^a1$MIC
    MICTemp=c(min(MIC1)-1,min(MIC1):max(MIC1),max(MIC1)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC1):max(MIC1))
    
    fit=ggplot(a1,aes(DIA,MIC))+geom_text(aes(label=Freq,color=factor(classification)),size=4,show_guide=FALSE)+
      geom_point(aes(group=factor(classification),color=factor(classification)),size=0)+
      geom_hline(yintercept=MICBrkpt1,lty=2,alpha=.4)+
      geom_hline(yintercept=MICBrkpt2,lty=2,alpha=.4)+
      geom_vline(xintercept=DIABrkpt1,lty=2,alpha=.4)+
      geom_vline(xintercept=DIABrkpt2,lty=2,alpha=.4)+
      labs(y='MIC (Dilution Test in ug/mL)',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(trans=log2_trans(),
           limits=c(min(MICTemp),max(MICTemp)),
           breaks=MICTemp,
           labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_color_manual(values=c('Correct'='Black','Minor'='blue','Major'='#CC9900','Very Major'='red'))+
      theme(
        legend.position='top',
        legend.title=element_blank(),
        legend.key=element_rect(fill="white",colour="white"),
        legend.text = element_text(size = 14))+
      guides(colour = guide_legend(override.aes = list(size=3,alpha = 1)))
    
    
  }
         
  return(fit)
  
}


bootStrapERB=function(MIC,DIA,MICBrkptL,MICBrkptU,VM1=10,M1=10,m1=40,VM2=2,M2=2,m2=5,
                      minWidth=4,maxWidth=20,session){
  
  MICBrkpt1=MICBrkptL
  MICBrkpt2=MICBrkptU
  
  #bootstrap
  DIABrkptL<<-rep(NA,12000); DIABrkptU<<-rep(NA,12000);
  n=length(MIC)

  withProgress(session, min=1, max=12000, {
    setProgress(message = 'Calculation in progress')
    for(i in 1:12000){
      #resample points
      if(i %% 1000==0) setProgress(value = i)
      idx=sample(seq(1,n,by=1),n,replace=T)
      MIC_sam=MIC[idx]
      DIA_sam=DIA[idx]
      
      #get weighted breakpoints
      parms=findBrkptsERBC(MIC_sam,DIA_sam,VM1,M1,m1,VM2,M2,m2,MICBrkpt1,MICBrkpt2,minWidth,maxWidth)
      DIABrkptL[i]<<-parms$D1;   DIABrkptU[i]<<-parms$D2
    }
  })

#   DIABrkptL=c(one[[1]]$D1,one[[2]]$D1,one[[3]]$D1,one[[4]]$D1)
#   DIABrkptU=c(one[[1]]$D2,one[[2]]$D2,one[[3]]$D2,one[[4]]$D2)
  
  
  #print results
  tab=table(DIABrkptL,DIABrkptU)
  a1=tab/margin.table(tab)*100
  a1=as.data.frame(a1)
  a1=a1[order(a1$Freq,decreasing=T),]
  a1=a1[a1$Freq!=0,]
  a2=a1
  a2$CumFreq=cumsum(a2$Freq)
  a2[,3]=round(a2[,3],2)
  a2[,4]=round(a2[,4],2)
  
  cat('Bootstrap samples = 12000 \n')
  cat('\n-------DIA Breakpoints by Confidence--------\n')
  temp=data.frame(DIABrkptL=a2[,1],DIABrkptU=a2[,2],Percent=a2[,3],Cumulative=a2[,4])
  temp[,1:4] = apply(temp[,1:4], 2, function(x) as.character(x));
  name.width <- max(sapply(names(temp), nchar))
  names(temp) <- format(names(temp), width = name.width, justify = "centre")
  print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
  return(a1)

}

plotBootDataERB=function(bootData){
  
  bootData$cumFreq=cumsum(bootData$Freq)
  bootData[,1]=as.numeric(as.character(bootData[,1]))
  bootData[,2]=as.numeric(as.character(bootData[,2]))
  D1=bootData[,1]
  D2=bootData[,2]
  idx=max(which(bootData$cumFreq<95))
  a1=bootData[1:(idx+1),]
  if((idx+1)!=nrow(bootData))
    a2=bootData[(idx+2):nrow(bootData),]
  a1$Freq=round(a1$Freq,2)
  
  fit=ggplot(data=a1,aes(x=DIABrkptL,y=DIABrkptU,label=Freq))+geom_text(size=7.5)+
    scale_x_continuous(breaks = seq(min(D1),max(D1),by=1),
                       limits = c(min(D1),max(D1)))+
    scale_y_continuous(breaks = seq(min(D2),max(D2),by=1),
                       limits = c(min(D2),max(D2)))+
    labs(x='Lower DIA Breakpoint',y='Upper DIA Breakpoint',
         title=expression(atop("DIA Breakpoint Distribution", atop("black points are outside 95% confidence", ""))))+
    theme(
      plot.title = element_text(size = 17, hjust=0.5,vjust=1),
      legend.position = "none")
  if((idx+1)!=nrow(bootData))
    fit=fit+geom_point(data=a2,(aes(x=DIABrkptL,y=DIABrkptU)),size=3)
  print(fit)
  
  invisible()
}






