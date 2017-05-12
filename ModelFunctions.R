
plotBrkptModel=function(modelBrkpt){
  
  modelBrkpt$cumFreq=cumsum(modelBrkpt$Freq)
  modelBrkpt[,1]=as.numeric(as.character(modelBrkpt[,1]))
  modelBrkpt[,2]=as.numeric(as.character(modelBrkpt[,2]))
  D1=modelBrkpt[,1]
  D2=modelBrkpt[,2]
  if(modelBrkpt$cumFreq[1]>=95){
    idx=0
  }else{
    idx=max(which(modelBrkpt$cumFreq<95))
  }
  
  a1=modelBrkpt[1:(idx+1),]
  if((idx+1)!=nrow(modelBrkpt))
    a2=modelBrkpt[(idx+2):nrow(modelBrkpt),]
  a1$Freq=round(a1$Freq,2)
  
  
  fit=ggplot(data=a1,aes(x=DIABrkptL,y=DIABrkptU,label=Freq))+geom_text(size=7.5)+
    scale_x_continuous(breaks = seq(min(D1),max(D1),by=1),
         limits = c(min(D1),max(D1)))+
    scale_y_continuous(breaks = seq(min(D2),max(D2),by=1),
         limits = c(min(D2),max(D2)))+
    labs(x='Lower DIA Breakpoint',y='Upper DIA Breakpoint',
         title=expression(atop("DIA Breakpoint Distribution", atop("black points are outside 95% probability", ""))))+
    theme(
      plot.title = element_text(size = 17, hjust=0.5,vjust=1),
      legend.position = "none")
  if((idx+1)!=nrow(modelBrkpt))
    fit=fit+geom_point(data=a2,(aes(x=DIABrkptL,y=DIABrkptU)),size=3)
  print(fit)
  
  invisible()
}


compareFitsPlot=function(dataSpline,dataLog,MIC,DIA,xcens,ycens,M1,M2,flipGraph,logConvert){
      
  logData=data.frame(grid=xgrid,fit=dataLog$fit)
  splineData=data.frame(grid=xgrid,fit=dataSpline$fit)
  
  M1=M1+.5
  M2=M2-.5
  xobs=MIC
  yobs=DIA
  xobs1=xobs
  DIA1=yobs
  xobs[xcens==1 & xobs==max(xobs)]=max(xobs)+1
  xobs[xcens==-1 & xobs==min(xobs)]=min(xobs)-1
  yobs[ycens==1 & yobs==max(yobs)]=max(yobs)+1
  yobs[ycens==-1 & yobs==min(yobs)]=min(yobs)-1
  a1=data.frame(xobs,yobs)
  a1=count(a1,c('xobs','yobs'))

  if(flipGraph=='Yes' && logConvert==TRUE){
    fit=ggplot(a1,aes(xobs,yobs))+geom_text(aes(label=freq),size=5)+
      geom_line(data=logData,aes(x=grid,y=fit))+
      geom_line(data=splineData,aes(x=grid,y=fit,color='red'))+
      geom_vline(xintercept=M1,lty=2,alpha=.5)+
      geom_vline(xintercept=M2,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in log(ug/mL))',y='DIA (Diffusion Test in mm)',title='Black = Logistic Fit, Red = Spline Fit')+
      scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
             labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
             limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
             labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
             limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(legend.position='NONE',
            plot.title = element_text(size = 17, hjust=0.5,vjust=1))
  }
  if(flipGraph=='Yes' && logConvert==FALSE){
    M1=2^M1
    M2=2^M2
    MIC2=xobs1
    a1$xobs=2^a1$xobs
    MICTemp=c(min(MIC2)-1,min(MIC2):max(MIC2),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    logData$grid=2^logData$grid
    splineData$grid=2^splineData$grid
    fit=ggplot(a1,aes(xobs,yobs))+geom_text(aes(label=freq),size=5)+
      geom_line(data=logData,aes(x=grid,y=fit))+
      geom_line(data=splineData,aes(x=grid,y=fit,color='red'))+
      geom_vline(xintercept=M1,lty=2,alpha=.5)+
      geom_vline(xintercept=M2,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in ug/mL)',y='DIA (Diffusion Test in mm)',title='Black = Logistic Fit, Red = Spline Fit')+
      scale_x_continuous(trans=log2_trans(),
           limits=c(min(MICTemp),max(MICTemp)),
           breaks=MICTemp,
           labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
           labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
           limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(legend.position='NONE',
            plot.title = element_text(size = 17, hjust=0.5,vjust=1))
  }
  if(flipGraph=='No' &&  logConvert==TRUE){
    fit=ggplot(a1,aes(yobs,xobs))+geom_text(aes(label=freq),size=5)+
      geom_line(data=logData,aes(x=fit,y=grid))+
      geom_line(data=splineData,aes(x=fit,y=grid,color='red'))+
      geom_hline(yintercept=M1,lty=2,alpha=.5)+
      geom_hline(yintercept=M2,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in log(ug/mL))',x='DIA (Diffusion Test in mm)',title='Black = Logistic Fit, Red = Spline Fit')+
      scale_y_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
           labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
           limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
           labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
           limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(legend.position='NONE',
            plot.title = element_text(size = 17, hjust=0.5,vjust=1))
  }
  if(flipGraph=='No' &&  logConvert==FALSE){
    
    theme_old <- theme_update(
      axis.text.y=element_text(size=15, color='black'),
      axis.text.x=element_text(size=14, color='black'),
      axis.title.x=element_text(size=17, margin=margin(t=5)),
      axis.title.y=element_text(size=17, vjust=0.2, angle=90),
      panel.background=element_rect(color='black', fill='grey93'),
      panel.grid.minor=element_blank()
    )
    
    
    M1=2^M1
    M2=2^M2
    MIC2=xobs1
    a1$xobs=2^a1$xobs
    MICTemp=c(min(MIC2)-1,min(MIC2):max(MIC2),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    logData$grid=2^logData$grid
    splineData$grid=2^splineData$grid
    fit=ggplot(a1,aes(yobs,xobs))+geom_text(aes(label=freq),size=5.2)+
      geom_line(data=logData,aes(x=fit,y=grid))+
      geom_line(data=splineData,aes(x=fit,y=grid,color='red'))+
      geom_hline(yintercept=M1,lty=2,alpha=.5)+
      geom_hline(yintercept=M2,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in ug/mL)',x='DIA (Diffusion Test in mm)',title='Black = Logistic Fit, Red = Spline Fit')+
      scale_y_continuous(trans=log2_trans(),
           limits=c(min(MICTemp),max(MICTemp)),
           breaks=MICTemp,
           labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
             labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
             limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(legend.position='NONE',
            plot.title = element_text(size = 17, hjust=0.5,vjust=1))
  }    
  
  return(fit)
}


probDIAClass=function(data,MICBrkptL,MICBrkptU,DIASet1,DIASet2){
    
  f=data$fit
  density=data$dens
  
  e=sqrt(.5^2+.5^2)
  t=sqrt(1.5^2+1.5^2)
  
  grid=xgrid[xgrid>=min(MIC) & xgrid<=max(MIC)]
  weights=density[xgrid>=min(MIC) & xgrid<=max(MIC)]
  f=f[xgrid>=min(MIC) & xgrid<=max(MIC)]
  weights=weights/sum(weights)
  
  
  M1=MICBrkptL-.5
  M2=MICBrkptU-.5
  D1=DIASet1[1];  D2=DIASet1[2]
  if(is.na(DIASet2[1])==TRUE){
    D11=0; D22=0
  }else{
    D11=DIASet2[1];  D22=DIASet2[2]    
  }
  
  probDIAC=rep(NA,length(grid))
  probDIAIC=rep(NA,length(grid))
  probDIAC1=rep(NA,length(grid))
  probDIAIC1=rep(NA,length(grid))
  probMIC=rep(NA,length(grid))
  
  print(c(D1,D2))
  
  ### DIA Probs
  for(i in 1:length(grid)){
    d=f[i]
    m=grid[i]
    
    ### MIC
    #Susceptible
    if (m<=M1)  probMIC[i]=pnorm(MICBrkptL,m,e)
    #intermediate
    if(m>M1 & m<M2) probMIC[i]=pnorm(MICBrkptU-1,m,e)-pnorm(MICBrkptL,m,e)
    #resistant
    if(m>=M2) probMIC[i]=1-pnorm(MICBrkptU-1,m,e)
    
    ###DIA
    #Susceptible
    if(m<=M1){
      probDIAC[i]=1-pnorm(D2-.5,d,t)
      probDIAC1[i]=1-pnorm(D22-.5,d,t)
    #intermediate
    }else if(m>M1 && m<M2){
      probDIAC[i]=pnorm(D2-.5,d,t)-pnorm(D1+.5,d,t)
      probDIAC1[i]=pnorm(D22-.5,d,t)-pnorm(D11+.5,d,t)
    #Resistant
    }else{
      probDIAC[i]=pnorm(D1+.5,d,t)
      probDIAC1[i]=pnorm(D11+.5,d,t)
    }
    probDIAIC[i]=1-probDIAC[i]
    probDIAIC1[i]=1-probDIAC1[i]
  }
  
  ### Print DIA
  if(is.na(DIASet2[1])==FALSE){
    cat('DIA Breakpoints \n')
    cat('Set 1: ',DIASet1[1],', ',DIASet1[2],'\n',sep='')
    cat('Set 2: ',DIASet2[1],', ',DIASet2[2],'\n',sep='')
	  CorrectM=sum(probMIC)/length(grid)
    Correct1=sum(probDIAC)/length(grid)
    Correct2=sum(probDIAC1)/length(grid)
    cat('\nProbability Correct Classification \n')
    temp=data.frame(CorrectM,Correct1,Correct2)
    temp[,1:3]=round(temp[,1:3],digits=3)
    temp[,1:3]=as.character(temp[,1:3])
    names(temp)=c('MIC Correct','DIA 1 Correct','DIA 2 Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
	  CorrectM=sum(probMIC*weights)
    Correct1=sum(probDIAC*weights)
    Correct2=sum(probDIAC1*weights)
    cat('\nProbability Correct Classification Weighted by Isolate Distribution \n')
    temp=data.frame(CorrectM,Correct1,Correct2)
    temp[,1:3]=round(temp[,1:3],digits=4)
    temp[,1:3]=as.character(temp[,1:3])
    names(temp)=c('MIC Correct','DIA 1 Correct','DIA 2 Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
    a1=data.frame(grid,probDIAC,probDIAC1,probMIC,weights)
  }else{
    cat('DIA Breakpoints: ',DIASet1[1],', ',DIASet1[2],'\n',sep='')
    Correct1=sum(probDIAC)/length(grid)
	  CorrectM=sum(probMIC)/length(grid)
    cat('\nProbability Correct Classification \n')
    temp=data.frame(CorrectM,Correct1)
    temp[,1:2]=round(temp[,1:2],digits=3)
    temp[,1:2]=as.character(temp[,1:2])
    names(temp)=c('MIC Correct','DIA Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
    Correct1=sum(probDIAC*weights)
	  CorrectM=sum(probMIC*weights)
    cat('\nProbability Correct Classification Weighted by Isolate Distribution \n')
    temp=data.frame(CorrectM,Correct1)
    temp[,1:2]=round(temp[,1:2],digits=4)
    temp[,1:2]=as.character(temp[,1:2])
    names(temp)=c('MIC Correct','DIA Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
    a1=data.frame(grid,probDIAC,probDIAC1,probMIC,weights)
  }
  
  return(a1)
  
}

plotUnderlyingDistibution=function(a1,logConvert,M1,M2){
  
  M1=M1-.5
  M2=M2-.5
  
  if(logConvert==TRUE){
    plt=ggplot(a1,aes(x=grid,y=weights))+geom_area(alpha=.3,fill='red',color='black')+
      labs(x='',y='',title='Isolate Distribution')+
      geom_vline(xintercept=M1,lty=2,alpha=1)+
      geom_vline(xintercept=M2,lty=2,alpha=1)+
      scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(MIC),max(MIC)))+
      theme(plot.title = element_text(size = 17, hjust=0.5,vjust=1))
  }else{
    
    M1=2^(M1)
    M2=2^(M2)
    a1$grid=2^a1$grid
    MICTemp=min(MIC):max(MIC)
    MICTemp=2^MICTemp
    
    plt=ggplot(a1,aes(x=grid,y=weights))+geom_area(alpha=.3,fill='red',color='black')+
      labs(x='',y='',title='Isolate Distribution')+
      geom_vline(xintercept=M1,lty=2,alpha=1)+
      geom_vline(xintercept=M2,lty=2,alpha=1)+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp)+
      theme(plot.title = element_text(size = 17, hjust=0.5,vjust=1))
  }
  
  return(plt)
  
}

plotProbDIAClass=function(a1,M1,M2,DIASet1,DIASet2,logConvert){
  
#   theme_old <- theme_update(
#     axis.text.y=element_text(size=15, color='black'),
#     axis.text.x=element_text(size=14, color='black'),
#     axis.title.x=element_text(size=17, vjust=0.2),
#     axis.title.y=element_text(size=17, margin=margin(l=5)),
#     panel.background=element_rect(color='black', fill='grey93'),
#     panel.grid.minor=element_blank(),
#     plot.title = element_text(size = 17, hjust=0.5,vjust=1),
#     legend.position='bottom',
#     legend.key=element_rect(fill="white",colour="white"),
#     legend.text = element_text(size = 16),
#     legend.title=element_text(size=15)
#   )
  
  
  M1=M1-.5
  M2=M2-.5
  
  xgrid=a1$grid
  probDIAC=a1$probDIAC
  if(is.na(DIASet2[1])==FALSE)
    probDIAC1=a1$probDIAC1
  probMIC=a1$probMIC
  
  
  ### Plot  
  grid1=xgrid[xgrid<M1-.02]
  grid2=xgrid[xgrid>M1+.02 & xgrid<M2-.02]
  grid3=xgrid[xgrid>M2+.02]
  
  probMIC1=probMIC[xgrid<M1-.02]
  probMIC2=probMIC[xgrid>M1+.02 & xgrid<M2-.02]
  probMIC3=probMIC[xgrid>M2+.02]
  
  ProbCorrectDIA1=probDIAC[xgrid<M1-.02]
  ProbCorrectDIA2=probDIAC[xgrid>M1+.02 & xgrid<M2-.02]
  ProbCorrectDIA3=probDIAC[xgrid>M2+.02]
  
  if(is.na(DIASet2[1])==FALSE){
    ProbCorrectDIA11=probDIAC1[xgrid<M1-.02]
    ProbCorrectDIA21=probDIAC1[xgrid>M1+.02 & xgrid<M2-.02]
    ProbCorrectDIA31=probDIAC1[xgrid>M2+.02]
  }
  
  
  if(is.na(DIASet2[1])==FALSE){
    a1=data.frame(grid=c(grid1,grid1,grid1),probC=c(probMIC1,ProbCorrectDIA1,ProbCorrectDIA11),
                  group=c(rep(0,length(grid1)),rep(1,length(grid1)),rep(2,length(grid1))))
    names(a1)=c('grid','Correct','Breakpoints')
    a2=data.frame(grid=c(grid2,grid2,grid2),probC=c(probMIC2,ProbCorrectDIA2,ProbCorrectDIA21),
                  group=c(rep(0,length(grid2)),rep(1,length(grid2)),rep(2,length(grid2))))
    names(a2)=c('grid','Correct','Breakpoints')
    a3=data.frame(grid=c(grid3,grid3,grid3),probC=c(probMIC3,ProbCorrectDIA3,ProbCorrectDIA31),
                  group=c(rep(0,length(grid3)),rep(1,length(grid3)),rep(2,length(grid3))))
    names(a3)=c('grid','Correct','Breakpoints')
    a1=melt(a1,id=c(1,3))
    a1[a1$Breakpoints==0,2]=paste('MIC ')
    a1[a1$Breakpoints==1,2]=paste('DIA Set 1: ',DIASet1[1],', ',DIASet1[2],sep='')
    a1[a1$Breakpoints==2,2]=paste('DIA Set 2: ',DIASet2[1],', ',DIASet2[2],sep='')
    a2=melt(a2,id=c(1,3))
    a2[a2$Breakpoints==0,2]=paste('MIC ')
    a2[a2$Breakpoints==1,2]=paste('DIA Set 1: ',DIASet1[1],', ',DIASet1[2],sep='')
    a2[a2$Breakpoints==2,2]=paste('DIA Set 2: ',DIASet2[1],', ',DIASet2[2],sep='')
    a3=melt(a3,id=c(1,3))
    a3[a3$Breakpoints==0,2]=paste('MIC ')
    a3[a3$Breakpoints==1,2]=paste('DIA Set 1: ',DIASet1[1],', ',DIASet1[2],sep='')
    a3[a3$Breakpoints==2,2]=paste('DIA Set 2: ',DIASet2[1],', ',DIASet2[2],sep='')
    
    a1$Breakpoints=factor(a1$Breakpoints,levels=c("MIC ",paste('DIA Set 1: ',DIASet1[1],', ',DIASet1[2],sep=''),
          paste('DIA Set 2: ',DIASet2[1],', ',DIASet2[2],sep='')))
    a2$Breakpoints=factor(a2$Breakpoints,levels=c("MIC ",paste('DIA Set 1: ',DIASet1[1],', ',DIASet1[2],sep=''),
          paste('DIA Set 2: ',DIASet2[1],', ',DIASet2[2],sep='')))
    a3$Breakpoints=factor(a3$Breakpoints,levels=c("MIC ",paste('DIA Set 1: ',DIASet1[1],', ',DIASet1[2],sep=''),
          paste('DIA Set 2: ',DIASet2[1],', ',DIASet2[2],sep='')))
    
    if(logConvert==TRUE){
      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_line(data=a3,aes(x=grid,y=value,color=Breakpoints))+  
        geom_vline(xintercept=M1,lty=2,alpha=.5)+
        geom_vline(xintercept=M2,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(MIC),max(MIC)))

    }
    if(logConvert==FALSE){
      M1=2^(M1)
      M2=2^(M2)
      a1$grid=2^a1$grid
      a2$grid=2^a2$grid
      a3$grid=2^a3$grid
      MICTemp=min(MIC):max(MIC)
      MICTemp=2^MICTemp
      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_line(data=a3,aes(x=grid,y=value,color=Breakpoints))+  
        geom_vline(xintercept=M1,lty=2,alpha=.5)+
        geom_vline(xintercept=M2,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(trans=log2_trans(),limits=c(min(MICTemp),max(MICTemp)),breaks=MICTemp)
    }
  }
  if(is.na(DIASet2[1])==TRUE){
    a1=data.frame(grid=c(grid1,grid1),probC=c(probMIC1,ProbCorrectDIA1),
                  group=c(rep(0,length(grid1)),rep(1,length(grid1))))
    names(a1)=c('grid','Correct','Breakpoints')
    a2=data.frame(grid=c(grid2,grid2),probC=c(probMIC2,ProbCorrectDIA2),
                  group=c(rep(0,length(grid2)),rep(1,length(grid2))))
    names(a2)=c('grid','Correct','Breakpoints')
    a3=data.frame(grid=c(grid3,grid3),probC=c(probMIC3,ProbCorrectDIA3),
                  group=c(rep(0,length(grid3)),rep(1,length(grid3))))
    names(a3)=c('grid','Correct','Breakpoints')
    a1=melt(a1,id=c(1,3))
    a1[a1$Breakpoints==1,2]=paste('DIA: ',DIASet1[1],', ',DIASet1[2],sep='')
    a1[a1$Breakpoints==0,2]=paste('MIC ')
    a2=melt(a2,id=c(1,3))
    a2[a2$Breakpoints==1,2]=paste('DIA: ',DIASet1[1],', ',DIASet1[2],sep='')
    a2[a2$Breakpoints==0,2]=paste('MIC ')
    a3=melt(a3,id=c(1,3))
    a3[a3$Breakpoints==1,2]=paste('DIA: ',DIASet1[1],', ',DIASet1[2],sep='')
    a3[a3$Breakpoints==0,2]=paste('MIC ')

    a1$Breakpoints=factor(a1$Breakpoints,levels=c("MIC ",paste('DIA: ',DIASet1[1],', ',DIASet1[2],sep='')))
    a2$Breakpoints=factor(a2$Breakpoints,levels=c("MIC ",paste('DIA: ',DIASet1[1],', ',DIASet1[2],sep='')))
    a3$Breakpoints=factor(a3$Breakpoints,levels=c("MIC ",paste('DIA: ',DIASet1[1],', ',DIASet1[2],sep='')))
    
    if(logConvert==TRUE){
      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_line(data=a3,aes(x=grid,y=value,color=Breakpoints))+  
        geom_vline(xintercept=M1,lty=2,alpha=.5)+
        geom_vline(xintercept=M2,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(MIC),max(MIC)))+
        theme(
          plot.title = element_text(size = 17, hjust=0.5,vjust=1),
          legend.position='bottom',
          legend.key=element_rect(fill="white",colour="white"),
          legend.text = element_text(size = 16),
          legend.title=element_text(size=15))
    }
    if(logConvert==FALSE){
      M1=2^(M1)
      M2=2^(M2)
      a1$grid=2^a1$grid
      a2$grid=2^a2$grid
      a3$grid=2^a3$grid
      MICTemp=min(MIC):max(MIC)
      MICTemp=2^MICTemp
      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_line(data=a3,aes(x=grid,y=value,color=Breakpoints))+  
        geom_vline(xintercept=M1,lty=2,alpha=.5)+
        geom_vline(xintercept=M2,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(trans=log2_trans(),
                           limits=c(min(MICTemp),max(MICTemp)),
                           breaks=MICTemp)+
        theme(
          plot.title = element_text(size = 17, hjust=0.5,vjust=1),
          legend.position='bottom',
          legend.key=element_rect(fill="white",colour="white"),
          legend.text = element_text(size = 16),
          legend.title=element_text(size=15))
    }
  }
  
  
  return(fit)
  
  
}

findDIAC=function(yobs,gridx,weights,fit,M1,M2,xsig,ysig,minWidth,maxWidth,minDIA,maxDIA){
  D1=0; D2=0; index=0
  lgrid=length(gridx)
  MICLowerObs=M1
  MICLowerTrue=M1-.5
  MICUpperObs=M2
  MICUpperTrue=M2-.5
  if(minDIA<min(yobs)) minDIA=min(yobs)
  if(maxDIA>max(yobs)) maxDIA=max(yobs)
  storage.mode(gridx) <- "double"
  storage.mode(weights) <- "double"
  storage.mode(fit) <- "double"
  storage.mode(MICLowerObs) <- "double"
  storage.mode(MICUpperObs) <- "double"
  storage.mode(MICLowerTrue) <- "double"
  storage.mode(MICUpperTrue) <- "double"
  storage.mode(xsig) <- "double"
  storage.mode(ysig) <- "double"
  storage.mode(minDIA) <- "double"
  storage.mode(maxDIA) <- "double"
  storage.mode(lgrid) <- "integer"
  storage.mode(D1) <- "double"
  storage.mode(D2) <- "double"
  storage.mode(index) <- "double"
  storage.mode(minWidth) <- "integer"
  storage.mode(maxWidth) <- "integer"
  temp=.C("findDIATrue",gridx,weights,fit,MICLowerObs,MICUpperObs,MICLowerTrue,MICUpperTrue,
          xsig,ysig,minDIA,maxDIA,lgrid,D1,D2,index,minWidth,maxWidth)
  D1=temp[[13]]
  D2=temp[[14]]
  index=temp[[15]]
  #   print(c(D1,D2,index))
  
  return(list(D1=D1,D2=D2,index=index))
}
