
basicPlotOne=function(MIC,DIA,xcens,ycens,MICBrkpt,flipGraph,logConvert){
    
  MIC1=MIC
  DIA1=DIA
  MIC[xcens==1 & MIC==max(MIC)]=max(MIC)+1
  MIC[xcens==-1 & MIC==min(MIC)]=min(MIC)-1
  DIA[ycens==1 & DIA==max(DIA)]=max(DIA)+1
  DIA[ycens==-1 & DIA==min(DIA)]=min(DIA)-1
  a1=data.frame(MIC,DIA)
  a1=count(a1,c('MIC','DIA'))
  Freq=a1$freq; MIC2=a1$MIC; DIA2=a1$DIA
  a1=data.frame(MIC2,DIA2,Freq)  
  
  if(flipGraph=='Yes' && logConvert==TRUE){
    fit=ggplot(a1,aes(MIC2,DIA2))+geom_text(aes(label=Freq),size=4)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in log(ug/mL))',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(breaks = seq(min(MIC1)-1,max(MIC1)+1,by=1),
         labels = c(paste("<",min(MIC1),sep=''),seq(min(MIC1),max(MIC1),by=1), paste(">",max(MIC1),sep='')),
           limits = c(min(MIC1)-1,max(MIC1)+1))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
           labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
           limits = c(min(DIA1)-1,max(DIA1)+1))
  }
  if(flipGraph=='Yes' && logConvert==FALSE){
    MICBrkpt=2^MICBrkpt
    a1$MIC2=2^a1$MIC2
    MICTemp=c(min(MIC1)-1,min(MIC1):max(MIC1),max(MIC1)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC1):max(MIC1))
    
    fit=ggplot(a1,aes(MIC2,DIA2))+geom_text(aes(label=Freq),size=4)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in ug/mL)',y='DIA (Diffusion Test in mm)')+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
           labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
           limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_x_continuous(trans=log2_trans(),
            limits=c(min(MICTemp),max(MICTemp)),
            breaks=MICTemp,
            labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))
  }
  if(flipGraph=='No' && logConvert==TRUE){
    fit=ggplot(a1,aes(DIA2,MIC2))+geom_text(aes(label=Freq),size=4)+
      geom_hline(yintercept=MICBrkpt,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in log(ug/mL))',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(breaks = seq(min(MIC1)-1,max(MIC1)+1,by=1),
               labels = c(paste("<",min(MIC1),sep=''),seq(min(MIC1),max(MIC1),by=1), paste(">",max(MIC1),sep='')),
               limits = c(min(MIC1)-1,max(MIC1)+1))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
               labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
               limits = c(min(DIA1)-1,max(DIA1)+1))
  }
  if(flipGraph=='No' && logConvert==FALSE){
    MICBrkpt=2^MICBrkpt
    a1$MIC2=2^a1$MIC2
    MICTemp=c(min(MIC1)-1,min(MIC1):max(MIC1),max(MIC1)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC1):max(MIC1))

    fit=ggplot(a1,aes(DIA2,MIC2))+geom_text(aes(label=Freq),size=4)+
      geom_hline(yintercept=MICBrkpt,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in ug/mL)',x='DIA (Diffusion Test in mm)')+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
           labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
           limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_y_continuous(trans=log2_trans(),
           limits=c(min(MICTemp),max(MICTemp)),
             breaks=MICTemp,
           labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))
    print(fit)
  }
  
  return(fit)
  
}


basicPlot=function(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,flipGraph,logConvert){
  
  MICBrkpt1=MICBrkptL+.5
  MICBrkpt2=MICBrkptU-.5
  
  MIC1=MIC
  DIA1=DIA
  MIC[xcens==1 & MIC==max(MIC)]=max(MIC)+1
  MIC[xcens==-1 & MIC==min(MIC)]=min(MIC)-1
  DIA[ycens==1 & DIA==max(DIA)]=max(DIA)+1
  DIA[ycens==-1 & DIA==min(DIA)]=min(DIA)-1
  a1=data.frame(MIC,DIA)
  a1=count(a1,c('MIC','DIA'))
  Freq=a1$freq; MIC2=a1$MIC; DIA2=a1$DIA
  a1=data.frame(MIC2,DIA2,Freq)  
  
  if(flipGraph=='Yes' && logConvert==TRUE){
    fit=ggplot(a1,aes(MIC2,DIA2))+geom_text(aes(label=Freq),size=4)+
      geom_vline(xintercept=MICBrkpt1,lty=2,alpha=.5)+
      geom_vline(xintercept=MICBrkpt2,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in log(ug/mL))',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(breaks = seq(min(MIC1)-1,max(MIC1)+1,by=1),
                         labels = c(paste("<",min(MIC1),sep=''),seq(min(MIC1),max(MIC1),by=1), paste(">",max(MIC1),sep='')),
                         limits = c(min(MIC1)-1,max(MIC1)+1))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))

  }
  if(flipGraph=='Yes' && logConvert==FALSE){
    MICBrkpt1=2^MICBrkpt1
    MICBrkpt2=2^MICBrkpt2
    a1$MIC2=2^a1$MIC2
    MICTemp=c(min(MIC1)-1,min(MIC1):max(MIC1),max(MIC1)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC1):max(MIC1))
    
    fit=ggplot(a1,aes(MIC2,DIA2))+geom_text(aes(label=Freq),size=4)+
      geom_vline(xintercept=MICBrkpt1,lty=2,alpha=.5)+
      geom_vline(xintercept=MICBrkpt2,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in ug/mL)',y='DIA (Diffusion Test in mm)')+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                   labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                   limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp,
                         labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))
  }
  if(flipGraph=='No' && logConvert==TRUE){
    fit=ggplot(a1,aes(DIA2,MIC2))+geom_text(aes(label=Freq),size=4)+
      geom_hline(yintercept=MICBrkpt1,lty=2,alpha=.5)+
      geom_hline(yintercept=MICBrkpt2,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in log(ug/mL))',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(breaks = seq(min(MIC1)-1,max(MIC1)+1,by=1),
                         labels = c(paste("<",min(MIC1),sep=''),seq(min(MIC1),max(MIC1),by=1), paste(">",max(MIC1),sep='')),
                         limits = c(min(MIC1)-1,max(MIC1)+1))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))
  }
  if(flipGraph=='No' && logConvert==FALSE){
    MICBrkpt1=2^MICBrkpt1
    MICBrkpt2=2^MICBrkpt2
    a1$MIC2=2^a1$MIC2
    MICTemp=c(min(MIC1)-1,min(MIC1):max(MIC1),max(MIC1)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC1):max(MIC1))
    
    fit=ggplot(a1,aes(DIA2,MIC2))+geom_text(aes(label=Freq),size=4)+
      geom_hline(yintercept=MICBrkpt1,lty=2,alpha=.5)+
      geom_hline(yintercept=MICBrkpt2,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in ug/mL)',x='DIA (Diffusion Test in mm)')+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      scale_y_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp,
                         labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))
    print(fit)
  }
  
  return(fit)
  
}

descriptiveStat=function(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU){
  
  
  #find information for plotting and display information
  N=length(MIC)
  withinOne=rep(0,length(MIC))
  withinOne[which(MIC>=MICBrkptL & MIC<=MICBrkptU)]=1
  numwithinOne=sum(withinOne)
  
  cat('--------------- Data Set Characteristics ---------------\n \n')
  cat('Number of Isolates: ',N,'\n')
  cat('Number of Observed Outside One of Intermediate Range: ',N-numwithinOne,'\n')
  cat('Number of Observed Within One of Intermediate Range: ',numwithinOne,'\n')
  cat('Number of MIC censored: ',sum(xcens!=0),'\n')
  cat('Number of DIA censored: ',sum(ycens!=0),'\n')
  invisible()
}

descriptiveStatOne=function(MIC,DIA,xcens,ycens,MICBrkpt){
  
  
  #find information for plotting and display information
  N=length(MIC)
  
  cat('--------------- Data Set Characteristics ---------------\n \n')
  cat('Number of Isolates: ',N,'\n')
  cat('Number Susceptible: ',sum(MIC<MICBrkpt),'\n')
  cat('Number Resistant: ',sum(MIC>MICBrkpt),'\n')
  cat('Number of MIC censored: ',sum(xcens!=0),'\n')
  cat('Number of DIA censored: ',sum(ycens!=0),'\n')
  invisible()
}

