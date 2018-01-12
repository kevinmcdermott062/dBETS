getDIABrkptsModel_twoMICShiny=function(MICDens,gx,xgrid,DIA,MICBrkptL,MICBrkptU,xsig=.707,ysig=2.121,minWidth=3,
                                  maxWidth=12,minDIA=min(DIA)+2,maxDIA=max(DIA)-2){
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Calculating DIA Breakpoints", value = 0)
  
  DIA_Brkpts=matrix(NA,nrow=nrow(MICDens),ncol=2)
  for(i in 1:nrow(MICDens)){
    if(i %% 200==0) progress$inc(200/nrow(MICDens))
    parms=findDIAC(DIA,xgrid,MICDens[i,],gx[i,],MICBrkptL,MICBrkptU,xsig,ysig,minWidth,maxWidth,minDIA,maxDIA)
    DIA_Brkpts[i,1]=parms$D1
    DIA_Brkpts[i,2]=parms$D2
  }
  a1=as.data.frame(table(DIA_Brkpts[,1],DIA_Brkpts[,2]))
  names(a1)=c('DIA_L','DIA_U','Freq')
  a1 = a1 %>% arrange(desc(Freq)) %>% mutate(Percent=format(round(Freq/sum(Freq)*100),nsmall=2),
                                             CumPerc=format(round(cumsum(Freq)/sum(Freq)*100),nsmall=2)) %>%
    dplyr::select(-Freq)
  
  return(a1)
}

getDIABrkptsModel_oneMICShiny=function(MICDens,gx,xgrid,DIA,MICBrkpt,xsig=.707,ysig=2.121){
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Calculating DIA Breakpoints", value = 0)
  
  DIA_Brkpt=rep(NA,nrow=nrow(MICDens))
  for(i in 1:nrow(MICDens)){
    if(i %% 200==0) progress$inc(200/nrow(MICDens))
    parms=findDIACOne(DIA,xgrid,MICDens[i,],gx[i,],MICBrkpt,xsig,ysig)
    DIA_Brkpt[i]=parms$DIABrkpt
  }
  a1=as.data.frame(table(DIA_Brkpt))
  names(a1)=c('DIA','Freq')
  a1 = a1 %>% arrange(desc(Freq)) %>% mutate(Percent=format(round(Freq/sum(Freq)*100),nsmall=2),
                                             CumPerc=format(round(cumsum(Freq)/sum(Freq)*100),nsmall=2)) %>%
    dplyr::select(-Freq)
  
  return(a1)
}