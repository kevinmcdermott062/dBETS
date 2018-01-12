bootStrapERBShiny=function(MIC,DIA,MICBrkptL,MICBrkptU,VM1=10,M1=10,m1=40,VM2=2,M2=2,m2=5,
                      minWidth=3,maxWidth=10){
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Calculating", value = 0)
  
  a1=data.frame(MIC=MIC,DIA=DIA)
  DIABrkptL=rep(NA,5000)
  DIABrkptU=rep(NA,5000)
  n=nrow(a1)
  
  for(i in 1:5000){
    if(i %% 200==0) progress$inc(200/5000)
    tmp=sample_n(a1,n,replace=TRUE)
    parms=findBrkptsERBC(tmp$MIC,tmp$DIA,VM1,M1,m1,VM2,M2,m2,MICBrkptL,MICBrkptU,minWidth,maxWidth)
    DIABrkptL[i]=parms$D1
    DIABrkptU[i]=parms$D2
  }
  
  
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
  
  cat('Bootstrap samples = 5000 \n')
  cat('\n-------DIA Breakpoints by Confidence--------\n')
  temp=data.frame(DIABrkptL=a2[,1],DIABrkptU=a2[,2],Percent=a2[,3],Cumulative=a2[,4])
  temp[,1:4] = apply(temp[,1:4], 2, function(x) as.character(x));
  name.width <- max(sapply(names(temp), nchar))
  names(temp) <- format(names(temp), width = name.width, justify = "centre")
  print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
  
  return(a1)
  
}


bootStrapERBOneShiny=function(MIC,DIA,MICBrkpt,VM=1,M=5){
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Calculating", value = 0)
  
  a1=data.frame(MIC=MIC,DIA=DIA)
  DIABrkpt=rep(NA,5000)
  n=nrow(a1)
  
  for(i in 1:5000){
    if(i %% 200==0) progress$inc(200/5000)
    tmp=sample_n(a1,n,replace=TRUE)
    parms=findBrkptsERBOneC(tmp$MIC,tmp$DIA,VM,M,MICBrkpt)
    DIABrkpt[i]=parms$DIABrkpt
  }
  
  
  #print results
  tab=table(DIABrkpt)
  a1=tab/sum(tab)*100
  a1=as.data.frame(a1)
  a1=a1[order(a1$Freq,decreasing=T),]
  a1=a1[a1$Freq!=0,]
  a2=a1
  a2$CumFreq=cumsum(a2$Freq)
  a2[,2]=round(a2[,2],2)
  a2[,3]=round(a2[,3],2)
  
  cat('Bootstrap samples = 5000 \n')
  cat('\n-------DIA Breakpoints by Confidence--------\n')
  temp=data.frame(DIABrkpt=a2[,1],Percent=a2[,2],Cumulative=a2[,3])
  temp[,1:3] = apply(temp[,1:3], 2, function(x) as.character(x));
  name.width <- max(sapply(names(temp), nchar))
  names(temp) <- format(names(temp), width = name.width, justify = "centre")
  print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
  
  invisible()
  
}
