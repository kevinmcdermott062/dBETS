

runSpline=function(MIC,DIA,xcens,ycens,MICBrkptL,MICBrkptU,minWidth,maxWidth,session){
  
  source('Spline/splineFunctions.R')
  
  ### Initialize
  xobs=MIC
  yobs=DIA
  M1=MICBrkptL
  M2=MICBrkptU
  xsig=.707; ysig=2.121
  nobs=length(xobs)
  xtrue1 = rnorm(length(xobs),xobs-.5,xsig)
  
  ### Initalize Isplines
  lowept=min(xobs)-.5
  upperept=max(xobs)+.5
  if(sum(xcens==1)>0 | sum(ycens==1)>0)
    upperept=max(xobs)+4
  if(sum(xcens==-1)>0 | sum(ycens==-1)>0)
    lowept=min(xobs)-4
  dist=1
  parms=Ispline(seq(min(xobs),max(xobs),by=dist),lowept,upperept)
  bases=parms$bases;
  knotseq=parms$knotseq
  
  #### Initial spline coefficients
  designMatrix=getIsplineC(xtrue1,knotseq,bases)
  icoefs1=findInitialSpline(xtrue1,bases,knotseq,yobs,designMatrix)
  ytrue1=icoefs1%*%t(designMatrix)
  
  ###Update Dens
  group=rep(0,length(xobs))
  parms=updateDens(xtrue1,group,nIter=500,burnin=400)
  dens=parms$post; group=parms$group
  
  numIter=10000
  burnin=7000
  coefMat=matrix(nrow=numIter,ncol=length(icoefs1))
  fitMat <<- matrix(nrow=numIter-burnin,ncol=length(xgrid))
  MICDens <<- matrix(nrow=numIter-burnin,ncol=length(xgrid))
#   sigmaM=rep(NA,numIter)
#   sigmaD=rep(NA,numIter)
  smoothSpline = .5
  
#   progressInit()
  withProgress(message = 'Running Model...', min=1, max=numIter, {
  
    for(iter in 1:numIter){
      
      if(iter%%1000==0) setProgress(value = iter)
          
      ### Update Density
      if(iter %% 500==0){
        parms=updateDens(xtrue1,group)
        dens=parms$post; group=parms$group
      }
      
      ### Update xtrue (and corresponding ytrue)
      if(iter %% 15==0 | iter==1){
        parms=updateSplineXtrue(xtrue1,ytrue1,knotseq,bases,lowept,upperept,icoefs1,xcens,ycens,xsig,ysig,xobs,yobs,nobs,dens)
        xtrue1=parms$xtrue1; ytrue1=as.numeric(parms$ytrue1)
      }
      
      
      ### Update Coefs
      parms=updateSplineCoefs(smoothSpline,icoefs1,scaleSig,iter,Sigma,bases,xtrue1,
                              ytrue1,knotseq,xobs,yobs,xcens,ycens,xsig,ysig,coefMat)
      icoefs1=parms$coefs; ytrue1=as.numeric(parms$ytrue1)
      
      ### update smoothness parameter for the spline coefs
      parms=updateSmoothParm(smoothSpline,icoefs1)
      smoothSpline=parms$smooth
      
      
      ### Update Measurement Errors
      #   if(iter %% 20==0 | iter==1){
      #     parms=updateMeasErrorM(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens)
      #     xsig=parms$xsig
      #     parms=updateMeasErrorD(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens)
      #     ysig=parms$ysig
      #   }
      
      
      ###save
      coefMat[iter,]=log(icoefs1)
  #     sigmaM[iter]=xsig
  #     sigmaD[iter]=ysig
      if(iter>burnin){    
        designMatrixGrid=getIsplineC(xgrid,knotseq,bases)
        fit=icoefs1%*%t(designMatrixGrid)
        fitMat[iter-burnin,] <<- fit
        MICDens[iter-burnin,] <<- dens
      }  
    }
  })
  
  ### Thin    
  MICDens=MICDens[seq(1,nrow(MICDens),by=5),]
  fitMat=fitMat[seq(1,nrow(fitMat),by=5),]
  
  ### Compute DIA Breakpoints
 withProgress(message = 'Calculating DIA Breakpoints...', min=1, max=nrow(MICDens), {
    DIABrkpts <<- matrix(nrow=nrow(MICDens),ncol=3)
    D1=min(yobs); D2=max(yobs)
    for(i in 1:nrow(MICDens)){
      if(i %% 50==0) setProgress(value = i)
      temp=findDIAC(yobs,xgrid,MICDens[i,],fitMat[i,],M1,M2,xsig,ysig,minWidth,maxWidth,minDIA=D1,maxDIA=D2)
      DIABrkpts[i,] <<- c(temp$D1,temp$D2,temp$index)
      D1=temp$D1-5; D2=temp$D2+5
    }
  })
  
  DIABrkptL=DIABrkpts[,1]
  DIABrkptU=DIABrkpts[,2]
  
  tab=table(DIABrkptL,DIABrkptU)*100
  a1=tab/margin.table(tab)*100
  a1=as.data.frame(a1)
  a1=a1[order(a1$Freq,decreasing=T),]
  a1=a1[a1$Freq!=0,]
  a2=a1
  a2$CumFreq=cumsum(a2$Freq)
  a2[,3]=round(a2[,3],2)
  a2[,4]=round(a2[,4],2)
  D1 = a2[1,1]; D2=a2[1,2]
  
#   cat('MIC Measurement Error: 0.707 \n')
#   cat('DIA Measurement Error: 2.121 \n \n')
  cat('-------Optimal DIA Breakpoints--------\n')
  temp=data.frame(DIABrkptL=a2[,1],DIABrkptU=a2[,2],Percent=a2[,3],Cumulative=a2[,4])
  temp[,1:4] = apply(temp[,1:4], 2, function(x) as.character(x));
  name.width <- max(sapply(names(temp), nchar))
  names(temp) <- format(names(temp), width = name.width, justify = "centre")
  print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
  
  medianDensity=apply(MICDens,2,median)
  medFit=apply(fitMat,2,median)
  lowerFit=apply(fitMat,2,function(x) quantile(x,probs=.025))
  upperFit=apply(fitMat,2,function(x) quantile(x,probs=.975))
  
  D1=as.numeric(as.character(a2[1,1]))
  D2=as.numeric(as.character(a2[1,2]))
  
  return(list(medianDensity=medianDensity,medFit=medFit,lowerFit=lowerFit,
              upperFit=upperFit,a1=a1,D1=D1,D2=D2))
}

  
