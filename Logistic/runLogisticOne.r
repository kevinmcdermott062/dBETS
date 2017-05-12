

runLogisticOne=function(MIC,DIA,xcens,ycens,MICBrkpt,session){
  
  source('Logistic/logisticFunctions.R')
    
  ### Initialize
  xobs=MIC
  yobs=DIA
  xsig=.707; ysig=2.121
  nobs=length(xobs)
  xtrue1 = rnorm(length(xobs),xobs-.5,xsig)
  
  #### Initial Logistic coefficients
  icoefs1 <- numeric(length=4)
  icoefs1[1] <- max(yobs)+1
  options(warn=-1)
  bc <- tryCatch(lm(logit(yobs/icoefs1[1])~xtrue1),
                 error=function(e){'e'})
  if(bc!='e'){
    icoefs1[2] <- -bc$coeff[1]/bc$coeff[2]
    icoefs1[3] <- -bc$coeff[2]
    icoefs1[4] <- -bc$coeff[2]
  }else{
    icoefs1[2:4]=c(1.5,.2,.4)
  }
  options(warn=1)
  ytrue1=getylogtrue1(icoefs1,xtrue1)

  
  ###Update Dens
  group=rep(0,length(xobs))
  parms=updateDens(xtrue1,group,nIter=500,burnin=400)
  dens=parms$post; group=parms$group

  numIter=10000
  burnin=7000
  coefMat=matrix(nrow=numIter,ncol=length(icoefs1))
  fitMat<<-matrix(nrow=numIter-burnin,ncol=length(xgrid))
  MICDens<<-matrix(nrow=numIter-burnin,ncol=length(xgrid))
#   sigmaM=rep(NA,numIter)
#   sigmaD=rep(NA,numIter)
  
#   progressInit()
  withProgress(message = 'Running Model...', min=1, max=numIter, {

    for(iter in 1:numIter){
          
      if(iter%%1000==0) setProgress(value = iter)
      
      ## Update Density
      if(iter %% 500==0){
        parms=updateDens(xtrue1,group)
        dens=parms$post; group=parms$group
      }
      
      ### Update xtrue (and corresponding ytrue)
      if(iter %% 25==0 | iter==1){
        parms=updateXtrueLog(xtrue1,ytrue1,icoefs1,xcens,ycens,xsig,ysig,xobs,yobs,nobs,dens)
        xtrue1=parms$xtrue1; ytrue1=as.numeric(parms$ytrue1)
      }
      
      ####update Coefs
      parms=updateLogCoef(iter,icoefs1,xobs,yobs,xtrue1,xcens,ycens,ytrue1,xsig,ysig,coefMat)
      icoefs1=parms$coefs; ytrue1=parms$ytrue1
      
      
      ### Update Measurement Errors
      #   parms=updateMeasErrorM(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens)
      #   xsig=parms$xsig
      #   parms=updateMeasErrorD(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens)
      #   ysig=parms$ysig
      
      
      ###save
      coefMat[iter,]=icoefs1
  #     sigmaM[iter]=xsig
  #     sigmaD[iter]=ysig
      if(iter>burnin){    
        fit=getylogtrue1(icoefs1,xgrid)
        fitMat[iter-burnin,]<<-fit
        MICDens[iter-burnin,]<<-dens
      }   
    }
  })
      
  
  ### Thin    
  MICDens=MICDens[seq(1,nrow(MICDens),by=5),]
  fitMat=fitMat[seq(1,nrow(fitMat),by=5),]
  
  ### Compute DIA Breakpoints
 withProgress(message = 'Calculating DIA Breakpoints...', min=1, max=nrow(MICDens), {
    DIABrkpt<<-rep(NA,nrow(MICDens))
    for(i in 1:nrow(MICDens)){
      if(i %% 50==0) setProgress(value = i)
      parms=findDIACOne(yobs,xgrid,MICDens[i,],fitMat[i,],MICBrkpt,xsig,ysig)
      DIABrkpt[i] <<- parms$DIABrkpt
    }
  })
  
  
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
  
  cat('-------Optimal DIA Breakpoints--------\n')
  temp=data.frame(DIABrkpt=a2[,1],Percent=a2[,2],Cumulative=a2[,3])
  temp[,1:3] = apply(temp[,1:3], 2, function(x) as.character(x));
  name.width <- max(sapply(names(temp), nchar))
  names(temp) <- format(names(temp), width = name.width, justify = "centre")
  print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
  medianDensity=apply(MICDens,2,median)
  medFit=apply(fitMat,2,median)
  lowerFit=apply(fitMat,2,function(x) quantile(x,probs=.025))
  upperFit=apply(fitMat,2,function(x) quantile(x,probs=.975))
  
  DIABrkpt=as.numeric(as.character(a1[1,1]))

  return(list(medianDensity=medianDensity,medFit=medFit,lowerFit=lowerFit,
              upperFit=upperFit,a1=a1,DIABrkpt=DIABrkpt))
}

