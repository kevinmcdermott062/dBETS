
getLLK=function(xobs,yobs,xtrue,ytrue,xcens,ycens,xsig,ysig){
  
  nobs=length(xobs)
  onevec <- c(rep(1,nobs))
  xcensl=rep(0,nobs)
  xcensl[xcens==-1] = 1
  xcensu=rep(0,nobs)
  xcensu[xcens==1] = 1
  
  pmic <- pnorm(xobs,xtrue,xsig)*(1-xcensu) + onevec*xcensu - 
    pnorm(xobs-1,xtrue,xsig)*(1-xcensl)
  pmic[pmic<.001]=.0001
  
  ycensl=rep(0,nobs)
  ycensl[ycens==-1] = 1
  ycensu=rep(0,nobs)
  ycensu[ycens==1] = 1
  
  pdia <- pnorm(yobs+0.5,ytrue,ysig)*(1-ycensu) + onevec*ycensu -
    pnorm(yobs-0.5,ytrue,ysig)*(1-ycensl)
  pdia[pdia<.001]=.0001
  
  llike <- log(pmic) + log(pdia)
  
  return(llike)
}


getPOPLLK=function(x,dens){
  ###interpolate  
  pxden=log(approx(xgrid,dens,xout=x,rule=2)$y)
  return(pxden)
}

getPRIORLLK=function(coefs,smoothSpline){
  
  priorPcoef=dunif(smoothSpline,0,5,log=T)
  
  lenCoef=length(coefs)
  lpriorcoef = dnorm(coefs[lenCoef],1,10,log=T)
  for(i in (lenCoef-1):1)
    lpriorcoef=lpriorcoef+dnorm(coefs[i],coefs[i+1],smoothSpline,log=T)
  
  return(priorPcoef+lpriorcoef)
}


findInitialSpline=function(xtrue,bases,knotseq,yobs,designMatrix){
  
  min=999999
  for(i in seq(.5,7,by=.1)){
    icoefs1=seq(.1,i,length=ncol(designMatrix))
    ytrue1=icoefs1%*%t(designMatrix)
    if(sum(abs(yobs-ytrue1))<min){ min=sum(abs(yobs-ytrue1)); save=i}
  }
  coef=seq(.1,save,length=ncol(designMatrix))
  return(coef)
}


Ispline=function(intknots,lowept,upperept){
  k=3
  #determine knot sequence
  knotseq=c(rep(lowept,k+1),intknots,rep(upperept,k+1))
  numbases=length(knotseq)-k-2
  
  #create matrix of bases
  bases=matrix(NA,nrow=numbases,ncol=2)
  for(i in 1:numbases) bases[i,]=c(knotseq[i+1],knotseq[i+k+1])
  
  return(list(bases=bases,knotseq=knotseq))
}


getIsplineC=function(xtrue,knotseq,bases){
  
  numBases=nrow(bases)
  lxtrue=length(xtrue)
  
  mat=rep(0,numBases*lxtrue)
  
  storage.mode(mat) <- "double"
  storage.mode(knotseq) <- "double"
  storage.mode(xtrue) <- "double"
  storage.mode(numBases) <- "integer"
  storage.mode(lxtrue) <- "integer"
  temp=.C("getIspline",xtrue,lxtrue,knotseq,mat,numBases)
  designMatrix=matrix(temp[[4]],ncol=numBases,nrow=lxtrue)
  return(designMatrix)
  
}


updateDens=function(x,group,nIter=150,burnin=100){
  alpha=1
  tau0=.1
  beta0=.1
  mu0=0
  kappa0=.1
  C=group
  N=length(x)
  NM=500
  tab=table(C)
  m=rep(0,500)
  for(i in 1:length(tab))
    m[i]=tab[i]
  densPosterior=rep(0,length(xgrid)*(nIter-burnin))
  
  storage.mode(x) <- "double"
  storage.mode(alpha) <- "double"
  storage.mode(tau0) <- "double"
  storage.mode(beta0) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(kappa0) <- "double"
  storage.mode(C) <- "integer"
  storage.mode(N) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(nIter) <- "integer"
  storage.mode(NM) <- "integer"
  storage.mode(xgrid) <- "double"
  storage.mode(densPosterior) <- "double"
  storage.mode(burnin) <- "integer"
  
  
  temp=.C("runDP",x,alpha,tau0,beta0,mu0,kappa0,C,N,m,nIter,NM,xgrid,densPosterior,burnin)
  posterior=temp[[13]]
  C=temp[[7]]
  ###Rescale C
  C=match(C, sort(unique(C)))-1
  
  post=matrix(nrow=(nIter-burnin),ncol=1200)
  idxSeq=seq(1,(nIter-burnin)*1200,by=1200)
  idx=1
  for(i in 1:(nIter-burnin)){
    post[i,]=posterior[idx:(idx+1199)]
    idx=idx+1200
  }
  
  post=apply(post,2,median)
  post=post/sum(post)
  
  return(list(post=post,group=C))
}



updateSplineXtrue=function(xtrue1,ytrue1,knotseq,bases,lowept,upperept,icoefs1,xcens,ycens,
         xsig,ysig,xobs,yobs,nobs,dens){
  
  pxtrue1 = xtrue1 + rnorm(nobs,0,.5)
  pxtrue1[pxtrue1<lowept] = lowept+.01
  pxtrue1[pxtrue1>upperept] = upperept-.01
  designMatrix=getIsplineC(pxtrue1,knotseq,bases)
  pytrue1=as.numeric(icoefs1%*%t(designMatrix))
  llike = getPOPLLK(pxtrue1,dens) +
    getLLK(xobs,yobs,pxtrue1,pytrue1,xcens,ycens,xsig,ysig) -
    getPOPLLK(xtrue1,dens) -
    getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig)
  llike[is.na(llike)]=-99
  mhrat = runif(nobs)
  cond1=log(mhrat)<llike
  xtrue1[cond1]=pxtrue1[cond1]
  designMatrix=getIsplineC(xtrue1,knotseq,bases)
  ytrue1=icoefs1%*%t(designMatrix)  
  
  return(list(xtrue1=xtrue1,ytrue1=ytrue1))
  
}


updateSplineCoefs=function(smoothSpline,icoefs1,scaleSig,iter,Sigma,bases,xtrue1,ytrue1,
        knotseq,xobs,yobs,xcens,ycens,xsig,ysig,coefMat){
  
  if(iter<1000){
    Sigma=diag(.001,nrow=ncol(coefMat),ncol=ncol(coefMat))
  }else{
    Sigma=cov(coefMat[(iter-900):(iter-1),])
  }
  
  propCoef=tryCatch(as.numeric(exp(rmvnorm(n=1,as.numeric(log(icoefs1)),sigma=Sigma))),
        error=function(e){as.numeric(exp(rmvnorm(n=1,as.numeric(log(icoefs1)),
        sigma=diag(.001,nrow=nrow(bases),ncol=nrow(bases)))))})
  
  accept=0
  designMatrix=getIsplineC(xtrue1,knotseq,bases)
  newY=propCoef%*%t(designMatrix)
  
  #log likelihood and priors
  oldLLK=sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig))
  oldPRIORLLK=sum(getPRIORLLK(icoefs1,smoothSpline))
  newLLK=sum(getLLK(xobs,yobs,xtrue1,newY,xcens,ycens,xsig,ysig))
  newPRIORLLK=sum(getPRIORLLK(propCoef,smoothSpline))
  
  #accept/reject
  A=newLLK+newPRIORLLK-oldLLK-oldPRIORLLK
  if(is.nan(A)==FALSE & log(runif(1)) < A){
    ytrue1=newY
    icoefs1=propCoef
    accept=1
  }
  
  return(list(acceptCoef=accept,coefs=icoefs1,ytrue1=ytrue1))
  
}

updateSmoothParm=function(smoothSpline,icoefs1){
  
  accept=0
  propSmooth=rnorm(1,smoothSpline,.15)
  if(propSmooth<0 | propSmooth>5)
    return(list(smooth=smoothSpline,accept=accept))
  oldPRIORLLK=sum(getPRIORLLK(icoefs1,smoothSpline))
  newPRIORLLK=sum(getPRIORLLK(icoefs1,propSmooth))
  #accept/reject
  if(log(runif(1)) < newPRIORLLK-oldPRIORLLK){
    smoothSpline=propSmooth
    accept=1
  }
  
  return(list(smooth=smoothSpline,accept=accept))
}

# updateMeasErrorM=function(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens){
#   
#   propXsig=rnorm(1,xsig,.07)
#   accept=0
#   if(propXsig<0)
#     return(list(xsig=xsig))
#   
#   oldLLK=sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig))
#   newLLK=sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,propXsig,ysig))
#   
#   if(log(runif(1)) < newLLK-oldLLK)
#     xsig=propXsig
#   
#   return(list(xsig=xsig))
# }
# 
# updateMeasErrorD=function(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens){
#   
#   propYsig=rnorm(1,ysig,.25)
#   if(propYsig<0)
#     return(list(ysig=ysig))
#   
#   oldLLK=sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig))
#   newLLK=sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,propYsig))
#   
#   if(log(runif(1)) < newLLK-oldLLK)
#     ysig=propYsig
#   
#   return(list(ysig=ysig))
# }

