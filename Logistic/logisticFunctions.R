
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

updateXtrueLog=function(xtrue1,ytrue1,icoefs1,xcens,ycens,xsig,ysig,xobs,yobs,nobs,dens){
  
  pxtrue1 = xtrue1 + rnorm(nobs,0,0.5)
  pytrue1=getylogtrue1(icoefs1,pxtrue1)
  llike = getPOPLLK(pxtrue1,dens) +
    getLLK(xobs,yobs,pxtrue1,pytrue1,xcens,ycens,xsig,ysig) -
    getPOPLLK(xtrue1,dens) -
    getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig)
  llike[is.na(llike)]=-99
  mhrat = runif(nobs)
  cond1=log(mhrat)<llike
  xtrue1[cond1]=pxtrue1[cond1]
  ytrue1=getylogtrue1(icoefs1,xtrue1)
  
  return(list(xtrue1=xtrue1,ytrue1=ytrue1))
  
}

getylogtrue1=function(coef,xtrue){
  
  mb = (2*coef[3]*coef[4])/(coef[3]+coef[4])
  fx = 1/(1+exp(-mb*(coef[2]-xtrue)))
  ytrue=coef[1]*(fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*exp(coef[4]*
        (coef[2]-xtrue)))/(1+fx*exp(coef[3]*(coef[2]-xtrue))+(1-fx)*
         exp(coef[4]*(coef[2]-xtrue)))
  
  return(ytrue)  
}

updateLogCoef=function(iter,icoefs1,xobs,yobs,xtrue1,xcens,ycens,ytrue1,xsig,ysig,coefMat){
  
  
  if(iter<1000){
    Sigma=diag(.001,nrow=ncol(coefMat),ncol=ncol(coefMat))
  }else{
    Sigma=cov(coefMat[(iter-900):(iter-1),])
  }
  
  icoefs1a=mvrnorm(n=1,mu=icoefs1,Sigma=Sigma)
  pytrue1 = getylogtrue1(icoefs1a,xtrue1)
  llike = sum(getLLK(xobs,yobs,xtrue1,pytrue1,xcens,ycens,xsig,ysig)) -
    sum(getLLK(xobs,yobs,xtrue1,ytrue1,xcens,ycens,xsig,ysig))
  accept=0
  if(is.na(llike)) llike=-99
  #Acceptance and make sure coeffecients are positive
  if(log(runif(1)) < llike & length(which(icoefs1a>0))==length(icoefs1a))
  {
    accept=1
    icoefs1 = icoefs1a
    ytrue1 = pytrue1
  }
  
  return(list(accept=accept,coefs=icoefs1,ytrue1=ytrue1))
}


# updateMeasErrorM=function(xsig,ysig,xobs,yobs,xtrue1,ytrue1,xcens,ycens){
#   
#   propXsig=rnorm(1,xsig,.2)
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



