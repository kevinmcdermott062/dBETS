
findBrkptsERBC=function(MIC,DIA,VM1,M1,m1,VM2,M2,m2,MICBrkpt1,MICBrkpt2,minWidth,maxWidth){
    
  #create weights
  temp=max(VM1,M1,m1,VM2,M2,m2)
  VM1=temp/VM1; M1=temp/M1; m1=temp/m1
  VM2=temp/VM2; M2=temp/M2; m2=temp/m2
  
  D1=0; D2=0
  
  #observations within different ranges
  MICwithinOne=MIC[which(MIC>=(MICBrkpt1) & MIC<=(MICBrkpt2))]
  DIAwithinOne=DIA[which(MIC>=(MICBrkpt1) & MIC<=(MICBrkpt2))]
  N1=length(MICwithinOne)
  MICOutsideOne=MIC[which(MIC<(MICBrkpt1) | MIC>(MICBrkpt2))]
  DIAOutsideOne=DIA[which(MIC<(MICBrkpt1) | MIC>(MICBrkpt2))]
  N2=length(MICOutsideOne)

  minDIA=min(DIA); maxDIA=max(DIA)
  
  storage.mode(minDIA) <- "integer"
  storage.mode(maxDIA) <- "integer"
  storage.mode(MICwithinOne) <- "double"
  storage.mode(MICOutsideOne) <- "double"
  storage.mode(DIAwithinOne) <- "double"
  storage.mode(DIAOutsideOne) <- "double"
  storage.mode(MICBrkpt1) <- "double"
  storage.mode(MICBrkpt2) <- "double"
  storage.mode(N1) <- "integer"
  storage.mode(N2) <- "integer"
  storage.mode(VM1) <- "double"
  storage.mode(M1) <- "double"
  storage.mode(m1) <- "double"
  storage.mode(VM2) <- "double"
  storage.mode(M2) <- "double"
  storage.mode(m2) <- "double"
  storage.mode(minWidth) <- "integer"
  storage.mode(maxWidth) <- "integer"
  storage.mode(D1) <- "integer"
  storage.mode(D2) <- "integer"
  
  temp=.C("ERB",minDIA,maxDIA,MICwithinOne,MICOutsideOne,DIAwithinOne,DIAOutsideOne,MICBrkpt1,MICBrkpt2,N1,N2,
          VM1,M1,m1,VM2,M2,m2,minWidth,maxWidth,D1,D2)
  return(list(D1=temp[[19]],D2=temp[[20]]))
  
}

findBrkptsERBOneC=function(MIC,DIA,VM,M,MICBrkpt){
  
  #create weights
  temp=max(VM,M)
  VM=temp/VM; M=temp/M
  
  DIABrkpt=0; 
  index=0
  N=length(MIC)
  
  minDIA=min(DIA)+.5; maxDIA=max(DIA)-.5
  
  storage.mode(minDIA) <- "double"
  storage.mode(maxDIA) <- "double"
  storage.mode(MIC) <- "double"
  storage.mode(DIA) <- "double"
  storage.mode(MICBrkpt) <- "double"
  storage.mode(N) <- "integer"
  storage.mode(VM) <- "double"
  storage.mode(M) <- "double"
  storage.mode(DIABrkpt) <- "double"
  storage.mode(index) <- "double"
  
  temp=.C("ERBOne",MIC,DIA,MICBrkpt,N,VM,M,DIABrkpt,minDIA,maxDIA,index)
  
  return(list(DIABrkpt=temp[[7]],index=temp[[10]]))
  
}

