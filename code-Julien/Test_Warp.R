Fonction_WarpT = function(w,a,b,n){
  L=rep(0,n)
  for(i in 1:length(a)){
    L=L+n*w[i]*pbeta(seq(1/n,1,length.out=n),a[i],b[i])
  }
  return(L)
}

BetaT = function(x,r,ref){
  n=(length(x)+1)/3
  w=x[paste0("w",1:(n-1))]
  a=x[paste0("a",1:n)]^2
  b=x[paste0("b",1:n)]^2
  
  L=0
  
  w=exp(w)/(1+sum(exp(w)))
  w=c(w, 1-w)
  
  #print(c(a,b))
  #print(Fonction_Warp(w,a,b,length(r)))
  
  f = approx(1:length(r),r,Fonction_WarpT(w,a,b,length(r)))$y
  
  L= sum((f-ref)^2)
  
  return(L)
}


Julien = function(v,Ref){
  init = c(a1=1,a2=1,b1=1,b2=1,w1=0.5)
  nb_melange = 2
  
  param = optim(init, BetaT, r=v, ref=Ref)$par
  
  w = exp(param["w1"])/(1+sum(exp(param["w1"])))
  
  w=c(w, 1-sum(w))
  
  a=param[paste0("a",1:nb_melange)]
  b=param[paste0("b",1:nb_melange)]
  
  Transforme = fonction_spline(Fonction_WarpT(w,a,b,length(v)),1:length(v),v)[,1]
  
  return(Transforme)
  
  
}

v=c(0.9, 0.7, 0.85, 0.6, 0.7, 0.4, 0.1)
ref  =c(0.8, 0.8, 0.6, 0.7, 0.55, 0.3, 0.05)

plot(ref,type = "l")
  
lines(v,col="orange")

lines(Julien(v,ref),col="blue")
  