library(MASS)

set.seed(2024)

# Calcul de la longueur de chaque intervalles [xi, xi+1]
Calcul_hi = function(x){
  h=rep(0,(length(x)-1))
  for(i in 1:(length(x)-1)){
    h[i]=x[i+1]-x[i]
  }
  return(h)
}

# S'assurer de mettre seulement les x et y (angles et x ou y) appropriés
Difference_divisee = function(x,y){
  n = length(x)
  if(n>2){
    ((Difference_divisee(x[2:n],y[2:n])-Difference_divisee(x[1:(n-1)],y[1:(n-1)]))
    /(x[n]-x[1]))
  }else{
    (y[2]-y[1])/(x[2]-x[1])
  }
}

Calcul_coefs_spline = function(x, y){
  n = length(x)
  
  h = Calcul_hi(x)
  
  f=y[-length(y)] #fi = f(xi)
  
  # M est la matrice des coefficients des fpp
  M = matrix(data = rep(0,n*n),nrow = n, ncol = n)
  M[1,1]=M[n,n] = 1
  
  # B est le vecteur des solutions de M*fpp 
  B = rep(0,n)
  B[1]=B[n] = 0
  
  for(i in seq(2,(n-1))){
    M[i,i-1] = h[i-1]/(h[i]+h[i-1])
    M[i,i] = 2
    M[i,i+1] = h[i]/(h[i]+h[i-1])
    
    B[i] = 6*Difference_divisee(x[(i-1):(i+1)],y[(i-1):(i+1)])
  }
  
  fpp = solve(M)%*%B
  
  fp=c()
  fppp=c()
  for(i in 1:(n-1)){
    fp = c(fp, Difference_divisee(x[i:(i+1)],y[i:(i+1)]) - (h[i]/3)*fpp[i] - (h[i]/6)*fpp[i+1])
    fppp = c(fppp, (fpp[i+1]-fpp[i])/h[i])
  }
  names(f)=paste0("f",1:length(f))
  names(fp)=paste0("fp",1:length(fp))
  names(fpp)=paste0("fpp",1:length(fpp))
  names(fppp)=paste0("fppp",1:length(fppp))
  return(c(f, fp, fpp, fppp))
}


fonction_spline = function(donneesX,donneesY,x){
  n = length(donneesX)
  
  p=c()
  derive_p=c()
  
  coefs = Calcul_coefs_spline(donneesX, donneesY)
  
  f=coefs[paste0("f",1:(n-1))]
  fp=coefs[paste0("fp",1:(n-1))]
  fpp=coefs[paste0("fpp",1:(n))]
  fppp=coefs[paste0("fppp",1:(n-1))]
  
  
  for(j in 1:length(x)){
    
    i = max(which(donneesX<=x[j]))
    if(i == n) {
      i=i-1
    }
    
    p = c(p,f[i] + fp[i]*(x[j]-donneesX[i]) + fpp[i]/2*(x[j]-donneesX[i])^2 + fppp[i]/6*(x[j]-donneesX[i])^3)
      
    derive_p = c(derive_p, fp[i] + fpp[i]*(x[j]-donneesX[i]) + fppp[i]/2*(x[j]-donneesX[i])^2)
    
  }
  return(matrix(data = c(p, derive_p), nrow = length(p), ncol = 2))
}




