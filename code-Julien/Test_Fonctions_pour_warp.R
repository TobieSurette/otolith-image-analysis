
# Nous donne l'axe principal (plus grand rayon)
Trouver_axe_principal = function(noyau, Contour){
  r=0
  axe_principal = matrix(data = c(noyau,0,0),2,2)
  for(i in 1:length(Contour[,1])){
    if(euclidien(n,Contour[i,]) > r){
      r = euclidien(n,Contour[i,])
      axe_principal[,2] = Contour[i,]
      indice_Contour = i
    }
  }
  # Retourne un vecteur contenant en ordre les coordonnées du noyau,
  # les coordonnées du point sur les contour, ainsi que l'indice du
  # point sur le contour
  return(c(axe_principal,indice_Contour))
}

# Fonction linéaire qui rajuste les intensités de gris des rayons pour
# ressembler au rayon de référence.
Intensite = function(r,ref){
  #for(i in 1:length(Transforme[,1])){
  d = data.frame(x = r, y = ref)
  I = lm(y~x, d)$coefficients
  
  R = I[2]*(r)+I[1]
  #}
  return(R)
}


# Fonction qui définit la distortion des coordonnées pour qu'elles nous
# rapportent au rayon de référence
Fonction_Warp = function(a,b,n,w=1){
  L=rep(0,n)
  for(i in 1:length(a)){
    L=L+w[i]*pbeta(seq(0,1,length.out=n),a[i],b[i])
  }
  return((n-1)*L+1)
}


# Fonction à optimiser pour obtenir les paramètres de la fonction de 
# distortion
Beta = function(x,r,ref){
  n=(length(x)+1)/3
  w=x[paste0("w",1:(n-1))]
  a=x[paste0("a",1:n)]^2
  b=x[paste0("b",1:n)]^2
  #m=x[paste0("m",1:n)]
  #k=x[paste0("m",1:n)]
  
  L=0
  
  w=exp(w)/(1+sum(exp(w)))
  w=c(w, 1-sum(w))
  
  f = fonction_spline(1:length(r),r,Fonction_Warp(a,b,length(r),w))[,1]
  
  #Intensite est une fonction linéaire qui rapporte les intensitees de gris
  # de f vers celles de ref
  L= sum((Intensite(f,ref)-ref)^2)
  
  return(L)
}

Transformation = function(noyau, Contour, Sweep, largeur=length(Sweep[,1])){
  axe_p = Trouver_axe_principal(noyau,Contour)[5] #On trouve l'axe principal
  
  Ref = Sweep[axe_p,] #Ref contient les valeurs de l'axe principal
  
  nb_melange = 2
  initG = c(a1=1,a2=1,b1=1,b2=1,w1=0.5)
  initD = c(a1=1,a2=1,b1=1,b2=1,w1=0.5)
  
  paramG = matrix(,length(Sweep[,1])/2,length(initG)+1)
  
  paramD = matrix(,length(Sweep[,1])/2,length(initD)+1)
  
  Transforme = matrix(,length(Sweep[,1]),length(Sweep[1,]))
  Transforme[axe_p,] = Ref
  
  avant = axe_p - 1
  apres = axe_p + 1
  i = 0
  k=1
  
  while(k*avant <= k*apres && i < largeur){
    i=i+1
    
    print(c(i,initG))
    print(c(i,initD))
    
    # Ici on calcul les rayons à gauches du rayon de référence en commençant
    # par celui directement à gauche
    paramG[i,] = c(optim(initG, Beta, r=Sweep[avant,], ref=Ref)$par, avant)
    
    paramD[i,] = c(optim(initD, Beta, r=Sweep[apres,], ref=Ref)$par, apres)
    
    initG[paste0("a",1:nb_melange)] = paramG[i,1:nb_melange]
    initG[paste0("b",1:nb_melange)] = paramG[i,(nb_melange+1):(2*nb_melange)]
    initG[paste0("w",1:(nb_melange-1))] = paramG[i,(2*nb_melange+1):length(initG)]
    
    initD[paste0("a",1:nb_melange)] = paramD[i,1:nb_melange]
    initD[paste0("b",1:nb_melange)] = paramD[i,(nb_melange+1):(2*nb_melange)]
    initD[paste0("w",1:(nb_melange-1))] = paramD[i,(2*nb_melange+1):length(initG)]
    
    w = exp(paramG[i,(2*nb_melange+1):length(initG)])/(1+sum(exp(paramG[i,(2*nb_melange+1):length(initG)])))
    
    w=c(w, 1-sum(w))
    
    Transforme[avant,] = fonction_spline(1:length(Sweep[avant,]),Sweep[avant,],Fonction_Warp(paramG[i,1:nb_melange]^2,paramG[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,]),w))[,1]
    Transforme[apres,] = fonction_spline(1:length(Sweep[apres,]),Sweep[apres,],Fonction_Warp(paramD[i,1:nb_melange]^2,paramD[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,]),w))[,1]
    
    if(avant-1 == 0){
      avant = length(Sweep[,1])
      apres = apres+1
      k=-1
    }else if(apres == length(Sweep[,1])){
      apres = 1
      avant = avant-1
      k=-1
    }else{
      avant = (avant-1)
      apres = (apres+1)
    }
  }
  
  return(list("Transforme" = Transforme, "param" = rbind(paramG, paramD)))
}

TransformationVieux = function(noyau, Contour, Sweep, largeur=length(Sweep[,1])){
  axe_p = Trouver_axe_principal(noyau,Contour)[5] #On trouve l'axe principal
  
  Ref = Sweep[axe_p,] #Ref contient les valeurs de l'axe principal
  
  init = c(a1=1,a2=1,a3=1,a4=1,a5=1,a6=1,a7=1,a8=1,a9=1,a10=1,
           b1=1,b2=1,b3=1,b4=1,b5=1,b6=1,b7=1,b8=1,b9=1,b10=1,
           w1=1,w2=1,w3=1,w4=1,w5=1,w6=1,w7=1,w8=1,w9=1)
  #init = c(a1=1,a2=1,b1=1,b2=1,w1=0.5)
  #init = c(a1=1, b1=1)
  nb_melange = 2
  
  paramG = matrix(,length(Sweep[,1])/2,length(init)+1)
  
  paramD = matrix(,length(Sweep[,1])/2,length(init)+1)
  
  Transforme = matrix(,length(Sweep[,1]),length(Sweep[1,]))
  Transforme[axe_p,] = Ref
  
  refG = avant = axe_p - 1
  refD = apres = axe_p + 1
  i = 0
  k=1
  
  while(k*avant <= k*apres && i < largeur){
    i=i+1
    # Ici on calcul les rayons à gauches du rayon de référence en commençant
    # par celui directement à gauche
    paramG[i,] = c(optim(init, Beta, r=Sweep[avant,], ref=Sweep[refG,])$par, avant)
    
    paramD[i,] = c(optim(init, Beta, r=Sweep[apres,], ref=Sweep[refD,])$par, apres)
    
    w = exp(paramG[i,(2*nb_melange+1):length(init)])/(1+sum(exp(paramG[i,(2*nb_melange+1):length(init)])))
    
    w=c(w, 1-sum(w))
    
    Transforme[avant,] = fonction_spline(1:length(Sweep[avant,]),Sweep[avant,],Fonction_Warp(paramG[i,1:nb_melange]^2,paramG[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,]),w))[,1]
    Transforme[apres,] = fonction_spline(1:length(Sweep[apres,]),Sweep[apres,],Fonction_Warp(paramD[i,1:nb_melange]^2,paramD[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,]),w))[,1]
    if(i>1){
      for(j in 1:(i-1)){
        # f(f(f(f...(transforme)))) = ref
        Transforme[avant,] = fonction_spline(1:length(Transforme[avant,]),Transforme[avant,],Fonction_Warp(paramG[i-j,1:nb_melange]^2,paramG[i-j,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,]),w))[,1]
        if(avant != apres)
          Transforme[apres,] = fonction_spline(1:length(Transforme[apres,]),Transforme[apres,],Fonction_Warp(paramD[i-j,1:nb_melange]^2,paramD[i-j,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,]),w))[,1]
      }
    }
    
    refG = avant
    refD = apres
    if(avant-1 == 0){
      avant = length(Sweep[,1])
      apres = apres+1
      k=-1
    }else if(apres == length(Sweep[,1])){
      apres = 1
      avant = avant-1
      k=-1
    }else{
      avant = (avant-1)
      apres = (apres+1)
    }
  }
  
  return(list("Transforme" = Transforme, "param" = rbind(paramG, paramD)))
}


