
# Nous donne l'axe principal (plus grand rayon)
Trouver_axe_principal = function(noyau, Contour){
  r=0
  axe_principal = matrix(data = c(noyau,0,0),2,2)
  for(i in 1:length(Contour[,1])){
    if(euclidien(noyau,Contour[i,]) > r){
      r = euclidien(noyau,Contour[i,])
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
    L=L+w[i]*((n-1)*pbeta(seq(0,1,length.out=n),a[i],b[i])+1)
  }
  return(L)
}

#########################################################################
Inverse_Warp = function (r,a,b,w=1) {
  L = 0
  if(length(a) == 1){
    L = qbeta(r,a,b)
  }else{
    for(i in 1:length(a)){
      L = L + (1/w[i])*qbeta(r,a[i],b[i])
    }
  }
  return(L)
}

##########################################################################

# Fonction à optimiser pour obtenir les paramètres de la fonction de 
# distortion
Beta = function(x,r,ref,n,param){
  w=x[paste0("w",1:(n-1))]
  a=x[paste0("a",1:n)]^2
  b=x[paste0("b",1:n)]^2
  m=x["m1"]
  k=x["k1"]
  
  L=0
  
  if(is.na(sum(w))){
    w=1
  }else{
    w=exp(w)/(1+sum(exp(w)))
    w=c(w, 1-sum(w))
  }

  S = spline(1:length(r),r,xout=Fonction_Warp(a,b,length(r),w), method = "natural")$y
  
  #Intensite est une fonction linéaire qui rapporte les intensitees de gris
  # de f vers celles de ref
  L= sum((m*S+k-ref)^2 + 0*(param[1:n]^2-a)^2 + 0*(param[(n+1):(2*n)]-b)^2 + 
          100*(sum(w - 1/n))^2)

  return(L)
}

Transformation = function(noyau, Contour, Sweep, largeur=length(Sweep[,1])){
  start.time = Sys.time()
  axe_p = Trouver_axe_principal(noyau,Contour)[5] #On trouve l'axe principal
  
  Ref = Sweep[axe_p,] #Ref contient les valeurs de l'axe principal
  
  nb_melange = 5
  initG = c(a1=1, a2=1, a3=1, a4=1, a5=1,
            b1=1, b2=1, b3=1, b4=1, b5=1,
            w1=0.2, w2=0.2, w3=0.2, w4=0.2,
            m1=1,k1=0)
  initD = c(a1=1, a2=1, a3=1, a4=1, a5=1,
            b1=1, b2=1, b3=1, b4=1, b5=1,
            w1=0.2, w2=0.2, w3=0.2, w4=0.2,
            m1=1,k1=0)
  
  # initG = c(a1=1, a2=1, b1=1, b2=1, w1=0.5, m1=1,k1=0)
  # initD = c(a1=1, a2=1, b1=1, b2=1, w1=0.5, m1=1,k1=0)
  
  # initG = c(a1=1, b1=1, m1=1, k1=0)
  # initD = c(a1=1, b1=1, m1=1, k1=0)
  
  
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
    
    # Ici on calcul les rayons à gauches du rayon de référence en commençant
    # par celui directement à gauche
    paramG[i,] = c(optim(initG, Beta, r=Sweep[avant,], ref=Ref, n=nb_melange, param=initG, control = list(reltol = 1e-10))$par, avant)
    
    paramD[i,] = c(optim(initD, Beta, r=Sweep[apres,], ref=Ref, n=nb_melange, param=initD, control = list(reltol = 1e-10))$par, apres)
    
    
    if(nb_melange == 1){
      nom = c("a1", "b1", "m1", "k1")
      initG = c(a1 = paramG[i,1], b1 = paramG[i,2], m1 = paramG[i,3], k1 = paramG[i,4])
      initD = c(a1 = paramD[i,1], b1 = paramD[i,2], m1 = paramD[i,3], k1 = paramD[i,4])
      wG=wD=1
    } else{
      nom = c(paste0("a", 1:nb_melange), paste0("b", 1:nb_melange), 
              paste0("w", 1:(nb_melange-1)), "m1", "k1")
      initG = setNames(paramG[i,-length(paramG[i,])], nom)
      initD = setNames(paramD[i,-length(paramD[i,])], nom)
      
      wG = exp(initG[paste0("w", 1:(nb_melange-1))])/(1+sum(exp(initG[paste0("w", 1:(nb_melange-1))])))
      wG = c(wG, 1-sum(wG))
      
      wD = exp(initD[paste0("w", 1:(nb_melange-1))])/(1+sum(exp(initD[paste0("w", 1:(nb_melange-1))])))
      wD = c(wD, 1-sum(wD))
    }
    
    Transforme[avant,] = paramG[i,ncol(paramG)-2]*spline(1:length(Sweep[avant,]), Sweep[avant,], xout = Fonction_Warp(paramG[i,1:nb_melange]^2,paramG[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,]), wG), method = "natural")$y+paramG[i,ncol(paramG)-1]
    Transforme[apres,] = paramD[i,ncol(paramD)-2]*spline(1:length(Sweep[apres,]), Sweep[apres,], xout = Fonction_Warp(paramD[i,1:nb_melange]^2,paramD[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,]), wD), method = "natural")$y+paramD[i,ncol(paramD)-1]
    
    if(avant-1 == 0){
      avant = length(Sweep[,1])
      apres = apres+1
      k=-1*k
    }else if(apres == length(Sweep[,1])){
      apres = 1
      avant = avant-1
      k=-1*k
    }else{
      avant = (avant-1)
      apres = (apres+1)
    }
  }
  end.time = Sys.time()
  print(end.time-start.time)
  
  P = rbind(paramG, c(rep(1,ncol(paramD)-1), axe_p), paramD)
  colnames(P) = c(nom, "indice")
  P[,paste0("a",1:nb_melange)] = P[,paste0("a",1:nb_melange)]^2
  P[,paste0("b",1:nb_melange)] = P[,paste0("b",1:nb_melange)]^2
  if(nb_melange > 1){
    P[,paste0("w", 1:(nb_melange-1))] = exp(P[,paste0("w", 1:(nb_melange-1))])/(1+sum(exp(P[,paste0("w", 1:(nb_melange-1))])))
  }
  P = P[order(P[,"indice"]),]
  return(list("Transforme" = Transforme, "param" = P[match(unique(P[,"indice"]),P[,"indice"]), ]))
}


# À mettre devant paramG[i,] = c(optim(initG, Beta, r=Sweep[avant,], ref=Ref, n=nb_melange, param=initG, control = list(reltol = 1e-10))$par, avant)
# pour utiliser
#
# eps_init = c(-0.01, 0, 0.01)
# valeur = c(0, 0, 0)
# for(j in 1:3){
#   initG_temp = initG
#   initG_temp[1:2] = initG[1:2] + eps_init[j]
#   valeur[j] = optim(initG_temp, Beta, r=Sweep[avant,], ref=Ref, n=nb_melange, param=initG, control = list(reltol = 1e-10))$value
# }
# initG = initG + eps_init[which(valeur == min(valeur))]
# 
# valeur = c(0, 0, 0)
# for(j in 1:3){
#   initD_temp = initD
#   initD_temp[1:2] = initD[1:2] + eps_init[j]
#   valeur[j] = optim(initD_temp, Beta, r=Sweep[avant,], ref=Ref, n=nb_melange, param=initG, control = list(reltol = 1e-10))$value
# }
# initD = initD + eps_init[which(valeur == min(valeur))]








Warp_spline = function(theta,image){
  x = theta[paste0("x", 1:nrow(image))]
  y = theta[paste0("y", 1:nrow(image))]
  
  S = spline(x, y, xout = x, method = "natural")
  
  dy = as.matrix(imgradient(as.cimg(image),"y")[,,1,1])
  
  print(cbind(floor(S$x),floor(S$y)))
  
  return(sum(dy%*%dy[cbind(floor(S$x),floor(S$y))])^2)
}














TransformationVieux = function(noyau, Contour, Sweep, largeur=length(Sweep[,1])){
  axe_p = Trouver_axe_principal(noyau,Contour)[5] #On trouve l'axe principal

  Ref = Sweep[axe_p,] #Ref contient les valeurs de l'axe principal

  #init = c(a1=1,a2=1,b1=1,b2=1,w1=0.5,m1=1,k1=0)
  init = c(a1=1, b1=1, m1=1, k1=0)
  nb_melange = 1

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
    paramG[i,] = c(optim(init, Beta, r=Sweep[avant,], ref=Sweep[refG,], n = nb_melange)$par, avant)

    paramD[i,] = c(optim(init, Beta, r=Sweep[apres,], ref=Sweep[refD,], n = nb_melange)$par, apres)

    # w = exp(paramG[i,(2*nb_melange+1):length(init)])/(1+sum(exp(paramG[i,(2*nb_melange+1):length(init)])))
    #   
    # w=c(w, 1-sum(w))

    Transforme[avant,] = spline(1:length(Sweep[avant,]), Sweep[avant,], xout = Fonction_Warp(paramG[i,1:nb_melange]^2,paramG[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,])), method = "natural")$y
    Transforme[apres,] = spline(1:length(Sweep[apres,]), Sweep[apres,], xout = Fonction_Warp(paramD[i,1:nb_melange]^2,paramD[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,])), method = "natural")$y
    if(i>1){
      for(j in 1:(i-1)){
        # f(f(f(f...(transforme)))) = ref
        Transforme[avant,] = paramG[i-j,ncol(paramG)-2]*spline(1:length(Transforme[avant,]), Transforme[avant,], xout = Fonction_Warp(paramG[i-j,1:nb_melange]^2,paramG[i-j,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,])), method = "natural")$y + paramG[i-j,ncol(paramG)-1]
        if(avant != apres)
          Transforme[apres,] = paramD[i-j,ncol(paramD)-2]*spline(1:length(Transforme[apres,]), Transforme[apres,], xout = Fonction_Warp(paramD[i-j,1:nb_melange]^2,paramD[i-j,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,])), method = "natural")$y + paramD[i-j,ncol(paramD)-1]
      }
    }

    refG = avant
    refD = apres
    if(avant-1 == 0){
      avant = length(Sweep[,1])
      apres = apres+1
      k=-1*k
    }else if(apres == length(Sweep[,1])){
      apres = 1
      avant = avant-1
      k=-1*k
    }else{
      avant = (avant-1)
      apres = (apres+1)
    }
  }  

  return(list("Transforme" = Transforme, "param" = rbind(paramG, paramD)))
}


