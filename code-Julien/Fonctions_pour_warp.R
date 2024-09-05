
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

# Fonction linéaire qui réajuste les intensités de gris des rayons pour
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
# Fonction qui inverse l'effet de la fonction qui warp les coordonnées
# (ne fonctionne pas encore pour les mélanges)
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

derive = function(vec){
  dv = 0
  h=1
  #V = splinefun(1:length(vec), vec, method = "natural")
  for(i in 1:(length(vec)-1)){
    dv = c(dv, (vec[i+h]-vec[i])/h)
  }
  return(dv)
}

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
  
  S = m*spline(1:length(r),r,xout=Fonction_Warp(a,b,length(r),w), method = "natural")$y+k
  
  S_temp = (S-mean(c(S,ref)))/sd(c(S,ref))
  ref_temp = (ref-mean(c(S,ref)))/sd(c(S,ref))

  S = S_temp
  ref = ref_temp
  
  # dref = derive(ref)
  # dS = derive(S)
  # 
  # dS_temp = (dS-mean(c(dS,dref)))/sd(c(dS,dref))
  # dref_temp = (dref-mean(c(dS,dref)))/sd(c(dS,dref))
  # 
  # dS = dS_temp
  # dref = dref_temp
  
  # Intensite est une fonction linéaire qui rapporte les intensitees de gris
  # de f vers celles de ref
  
  # for(i in 1:length(S)){
  #   L = L + ((1+exp(-(S[i]-ref[i])^2))*(S[i] - ref[i]))^2#+k^2
  # }
  L= sum((S-ref)^2) #+ sum((dS - dref)^2)# + 0.1*sum((param[1:n]^2-a)^2) + 0.1*sum((param[(n+1):(2*n)]-b)^2)
  
  return(L)
}



# # Fonction à optimiser pour obtenir les paramètres de la fonction de
# # distortion
# Beta = function(x,r,ref,n,param){
#   w=x[paste0("w",1:(n-1))]
#   a=x[paste0("a",1:n)]^2
#   b=x[paste0("b",1:n)]^2
#   m=x["m1"]
#   k=x["k1"]
# 
#   L=0
# 
#   if(is.na(sum(w))){
#     w=1
#   }else{
#     w=exp(w)/(1+sum(exp(w)))
#     w=c(w, 1-sum(w))
#   }
# 
#   S = m*spline(1:length(r),r,xout=Fonction_Warp(a,b,length(r),w), method = "natural")$y+k
#   
#   dref = derive(ref)
#   dS = derive(S)
# 
#   #Intensite est une fonction linéaire qui rapporte les intensitees de gris
#   # de f vers celles de ref
#   L= sum((S-ref)^2) + sum((dS - dref)^2)# + 0.1*sum((param[1:n]^2-a)^2) + 0.1*sum((param[(n+1):(2*n)]-b)^2)
# 
#   return(L)
# }

# Fonction qui fait la transformation comme tel pour aligner les anneaux. Elle retourne l'image
# transformée ainsi que les paramètres
Transformation = function(noyau, Contour, Sweep, Sweep_smooth = Sweep, largeur=length(Sweep[,1]), nb_melange = 1){
  start.time = Sys.time()
  axe_p = Trouver_axe_principal(noyau,Contour)[5] #On trouve l'axe principal

  Ref = Sweep[axe_p,] #Ref contient les valeurs de l'axe principal

  if(nb_melange == 1){
    initG = c(a1=1, b1=1, m1=1, k1=0)
    initD = c(a1=1, b1=1, m1=1, k1=0)

  } else{
    nom = c(paste0("a", 1:nb_melange), paste0("b", 1:nb_melange),
            paste0("w", 1:(nb_melange-1)), "m1", "k1")
    initG = c(rep(1,nb_melange), rep(1,nb_melange),
              rep(1/nb_melange,(nb_melange-1)), 1, 0)
    initD = c(rep(1,nb_melange), rep(1,nb_melange),
              rep(1/nb_melange,(nb_melange-1)), 1, 0)

    initG = setNames(initG, nom)
    initD = setNames(initD, nom)
  }

  paramG = matrix(,length(Sweep[,1])/2,length(initG)+1)

  paramD = matrix(,length(Sweep[,1])/2,length(initD)+1)

  # Transforme = matrix(,length(Sweep[,1]),length(Sweep[1,]))
  # Transforme[axe_p,] = Ref
  Transforme = Sweep

  avant = axe_p - 1
  apres = axe_p + 1
  i = 0
  k=1
  
  while(k*avant <= k*apres && i < largeur){
    # if(mod(i,40) == 0){
    #   plot(as.cimg(Transforme), ylim=c(0,400))
    # }
    i=i+1

    # Ici on calcul les rayons à gauches du rayon de référence en commençant
    # par celui directement à gauche

    param = optim(initG, Beta, r=Sweep_smooth[avant,], ref=Ref, n=nb_melange, param=initG)#, control = list(maxit = 10000))
    if(param$convergence != 0){
      warning(paste0("Pas convergence : erreur de convergence ", param$convergence,
                     " au rayon ", avant, "."))
    }
    paramG[i,] = c(param$par, avant)

    param = optim(initD, Beta, r=Sweep_smooth[apres,], ref=Ref, n=nb_melange, param=initD)#, control = list(maxit = 10000))
    if(param$convergence != 0){
      warning(paste0("Pas convergence : erreur de convergence ", param$convergence,
                     " au rayon ", apres, "."))
    }
    paramD[i,] = c(param$par, apres)


    if(nb_melange == 1){
      nom = c("a1", "b1", "m1", "k1")
      initG = c(a1 = paramG[i,1], b1 = paramG[i,2], m1 = paramG[i,3], k1 = paramG[i,4])
      initD = c(a1 = paramD[i,1], b1 = paramD[i,2], m1 = paramD[i,3], k1 = paramD[i,4])
      wG=wD=1
    } else{
      nom = c(paste0("a", 1:nb_melange), paste0("b", 1:nb_melange),
              paste0("w", 1:(nb_melange-1)), "m1", "k1")
      #if (i==1){
      initG = setNames(paramG[i,-length(paramG[i,])], nom)
      initD = setNames(paramD[i,-length(paramD[i,])], nom)
      #}else if(i<=10){
       # initG = setNames(colMeans(paramG[1:i,-length(paramG[i,])]), nom)
      #  initD = setNames(colMeans(paramD[1:i,-length(paramD[i,])]), nom)
      #}else{
      #  initG = setNames(colMeans(paramG[(i-10):i,-length(paramG[i,])]), nom)
      #  initD = setNames(colMeans(paramD[(i-10):i,-length(paramD[i,])]), nom)
      #}

      wG = (exp(initG[paste0("w", 1:(nb_melange-1))]))/(1+sum(exp(initG[paste0("w", 1:(nb_melange-1))])))
      wG = c(wG, 1-sum(wG))

      wD = (exp(initD[paste0("w", 1:(nb_melange-1))]))/(1+sum(exp(initD[paste0("w", 1:(nb_melange-1))])))
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
  P = P[match(unique(P[,"indice"]),P[,"indice"]), ]
  P[,paste0("a",1:nb_melange)] = P[,paste0("a",1:nb_melange)]^2
  P[,paste0("b",1:nb_melange)] = P[,paste0("b",1:nb_melange)]^2
  if(nb_melange > 1){
    P[,paste0("w", 1:(nb_melange-1))] = (exp(P[,paste0("w", 1:(nb_melange-1))]))/(1+rowSums(exp(P[,paste0("w", 1:(nb_melange-1))])))
  }
  P = P[order(P[,"indice"]),]
  return(list("Transforme" = Transforme, "param" = P[match(unique(P[,"indice"]),P[,"indice"]), ]))
}

# Fonction qui fait la transformation comme tel pour aligner les anneaux. Elle retourne l'image 
# transformée ainsi que les paramètres
# Transformation = function(noyau, Contour, Sweep, largeur=length(Sweep[,1]), nb_melange = 1){
#   start.time = Sys.time()
#   axe_p = Trouver_axe_principal(noyau,Contour)[5] #On trouve l'axe principal
#   
#   Ref = Sweep[axe_p,] #Ref contient les valeurs de l'axe principal
#   
#   if(nb_melange == 1){
#     initG = c(a1=1, b1=1, m1=1, k1=0)
#     initD = c(a1=1, b1=1, m1=1, k1=0)
#     
#   } else{
#     nom = c(paste0("a", 1:nb_melange), paste0("b", 1:nb_melange),
#             paste0("w", 1:(nb_melange-1)), "m1", "k1")
#     initG = c(rep(1,nb_melange), rep(1,nb_melange),
#               rep(1/nb_melange,(nb_melange-1)), 1, 0)
#     initD = c(rep(1,nb_melange), rep(1,nb_melange),
#               rep(1/nb_melange,(nb_melange-1)), 1, 0)
#     
#     initG = setNames(initG, nom)
#     initD = setNames(initD, nom)
#   }
#   
#   paramG = matrix(,length(Sweep[,1])/2,length(initG)+1)
#   
#   paramD = matrix(,length(Sweep[,1])/2,length(initD)+1)
#   
#   Transforme = matrix(,length(Sweep[,1]),length(Sweep[1,]))
#   Transforme[axe_p,] = Ref
#   
#   avant = axe_p - 1
#   apres = axe_p + 1
#   i = 0
#   k=1
#   
#   while(k*avant <= k*apres && i < largeur){
#     i=i+1
#     
#     # Ici on calcul les rayons à gauches du rayon de référence en commençant
#     # par celui directement à gauche
#     valeurG = 10e10
#     valeurD = 10e10
#     Gauss_parG = matrix(, 10, 2*nb_melange)
#     Gauss_parG = set_colnames(Gauss_parG,nom[1:(2*nb_melange)])
#     Gauss_parG[1,] = initG[1:(2*nb_melange)]
#     Gauss_parD = matrix(, 10, 2*nb_melange)
#     Gauss_parD = set_colnames(Gauss_parD,nom[1:(2*nb_melange)])
#     Gauss_parD[1,] = initD[1:(2*nb_melange)]
#     for(j in 1:(2*nb_melange)){
#       Gauss_parG[-1,j] = rnorm(9, mean = initG[j], sd = 1)
#       Gauss_parD[-1,j] = rnorm(9, mean = initD[j], sd = 1)
#     }
#     for(j in 1:10){
#       opG = optim(c(Gauss_parG[j,], initG[(2*nb_melange+1):length(initG)]), Beta, r=Sweep[avant,], ref=Ref, n=nb_melange, param=initG)
#       opD = optim(c(Gauss_parD[j,], initD[(2*nb_melange+1):length(initD)]), Beta, r=Sweep[avant,], ref=Ref, n=nb_melange, param=initG)
#       
#       if(opG$value < valeurG){
#         valeurG = opG$value
#         paramG[i,] = c(opG$par, avant)
#       }
#       if(opD$value < valeurD){
#         valeurD = opD$value
#         paramD[i,] = c(opD$par, apres)
#       }
#     }
#     # param = optim(initG, Beta, r=Sweep[avant,], ref=Ref, n=nb_melange, param=initG)#, control = list(maxit = 10000))
#     # if(param$convergence != 0){
#     #   warning(paste0("Pas convergence : erreur de convergence ", param$convergence,
#     #                  " au rayon ", avant, "."))
#     # }
#     # paramG[i,] = c(param$par, avant)
#     # 
#     # param = optim(initD, Beta, r=Sweep[apres,], ref=Ref, n=nb_melange, param=initD)#, control = list(maxit = 10000))
#     # if(param$convergence != 0){
#     #   warning(paste0("Pas convergence : erreur de convergence ", param$convergence,
#     #                  " au rayon ", apres, "."))
#     # }
#     # paramD[i,] = c(param$par, apres)
#     
#     
#     if(nb_melange == 1){
#       nom = c("a1", "b1", "m1", "k1")
#       initG = c(a1 = paramG[i,1], b1 = paramG[i,2], m1 = paramG[i,3], k1 = paramG[i,4])
#       initD = c(a1 = paramD[i,1], b1 = paramD[i,2], m1 = paramD[i,3], k1 = paramD[i,4])
#       wG=wD=1
#     } else{
#       nom = c(paste0("a", 1:nb_melange), paste0("b", 1:nb_melange),
#               paste0("w", 1:(nb_melange-1)), "m1", "k1")
#       initG = setNames(paramG[i,-length(paramG[i,])], nom)
#       initD = setNames(paramD[i,-length(paramD[i,])], nom)
#       
#       wG = (exp(initG[paste0("w", 1:(nb_melange-1))]))/(1+sum(exp(initG[paste0("w", 1:(nb_melange-1))])))
#       wG = c(wG, 1-sum(wG))
#       
#       wD = (exp(initD[paste0("w", 1:(nb_melange-1))]))/(1+sum(exp(initD[paste0("w", 1:(nb_melange-1))])))
#       wD = c(wD, 1-sum(wD))
#     }
#     
#     Transforme[avant,] = paramG[i,ncol(paramG)-2]*spline(1:length(Sweep[avant,]), Sweep[avant,], xout = Fonction_Warp(paramG[i,1:nb_melange]^2,paramG[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,]), wG), method = "natural")$y+paramG[i,ncol(paramG)-1]
#     Transforme[apres,] = paramD[i,ncol(paramD)-2]*spline(1:length(Sweep[apres,]), Sweep[apres,], xout = Fonction_Warp(paramD[i,1:nb_melange]^2,paramD[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,]), wD), method = "natural")$y+paramD[i,ncol(paramD)-1]
#     
#     if(avant-1 == 0){
#       avant = length(Sweep[,1])
#       apres = apres+1
#       k=-1*k
#     }else if(apres == length(Sweep[,1])){
#       apres = 1
#       avant = avant-1
#       k=-1*k
#     }else{
#       avant = (avant-1)
#       apres = (apres+1)
#     }
#   }
#   end.time = Sys.time()
#   print(end.time-start.time)
#   
#   P = rbind(paramG, c(rep(1,ncol(paramD)-1), axe_p), paramD)
#   colnames(P) = c(nom, "indice")
#   P = P[match(unique(P[,"indice"]),P[,"indice"]), ]
#   P[,paste0("a",1:nb_melange)] = P[,paste0("a",1:nb_melange)]^2
#   P[,paste0("b",1:nb_melange)] = P[,paste0("b",1:nb_melange)]^2
#   if(nb_melange > 1){
#     P[,paste0("w", 1:(nb_melange-1))] = (exp(P[,paste0("w", 1:(nb_melange-1))]))/(1+rowSums(exp(P[,paste0("w", 1:(nb_melange-1))])))
#   }
#   P = P[order(P[,"indice"]),]
#   return(list("Transforme" = Transforme, "param" = P[match(unique(P[,"indice"]),P[,"indice"]), ]))
# }


Warp_spline = function(theta,image,Contour,noyau){
  x = theta[paste0("x", 1:400)]
  y = theta[paste0("y", 1:400)]

  u = splinefun(x, y, method = "natural")
  du = u(x, deriv = 1)
  
  dx = as.matrix(imgradient(as.cimg(image),"x")[,,1,1])
  dy = as.matrix(imgradient(as.cimg(image),"y")[,,1,1])
  d=dx+dy

  #print(cbind(floor(S$x),floor(S$y)))

  return(((d[cbind(seq(noyau[1],Contour[100,1],len=400),seq(noyau[2],Contour[100,2],len=400))])[1:399]%*%du[1:399])^2)
}













# La vielle façon de warper les coordonnées
TransformationVieux = function(noyau, Contour, Sweep, largeur=length(Sweep[,1]), nb_melange = 1){
  start.time = Sys.time()
  axe_p = Trouver_axe_principal(noyau,Contour)[5] #On trouve l'axe principal

  Ref = Sweep[axe_p,] #Ref contient les valeurs de l'axe principal

  if(nb_melange == 1){
    init = c(a1=1, b1=1, m1=1, k1=0)
    
  } else{
    nom = c(paste0("a", 1:nb_melange), paste0("b", 1:nb_melange),
            paste0("w", 1:(nb_melange-1)), "m1", "k1")
    init = c(rep(1,nb_melange), rep(1,nb_melange),
              rep(1/nb_melange,(nb_melange-1)), 1, 0)
    
    init = setNames(init, nom)
  }

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

    if(nb_melange == 1){
      wG=wD=1
    } else{
      wG = (exp(paramG[i,(2*nb_melange+1):(ncol(paramG)-3)]))/(1+sum(exp(paramG[i,(2*nb_melange+1):(ncol(paramG)-3)])))
      wG = c(wG, 1-sum(wG))
      
      wD = (exp(paramD[i,(2*nb_melange+1):(ncol(paramD)-3)]))/(1+sum(exp(paramD[i,(2*nb_melange+1):(ncol(paramD)-3)])))
      wD = c(wD, 1-sum(wD))
    }

    Transforme[avant,] = paramG[i,ncol(paramG)-2]*spline(1:length(Sweep[avant,]), Sweep[avant,], xout = Fonction_Warp(paramG[i,1:nb_melange]^2,paramG[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,]),wG), method = "natural")$y + paramG[i,ncol(paramG)-1]
    Transforme[apres,] = paramD[i,ncol(paramD)-2]*spline(1:length(Sweep[apres,]), Sweep[apres,], xout = Fonction_Warp(paramD[i,1:nb_melange]^2,paramD[i,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,]),wD), method = "natural")$y + paramD[i,ncol(paramD)-1]
    if(i>1){
      for(j in 1:(i-1)){
        if(nb_melange == 1){
          wG=wD=1
        } else{
          wG = (exp(paramG[i-j,(2*nb_melange+1):(ncol(paramG)-3)]))/(1+sum(exp(paramG[i,(2*nb_melange+1):(ncol(paramG)-3)])))
          wG = c(wG, 1-sum(wG))
          
          wD = (exp(paramD[i-j,(2*nb_melange+1):(ncol(paramD)-3)]))/(1+sum(exp(paramD[i,(2*nb_melange+1):(ncol(paramD)-3)])))
          wD = c(wD, 1-sum(wD))
        }
        # f(f(f(f...(transforme)))) = ref
        Transforme[avant,] = paramG[i-j,ncol(paramG)-2]*spline(1:length(Transforme[avant,]), Transforme[avant,], xout = Fonction_Warp(paramG[i-j,1:nb_melange]^2,paramG[i-j,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[avant,]),wG), method = "natural")$y + paramG[i-j,ncol(paramG)-1]
        if(avant != apres)
          Transforme[apres,] = paramD[i-j,ncol(paramD)-2]*spline(1:length(Transforme[apres,]), Transforme[apres,], xout = Fonction_Warp(paramD[i-j,1:nb_melange]^2,paramD[i-j,(nb_melange+1):(2*nb_melange)]^2,length(Sweep[apres,]),wD), method = "natural")$y + paramD[i-j,ncol(paramD)-1]
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
  end.time = Sys.time()
  print(end.time-start.time)

  P = rbind(paramG, c(rep(1,ncol(paramD)-1), axe_p), paramD)
  colnames(P) = c(nom, "indice")
  P = P[match(unique(P[,"indice"]),P[,"indice"]), ]
  P[,paste0("a",1:nb_melange)] = P[,paste0("a",1:nb_melange)]^2
  P[,paste0("b",1:nb_melange)] = P[,paste0("b",1:nb_melange)]^2
  if(nb_melange > 1){
    P[,paste0("w", 1:(nb_melange-1))] = (exp(P[,paste0("w", 1:(nb_melange-1))]))/(1+rowSums(exp(P[,paste0("w", 1:(nb_melange-1))])))
  }
  P = P[order(P[,"indice"]),]
  return(list("Transforme" = Transforme, "param" = P[match(unique(P[,"indice"]),P[,"indice"]), ]))
}


