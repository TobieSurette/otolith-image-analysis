
library(zoo)

euclidien = function(a, b){
  if(length(a) == 2 && length(b) == 2){
    return(sqrt(sum(a - b)^2))
  } else if(length(a) > 2 && length(b) == 2){
    return(sqrt((a[,1] - b[1])^2 + (a[,2] - b[2])^2))
  } else if(length(a) == 2 && length(b) > 2){
    return(sqrt((a[1] - b[,1])^2 + (a[2] - b[,2])^2))
  } else
    return("Seulement une matrice !")
} 

# Fait un sweep de tout le tour de l'otolithe à partir du noyau et 
# l'étend 
Sweep_anneaux = function(point, image, Contour, n_points){
  L = matrix(,nrow = n_points, ncol = length(Contour[,1]))
  
  for(j in 1:length(Contour[,1])){
    indiceX = floor(seq(point[1],Contour[j,1],length.out=n_points))
    indiceY = floor(seq(point[2],Contour[j,2],length.out=n_points))
    L[,j] = image[cbind(indiceX, indiceY)]
  }
  
  return(L)
}

Tranche_Contour = function(image,noyau,Contour,angle = pi/4, n_points = 100){
  axe_p = Trouver_axe_principal(noyau, Contour)[5]
  angle_p = seq(0,2*pi,len=length(Contour[,1]))[axe_p]
  
  angles = seq(0,2*pi, length.out = length(Contour[,1]))
  
  angles_tranche = seq(angle_p - angle/2, angle_p + angle/2, length.out = n_points)
  for(i in 1:n_points){
    if(angles_tranche[i] < 0){
      angles_tranche[i] = angles_tranche[i]+2*pi
    }
  }
  
  TrancheX = floor(spline(angles, Contour[,1], xout = angles_tranche)$y)
  TrancheY = floor(spline(angles, Contour[,2], xout = angles_tranche)$y)
  
  return(matrix(c(TrancheX,TrancheY), length(TrancheX), 2))
}

Enrouler = function(Sweep, Contour, point){
  n_points=length(Sweep[1,])
  
  L = matrix(,nrow = max(Contour[,1]), ncol = max(Contour[,2]))
  
  for(j in 1:length(Contour[,1])){
    indiceX = floor(seq(point[1],Contour[j,1],length.out=n_points))
    indiceY = floor(seq(point[2],Contour[j,2],length.out=n_points))
    L[cbind(indiceX,indiceY)] = Sweep[j,]
  }
  return(L)
}

Similitude_des_rayons = function(point, image, Contour, n_points = 400){
  L = t(Sweep_anneaux(point,image,Contour,n_points))
  
  # Prend la moyenne de gris pour chaques anneaux.
  moy = colMeans(L)
  
  d = derive(moy)
  
  S=0
  
  for(i in 1:length(L[1,])){
    S=(S+sum((L[,i] - moy)^2 + (derive(L[,i]) - d)^2))
  }
  return(S)
}

Nombre_Formes = function(image){
  a = as.matrix(label(as.cimg(image)))
  a = unique(as.vector(a))
  return(length(a)-1)
}

Trouver_noyau = function(image, i=0, g=30){
  k = 0.01
  while (TRUE) {
    while(sum(clean(as.cimg(image)>i+k, g))>1){
      i=i+k
    }
    if (Nombre_Formes(matrix(clean(as.cimg(image)>i, g),nrow(image),ncol(image))) > 1){
      k = k*0.5
    } else{
      break
    }
  }
  
  noyau = which(matrix(clean(as.cimg(image)>i, g),nrow(image),ncol(image)) == 1, arr.ind = TRUE)
  
  centre = c(mean(noyau[,1]), mean(noyau[,2]))
  
  
  return(round(centre))
}

Distance_Max = function(Coord){
  distmax = 0
  for(i in 1:(nrow(Coord)-1)){
    for(j in (i+1):nrow(Coord)){
      d = euclidien(Coord[i,], Coord[j,])
      if(d > distmax){
        distmax = d
        C1 = Coord[i,]
        C2 = Coord[j,]
      }
    }
  }
  return(matrix(c(C1[1], C2[1], C1[2], C2[2]),2,2))
}


TrouverFormes = function(x,image){
  IM = as.cimg(clean(image>x[1],x[2]))
  return((mean(image[IM==1])-1)^2+100*(mean(image[IM==0])-0)^2)
}


TrouverContour = function(image_grise,n_points, nb_centre = 3){
  #tc = optim(c(0.2,30),TrouverFormes,image = as.cimg(image_grise))$par
  image = as.vector(clean(as.cimg(image_grise)>0.2,30))
  image = matrix(image, nrow(image_grise), ncol(image_grise))
  
  NP = which(image == 1, arr.ind = TRUE)
  rand_i = runif(nb_centre,1,nrow(NP))
  
  C = c()
  x=c()
  y=c()
  angles = seq(0,2*pi, length.out = n_points)
  
  for(j in 1:nb_centre){
    if(nb_centre == 1){
      centre = c(nrow(image_grise)/2, ncol(image_grise)/2)
    } else{
      centre = c(NP[rand_i[j], 1], NP[rand_i[j], 2])
    }
    for(theta in angles){
      r=0
      k=100
      while(TRUE){
        indice = c(floor(centre[1]+(r+k)*cos(theta)),floor(centre[2]+(r+k)*sin(theta)))
        if(indice[1] > nrow(image) || indice[1] < 1 || indice[2] > ncol(image) || indice[2] < 1){
          k=k*0.1
        }
        else {
          if(image[indice[1], indice[2]] == 1){
            r=r+k
          }else{
            k=k*0.1
          }
        }
        if(k<1e-10) break
      }
      x=c(x,floor(centre[1]+r*cos(theta)))
      y=c(y,floor(centre[2]+r*sin(theta)))
    }
    C = rbind(C, centre)
  }
  
  Contour=matrix(data=c(x,y), nrow = length(x), ncol = 2)
  
  decoupe = image
  decoupe[which(image == 1, arr.ind = TRUE)] = image_grise[which(image == 1, arr.ind = TRUE)]
  decoupe[which(decoupe == 0, arr.ind = TRUE)] = 0
  
  return(list(coord = Contour, param = c(0.2,30), decoupe = decoupe, centres = C))
}

# Step est un pourcentage du montant de scaling de 0 à 1, 0 étant le contour ne 
# change pas et 1 est le noyau
Contour_scaling = function(Contour, noyau, step = 0.8){
  new_Contour = Contour
  for(i in 1:nrow(Contour)){
    new_Contour[i,1] = seq(Contour[i,1], noyau[1], len=1000)[round(1000*step)]
    new_Contour[i,2] = seq(Contour[i,2], noyau[2], len=1000)[round(1000*step)]
  }
  return(new_Contour)
}

Level_set = function(img, lambda){
  return(list(lower = matrix(fill(clean(img <= lambda,5),5),nrow(img),ncol(img)), upper = matrix(fill(clean(img >= lambda,5),5),nrow(img),ncol(img))))
  #return(list(lower = medianblur(img,10) <= lambda, upper = medianblur(img,10) >= lambda))
}

Level_lines = function(lvl_set, alpha = 5){
  return(which(as.matrix(dilate_square(as.cimg(lvl_set),alpha) - erode_square(as.cimg(lvl_set),alpha)) == TRUE, arr.ind=TRUE))
}

# voisinage = function(i,j,img,taille = 1){
#   haut = i + taille
#   bas = i - taille
#   gauche = j - taille
#   droite = j + taille
#   if(haut > nrow(img)) haut = nrow(img)
#   if(bas < 1) bas = 1
#   if(gauche < 1) gauche = 1
#   if(droite > ncol(img)) droite = ncol(img)
#   
#   return(list(gauche = gauche, droite = droite, bas = bas, haut = haut))
# }





shape_Contour = function(lvl_set){
  shape = c()
  longueur = c()
  for(i in 1:nrow(lvl_set)){
    for(j in 1:ncol(lvl_set)){
      ligne = c()
      vois = voisinage(i,j,lvl_set)
      gauche = vois$gauche
      droite = vois$droite
      bas = vois$bas
      haut = vois$haut
      vois = lvl_set[bas:haut,gauche:droite]
      #vois[1,1] = vois[1,length(gauche:droite)] = vois[length(bas:haut),1] = vois[length(bas:haut),length(gauche:droite)] = 0

      if(0 < mean(lvl_set[bas:haut,gauche:droite]) && mean(lvl_set[bas:haut,gauche:droite]) < 1 && lvl_set[i,j] == 1 && !coordinate_in_matrix(shape,c(i,j))){
        ligne = rbind(ligne, c(i,j))
        i_avant = i
        j_avant = j
        while(TRUE){
          new_shape=TRUE
          duplicate = FALSE
          ind = which(vois == 1, arr.ind = TRUE)
          
          # On met les coins en dernier, car on veut prioriser ne pas bouger en diagonal sauf nécessaire
          ind = ind[unique(c(which(ind[,1]==2), which(ind[,2]==2),which(ind[,1]==1), which(ind[,2]==1),which(ind[,1]==3), which(ind[,2]==3))), ]
          ind = ind-2
          ind[,1] = i_avant+ind[,1] 
          ind[,2] = j_avant+ind[,2]
          for(v in 1:nrow(ind)){
            if(coordinate_in_matrix(shape,ind[v,])){
              new_shape = FALSE
              duplicate = TRUE
              break
            }
            if(!coordinate_in_matrix(ligne,ind[v,])){
              vois = voisinage(ind[v,1],ind[v,2],lvl_set)
              gauche = vois$gauche
              droite = vois$droite
              bas = vois$bas
              haut = vois$haut
              vois = lvl_set[bas:haut,gauche:droite]
              #vois[1,1] = vois[1,length(gauche:droite)] = vois[length(bas:haut),1] = vois[length(bas:haut),length(gauche:droite)] = 0
              if(0 < mean(lvl_set[bas:haut,gauche:droite]) && mean(lvl_set[bas:haut,gauche:droite]) < 1){
                #print(c(ind[v,1],ind[v,2]))
                ligne = rbind(ligne, c(ind[v,1],ind[v,2]))
                i_avant = ind[v,1]
                j_avant = ind[v,2]
                new_shape = FALSE
                break
              }
            }
          }
          if(new_shape){
            longueur = c(longueur, nrow(ligne))
            shape = rbind(shape, ligne)
            break
          }
          if(duplicate){
            break
          }
        }
      }
    }
  }
  return(list(shape = shape, longueur = longueur))
}



build_shape_tree = function(lvl_set){
  
  Arbre = c()
  for(i in length(lvl_set[1,1,]):1){
    # On diminue la résolution de l'image pour accélérer le calcul
    resX = round(nrow(lvl_set[,,i])/5)
    resY = round(ncol(lvl_set[,,i])/5)
    lvl = as.matrix(resize(as.cimg(lvl_set[,,i]), resX, resY))
    num_shape = 0
    
    # if qui s'assure que l'image n'est pas tout blanc ou tout noir
    if(mean(lvl) > 0 && mean(lvl) < 1){
      start.time = Sys.time()
      
      # blanc contient tous les pixels blancs de l'image
      blanc = which(lvl == 1, arr.ind = TRUE)
      
      # Ici on trouve tous les pixels qui compose la bordure des formes
      bordure = blanc[which(apply(blanc, 1, function(point) mean(lvl[voisinage(point)]) != 1) == TRUE), ]

      # On initialise la pile dans laquelle on viendra ajouter les pixels de bordure de chaque formes
      pile = matrix(bordure[1,],1,2)
      bordure = bordure[-1, ]
      shape=c()
      
      while(TRUE){
        
        # On ajoute à la forme courante le pixel au top de la pile et on le supprime de cette pile
        if(is.null(nrow(pile))){
          shape = rbind(shape, pile)
          direct = direction(shape)
          vois = voisinage(pile,direct)
          pile=c()
        }else{
          shape = rbind(shape, pile[1,])
          direct = direction(shape)
          vois = voisinage(pile[1,],direct)
          pile = pile[-1,]
        }
        
        # On passe au travers de tous les pixels du voisinage du nouveau pixel qui était au top de la pile
        vois_bordure = coordinate_in_matrix(bordure,vois)
        vois = vois[which(vois_bordure$logic == TRUE), ]
        pile = rbind(pile, vois)
        # if(length(vois) > 2)
        #   pile = rbind(pile, vois[1,])
        # else
        #   pile = rbind(pile, vois)
        if(!is.null(vois_bordure$row))
          bordure = bordure[-vois_bordure$row, ]
        
        # Une fois la pile vidée, on regarde s'il reste plus d'une coordonnée dans bordure,
        # si oui, cela veut dire qu'il reste encore des formes dans l'image. Dans ce cas,
        # on ajoute dans la pile le premier élément de bordure, puisqu'il doit sûrement se 
        # trouver sur une autre forme et on l'enlève de bordure. Sinon, ça veut dire qu'on a
        # trouver toutes les formes de l'image et on arrête
        if(length(pile) == 0){
          if(length(bordure)>2){
            pile = bordure[1,]
            bordure = bordure[-1,]
            num_shape = num_shape+1
            # On ne prends pas les formes trop petites pour être le noyau
            if(area(shape, lvl)$aire > 40){
              # e = enfant(Arbre, shape, i)
              # Arbre = e$arbre
              # e = e$enfant
              Arbre = rbind(Arbre, list(shape = shape, arbre = i))
            }
            shape = c()
          } else{
            num_shape = num_shape+1
            # On ne prends pas les formes trop petites pour être le noyau
            if(area(shape, lvl)$aire > 40){
              # e = enfant(Arbre, shape, i)
              # Arbre = e$arbre
              # e = e$enfant
              Arbre = rbind(Arbre, list(shape = shape, arbre = i))
            }
            shape = c()
            break
          }
        }
        
      }
      
      end.time = Sys.time()
      print(end.time - start.time)
    }
  }
  # On essaye de placer les points du contour en ordre, en assumant que c'est un cercle (ce qui n'est pas souvent le cas)
  for(i in 1:nrow(Arbre)){
    Arbre[i,]$shape = Ordre_cercle(cbind(round((nrow(lvl_set[,,1])/resX)*Arbre[i,]$shape[,1]), round((ncol(lvl_set[,,1])/resY)*Arbre[i, ]$shape[,2])))
  }
  
  e=c()
  for(i in nrow(Arbre):1){
    child = enfant(Arbre, Arbre[i,1]$shape, Arbre[i,2]$arbre)
    Arbre = child$arbre
    
    e = rbind(e, list(enfant = child$enfant))
    
  }
  Arbre = cbind(Arbre, e[length(e):1])
  return(Arbre)
}


Maxima_elimination = function(Arbre){
  E = c()
  for(i in 1:length(Arbre)){
    enfant = FALSE
    print(i)
    for(j in 1:length(Arbre)){
      if(Arbre[[j,2]] > Arbre[[i,2]]){
        enfant = enfant | polygon_in_polygon(Arbre[[j,1]], Arbre[[i,1]])
        if(enfant){
          break
        }
      }
    }
    E = c(E, enfant)
  }
  Arbre = Arbre[-which(E == TRUE)]
  return(Arbre)
}

distance_point_to_line <- function(m, b, x0, y0) {
  numerator <- abs(m * x0 - y0 + b)
  denominator <- sqrt(m^2 + 1)
  distance <- numerator / denominator
  return(distance)
}


spline_Contour = function(Contour, point){
  angles = seq(0,2*pi, length.out = length(Contour[,1])+1)
  SX = splinefun(angles, c(Contour[,1], Contour[1,1]), method = "periodic")
  SY = splinefun(angles, c(Contour[,2], Contour[1,2]), method = "periodic")
  
  return(data.frame(x = SX$y, y = SY$y))
}


Trouver_niveau = function(Arbre, num_node){
  compteur = 0
  niveau = c()
  
  children = Arbre[num_node, ]$enfant
  
  if(is.null(children)){
    return(0)
  }
  
  depths = c()
  for(child in children){
    depths = c(depths, 1 + Trouver_niveau(Arbre, child))
  }
  
  return(max(depths))
}







minima_elimination = function(Arbre){
  for(i in 1:nrow(Arbre)){
    e = Trouver_enfants(Arbre, i)
    if(!is.null(e)){
      
    }
  }
}



