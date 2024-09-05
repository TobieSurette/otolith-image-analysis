
#########################################################################
# 
# Construction de l'arbre des formes de l'otolithe
#
#########################################################################

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
              Arbre = rbind(Arbre, list(shape = shape, arbre = i))
            }
            shape = c()
          } else{
            num_shape = num_shape+1
            # On ne prends pas les formes trop petites pour être le noyau
            if(area(shape, lvl)$aire > 40){
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


#########################################################################
#
# Trouve le voisinage d'un point
#
#########################################################################

voisinage = function(point, d = -1){
  i=point[1]
  j=point[2]
  
  # Nous donne le voisinage immédiat du point
  vois = cbind(rep(c(i-1, i, i+1),3), c(rep(j-1, 3), rep(j, 3), rep(j+1, 3)))
  vois = vois[-5,] # On enlève le point comme tel
  
  vois = set_rownames(vois, c("SO", "S", "SE", "O", "E", "NO", "N", "NE"))
  
  # On classe les points du voisinnage par ordre d'importance selon la direction
  # donnée. On ne fait rien si aucune direction donnée.
  switch (d,
          NE = {vois = vois[c("E", "N", "S", "O", "NE","SE", "NO", "SO"), ]},
          E = {vois = vois[c("E", "N", "S", "O", "NE","SE", "NO", "SO"), ]},
          SE = {vois = vois[c("E", "S", "N", "O", "SE","NE", "SO", "NO"), ]},
          N = {vois = vois[c("N", "E", "O", "S", "NE", "NO", "SO", "SE"), ]},
          S = {vois = vois[c("S", "E", "O", "N", "SO", "SE", "NO", "NE"), ]},
          NO = {vois = vois[c("O", "N", "S", "E", "NO", "NE", "SO", "SE"), ]},
          O = {vois = vois[c("O", "N", "S", "E", "NO", "NE", "SE", "SO"), ]},
          SO = {vois = vois[c("O", "S", "E", "N", "SO", "SE", "NO", "NE"), ]},
  )
  return(vois)
}


#########################################################################
#
# Donne la direction où l'on se dirige dans la définition du contour d'une forme
#
#########################################################################

direction = function(shape){
  # On va chercher les deux derniers points dans shape
  if(!is.null(nrow(shape)) && nrow(shape) > 1){
    dernier = shape[nrow(shape), ]
    avant_dernier = shape[nrow(shape)-1, ]
    
    # On fait leur soustraction et selon si x et/ou y sont négatifs ou positif,
    # on obtient la direction
    d = dernier - avant_dernier
    
    if(d[1] > 0){
      if(d[2] > 0) direc = "NE"
      else if(d[2] == 0) direc = "E"
      else if(d[2] < 0) direc = "SE"
      
    } else if(d[1] == 0){
      if(d[2] > 0) direc = "N"
      else if(d[2] < 0) direc = "S"
      
    } else if(d[1] < 0){
      if(d[2] > 0) direc = "NO"
      else if(d[2] == 0) direc = "O"
      else if(d[2] < 0) direc = "SO"
    }
  } else direc = "E"
  return(direc)
}


#########################################################################
#
# Nous dis si un point fait partie d'un ensemble de points
#
#########################################################################

coordinate_in_matrix = function(mat, coord) {
  # Check si coordonnées sont dans la matrice
  
  logique = c()
  rangee = c()
  
  if(is.null(nrow(coord))){
    coord = matrix(coord,1,2)
  }
  
  if(is.null(nrow(mat))){
    for(i in 1:nrow(coord)){
      logique = identical(mat, coord[i,])
    }
    return(list(logic = logique, row = NULL))
  } else{
    for(i in 1:nrow(coord)){
      logique = c(logique, any(apply(mat, 1, function(row) all(row == coord[i,]))))
      if(logique[length(logique)])
        rangee = c(rangee, which(apply(mat, 1, function(row) all(row == coord[i,])) == TRUE))
      else
        rangee = c(rangee, NULL)
    }
    
    return(list(logic = logique, row = rangee))
  }
}


##########################################################################
#
# Nous donne les enfants direct d'un noeud
#
##########################################################################

enfant = function(Arbre, shape, num_arbre){
  e = c()
  # n contient les indices des Arbres au niveau de threshold supérieur
  n = which(Arbre[,2] == num_arbre+1)
  
  # Si n est vide, ça veut dire qu'il n'y a pas d'enfants.
  if(length(n) == 0){
    return(list(enfant = NULL, arbre = Arbre))
  }
  ArbreN = Arbre[n,]
  
  # On place dans e les n qui sont les enfants de la forme
  if(length(ArbreN) == 0){
    return(NULL)
    
  } else if(is.null(nrow(ArbreN))){
    if(polygon_in_polygon(ArbreN$shape, shape)){
      e = n
    }
  } else{
    for(i in 1:nrow(ArbreN)){
      if(polygon_in_polygon(ArbreN[[i,1]], shape)){
        e = c(e, n[i])
      }
    }
  }
  
  # Parmi les enfants trouvés on enlève ceux qui sont enfants des enfants trouvés
  # on augmente aussi leur niveau d'arbre de 1 pour être traité plus tard
  for(i in e){
    for(j in e){
      if(i != j){
        if(polygon_in_polygon(Arbre[i,1]$shape, Arbre[j,1]$shape)){
          e = e[-which(e == i)]
          Arbre[i,2]$arbre = Arbre[i,2]$arbre + 1
        }
      }
    }
  }
  # On retourne ensuite le vecteur d'indices des enfants ainsi que la nouvelle matrice Arbre
  # qui a peut-être été changé
  return(list(enfant = e, arbre = Arbre))
}


##########################################################################
#
# Nous donne tous les enfants d'un noeud
#
##########################################################################

Trouver_enfants <- function(Arbre, num_node) {
  # Initialize an empty vector to hold all descendants
  descendants = c()
  
  # Helper function to recursively find descendants
  find_recursive <- function(current_node) {
    children <- Arbre[current_node,]$enfant
    for (child in children) {
      descendants <<- c(descendants, child) 
      find_recursive(child)  
    }
  }
  
  # Start the recursive search from the given node
  find_recursive(num_node)
  
  return(descendants[order(unlist(Arbre[descendants,2]))])
}


##########################################################################
#
# Nous donne si un point ou une forme fait partie d'une autre forme
#
##########################################################################

point_in_polygon <- function(point, poly) {
  
  if(coordinate_in_matrix(poly,point)$logic){
    return(TRUE)
  }
  x <- point[1]
  y <- point[2]
  n <- nrow(poly)
  
  inside <- FALSE
  j <- n
  for (i in 1:n) {
    xi <- poly[i, 1]
    yi <- poly[i, 2]
    xj <- poly[j, 1]
    yj <- poly[j, 2]
    
    intersect <- ((yi > y) != (yj > y)) &&
      (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
    if (intersect) inside <- !inside
    j <- i
  }
  
  return(inside)
}

polygon_in_polygon <- function(polygon1, polygon2, tol = 0.9) {
  polygon1 = Ordre_cercle(polygon1)
  polygon2 = Ordre_cercle(polygon2)
  A = sapply(1:nrow(polygon1), function(i) point_in_polygon(polygon1[i, ], polygon2))
  A = length(which(A == TRUE))/length(A) >= tol
  return(A)
}


##########################################################################
#
# Réarrange les points en ordre circulaire
#
##########################################################################

Ordre_cercle = function(coords){
  centroid = colMeans(coords)
  
  # Calcul les angles
  angles = atan2(coords[,2] - centroid[2], coords[,1] - centroid[1])
  
  # Ordonne selonles angles
  ordered_indices = order(angles)
  return(coords[ordered_indices, ])
}


##########################################################################
#
# Donne l'axe principal (plus long) de tout l'otolithe
#
##########################################################################

Principal_Axis = function(Contour){
  dmax = 0
  # On passe tous les points du contour et garde les deux qui sont les plus loin
  for(i in 1:(nrow(Contour)-1)){
    for(j in (i+1):nrow(Contour)){
      distance = euclidien(Contour[i,],Contour[j,])
      if(distance > dmax){
        dmax = distance
        axe = rbind(Contour[i,], Contour[j,])
      }
    }
  }
  return(axe)
}


##########################################################################
#
# Étape 2.2.2 de l'article (Principal axis proximity) ne fonctionne pas encore
#
##########################################################################

# Principal_Axis_Proximity = function(Arbre, Contour, decoupe){
#   axe_p = Principal_Axis(Contour)
#   axe_p = cbind(seq(axe_p[1,1], axe_p[2,1], len = 100), seq(axe_p[1,2], axe_p[2,2], len = 100))
#   reg = lm(y~x, data.frame(x = axe_p[,1], y = axe_p[,2]))
#   mp = -1/reg$coefficients[2]
#   m = -1/reg$coefficients[2]
#   #b = reg$coefficients[1]
#   
#   for(i in 1:100){
#     x = axe_p[i,1]
#     b = axe_p[i,2] - m*x
#     perp = c()
#     print(i)
#     while(decoupe[x, round(m*x+b)] != 0){
#       perp = rbind(perp, cbind(x, m*x+b))
# 
#       for(arbre in 1:nrow(Arbre)){
# 
#         if(point_in_polygon(c(x, m*x+b), Arbre[[arbre, 1]])){
# 
#           distance = euclidien(axe_p[i,], c(x, m*x+b))
# 
#           D = rbind(D, c(distance, i, arbre))
#         }
#       }
#       x=x+1
#     }
#     x = axe_p[i,1]
#     while(decoupe[x, round(m*x+b)] != 0){
#       perp = rbind(perp, cbind(x, m*x+b))
# 
#       for(arbre in 1:nrow(Arbre)){
# 
#         if(point_in_polygon(c(x, m*x+b),Arbre[[arbre, 1]])){
# 
#           distance = euclidien(axe_p[i,], c(x, m*x+b))
# 
#           D = rbind(D, c(distance, i, arbre))
#         }
#       }
#       x=x-1
#     }
#   }
# }


Principal_Axis_Proximity = function(Arbre, Contour, decoupe){
  D = c()
  B = c()
  
  axe_p = Principal_Axis(Contour)
  axe_p = cbind(seq(axe_p[1,1], axe_p[2,1], len = 100), seq(axe_p[1,2], axe_p[2,2], len = 100))
  reg = lm(y~x, data.frame(x = axe_p[,1], y = axe_p[,2]))
  mp = -1/reg$coefficients[2]
  m = reg$coefficients[2]
  b = reg$coefficients[1]
  
  #axe_p = splinefun(axe_p, method = "natural")
  
  
  for(i in 1:nrow(Arbre)){
    dmin = 100000000
    #b = c()
    proj = c()
    
    proj = Projection_axe_principal(Arbre[[i,1]], m, b)
    
    for(j in 1:nrow(Arbre[[i,1]])){
      
      distance = distance_point_to_line(m,b,Arbre[[i,1]][j,1], Arbre[[i,1]][j,2])
      
      if(distance < dmin){
        dmin = distance
      }
      
      D = rbind(D, c(proj[j,], distance, i))
      
    }
  }
  
  enleve = c()
  
  for(i in nrow(Arbre) : 1){
    children = Arbre[[i,3]]
    
    if(!is.null(children)){
      D_child = c()
      for(child in children){
        D_child = rbind(D_child, D[which(D[,4] == child), ])
      }
      
      #D_child = D_child[order(D_child[,1]), ]
      
      # minmax = c()
      # for(child in children){
      #   minmax = rbind(minmax, c(child, min(which(D_child[,4] == child)), max(which(D_child[,4] == child))))
      # }
      
      proche = c()
      if(length(children) > 1){
        for(j in min(D_child[,1]):max(D_child[,1])){
          di = 0
          reste = TRUE
          while(reste){
            for(child in children){
              
              temp = D_child[order(D_child[which(D_child[,4] == child),1] - j),]
              temp = temp[1:2,]
              
              if(!between(j, temp[1,1], temp[2,1])){
                break
              }
              
              if(any(di >= temp[,3])){
                proche = c(proche, D_child[which(D_child[,4] == child), 4])
                reste = FALSE
                print(c(j,di))
              }
            }
            di = di+1
            
          }
          
        }
      }
    } else{
      proche = D_child[,4]
    }
    
    browser()
  }
  return(proche)
  
  # for(i in 1:(nrow(Arbre)-1)){
  #   Di = D[which(D[,4] == i), ]
  #   for(j in (i+1):nrow(Arbre)){
  #     Dj = D[which(D[,4] == j), ]
  # 
  #     skip = FALSE
  #     if(any(j == Trouver_enfants(Arbre, i)) ){
  #       skip = TRUE
  #       enleve = c(enleve, -i)
  #       print(c("a",i))
  #     }
  #     # if(any(i == Trouver_enfants(Arbre, j))){
  #     #   skip = TRUE
  #     #   enleve = c(enleve, -j)
  #     # }
  # 
  # 
  #     if(!skip && i!=j && ((between(min(Di[,1]), min(Dj[,1]), max(Dj[,1]))) || (between(max(Di[,1]), min(Dj[,1]), max(Dj[,1])))) ){
  # 
  # 
  #       if(min(Dj[,1]) <= min(Di[,1])){
  #         borne_inf = min(Di[,1])
  #       } else{
  #         borne_inf = min(Dj[,1])
  #       }
  # 
  #       if(max(Dj[,1]) <= max(Di[,1])){
  #         borne_sup = max(Dj[,1])
  #       } else{
  #         borne_sup = max(Di[,1])
  #       }
  # 
  #       dmin_i = min(Di[which(between(Di[,1], borne_inf, borne_sup) == TRUE), 3])
  #       dmin_j = min(Dj[which(between(Dj[,1], borne_inf, borne_sup) == TRUE), 3])
  # 
  #       if(dmin_i < dmin_j){
  #         enleve = c(enleve, -j)
  #         print(c("b",j))
  #       } else{
  #         enleve = c(enleve, -i)
  #         print(c("c",i))
  #       }
  # 
  #     }
  #   }
  # }
  # 
  # # for(i in 1:nrow(Arbre)){
  # #   for(j in 1:nrow(Arbre)){
  # #
  # #     if( (i != j) && ((B[i,1] <= B[j,1] && B[i,2] >= B[j,1]) || (B[i,1] <= B[j,2] && B[i,2] >= B[j,2])) ){
  # #       if(D[j] >= D[i]){
  # #         enleve = c(enleve, -j)
  # #       }
  # #     }
  # #   }
  # # }
  # 
  # # D = set_colnames(D,c("distance", "ligne", "arbre"))
  # #
  # # enleve = c()
  # # for(i in unique(D[,"ligne"])){
  # #   J = D[which(D[,"ligne"] == i),]
  # #   U = c()
  # #   for(j in unique(J[,"arbre"])){
  # #     U = c(U,min(J[which(J[,"arbre"] == j), "distance"]))
  # #     if(min(U) != U[length(U)]){
  # #       enleve = c(enleve, -j)
  # #     }
  # #   }
  # # }
  # 
  # 
  # return(list(a = Arbre[unique(enleve),], e = unique(enleve)))
}





##########################################################################
#
# Donne l'aire d'une forme
#
##########################################################################

area = function(shape, lvl_set){
  aire = 0
  point = c()
  
  centre = colMeans(shape)
  
  if(lvl_set[centre[1],centre[2]] == 1){
    couleur = 1
  } else{
    couleur = 0
  }
  # L'aire est le nombre de points de la même que le centre de la forme
  for(i in min(shape[,1]):max(shape[,1])){
    for(j in min(shape[,2]):max(shape[,2])){
      if(lvl_set[i,j] == couleur){
        aire = aire + 1
        point = rbind(point, c(i,j))
      }
    }
  }
  return(list(aire = aire, points = point))
}


##########################################################################
#
# Dis si un nombre se trouve entre deux autres
#
##########################################################################

between = function(point, inf, sup){
  entre = c()
  for(i in 1:length(point)){
    entre = c(entre, point[i] >= inf & point[i] <= sup)
  }
  return(entre)
}


##########################################################################
#
# Fait la projection d'une forme sur une droite
#
##########################################################################

Projection_axe_principal = function(shape, m, b){
  shape_proj = c()
  if(length(shape) == 2){
    x_proj = (shape[1] + m * (shape[2] - b)) / (1 + m^2)
    y_proj = m * x_proj + b
    
    shape_proj = rbind(shape_proj, c(x_proj, y_proj))
    
  } else{
    
    for(i in 1:nrow(shape)){
      x_proj = (shape[i,1] + m * (shape[i,2] - b)) / (1 + m^2)
      y_proj = m * x_proj + b
      
      shape_proj = rbind(shape_proj, c(x_proj, y_proj))
    }
  }
  
  return(shape_proj)
}


##########################################################################
#
# La longueur de l'étendu des points de la forme lorsque projetés
# sur l'axe principal 
#
##########################################################################

Diametre_Projection_axe_principal = function(shape, Contour){
  axe_p = Principal_Axis(Contour)
  reg = lm(y~x, data.frame(x = axe_p[,1], y = axe_p[,2]))
  m = -1 / reg$coefficients[2]
  b = reg$coefficients[1]
  
  P = Projection_axe_principal(shape, m, b)
  
  D = sqrt((max(P[,1]) - min(P[,1]))^2 + (max(P[,2]) - min(P[,2]))^2)
  
  return(D)
}


##########################################################################
#
# Fi de la section 2.2.3. Représentent 3 qualités géométriques d'un bon noyau
#
##########################################################################

Fi = function(S, lvl_set, Contour){
  
  # Calcul l'aire de la forme
  aire = area(S, lvl_set)$aire
  
  # Calcul la distance entre l'axe principal et la forme
  d = dist_axe_principal(S, Contour, minimum = FALSE)
  
  # La longueur de l'étendu des points de la forme lorsque projetés
  # sur l'axe principal
  D = Diametre_Projection_axe_principal(S, Contour)
  
  f1 = aire
  f2 = aire / (d^2)
  f3 = aire / (D^2)
  #f4 = Similitude_des_rayons(colMeans(S), image, Contour)
  
  return(c(f1 = f1, f2 = f2, f3 = f3))
}


##########################################################################
#
# Hi et E1 de la section 2.2.3
#
##########################################################################

Hi = function(mu, fi, i){
  nb = 0
  for(j in 1:nrow(fi)){
    if(fi[j, i] > mu){
      nb = nb+1
    }
  }
  
  Hi = nb / nrow(fi) # nrow(fi) = # de shapes
  
  return(Hi)
}

E1 = function(fi_S, fi_total){
  H = c()
  for(i in 1:3){
    H = c(H, Hi(mu = fi_S[i], fi_total, i))
  }
  
  E1 = max(H)^3
}


##########################################################################
#
# Distance entre l'axe principal et la forme
#
##########################################################################

dist_axe_principal = function(shape, Contour, minimum = TRUE){
  D = c()
  
  axe_p = Principal_Axis(Contour)
  
  reg = lm(y~x, data.frame(x = axe_p[,1], y = axe_p[,2]))
  m = reg$coefficients[2]
  b = reg$coefficients[1]
  
  if(minimum){
    dmin = 100000000
    y = c()
    x = c()
    for(j in 1:nrow(shape)){
      distance = distance_point_to_line(m,b,shape[j,1], shape[j,2])
      
      if(distance < dmin){
        dmin = distance
      }
    }
    return(dmin)
  } else{
    dmax = 0
    y = c()
    x = c()
    for(j in 1:nrow(shape)){
      distance = distance_point_to_line(m,b,shape[j,1], shape[j,2])
      
      if(distance > dmax){
        dmax = distance
      }
    }
    return(dmax)
  }
  
}


cercle = list()
centre = c()
for(i in 1:nrow(Arbre)){
  point_pattern <- ppp(x = Arbre[[i,1]][, 1], y = Arbre[[i,1]][, 2], window = owin(range(Arbre[[i,1]][, 1]), range(Arbre[[i,1]][, 2])))
  ch <- convexhull(point_pattern)
  enclosing_circle <- boundingcircle(ch)
  
  cercle = c(cercle, list(enclosing_circle$bdry[[1]]))
  centre = rbind(centre, c(mean(cercle[[i]]$x), mean(cercle[[i]]$y)))
}

for(i in 1:nrow(centre)){
  lines(cercle[[i]]$x, cercle[[i]]$y, col = "blue")
  Sys.sleep(1)
}

