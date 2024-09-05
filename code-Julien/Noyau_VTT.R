# VTT = which(lvl_set_up[,,2] == TRUE, arr.ind = TRUE)
# VTT = VTT[runif(1000, 1, nrow(VTT)), ]
# # 
# # plot(as.cimg(decoupe))
# # points(VTT[1], VTT[2])
# # 
# # dxy = as.matrix(imgradient(as.cimg(dx), "y"))
# # dxy = sqrt(dx^2+dy^2)
# # 
# # dxy = as.matrix(medianblur(as.cimg(Contraste2(dx+dy)), 10))
# # 
# # plot(medianblur(as.cimg(Contraste2(dx+dy)), 10))
# 
# 
# 
# Trouver_vallee = function(VTT, img, step = 1, tol = 1e-6, max_iter = 1e4){
#   plot(as.cimg(img))
#   dx = as.matrix(imgradient(as.cimg(img), "x"))
#   dx = (dx-mean(dx[which(!is.na(dx), arr.ind = TRUE)]))/sd(dx[which(!is.na(dx), arr.ind = TRUE)])
#   dy = as.matrix(imgradient(as.cimg(img), "y"))
#   dx = (dy-mean(dy[which(!is.na(dy), arr.ind = TRUE)]))/sd(dy[which(!is.na(dy), arr.ind = TRUE)])
#   
#   iter = 0
#   converge = FALSE
#   convergence = VTT
#   
#   while(converge == FALSE){
#     direction = c(-dx[VTT[1], VTT[2]], -dy[VTT[1], VTT[2]])
#     
#     VTT = c(VTT[1] + step * direction[1], VTT[2] + step * direction[2])
#     
#     convergence = rbind(convergence, VTT)
#     
#     points(VTT[1], VTT[2])
#     
#     iter = iter + 1
#     if(sqrt(dx[VTT[1], VTT[2]]^2 + dy[VTT[1], VTT[2]]^2) < tol || iter >= max_iter){
#       converge = TRUE
#       if(iter >= max_iter)
#         warning("Max iter atteint")
#     }
#   }
#   
#   return(list(point = VTT, historique = convergence))
# }


# Trouver_vallee = function(X){
#   if(is.null(nrow(X))){
#     return(decoupe[X[1], X[2]])
#   }
#   return(sum(decoupe[X]))
# }

Trouver_vallee = function(X, img){
  img[which(is.na(img) == TRUE, arr.ind = TRUE)] = 1
  #print(X)
  if(X[1] > nrow(img) || X[2] > ncol(img)){
    return(img[min(c(X[1], nrow(img))), min(c(X[2], ncol(img)))])
  }
  return(img[X[1], X[2]])
}

points_arc_anneaux = function(ridges, k=0.15){
  shape = c()
  if(length(which(is.na(ridges))) > 0)
    ridges[which(is.na(ridges), arr.ind = TRUE)] = 0
  
  lab = unique(as.vector(label(as.cimg(ridges < k))))
  
  for(l in lab[c(-1, -2)]){
    if(nrow(which(as.matrix(label(as.cimg(ridges < k))) == l, arr.ind = TRUE)) > 1000)
      shape = c(shape, list(which(as.matrix(label(as.cimg(ridges < k))) == l, arr.ind = TRUE)))
  }
  
  return(shape)
}

cost = function(point, courbe){
  cout = 0
  for(i in 1:nrow(point)){
    cout = cout + min(euclidien(point[i,], courbe))
  }
  return(cout)
}

arc_anneaux = function(ridges, k = 0.15, graph = FALSE){
  points_arcs = points_arc_anneaux(ridges, k)
  arc = c()
  for(i in 1:length(points_arcs)){
    model = lm(y ~ poly(x, 3), data = data.frame(x = points_arcs[[i]][,1], y = points_arcs[[i]][,2]))
    
    model_inv = lm(y ~ poly(x, 3), data = data.frame(x = points_arcs[[i]][,2], y = points_arcs[[i]][,1]))
    
    arc_xy = cbind(seq(min(points_arcs[[i]][,1]), max(points_arcs[[i]][,1]), len = 1000), predict(model, data.frame(x = seq(min(points_arcs[[i]][,1]), max(points_arcs[[i]][,1]), len = 1000))))
    
    arc_yx = cbind(predict(model_inv, data.frame(x = seq(min(points_arcs[[i]][,2]), max(points_arcs[[i]][,2]), len = 1000))), seq(min(points_arcs[[i]][,2]), max(points_arcs[[i]][,2]), len = 1000))
    
    if(cost(points_arcs[[i]], arc_xy) < cost(points_arcs[[i]], arc_yx)){
      arc = c(arc, list(arc_xy))
    } else{
      arc = c(arc, list(arc_yx))
    }
  }
  
  if(graph){
    plot(as.cimg(decoupe))
    for(i in 1:length(a)){
      points(points_arcs[[i]], col = "red", pch = ".")
      lines(arc[[i]], col = "blue", lwd = 2)
    }
  }
  
  return(arc)
}

perpendiculaire = function(courbe){
  #x = min(courbe[,1]):max(courbe[,1])
  perp = c()
  for(i in 1:(nrow(courbe)-1)){
    x = (courbe[i,1]-500):(courbe[i,1]+500)
    m = (courbe[i+1, 2] - courbe[i, 2])/(courbe[i+1, 1] - courbe[i, 1])
    m = -1/m
    
    perp = rbind(perp, list(x, m*(x - courbe[i, 1]) + courbe[i, 2]))
  }
  return(perp)
}

# library(viridis)
# 
# couleur = viridis(length(arc))
# 
# for(i in 1:length(arc)){
#   perp = perpendiculaire(arc[[i]])
#   
#   for(indice in unique(round(seq(1, nrow(perp), len = 100)))){
#     points(arc[[i]][indice,1], arc[[i]][indice,2], col = couleur[i])
#     lines(perp[[indice,1]], perp[[indice,2]], col = couleur[i])
#   }
#   
# }

dx = as.matrix(imgradient(as.cimg(decoupe), "x"))
dy = as.matrix(imgradient(as.cimg(decoupe), "y"))

dx = (dx-mean(dx))/sd(dx)
dy = (dy-mean(dy))/sd(dy)

plot(as.cimg(decoupe))
point = matrix(Contour[1,],1,2)
iter = 0
while(iter < 1000){
  iter = iter+1
  
  point = rbind(point, c(point[iter,1] + dx[point[iter,1], point[iter,2]], point[iter,2] + dy[point[iter,1], point[iter,2]]))
  # if(iter <= 1)
  #   vois = voisinage(point)
  # else
  #   vois = voisinage(point[nrow(point), ])
  # 
  # if(iter > 1){
  #   #print(point)
  #   if(iter == 815){
  #     browser()
  #   }
  #   while(coordinate_in_matrix(point, vois[which.max(abs(decoupe[vois] - decoupe[rbind(point[nrow(point), ])])), ])$logic){
  #     vois = vois[-which.max(abs(decoupe[vois] - decoupe[rbind(point[nrow(point), ])])), ]
  #   }
  #   point = rbind(point, vois[which.max(abs(decoupe[vois] - decoupe[rbind(point[nrow(point), ])])), ])
  # } else{
  #   point = rbind(point, vois[which.max(abs(decoupe[vois] - decoupe[point[1], point[2]])), ])
  # }

  
  
  
  #points(point[iter, 1], point[iter, 2], col = "red", pch = ".", cex = 2)
}

# a = c()
# blur = as.matrix(medianblur(as.cimg(decoupe), 15))
# for(i in 1:nrow(VTT)){
#   a = rbind(a, optim(VTT[i,], Trouver_vallee, img = decoupe, control = c("maxit" = 1e6))$par)
# }
# 
# ridges = load.image("C:/Users/THIBODEAUJU/Documents/Otolithe/ridge.jpg")
# 
# ridges = as.matrix(grayscale(ridges))
# ridges = resize(as.cimg(ridges), nrow(decoupe), ncol(decoupe))
# 
# for(i in 1:nrow(decoupe)){
#   ridges[i, which(is.na(decoupe[i,]) == TRUE, arr.ind = TRUE)] = NA
# }


