# Seulement testé pour l'image plaice.jpg

rangee = c()

temp = Contour

for(j in (min(Contour[,2])) : (max(Contour[,2]))){
  ligne = c()
  for(i in min(Contour[,1]) : max(Contour[,1])){
    if(decoupe[i,j] != 0){
      ligne = rbind(ligne, c(i, decoupe[i,j]))
    }
  }
  rangee = rbind(rangee,c(list(ligne), j))
}

colonne = c()
for(i in min(Contour[,1]) : max(Contour[,1])){
  ligne = c()
  for(j in (min(Contour[,2])) : (max(Contour[,2]))){
    if(decoupe[i,j] != 0){
      ligne = rbind(ligne, c(j, decoupe[i,j]))
    }
  }
  colonne = rbind(colonne,c(list(ligne), i))
}

Transformation_0_1 = function(v){
  a = min(v)
  b = max(v)
  return((v-a)/(b-a))
}

norm_gauche_droite = function(point, v){
  n = nrow(v)
  pi = which(v[,1] == point)
  at = min(v[,1])
  bt = max(v[,1])
  a1 = at
  a2 = point
  b1 = point
  b2 = bt
  temp = v
  temp[1:pi, 1] = (v[1:pi,1] - a1)/(2*(b1 - a1))
  temp[pi:n, 1] = temp[pi,1] + (v[pi:n,1] - a2)/(2*(b2 - a2))
  
  v[,1] = temp[,1]#*(bt-at) + at
  return(v)
}

sym_warp = function(a, b, x){
  c = min(x[,1])
  k = max(x[,1])
  
  x[,1] = pbeta((x[,1] - c) / (k - c), a, b)
  x[,1] = x[,1]*(k-c) + c
  
  return(x)
}

sym_warp_objective = function(rangee){
  cout = c()
  for(point in rangee[c(-1,-nrow(rangee)),1]){
    r = norm_gauche_droite(point, rangee)
    
    milieu = which(r[,1] == 0.5)
    gauche = spline(x = r[1:milieu, 1], y = r[1:milieu, 2], n = 100)
    droite = spline(x = r[milieu:nrow(rangee), 1], y = r[milieu:nrow(rangee), 2], n = 100)
    
    cout = c(cout, sum(gauche$y - droite$y[length(droite$y):1])^2)
  }
  
  return(rangee[which.min(cout)+1,1])
}

symetrie = function(rangee, is.colonne = FALSE){
  p = c()
  if(is.colonne){
    for(i in 1:nrow(rangee)){
      p = rbind(p, c(rangee[[i,2]], sym_warp_objective(rangee[[i,1]])))
    }
  } else{
    for(i in 1:nrow(rangee)){
      p = rbind(p, c(sym_warp_objective(rangee[[i,1]]), rangee[[i,2]]))
    }
  }
  
  return(p)
}

Trouver_noyau_sym = function(rangee, colonne, img){
  r = symetrie(rangee)
  c = symetrie(colonne, is.colonne = TRUE)
  
  plot(as.cimg(img))
  lines(r[,1], r[,2], col = "red")
  lines(c[,1], c[,2], col = "blue")
  
  m1 = lm(y~x, data.frame(x=r[,1], y=r[,2]))
  abline(m1, col = "red3", lwd = 2)
  m2 = lm(y~x, data.frame(x=c[,1], y=c[,2]))
  abline(m2, col = "blue3", lwd = 2)
  
  x = (m1$coefficients[1] - m2$coefficients[1])/(m2$coefficients[2] - m1$coefficients[2])
  y = m1$coefficients[2]*x + m1$coefficients[1]
  
  n = c(x, y)
  points(n[1], n[2], lwd = 3)
  
  return(list(n, r, c))
}


