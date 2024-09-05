
i=0
decoupe_blur = medianblur(as.cimg(decoupe),20)
lvl_set_up = lvl_set_low = array(,c(nrow(decoupe), ncol(decoupe), 51))
for(lambda in seq(0,1,len=51)){
  i=i+1
  lvl_set = Level_set(as.cimg(decoupe_blur), lambda)
  lvl_set_up[,,i] = lvl_set$upper
  lvl_set_low[,,i] = lvl_set$lower
}

a = build_shape_tree(lvl_set_up)

axe_p = Principal_Axis(Contour)

plot(as.cimg(decoupe))
couleur = c("red", "blue", "orange3", "green", "purple")
for(i in 1:nrow(Arbre)){
  points(Arbre[i,]$shape, col=couleur[mod(i, 5)+1], pch = ".", cex = 1.5)
}
lines(axe_p, lwd = 2)

indice = c()
for(i in 1:nrow(Arbre)){
  if(is.null(Arbre[[i,3]])){
    indice = c(indice, i)
  }
}
aMax = Arbre[indice, ]

plot(as.cimg(decoupe))
couleur = c("red", "blue", "orange3", "green", "purple")
for(i in 1:nrow(aMax)){
  points(aMax[i,]$shape, col=couleur[mod(i, 5)+1], pch = ".")
}
lines(axe_p, lwd = 2)

fi = c()
for(i in 1:nrow(aMax)){
  fi = rbind(fi, c(Fi(aMax[[i,1]], lvl_set_up[,,aMax[[i,2]]], Contour), shape = i))
}


e1 = c()
for(i in 1:nrow(aMax)){
  e1 = c(e1, E1(fi[i,], fi))
}

fi = fi[order(cbind(fi,e1)[,5], decreasing = FALSE), ]

plot(as.cimg(decoupe))
for(i in 1:nrow(aMax)){
  points(aMax[[which(fi[,4] == i),1]], col="red", pch=".")
  Sys.sleep(2)
}

for(i in 1:length(AMax2)){
  print(Similitude_des_rayons(colMeans(AMax2[[paste0("shape",fi[i,"shape"])]]$shape), decoupe, Contour))
}

for(i in 1:nrow(a[[2,1]])){
  points(a[[2,1]][i,1]/5, a[[2,1]][i,2]/5, col="red", pch = ".")
  Sys.sleep(0.5)
}

