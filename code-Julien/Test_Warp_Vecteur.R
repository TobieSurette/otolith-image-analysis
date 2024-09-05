#R = colMeans(Sweep2[(axe_principal[5]-20):(axe_principal[5]-14),])
#Ref = colMeans(Sweep2[(axe_principal[5]-3):(axe_principal[5]+3),])

Ref = Sweep2[axe_principal[5],]
R = Sweep2[axe_principal[5]-50,]

plot(Ref, type="l")
lines(R, col="orange")

init = c(a1 = 0.7, b1 = 0.1)

meilleurS=10000000000000

for(i in seq(0.1,10,0.1)){
  for(j in seq(0.1,10,0.1)){
    init = c(a1=i,b1=j)
    
    optimisation = optim(init, Beta, r=R, ref=Ref)
    S = optimisation$value
    
    if(S<meilleurS){
      meilleurS = S
      param = optimisation$par
    }
  }
}

a1=1.61#param["a1"]
b1=1.39#param["b1"]

new_R = approx(1:length(R),R,Fonction_Warp(a1,a2,length(R)))$y


lines(new_R, col="blue")

legend(x="topright",legend = c("Référence", "Vecteur non-transformé", "Vecteur transformé"),
       fill = c("black", "orange", "blue"), cex = 0.7)

