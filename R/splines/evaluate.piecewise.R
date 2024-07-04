evaluate.piecewise <- function(p, x){
   # EVALUATE.PIECEWISE - Evaluate a 'piecewise' object using Horner's method.

   # Piecewise polynomial knots:
   k <- p$knots
   
   # Initialize result matrix:
   y <- NA*x
   
   # Loop over piecewise intervals:
   for (i in 1:(length(k)-1)){
      if ((k[i] - k[i+1]) != 0){
         if (k[i] < k[i+1]){
            index <- (x >= k[i]) & (x <= k[i+1])
         }else{
            index <- (x <= k[i]) & (x >= k[i+1])
         }
         a <- p$coefficients[i, ]
         if (!all(a == 0)){
            if ((a[1] > 0) & all(a[2:length(a)] == 0)){
               y[index] <- a[1]
            }else{
               b <- rep(a[length(a)], sum(index))
               if (length(a) > 1){
                  for (j in 2:length(a)){
                     b <- a[length(a)-j+1] + b * x[index]
                  }
               }
               y[index] <- b
            }
         }else{
            y[index] <- 0
         }
      }
   }

   return(y)
}
