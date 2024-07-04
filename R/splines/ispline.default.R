ispline.default <- function(...){
   # ISPLINE.DEFAULT - Create an 'ispline' object.

   # Integrate B-splien basis function:
   r <- integrate(bspline(...))

   # Normalize B-spline basis function:
   
   coef(r) <- (coef(r) / sum(coef(r)))
   for (i in 1:length(r$basis)){
      #v <- evaluate(r$basis[[i]], r$basis[[i]]$knots[length(r$basis[[i]]$knots)])
      #if (v != 0) r$basis[[i]] <- (1 / v) * r$basis[[i]]
      r$basis[[i]] <- r$dim * r$basis[[i]]
   }

   class(r) <- c("ispline", class(r))
   
   return(r)
}
