cluster.image <- function(I){
   L <- I * NA
   label <- 1
   
   # One row at a time:
   for (i in 1:nrow(I)){
      if (I[i,1]) L[i,1] <- label
      
      for (j in 2:ncol(I)){
         if (I[i,j] & !I[i,j-1]){
            label <- label + 1
            L[i,j] <- label
         }
         
         if (I[i,j] & I[i,j-1]) L[i,j] <- L[i,j-1]
      }
   }
   
   # Collapse labels when blobs touch:
   for (i in 2:nrow(I)){
      print(i)
      for (j in 2:ncol(I)){
         if (I[i,j] & I[i,j-1]){
            if (L[i,j] > L[i,j-1]){
               ix <- L == L[i,j]
               L[ix] <- L[i,j-1]
            }
            if (L[i,j] < L[i,j-1]){
               ix <- L == L[i,j-1]
               L[ix] <- L[i,j]
            }
         }
         if (I[i,j] & I[i-1,j]){
            if (L[i,j] > L[i-1,j]){
               ix <- L == L[i,j]
               L[ix] <- L[i-1,j] 
            }
            if (L[i,j] < L[i-1,j]){
               ix <- L == L[i-1,j]
               L[ix] <- L[i,j] 
            }         
         }
      }
   }
   
   # Make sure the labels are sequential:
   t <- rev(sort(table(L)))
   V <- match(L, as.numeric(names(t)))
   dim(V) <- dim(L)
   L <- V
   
   return(L)
}


