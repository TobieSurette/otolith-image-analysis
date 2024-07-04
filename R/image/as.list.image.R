as.list.image <- function(I){
   # AS.LIST.IMAGE - Convert 'image' object to list.

   class(I) <- "list"
   return(I)
}
