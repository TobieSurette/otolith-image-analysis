print.image <- function(I, indent = 0, cr = TRUE){
   # PRINT.IMAGE - Print an 'image' object to the R console.
   
   # Build indentation str:
   istr <- paste(rep(" ", indent), collapse = "")

   # Carriage return:
   if (cr) cat("\n")
   
   cat(paste(istr, "'", class(I)[1], "' object:\n", sep = ""))

   str <- gsub(" ", "", formatC(c(I$x[1], I$x[length(I$x)]), digits = 5), fixed = FALSE)
   cat(paste(istr, "   x : [", paste(str, collapse = ", ..., ") , "], step = ",  mean(unique(diff(I$x))), " (n = ",  length(I$x), ")\n", sep = ""))
   
   str <- gsub(" ", "", formatC(c(I$y[1], I$y[length(I$y)]), digits = 5), fixed = FALSE)
   cat(paste(istr, "   y : [", paste(str, collapse = ", ..., ") , "], step = ",  mean(unique(diff(I$y))), " (n = ",  length(I$y), ")\n", sep = ""))

   str <- gsub(" ", "", formatC(dim(I$z), digits = 5), fixed = FALSE)
   cat(paste(istr, "   z : [", paste(str, collapse = " x ") , "]\n", sep = ""))
}
