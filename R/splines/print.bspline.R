print.bspline <- function(b){
   # PRINT.BSPLINE - Display a 'bspline' object on the R console.
   
   cat("\n")
   cat(paste("'", class(b)[1], "' object:\n", sep = ""))
   str <- gsub(" ", "", formatC(b$knots, digits = 5), fixed = FALSE)
   cat(paste("   knots        : [", paste(str, collapse = ", "), "]\n", sep = ""))
   str <- gsub(" ", "", formatC(b$coefficients, digits = 5), fixed = FALSE)
   cat(paste("   coefficients : [", paste(str, collapse = ", "), "]\n", sep = ""))
   cat(paste("   dim          : ", b$dim, "\n", sep = ""))
   cat(paste("   degree       : ", b$degree, "\n", sep = ""))
   cat(paste("   type         : '", b$type, "'\n", sep = ""))
   str <- gsub(" ", "", format(b$is.fixed), fixed = FALSE)
   cat(paste("   is.fixed     : [", paste(str, collapse = ", "), "]\n", sep = ""))
}
