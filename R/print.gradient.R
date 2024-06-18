print.gradient <- function(I, indent = 0, cr = TRUE){
   # PRINT.GRADIENT - Print a 'gradient' object to the R console.
   
   # Build indentation str:
   istr <- paste(rep(" ", indent), collapse = "")

   # Carriage return:
   if (cr) cat("\n")

   vars <- names(I)
   
   cat(paste(istr, "'", vars[1], "'\n", sep = ""))
   print(I[[vars[1]]], indent + 3, FALSE)
   cat(paste(istr, "'", vars[2], "'\n", sep = ""))
   print(I[[vars[2]]], indent + 3, FALSE)
   cat(paste(istr, "polar : ", G$polar, "\n", sep = ""))
}
