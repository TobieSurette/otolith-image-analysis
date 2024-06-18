as.character.polynomial <- function(p){
   # AS.CHARACTER - Convert a 'polynomial' object to a character string.
   
   # Extract polynomial coefficients:
   beta <- p$coefficients
   
   # Build monomial string:
   str <- c("", "x")
   if (length(beta) >= 2) str <- c(str, paste("x^", as.character(2:length(beta)), sep = ""))

   # Remove zero coefficients:
   str <- str[beta != 0]
   beta <- beta[beta != 0]
   
   # Append coefficients:
   str <- paste(beta, str, sep = "")
   
   # Add plus signs:
   str <- paste(str, collapse = " + ")
   
   return(str)
}
