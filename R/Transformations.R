#' Adstock transformation function.
#' 
#' Computes an adstock transformed variable
#' 
#' @param x An ordered numeric none-na vector
#' @param rate A decimal between 0-1
#' 
#' @author Gabriel Mohanna
#' 
#' @details
#' Transforms variable 'x' into adstock version with rate 'rate'.
#' 
#' @seealso
#' \url{http://analyticsartist.wordpress.com/2013/11/02/calculating-adstock-effect/}
#' 
#' @examples
#' adstock(x=c(10, 14, 30, 32, 0, 0, 0, 10), rate=.35)
#' 
#' @export
adstock <- function(x, rate=0){
    return(as.numeric(filter(x=x, filter=rate, method="recursive")))
}

#' MMR transformation function.
#' 
#' Function Not Yet Done. Transforms advertising into adstock, power and lag version.
#' 
#' @param x An ordered numeric vector
#' @param adstock A decimal between 0-1
#' @param power A decimal between 0-1
#' @param lag An positive integer
#' 
#' @author Gabriel Mohanna
#' 
#' @details
#' Transforms advertising into adstock version
MMR.Transform <- function(x, adstock=0, power=1, lag=0){
    return(as.numeric(filter(x=x, filter=adstock, method="recursive"))^power)
}
