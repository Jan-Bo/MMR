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


#' Lag Time-series Function
#' 
#' Compute a lagged version of a time series, shifting the time forward or backward by a given number of observations.
#' 
#' @param x An ordered numeric vector
#' @param k The number of lags (in units of observations). Positive \code{k} shifts the time-series forward while a negative \code{k} shifts the time-series backward.
#' @param pad A value to added
#' 
#' @author Gabriel Mohanna
#' 
#' @examples
#' lagpad(x=1:10, k=2)
#' lagpad(x=1:10, k=-2)
#' 
#' @export
lagpad <- function(x, k=0, pad=0) {
    if(k>=0) {
        c(rep(pad, k), x[1:(length(x)-min(k,length(x)))])[1:length(x)]
    } else {
        if(-k>=length(x)) {
            rep(pad, length(x))
        } else {
            c(x[(abs(k)+1):length(x)], rep(pad, -k))
        }
    }
}
