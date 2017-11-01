#' @title Tridiagonal matrix
#' @description Generate tridiagonal matrix Q based on three input vectors.
#'
#' @param upper The upper diagonal vector.
#' @param lower The lower diagonal vector.
#' @param main The main diagonal vector.
#' @return A tridiagonal matrix is returned.
#' @examples a = c(1:7)^2
#' b = c(1:7)^2
#' c = -c(1:8)^2
#' tridiag(b, a, c)

#' @export
branch <- function(alpha, N) {
    p0 = alpha/2
    
    p = c(0, (2 - alpha)/2^(2:(N - 1)))
    
    Q11 = upper.tri(toeplitz(c(-1, p[-1])), diag = T) * toeplitz(c(-1, 
        p[-1]))
    Q11[cbind(2:(N - 1), 1:(N - 2))] = p0
    
    Q12 = c(rep(0, N - 2), p0)
    Q21 = (2 - alpha)/2 - c(rev(cumsum(p[-1])), 0)
    Q22 = -p0
    
    Q = rbind(cbind(Q11, Q21), c(Q12, Q22)) * (1:N)
    
    return(Q)
}
