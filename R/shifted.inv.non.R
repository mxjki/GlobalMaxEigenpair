#' @title Shift inverse iteration
#' @description Rayleigh quotient iteration algorithm to computing the maximal eigenpair of
#' matrix Q.
#'
#' @param Q The input matrix to find the maximal eigenpair.
#' @param v0_tilde The unnormalized initial vector \eqn{\tilde{v0}}.
#' @param zstart The initial \eqn{z_0} as an approximation of \eqn{\rho(Q)}.
#' @param digit.thresh The precise level of output results.
#' @return A list of eigenpair object are returned, with components '\eqn{z}' and '\eqn{v}'.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating sequence of the corresponding eigenvector.}
#'
#' @examples
#' Q = matrix(c(1, 1, 3, 2, 2, 2, 3, 1, 1), 3, 3)
#' ray.quot(Q, mu=rep(1,dim(Q)[1]), v0_tilde=rep(1,dim(Q)[1]), zstart=6,
#'  digit.thresh = 6)

#' @export
shifted.inv.non = function(Q, w0, zstart, digit.thresh = 6) {
    
    z = list()
    rz = list()
    v = list()
    w = list()
    
    v[[1]] = w0/sqrt(sum(w0^2))
    
    z[[1]] = zstart
    rz[[1]] = round(zstart, digit.thresh)
    
    ratio = 1
    iter = 0
    while (ratio >= 10^(-digit.thresh)) {
        iter = iter + 1
        w = append(w, list(solve(-Q - z[[iter]] * diag(1, length(w0)), 
            v[[iter]], tol = 1e-100)))
        v = append(v, list(w[[iter]]/sqrt(sum(w[[iter]]^2))))
        
        z = append(z, list(min((-Q %*% w[[iter]])/w[[iter]])))
        
        ratio = abs(round(z[[iter + 1]], digit.thresh) - round(z[[iter]], 
            digit.thresh))
        
        rz[[iter + 1]] = round(z[[iter + 1]], digit.thresh)
    }
    
    if (ratio == 0) {
        v = v[-(iter + 1)]
        rz = rz[-(iter + 1)]
        
        iter = iter - 1
    }
    
    return(list(v = v, z = rz, iter = iter))
}
