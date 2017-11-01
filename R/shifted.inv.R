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
shifted.inv = function(A, complex = F, w0, zstart, digit.thresh = 6, 
    digit.thresh2 = 8) {
    z = list()
    rz = list()
    y = list()
    v = list()
    w = list()
    
    v[[1]] = w0/sqrt(sum(w0^2))
    
    z[[1]] = zstart
    rz[[1]] = round(zstart, digit.thresh)
    y[[1]] = zstart
    
    ratio = 1
    iter = 0
    if (!complex) {
        while (ratio >= 10^(-digit.thresh)) {
            iter = iter + 1
            w = append(w, list(solve(z[[iter]] * diag(1, length(w0)) - 
                A, v[[iter]], tol = 1e-100)))
            v = append(v, list(w[[iter]]/sqrt(sum(w[[iter]]^2))))
            
            z = append(z, list(max((A %*% w[[iter]])/w[[iter]])))
            
            ratio = abs(round(z[[iter + 1]], digit.thresh) - round(z[[iter]], 
                digit.thresh))
            
            rz[[iter + 1]] = round(z[[iter + 1]], digit.thresh)
        }
    } else {
        while (ratio >= 10^(-digit.thresh)) {
            iter = iter + 1
            w = append(w, list(solve(z[[iter]] * diag(1, length(w0)) - 
                A, v[[iter]], tol = 1e-100)))
            v = append(v, list(w[[iter]]/sqrt(sum(Conj(w[[iter]]) * 
                w[[iter]]))))
            
            z = append(z, list(max((Re(A) %*% Re(w[[iter]]))/Re(w[[iter]]))))
            
            y = append(y, list(sum(Conj(v[[iter + 1]]) * (A %*% v[[iter + 
                1]]))))
            
            ratio = round(sqrt((Re(y[[iter + 1]]) - Re(y[[iter]]))^2 + 
                (Im(y[[iter + 1]]) - Im(y[[iter]]))^2), digit.thresh)
            
            rz[[iter + 1]] = round(Re(y[[iter + 1]]), digit.thresh) + 
                round(Im(y[[iter + 1]]), digit.thresh2) * (0 + (0 + 
                  (0 + (0 + (0 + (0 + (0 + (0 + (0 + (0 + (0 + (0 + 
                    (0+1i)))))))))))))
        }
    }
    
    if (ratio == 0) {
        v = v[-(iter + 1)]
        rz = rz[-(iter + 1)]
        
        iter = iter - 1
    }
    
    return(list(v = v, z = rz, iter = iter))
}
