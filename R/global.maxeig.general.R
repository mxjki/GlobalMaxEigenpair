#' @title The maximal eigenpair of the general matrices.
#'
#' @description Use global algorithms to calculate the maximal eigenpair for the general matrices.
#'
#' @param A The input general matrix.
#' @param complex Whether the input matrix A is complex or not.
#' @param alg The detailed global algorithm used to calculate the maximal eigenpair.
#' @param w0 The unnormalized initial vector \eqn{w0}.
#' @param z0 The type of initial \eqn{z_0} used to calculate the approximation of \eqn{\rho(A)}.
#'        There are three types: 'fixed', 'convex' and 'numeric' corresponding to three choices
#'        of \eqn{z_0} in paper.
#' @param z0numeric The numerical value assigned to initial \eqn{z_0} as an approximation of
#'        \eqn{\rho(A)} when z_0='numeric'.
#' @param xi The coefficient used to form the convex combination of \eqn{\min(A*w0)} and
#'        \eqn{(v_0,A*v_0)_\mu}, it should between 0 and 1.
#' @param digit.thresh The precise level of output results.
#'
#' @return A list of eigenpair object are returned, with components \eqn{z}, \eqn{v} and \eqn{iter}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating sequence of the corresponding eigenvector.}
#' \item{iter}{The number of iterations.}
#' 
#' @seealso \code{\link{global.maxeig.non}} for maximal eigenpair of the matrices with non-negative 
#'          off-diagonal elements.

#' @examples 
#' b4 = c(0.01, 1, 100, 10^4)
#' digits = c(9, 7, 6, 6)

#' for (i in 1:4) {
#'   A = matrix(c(-3, 4, 0, 10, 0, 2, -7, 5, 0, 0, 0, 3, -5, 0, 0,
#'                1, 0, 0, -16, 11, 0, 0, 0, 6, -11 - b4[i]), 5, 5)
#'                
#'   print(-global.maxeig.general(A, alg = 'ray.quot', w0 = rep(1, dim(A)[1]), z0 = 'fixed', digit.thresh = digits[i])$z[-1])
#'   
#'   print(-global.maxeig.general(A, alg = 'shifted.inv', w0 = rep(1, dim(A)[1]), z0 = 'fixed', digit.thresh = digits[i])$z[-1])
#' }

#' @export
global.maxeig.general = function(A, complex = F, alg = "ray.quot", 
    w0 = NULL, z0 = NULL, z0numeric, xi = 1, digit.thresh = 6) {
    
    if (xi < 0 | xi > 1) 
        stop("The coefficient xi should between 0 and 1!")
    
    N = dim(A)[1]
    
    v0 = w0/sqrt(sum(w0^2))
    if (!complex) {
        zstart = switch(z0, fixed = max(A %*% w0), convex = xi * 
            sum(v0 * (A %*% v0)) + (1 - xi) * min(A %*% w0), numeric = z0numeric)
        
        maxeig = switch(alg, ray.quot = ray.quot(A = A, w0 = w0, 
            zstart = zstart, digit.thresh = digit.thresh), shifted.inv = shifted.inv(A = A, 
            w0 = w0, zstart = zstart, digit.thresh = digit.thresh))
        
        
    } else {
        zstart = max(Re(A) %*% w0)
        
        maxeig = switch(alg, ray.quot = ray.quot(A = A, w0 = w0, 
            zstart = zstart, digit.thresh = digit.thresh), shifted.inv = shifted.inv(A = A, 
            complex = T, w0 = w0, zstart = zstart, digit.thresh = digit.thresh))
    }
    
    
    return(list(z = unlist(maxeig$z), v = maxeig$v, iter = maxeig$iter))
}
