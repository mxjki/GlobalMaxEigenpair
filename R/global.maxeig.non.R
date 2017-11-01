#' @title The maximal eigenpair of the matrices with non-negative off-diagonal elements.
#'
#' @description Use global algorithm to calculate the maximal eigenpair for the matrices 
#'              with non-negative off-diagonal elements.
#'
#' @param Q The input general matrix.
#' @param complex Whether the input matrix A is complex or not.
#' @param w0 The unnormalized initial vector \eqn{w0}.
#' @param z0 The type of initial \eqn{z_0} used to calculate the approximation of \eqn{\rho(A)}.
#'        There are four types: 'fixed', 'convex_1', 'convex_2' and 'numeric' corresponding to
#'        four choices of \eqn{z_0} in paper.
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
#' @seealso \code{\link{global.maxeig.general}} for the maximal eigenpair of the 
#'           general matrices.

#' @examples nn = c(8, 16, 32, 50, 100, 500, 1000, 5000, 10^4)-1
#'
#' iter.ray.quot.non = rep(0, 9)
#' iter.mod.ray.quot.non = rep(0, 5)
#' iter.conv.comb.ray.quot.non = rep(0, 9)
#'
#' for (i in 1:9) {
#'   
#'  Q=matrix(0,nn[i]+1,nn[i]+1)
#'  a = 1/seq(2,nn[i]+1)
#'  b = c(1:nn[i])
#'  
#'  Q[,1]=c(-1,a)
#'  diag(Q)=c(-1,-(a[1:(nn[i] - 1)]+b[2:nn[i]]),-a[nn[i]]-nn[i]-1)
#'  Q[cbind(seq(1,nn[i]),seq(2,nn[i]+1))]=b
#'  
#'  global.maxeig.ray.quot.non = global.maxeig.non(Q, alg = 'ray.quot.non', w0 = rep(1, dim(Q)[1]), z0 = 'fixed')
#'  
#'  iter.ray.quot.non[i]=global.maxeig.ray.quot.non$iter
#'  
#'  if (i >4){
#'    global.maxeig.mod.ray.quot.non = global.maxeig.non(Q, alg = 'mod.ray.quot.non', w0 = rep(1, dim(Q)[1]), z0 = 'fixed')
#'    
#'    iter.mod.ray.quot.non[i-4]=global.maxeig.mod.ray.quot.non$iter
#'    
#'  }
#'  
#'  global.maxeig.conv.comb.non = global.maxeig.non(Q, alg = 'ray.quot.non', w0 = rep(1, dim(Q)[1]), z0 = 'convex_1', xi = 0.34189)
#'  
#'}
#'
#' print(iter.ray.quot.non)
#' print(iter.mod.ray.quot.non)
#' print(iter.conv.comb.ray.quot.non)

#' @export
global.maxeig.non = function(Q, alg = "ray.quot.non", w0 = NULL, 
    z0 = NULL, z0numeric, xi = 1, digit.thresh = 6) {
    
    if (xi < 0 | xi > 1) 
        stop("The coefficient xi should between 0 and 1!")
    
    N = dim(Q)[1]
    
    v0 = w0/sqrt(sum(w0^2))
    zstart = switch(z0, fixed = 0, convex_1 = xi * sum(v0 * (-Q %*% 
        v0)), convex_2 = xi * sum(v0 * (-Q %*% v0)) + (1 - xi) * 
        max(-Q %*% w0), numeric = z0numeric)
    
    maxeig = switch(alg, ray.quot.non = ray.quot.non(Q = Q, w0 = w0, 
        zstart = zstart, digit.thresh = digit.thresh), shifted.inv.non = shifted.inv.non(Q = Q, 
        w0 = w0, zstart = zstart, digit.thresh = digit.thresh), mod.ray.quot.non = mod.ray.quot.non(Q = Q, 
        w0 = w0, zstart = zstart, digit.thresh = digit.thresh))
    
    return(list(z = unlist(maxeig$z), v = maxeig$v, iter = maxeig$iter))
}
