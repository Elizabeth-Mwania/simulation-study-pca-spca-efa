# # Minimal Lp (L1) rotation function for simulation
# LpRotation <- function(L, p = 1, maxit = 500, tol = 1e-5) {
#   # L: loadings (variables x factors)
#   # p: Lp norm (set p=1 for L1)
#   # maxit: maximum iterations
#   # tol: convergence tolerance
  
#   # Initialization
#   d <- ncol(L)
#   p_vars <- nrow(L)
#   R <- diag(d)
#   obj_old <- sum(abs(L %*% R)^p)
  
#   for (iter in 1:maxit) {
#     # Compute weights for IRLS
#     L_rot <- L %*% R
#     W <- (abs(L_rot) + 1e-8)^(p - 2)
#     W[is.infinite(W)] <- 1e8  # Avoid Inf weights
    
#     # Weighted SVD update
#     for (j in 1:d) {
#       w <- W[, j]
#       z <- L[, j]
#       # Weighted least squares
#       R[, j] <- z / sqrt(sum((z^2) * w))
#     }
#     # Orthonormalize R (QR decomposition)
#     qrR <- qr(R)
#     R <- qr.Q(qrR)
    
#     # Check convergence
#     obj_new <- sum(abs(L %*% R)^p)
#     if (abs(obj_new - obj_old) < tol) break
#     obj_old <- obj_new
#   }
#   # Return rotated loadings
#   return(L %*% R)
# }


# Minimal Lp (L1) rotation function for simulation
LpRotation <- function(L, lp_norm = 1, maxit = 500, tol = 1e-5) {
  # L: loadings (variables x factors)
  # lp_norm: Lp norm (set lp_norm=1 for L1)
  # maxit: maximum iterations
  # tol: convergence tolerance

  d <- ncol(L)
  R <- diag(d)
  obj_old <- sum(abs(L %*% R)^lp_norm)

  for (iter in 1:maxit) {
    L_rot <- L %*% R
    W <- (abs(L_rot) + 1e-8)^(lp_norm - 2)
    W[is.infinite(W)] <- 1e8

    for (j in 1:d) {
      w <- W[, j]
      z <- L[, j]
      # Weighted least squares
      R[, j] <- z / sqrt(sum((z^2) * w))
    }
    # Orthonormalize R (QR decomposition)
    qrR <- qr(R)
    R <- qr.Q(qrR)

    obj_new <- sum(abs(L %*% R)^lp_norm)
    if (abs(obj_new - obj_old) < tol) break
    obj_old <- obj_new
  }
  return(L %*% R)
}
