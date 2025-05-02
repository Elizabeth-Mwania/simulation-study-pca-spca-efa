# --- Load Libraries ---
library(parallel)
library(doParallel)
library(foreach)
library(psych)
library(sparsepca)
library(GPArotation)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(corrplot)
library(vegan)

# --- Source Lp Rotation Function ---
source("cl1121_Lp_fun.R") # Ensure this path is correct

# --- Output Directory ---
output_dir <- "simulation_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --- Parallel Processing Setup ---
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# --- Data Generation Functions ---
create_loading_matrix <- function(p_vars, m_factors, alpha, gamma) {
  L <- matrix(0, nrow = p_vars, ncol = m_factors)
  alpha_p <- floor(alpha * p_vars)
  L[1:alpha_p, 1] <- runif(alpha_p, 0.7, 0.9)
  if (2 * alpha_p <= p_vars) {
    L[(alpha_p + 1):(2 * alpha_p), 2] <- runif(alpha_p, 0.7, 0.9)
  }
  if (2 * alpha_p < p_vars) {
    dense_indices <- (2 * alpha_p + 1):p_vars
    L[dense_indices, ] <- matrix(runif(length(dense_indices) * m_factors, 0.1, 0.4), 
                                nrow = length(dense_indices), ncol = m_factors)
  }
  if (gamma > 0) {
    for (j in 1:m_factors) {
      decay_factor <- j^(-gamma)
      L[, j] <- L[, j] * sqrt(decay_factor)
    }
  }
  return(L)
}

generate_data <- function(n_obs, p_vars, m_factors, alpha, gamma, sigma_sq) {
  L <- create_loading_matrix(p_vars, m_factors, alpha, gamma)
  F <- matrix(rnorm(n_obs * m_factors), nrow = n_obs, ncol = m_factors)
  E <- matrix(rnorm(n_obs * p_vars, 0, sqrt(sigma_sq)), nrow = n_obs, ncol = p_vars)
  X <- F %*% t(L) + E
  return(list(X = X, L = L, F = F, E = E))
}

# --- Lp Rotation Wrapper ---
apply_lp_rotation <- function(loadings, lp_p = 1, max_iter = 500, tol = 1e-5) {
  rotated <- tryCatch({
    LpRotation(loadings, lp_norm = lp_p, maxit = max_iter, tol = tol)
  }, error = function(e) {
    warning("LpRotation failed: ", e$message)
    matrix(NA, nrow = nrow(loadings), ncol = ncol(loadings))
  })
  return(rotated)
}

# --- PCA with Varimax ---
apply_pca <- function(X, m) {
  tryCatch({
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    pca_result <- prcomp(X_centered, center = FALSE)
    loadings <- pca_result$rotation[, 1:m, drop = FALSE]
    scores <- pca_result$x[, 1:m, drop = FALSE]
    rot_result <- GPArotation::Varimax(loadings, normalize = FALSE)
    loadings_rotated <- rot_result$loadings
    scores_rotated <- scores %*% rot_result$rotmat
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("PCA with Varimax failed: ", e$message)
    list(loadings = matrix(NA, nrow = ncol(X), ncol = m),
         scores = matrix(NA, nrow = nrow(X), ncol = m),
         X_reconstructed = matrix(NA, nrow = nrow(X), ncol = ncol(X)))
  })
}

# --- PCA with Lp Rotation ---
apply_pca_lp <- function(X, m, lp_p = 1) {
  tryCatch({
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    pca_result <- prcomp(X_centered, center = FALSE)
    loadings <- pca_result$rotation[, 1:m, drop = FALSE]
    scores <- pca_result$x[, 1:m, drop = FALSE]
    loadings_rotated <- apply_lp_rotation(loadings, lp_p = lp_p)
    X_reconstructed <- scores %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("PCA with Lp Rotation failed: ", e$message)
    list(loadings = matrix(NA, nrow = ncol(X), ncol = m),
         scores = matrix(NA, nrow = nrow(X), ncol = m),
         X_reconstructed = matrix(NA, nrow = nrow(X), ncol = ncol(X)))
  })
}

# --- EFA with Varimax ---
apply_efa <- function(X, m) {
  tryCatch({
    efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "ml")
    loadings <- as.matrix(efa_result_unrotated$loadings)
    scores <- efa_result_unrotated$scores
    rot_result <- GPArotation::Varimax(loadings, normalize = FALSE)
    loadings_rotated <- rot_result$loadings
    scores_rotated <- scores %*% rot_result$rotmat
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("EFA with Varimax failed: ", e$message)
    list(loadings = matrix(NA, nrow = ncol(X), ncol = m),
         scores = matrix(NA, nrow = nrow(X), ncol = m),
         X_reconstructed = matrix(NA, nrow = nrow(X), ncol = ncol(X)))
  })
}

# --- EFA with Lp Rotation ---
apply_efa_lp <- function(X, m, lp_p = 1) {
  tryCatch({
    efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "ml")
    loadings <- as.matrix(efa_result_unrotated$loadings)
    scores <- efa_result_unrotated$scores
    loadings_rotated <- apply_lp_rotation(loadings, lp_p = lp_p)
    X_reconstructed <- scores %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("EFA with Lp Rotation failed: ", e$message)
    list(loadings = matrix(NA, nrow = ncol(X), ncol = m),
         scores = matrix(NA, nrow = nrow(X), ncol = m),
         X_reconstructed = matrix(NA, nrow = nrow(X), ncol = ncol(X)))
  })
}

# --- Cross-validation for Sparse PCA Tuning ---
find_optimal_sparsity_cv <- function(X, m, n_folds = 5, alpha_grid = seq(0.1, 0.9, by = 0.1)) {
  tryCatch({
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    n <- nrow(X)
    fold_indices <- sample(rep(1:n_folds, length.out = n))
    cv_errors <- matrix(NA, nrow = n_folds, ncol = length(alpha_grid))
    colnames(cv_errors) <- alpha_grid

    for (fold in 1:n_folds) {
      test_indices <- which(fold_indices == fold)
      train_indices <- which(fold_indices != fold)
      X_train <- X_centered[train_indices, ]
      X_test <- X_centered[test_indices, ]
      for (i in seq_along(alpha_grid)) {
        alpha <- alpha_grid[i]
        sparse_pca_result <- sparsepca::spca(X_train, k = m, alpha = alpha)
        loadings <- sparse_pca_result$loadings
        scores_test <- X_test %*% loadings
        X_reconstructed <- scores_test %*% t(loadings)
        cv_errors[fold, i] <- mean((X_test - X_reconstructed)^2)
      }
    }
    mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
    if (all(is.na(mean_cv_errors))) {
      warning("All sparsity levels failed. Using default alpha = 0.5")
      return(0.5)
    }
    optimal_alpha <- alpha_grid[which.min(mean_cv_errors)]
    return(optimal_alpha)
  }, error = function(e) {
    warning("CV for Sparse PCA failed: ", e$message)
    return(0.5) # Default alpha if CV fails
  })
}

# --- Fit Sparse PCA with CV-tuned penalty ---
apply_sparse_pca <- function(X, m, sparsity_level) {
  tryCatch({
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    means <- colMeans(X)
    sparse_pca_result <- sparsepca::spca(X_centered, k = m, alpha = sparsity_level)
    loadings <- sparse_pca_result$loadings
    scores <- X_centered %*% loadings
    X_reconstructed <- scores %*% t(loadings) + matrix(rep(means, each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("Sparse PCA failed: ", e$message)
    list(loadings = matrix(NA, nrow = ncol(X), ncol = m),
         scores = matrix(NA, nrow = nrow(X), ncol = m),
         X_reconstructed = matrix(NA, nrow = nrow(X), ncol = ncol(X)))
  })
}

# --- Evaluation Metrics ---
evaluate_results <- function(true_X, reconstructed_X, true_L, estimated_L, true_F, estimated_scores, true_L_is_zero) {
  tryCatch({
    if(any(is.na(estimated_L)) || any(is.na(estimated_scores)) || any(is.na(reconstructed_X))) {
      return(list(
        mse = NA,
        loading_error = NA,
        false_positive_rate = NA,
        false_negative_rate = NA,
        mean_factor_score_cor = NA
      ))
    }
    mse <- mean((true_X - reconstructed_X)^2)
    loading_error <- tryCatch({
      procrustes_result <- vegan::procrustes(true_L, estimated_L)
      procrustes_result$ss
    }, error = function(e) {
      warning("Procrustes failed: ", e$message)
      NA
    })
    estimated_L_is_zero <- abs(estimated_L) < 0.05
    fp <- if(sum(true_L_is_zero == TRUE) == 0) NA else sum(estimated_L_is_zero == FALSE & true_L_is_zero == TRUE) / sum(true_L_is_zero == TRUE)
    fn <- if(sum(true_L_is_zero == FALSE) == 0) NA else sum(estimated_L_is_zero == TRUE & true_L_is_zero == FALSE) / sum(true_L_is_zero == FALSE)
    mean_factor_score_cor <- tryCatch({
      factor_score_cors <- diag(cor(true_F, estimated_scores))
      mean(factor_score_cors)
    }, error = function(e) {
      warning("Factor score correlation failed: ", e$message)
      NA
    })
    return(list(
      mse = mse,
      loading_error = loading_error,
      false_positive_rate = fp,
      false_negative_rate = fn,
      mean_factor_score_cor = mean_factor_score_cor
    ))
  }, error = function(e) {
    warning("Evaluation failed: ", e$message)
    list(mse = NA, loading_error = NA, false_positive_rate = NA, false_negative_rate = NA, mean_factor_score_cor = NA)
  })
}

get_true_zero_loadings <- function(L) {
  return(abs(L) < 1e-10)
}

# --- Simulation Function ---
run_simulation <- function(params) {
  tryCatch({
    set.seed(params$seed)
    n <- params$n
    p_vars <- params$p
    m <- params$m
    alpha <- params$alpha
    gamma <- params$gamma
    sigma_sq <- params$sigma_sq
    rep <- params$rep

    data <- generate_data(n, p_vars, m, alpha, gamma, sigma_sq)
    X <- data$X
    true_L <- data$L
    true_F <- data$F
    true_L_is_zero <- get_true_zero_loadings(true_L)
    optimal_sparsity <- find_optimal_sparsity_cv(X, m)
    
    pca_result <- apply_pca(X, m)
    pca_lp_result <- apply_pca_lp(X, m, lp_p = 1)
    efa_result <- apply_efa(X, m)
    efa_lp_result <- apply_efa_lp(X, m, lp_p = 1)
    sparse_pca_result <- apply_sparse_pca(X, m, optimal_sparsity)

    results <- list(
      PCA = evaluate_results(X, pca_result$X_reconstructed, true_L, pca_result$loadings, true_F, pca_result$scores, true_L_is_zero),
      PCA_L1 = evaluate_results(X, pca_lp_result$X_reconstructed, true_L, pca_lp_result$loadings, true_F, pca_lp_result$scores, true_L_is_zero),
      EFA = evaluate_results(X, efa_result$X_reconstructed, true_L, efa_result$loadings, true_F, efa_result$scores, true_L_is_zero),
      EFA_L1 = evaluate_results(X, efa_lp_result$X_reconstructed, true_L, efa_lp_result$loadings, true_F, efa_lp_result$scores, true_L_is_zero),
      SPCA = evaluate_results(X, sparse_pca_result$X_reconstructed, true_L, sparse_pca_result$loadings, true_F, sparse_pca_result$scores, true_L_is_zero)
    )
    return(list(params = params, results = results))
  }, error = function(e) {
    warning("Simulation run failed: ", e$message)
    return(list(params = params, results = NULL)) # Or a list of NAs
  })
}

# --- Example Simulation Run ---
# Define parameter grid
params_grid <- expand.grid(
  n = c(100, 500),
  p = c(50, 100),
  m = c(3, 5),
  alpha = c(0.2, 0.5),
  gamma = c(0, 1),
  sigma_sq = c(0.1, 1),
  rep = 1:2 # For demo, increase to 100-500 in real runs
)
params_list <- split(params_grid, seq(nrow(params_grid)))
for (i in seq_along(params_list)) params_list[[i]]$seed <- 1000 + i

# Run in parallel
simulation_results <- foreach(param = params_list, .packages = c('psych', 'sparsepca', 'GPArotation', 'vegan')) %dopar% {
  run_simulation(param)
}

# Save results
saveRDS(simulation_results, file = file.path(output_dir, "simulation_results.rds"))
stopCluster(cl)
