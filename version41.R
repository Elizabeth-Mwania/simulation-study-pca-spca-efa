# Load libraries
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

# ---- SOURCE Lp ROTATION FUNCTION ----
source("cl1121_Lp_fun.R") # <-- Make sure the path is correct

# Set up output directory
output_dir <- "simulation_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Setup parallel processing
n_cores <- detectCores() - 1  
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Function to create loading matrix
create_loading_matrix <- function(p, m, alpha, gamma) {
  L <- matrix(0, nrow = p, ncol = m)
  alpha_p <- floor(alpha * p)
  L[1:alpha_p, 1] <- runif(alpha_p, 0.7, 0.9)
  if (2*alpha_p <= p) {
    L[(alpha_p+1):(2*alpha_p), 2] <- runif(alpha_p, 0.7, 0.9)
  }
  if (2*alpha_p < p) {
    dense_indices <- (2*alpha_p+1):p
    L[dense_indices, ] <- matrix(runif(length(dense_indices) * m, 0.1, 0.4), 
                                nrow = length(dense_indices), ncol = m)
  }
  if (gamma > 0) {
    for (j in 1:m) {
      decay_factor <- j^(-gamma)
      L[, j] <- L[, j] * sqrt(decay_factor)
    }
  }
  return(L)
}

# Data generation
generate_data <- function(n, p, m, alpha, gamma, sigma_sq) {
  L <- create_loading_matrix(p, m, alpha, gamma)
  F <- matrix(rnorm(n * m), nrow = n, ncol = m)
  E <- matrix(rnorm(n * p, 0, sqrt(sigma_sq)), nrow = n, ncol = p)
  X <- F %*% t(L) + E
  return(list(X = X, L = L, F = F, E = E))
}

# --- L1 (Lp, p=1) ROTATION WRAPPER ---
apply_lp_rotation <- function(loadings, lp_p = 1, max_iter = 500, tol = 1e-5) {
  # LpRotation from cl1121_Lp_fun.R expects loadings as input
  rotated <- LpRotation(loadings, p = lp_p, maxit = max_iter, tol = tol)
  # If the function returns a list, get the rotated loadings
  if (is.list(rotated) && !is.null(rotated$L)) {
    return(rotated$L)
  } else {
    return(rotated)
  }
}

# --- PCA with Varimax ---
apply_pca <- function(X, m) {
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  tryCatch({
    pca_result <- prcomp(X_centered, center = FALSE)
    m_actual <- min(m, ncol(pca_result$rotation))
    loadings <- pca_result$rotation[, 1:m_actual, drop = FALSE]
    scores <- pca_result$x[, 1:m_actual, drop = FALSE]
    if(m_actual > 1) {
      rot_result <- GPArotation::Varimax(loadings, normalize = FALSE)
      loadings_rotated <- rot_result$loadings
      scores_rotated <- scores %*% rot_result$rotmat
    } else {
      loadings_rotated <- loadings
      scores_rotated <- scores
    }
    if(m_actual < m) {
      loadings_rotated <- cbind(loadings_rotated, matrix(NA, nrow = nrow(loadings_rotated), ncol = m - m_actual))
      scores_rotated <- cbind(scores_rotated, matrix(NA, nrow = nrow(scores_rotated), ncol = m - m_actual))
    }
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("PCA failed. Returning NA matrices. Error: ", e$message)
    p <- ncol(X)
    n <- nrow(X)
    loadings <- matrix(NA, nrow = p, ncol = m)
    scores <- matrix(NA, nrow = n, ncol = m)
    X_reconstructed <- matrix(NA, nrow = n, ncol = p)
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  })
}

# --- PCA with L1 Rotation ---
apply_pca_lp <- function(X, m, lp_p = 1) {
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  tryCatch({
    pca_result <- prcomp(X_centered, center = FALSE)
    m_actual <- min(m, ncol(pca_result$rotation))
    loadings <- pca_result$rotation[, 1:m_actual, drop = FALSE]
    scores <- pca_result$x[, 1:m_actual, drop = FALSE]
    if(m_actual > 1) {
      loadings_rotated <- apply_lp_rotation(loadings, lp_p = lp_p)
      scores_rotated <- scores # No rotation matrix from Lp
    } else {
      loadings_rotated <- loadings
      scores_rotated <- scores
    }
    if(m_actual < m) {
      loadings_rotated <- cbind(loadings_rotated, matrix(NA, nrow = nrow(loadings_rotated), ncol = m - m_actual))
      scores_rotated <- cbind(scores_rotated, matrix(NA, nrow = nrow(scores_rotated), ncol = m - m_actual))
    }
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("PCA-L1 failed. Returning NA matrices. Error: ", e$message)
    p <- ncol(X)
    n <- nrow(X)
    loadings <- matrix(NA, nrow = p, ncol = m)
    scores <- matrix(NA, nrow = n, ncol = m)
    X_reconstructed <- matrix(NA, nrow = n, ncol = p)
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  })
}

# --- EFA with Varimax ---
apply_efa <- function(X, m) {
  tryCatch({
    efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "ml")
    loadings <- as.matrix(efa_result_unrotated$loadings)
    scores <- efa_result_unrotated$scores
    if(m > 1) {
      rot_result <- GPArotation::Varimax(loadings, normalize = FALSE)
      loadings_rotated <- rot_result$loadings
      scores_rotated <- scores %*% rot_result$rotmat
    } else {
      loadings_rotated <- loadings
      scores_rotated <- scores
    }
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    tryCatch({
      efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "minres")
      loadings <- as.matrix(efa_result_unrotated$loadings)
      scores <- efa_result_unrotated$scores
      if(m > 1) {
        rot_result <- GPArotation::Varimax(loadings, normalize = FALSE)
        loadings_rotated <- rot_result$loadings
        scores_rotated <- scores %*% rot_result$rotmat
      } else {
        loadings_rotated <- loadings
        scores_rotated <- scores
      }
      X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
      return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
    }, error = function(e) {
      p <- ncol(X)
      n <- nrow(X)
      loadings <- matrix(NA, nrow = p, ncol = m)
      scores <- matrix(NA, nrow = n, ncol = m)
      X_reconstructed <- matrix(NA, nrow = n, ncol = p)
      warning("EFA failed to converge. Returning NA matrices.")
      return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
    })
  })
}

# --- EFA with L1 Rotation ---
apply_efa_lp <- function(X, m, lp_p = 1) {
  tryCatch({
    efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "ml")
    loadings <- as.matrix(efa_result_unrotated$loadings)
    scores <- efa_result_unrotated$scores
    if(m > 1) {
      loadings_rotated <- apply_lp_rotation(loadings, lp_p = lp_p)
      scores_rotated <- scores
    } else {
      loadings_rotated <- loadings
      scores_rotated <- scores
    }
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    tryCatch({
      efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "minres")
      loadings <- as.matrix(efa_result_unrotated$loadings)
      scores <- efa_result_unrotated$scores
      if(m > 1) {
        loadings_rotated <- apply_lp_rotation(loadings, lp_p = lp_p)
        scores_rotated <- scores
      } else {
        loadings_rotated <- loadings
        scores_rotated <- scores
      }
      X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
      return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
    }, error = function(e) {
      p <- ncol(X)
      n <- nrow(X)
      loadings <- matrix(NA, nrow = p, ncol = m)
      scores <- matrix(NA, nrow = n, ncol = m)
      X_reconstructed <- matrix(NA, nrow = n, ncol = p)
      warning("EFA-L1 failed to converge. Returning NA matrices.")
      return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
    })
  })
}


# Cross-validation for sparse PCA tuning
find_optimal_sparsity_cv <- function(X, m, n_folds = 5, alpha_grid = seq(0.1, 0.9, by = 0.1)) {
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  n <- nrow(X)
  
  fold_indices <- cut(seq(1, n), breaks = n_folds, labels = FALSE)
  fold_indices <- sample(fold_indices)
  
  cv_errors <- matrix(NA, nrow = n_folds, ncol = length(alpha_grid))
  colnames(cv_errors) <- alpha_grid
  
  for (fold in 1:n_folds) {
    test_indices <- which(fold_indices == fold)
    train_indices <- which(fold_indices != fold)
    
    X_train <- X_centered[train_indices, ]
    X_test <- X_centered[test_indices, ]
    
    for (i in seq_along(alpha_grid)) {
      alpha <- alpha_grid[i]
      tryCatch({
        sparse_pca_result <- sparsepca::spca(X_train, k = m, alpha = alpha)
        loadings <- sparse_pca_result$loadings
        scores_test <- X_test %*% loadings
        X_reconstructed <- scores_test %*% t(loadings)
        cv_errors[fold, i] <- mean((X_test - X_reconstructed)^2)
      }, error = function(e) {
        cv_errors[fold, i] <- NA
      })
    }
  }
  
  mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
  if (all(is.na(mean_cv_errors))) {
    warning("All sparsity levels failed. Using default alpha = 0.5")
    return(0.5)
  }
  optimal_alpha <- alpha_grid[which.min(mean_cv_errors)]
  return(optimal_alpha)
}

# Fit Sparse PCA with CV-tuned penalty
apply_sparse_pca <- function(X, m, sparsity_level) {
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  means <- colMeans(X)
  tryCatch({
    sparse_pca_result <- sparsepca::spca(X_centered, k = m, alpha = sparsity_level)
    loadings <- sparse_pca_result$loadings
    scores <- X_centered %*% loadings
    X_reconstructed <- scores %*% t(loadings) + matrix(rep(means, each = nrow(X)), nrow = nrow(X))
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("Sparse PCA failed. Returning NA matrices. Error: ", e$message)
    p_vars <- ncol(X)
    n <- nrow(X)
    loadings <- matrix(NA, nrow = p_vars, ncol = m)
    scores <- matrix(NA, nrow = n, ncol = m)
    X_reconstructed <- matrix(NA, nrow = n, ncol = p_vars)
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  })
}

# Evaluation metrics
evaluate_results <- function(true_X, reconstructed_X, true_L, estimated_L, true_F, estimated_scores, true_L_is_zero) {
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
  
  tryCatch({
    procrustes_result <- vegan::procrustes(true_L, estimated_L)
    loading_error <- procrustes_result$ss
  }, error = function(e) {
    warning("Procrustes analysis failed. Using NA for loading error.")
    loading_error <- NA
  })
  
  estimated_L_is_zero <- abs(estimated_L) < 0.05
  
  if(sum(true_L_is_zero == TRUE) == 0) {
    fp <- NA
  } else {
    fp <- sum(estimated_L_is_zero == FALSE & true_L_is_zero == TRUE) / sum(true_L_is_zero == TRUE)
  }
  
  if(sum(true_L_is_zero == FALSE) == 0) {
    fn <- NA
  } else {
    fn <- sum(estimated_L_is_zero == TRUE & true_L_is_zero == FALSE) / sum(true_L_is_zero == FALSE)
  }
  
  tryCatch({
    factor_score_cors <- diag(cor(true_F, estimated_scores))
    mean_factor_score_cor <- mean(factor_score_cors)
  }, error = function(e) {
    warning("Factor score correlation calculation failed. Using NA.")
    mean_factor_score_cor <- NA
  })
  
  return(list(
    mse = mse,
    loading_error = loading_error,
    false_positive_rate = fp,
    false_negative_rate = fn,
    mean_factor_score_cor = mean_factor_score_cor
  ))
}

# Function to determine true zero loadings
get_true_zero_loadings <- function(L) {
  return(abs(L) < 1e-10)
}

# --- Simulation function ---
run_simulation <- function(params) {
  tryCatch({
    set.seed(params$seed)
    n <- params$n
    p <- params$p
    m <- params$m
    alpha <- params$alpha
    gamma <- params$gamma
    sigma_sq <- params$sigma_sq
    rep <- params$rep

    data <- generate_data(n, p, m, alpha, gamma, sigma_sq)
    X <- data$X
    true_L <- data$L
    true_F <- data$F
    true_L_is_zero <- get_true_zero_loadings(true_L)
    optimal_sparsity <- find_optimal_sparsity_cv(X, m)

    # Fit each method
    pca_result <- apply_pca(X, m)
    sparse_pca_result <- apply_sparse_pca(X, m, optimal_sparsity)
    efa_result <- apply_efa(X, m)
    pca_lp_result <- apply_pca_lp(X, m, lp_p = 1)
    efa_lp_result <- apply_efa_lp(X, m, lp_p = 1)

    # Evaluate
    pca_eval <- evaluate_results(X, pca_result$X_reconstructed, true_L, pca_result$loadings, true_F, pca_result$scores, true_L_is_zero)
    sparse_pca_eval <- evaluate_results(X, sparse_pca_result$X_reconstructed, true_L, sparse_pca_result$loadings, true_F, sparse_pca_result$scores, true_L_is_zero)
    efa_eval <- evaluate_results(X, efa_result$X_reconstructed, true_L, efa_result$loadings, true_F, efa_result$scores, true_L_is_zero)
    pca_lp_eval <- evaluate_results(X, pca_lp_result$X_reconstructed, true_L, pca_lp_result$loadings, true_F, pca_lp_result$scores, true_L_is_zero)
    efa_lp_eval <- evaluate_results(X, efa_lp_result$X_reconstructed, true_L, efa_lp_result$loadings, true_F, efa_lp_result$scores, true_L_is_zero)

    results <- data.frame(
      method = c("PCA", "Sparse PCA", "EFA", "PCA_L1_Rotation", "EFA_L1_Rotation"),
      mse = c(pca_eval$mse, sparse_pca_eval$mse, efa_eval$mse, pca_lp_eval$mse, efa_lp_eval$mse),
      loading_error = c(pca_eval$loading_error, sparse_pca_eval$loading_error, efa_eval$loading_error, pca_lp_eval$loading_error, efa_lp_eval$loading_error),
      false_positive_rate = c(pca_eval$false_positive_rate, sparse_pca_eval$false_positive_rate, efa_eval$false_positive_rate, pca_lp_eval$false_positive_rate, efa_lp_eval$false_positive_rate),
      false_negative_rate = c(pca_eval$false_negative_rate, sparse_pca_eval$false_negative_rate, efa_eval$false_negative_rate, pca_lp_eval$false_negative_rate, efa_lp_eval$false_negative_rate),
      mean_factor_score_cor = c(pca_eval$mean_factor_score_cor, sparse_pca_eval$mean_factor_score_cor, efa_eval$mean_factor_score_cor, pca_lp_eval$mean_factor_score_cor, efa_lp_eval$mean_factor_score_cor),
      optimal_alpha = c(NA, optimal_sparsity, NA, NA, NA),
      n = n,
      p = p,
      m = m,
      alpha = alpha,
      gamma = gamma,
      sigma_sq = sigma_sq,
      rep = rep
    )

    # Save loadings
    loadings_file <- sprintf("%s/loadings_n%d_p%d_m%d_alpha%.1f_gamma%d_sigma%.1f_rep%d.rds", 
                            output_dir, n, p, m, alpha, gamma, sigma_sq, rep)
    loadings_data <- list(
      true_L = true_L,
      pca_L = pca_result$loadings,
      sparse_pca_L = sparse_pca_result$loadings,
      efa_L = efa_result$loadings,
      pca_lp_L = pca_lp_result$loadings,
      efa_lp_L = efa_lp_result$loadings,
      optimal_alpha = optimal_sparsity
    )
    saveRDS(loadings_data, loadings_file)
    return(results)
  }, error = function(e) {
    warning("Error in simulation run: ", e$message)
    data.frame(
      method = c("PCA", "Sparse PCA", "EFA", "PCA_L1_Rotation", "EFA_L1_Rotation"),
      mse = rep(NA, 5),
      loading_error = rep(NA, 5),
      false_positive_rate = rep(NA, 5),
      false_negative_rate = rep(NA, 5),
      mean_factor_score_cor = rep(NA, 5),
      optimal_alpha = rep(NA, 5),
      n = params$n,
      p = params$p,
      m = params$m,
      alpha = params$alpha,
      gamma = params$gamma,
      sigma_sq = params$sigma_sq,
      rep = params$rep,
      error_message = rep(e$message, 5)
    )
  })
}

# --- Parameter grid and simulation loop remain unchanged ---
params <- expand.grid(
  n = c(100, 500),
  p = c(50, 100, 200),
  m = c(3, 5),
  alpha = c(0.2, 0.5, 0.8),
  gamma = c(0, 1, 2),
  sigma_sq = c(0.1, 1, 5),
  rep = 1:100
)
params$seed <- 1000 + (1:nrow(params))

start_time <- Sys.time()
cat("Starting simulation at", as.character(start_time), "\n")
results <- foreach(i = 1:nrow(params), .combine = rbind, .packages = c("psych", "sparsepca", "GPArotation", "dplyr", "vegan")) %dopar% {
  cat("Running scenario", i, "of", nrow(params), "\n")
  run_simulation(params[i, ])
}
end_time <- Sys.time()
cat("Simulation completed at", as.character(end_time), "\n")
cat("Total time:", difftime(end_time, start_time, units = "hours"), "hours\n")
stopCluster(cl)
write_csv(results, file.path(output_dir, "simulated_results_cv_tuned_with_L1.csv"))
