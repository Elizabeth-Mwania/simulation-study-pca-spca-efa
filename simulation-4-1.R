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
  
  # Loading with sparsity structure
  alpha_p <- floor(alpha * p)
  
  # Variables to load only on factor 1 
  L[1:alpha_p, 1] <- runif(alpha_p, 0.7, 0.9)
  
  # Variables to load only on factor 2
  if (2*alpha_p <= p) {
    L[(alpha_p+1):(2*alpha_p), 2] <- runif(alpha_p, 0.7, 0.9)
  }
  
  # Dense set up
  if (2*alpha_p < p) {
    dense_indices <- (2*alpha_p+1):p
    L[dense_indices, ] <- matrix(runif(length(dense_indices) * m, 0.1, 0.4), 
                                nrow = length(dense_indices), ncol = m)
  }
  
  # Eigenvalue decay
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
  # Loading matrix
  L <- create_loading_matrix(p, m, alpha, gamma)
  
  # Factors
  F <- matrix(rnorm(n * m), nrow = n, ncol = m)
  
  # Error or noise
  E <- matrix(rnorm(n * p, 0, sqrt(sigma_sq)), nrow = n, ncol = p)
  
  # Data matrix
  X <- F %*% t(L) + E
  
  return(list(X = X, L = L, F = F, E = E))
}

# Fit PCA 
# Fit PCA with consistent rotation
apply_pca <- function(X, m) {
  # Center data 
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  tryCatch({
    pca_result <- prcomp(X_centered, center = FALSE)
    
    # Check for components
    if(ncol(pca_result$rotation) < m) {
      warning("PCA could not extract enough components. Using available components.")
      m_actual <- ncol(pca_result$rotation)
    } else {
      m_actual <- m
    }
    
    loadings <- pca_result$rotation[, 1:m_actual, drop = FALSE]
    scores <- pca_result$x[, 1:m_actual, drop = FALSE]
    
    # Apply consistent varimax rotation using GPArotation
    if(m_actual > 1) {
      rot_result <- GPArotation::Varimax(loadings, normalize = FALSE)
      loadings_rotated <- rot_result$loadings
      scores_rotated <- scores %*% rot_result$rotmat
    } else {
      loadings_rotated <- loadings
      scores_rotated <- scores
    }
    
    # Maintain dimension
    if(m_actual < m) {
      loadings_rotated <- cbind(loadings_rotated, matrix(NA, nrow = nrow(loadings_rotated), ncol = m - m_actual))
      scores_rotated <- cbind(scores_rotated, matrix(NA, nrow = nrow(scores_rotated), ncol = m - m_actual))
    }
    
    # Reconstruct data
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

# Fit EFA with consistent rotation
apply_efa <- function(X, m) {
  tryCatch({
    # First get unrotated solution
    efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "ml")
    loadings <- as.matrix(efa_result_unrotated$loadings)
    scores <- efa_result_unrotated$scores
    
    # Apply same varimax rotation as PCA
    if(m > 1) {
      rot_result <- GPArotation::Varimax(loadings, normalize = FALSE)
      loadings_rotated <- rot_result$loadings
      scores_rotated <- scores %*% rot_result$rotmat
    } else {
      loadings_rotated <- loadings
      scores_rotated <- scores
    }
    
    # Reconstruct data
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    
    return(list(loadings = loadings_rotated, scores = scores_rotated, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    # Fallback to minimum residual if ML fails
    tryCatch({
      efa_result_unrotated <- fa(X, nfactors = m, rotate = "none", fm = "minres")
      loadings <- as.matrix(efa_result_unrotated$loadings)
      scores <- efa_result_unrotated$scores
      
      # Apply same rotation
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
      # If both methods fail, return NA values
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
# Cross-validation for sparse PCA tuning
find_optimal_sparsity_cv <- function(X, m, n_folds = 5, alpha_grid = seq(0.1, 0.9, by = 0.1)) {
  # Center the data
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  
  # Create folds
  fold_indices <- cut(seq(1, n), breaks = n_folds, labels = FALSE)
  fold_indices <- sample(fold_indices)  # Randomize
  
  cv_errors <- matrix(NA, nrow = n_folds, ncol = length(alpha_grid))
  colnames(cv_errors) <- alpha_grid
  
  for (fold in 1:n_folds) {
    # Split data
    test_indices <- which(fold_indices == fold)
    train_indices <- which(fold_indices != fold)
    
    X_train <- X_centered[train_indices, ]
    X_test <- X_centered[test_indices, ]
    
    for (i in seq_along(alpha_grid)) {
      alpha <- alpha_grid[i]
      
      tryCatch({
        # Fit sparse PCA on training data
        sparse_pca_result <- sparsepca::spca(X_train, k = m, alpha = alpha)
        
        # Calculate reconstruction error on test data
        loadings <- sparse_pca_result$loadings
        scores_test <- X_test %*% loadings
        X_reconstructed <- scores_test %*% t(loadings)
        
        # Store MSE for this fold and alpha
        cv_errors[fold, i] <- mean((X_test - X_reconstructed)^2)
      }, error = function(e) {
        cv_errors[fold, i] <- NA
      })
    }
  }
  
  # Calculate mean CV error for each alpha
  mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
  
  # Return alpha with minimum CV error
  if (all(is.na(mean_cv_errors))) {
    warning("All sparsity levels failed. Using default alpha = 0.5")
    return(0.5)
  }
  
  optimal_alpha <- alpha_grid[which.min(mean_cv_errors)]
  return(optimal_alpha)
}

# Fit Sparse PCA with CV-tuned penalty
apply_sparse_pca <- function(X, m, sparsity_level) {
  # Center data
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  means <- colMeans(X)
  tryCatch({
    sparse_pca_result <- sparsepca::spca(X_centered, k = m, alpha = sparsity_level)
    
    loadings <- sparse_pca_result$loadings
    scores <- X_centered %*% loadings
    
    # Reconstruct data
    X_reconstructed <- scores %*% t(loadings) + matrix(rep(means, each = nrow(X)), nrow = nrow(X))
    
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    warning("Sparse PCA failed. Returning NA matrices. Error: ", e$message)
    
    p <- ncol(X)
    n <- nrow(X)

    loadings <- matrix(NA, nrow = p, ncol = m)
    scores <- matrix(NA, nrow = n, ncol = m)
    X_reconstructed <- matrix(NA, nrow = n, ncol = p)
    
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  })
}

# Evaluation metrics
evaluate_results <- function(true_X, reconstructed_X, true_L, estimated_L, true_F, estimated_scores, true_L_is_zero) {
  # Check for NA values
  if(any(is.na(estimated_L)) || any(is.na(estimated_scores)) || any(is.na(reconstructed_X))) {
    return(list(
      mse = NA,
      loading_error = NA,
      false_positive_rate = NA,
      false_negative_rate = NA,
      mean_factor_score_cor = NA
    ))
  }
  
  # Reconstruction error (MSE)
  mse <- mean((true_X - reconstructed_X)^2)
  
  # Align estimated loadings using procrustes analysis
  tryCatch({
    procrustes_result <- vegan::procrustes(true_L, estimated_L)
    loading_error <- procrustes_result$ss  
  }, error = function(e) {
    warning("Procrustes analysis failed. Using NA for loading error.")
    loading_error <- NA
  })
  
  # Sparsity measures
  estimated_L_is_zero <- abs(estimated_L) < 0.05

  if(sum(true_L_is_zero == TRUE) == 0) {
    fp <- NA
  } else {
    # False positive rate
    fp <- sum(estimated_L_is_zero == FALSE & true_L_is_zero == TRUE) / sum(true_L_is_zero == TRUE)
  }
  
  if(sum(true_L_is_zero == FALSE) == 0) {
    fn <- NA
  } else {
    # False negative rate
    fn <- sum(estimated_L_is_zero == TRUE & true_L_is_zero == FALSE) / sum(true_L_is_zero == FALSE)
  }
  
  # Factor score correlations
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

# Simulation function
run_simulation <- function(params) {
  tryCatch({
    set.seed(params$seed)
    
    # Parameters
    n <- params$n
    p <- params$p
    m <- params$m
    alpha <- params$alpha
    gamma <- params$gamma
    sigma_sq <- params$sigma_sq
    rep <- params$rep
    
    # Generate data
    data <- generate_data(n, p, m, alpha, gamma, sigma_sq)
    X <- data$X
    true_L <- data$L
    true_F <- data$F
    
    # Identify true zero loadings
    true_L_is_zero <- get_true_zero_loadings(true_L)
    
    # Find optimal sparsity using CV
    optimal_sparsity <- find_optimal_sparsity_cv(X, m)
    
    # Fit each method
    pca_result <- apply_pca(X, m)
    sparse_pca_result <- apply_sparse_pca(X, m, optimal_sparsity)
    efa_result <- apply_efa(X, m)
    
    # Evaluate results
    pca_eval <- evaluate_results(X, pca_result$X_reconstructed, true_L, pca_result$loadings, 
                                true_F, pca_result$scores, true_L_is_zero)
    
    sparse_pca_eval <- evaluate_results(X, sparse_pca_result$X_reconstructed, true_L, 
                                      sparse_pca_result$loadings, true_F, 
                                      sparse_pca_result$scores, true_L_is_zero)
    
    efa_eval <- evaluate_results(X, efa_result$X_reconstructed, true_L, efa_result$loadings, 
                                true_F, efa_result$scores, true_L_is_zero)
    
    # Create result dataframe 
    results <- data.frame(
      method = c("PCA", "Sparse PCA", "EFA"),
      mse = c(pca_eval$mse, sparse_pca_eval$mse, efa_eval$mse),
      loading_error = c(pca_eval$loading_error, sparse_pca_eval$loading_error, efa_eval$loading_error),
      false_positive_rate = c(pca_eval$false_positive_rate, sparse_pca_eval$false_positive_rate, efa_eval$false_positive_rate),
      false_negative_rate = c(pca_eval$false_negative_rate, sparse_pca_eval$false_negative_rate, efa_eval$false_negative_rate),
      mean_factor_score_cor = c(pca_eval$mean_factor_score_cor, sparse_pca_eval$mean_factor_score_cor, efa_eval$mean_factor_score_cor),
      optimal_alpha = c(NA, optimal_sparsity, NA),  # Record optimal alpha for Sparse PCA
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
      optimal_alpha = optimal_sparsity
    )
    
    saveRDS(loadings_data, loadings_file)
    
    return(results)
  }, error = function(e) {
    warning("Error in simulation run: ", e$message)
    
    # Return a data frame with NAs
    data.frame(
      method = c("PCA", "Sparse PCA", "EFA"),
      mse = rep(NA, 3),
      loading_error = rep(NA, 3),
      false_positive_rate = rep(NA, 3),
      false_negative_rate = rep(NA, 3),
      mean_factor_score_cor = rep(NA, 3),
      optimal_alpha = rep(NA, 3),
      n = params$n,
      p = params$p,
      m = params$m,
      alpha = params$alpha,
      gamma = params$gamma,
      sigma_sq = params$sigma_sq,
      rep = params$rep,
      error_message = rep(e$message, 3)
    )
  })
}

# Design parameter grid
params <- expand.grid(
  n = c(100, 500),
  p = c(50, 100, 200),
  m = c(3, 5),
  alpha = c(0.2, 0.5, 0.8),
  gamma = c(0, 1, 2),
  sigma_sq = c(0.1, 1, 5),
  rep = 1:1  
)

# Add seed for reproducibility
params$seed <- 1000 + (1:nrow(params))

# Parallel simulation
start_time <- Sys.time()
cat("Starting simulation at", as.character(start_time), "\n")

results <- foreach(i = 1:nrow(params), .combine = rbind, .packages = c("psych", "sparsepca", "GPArotation", "dplyr", "vegan")
) %dopar% {
  cat("Running scenario", i, "of", nrow(params), "\n")
  run_simulation(params[i, ])
}

end_time <- Sys.time()
cat("Simulation completed at", as.character(end_time), "\n")
cat("Total time:", difftime(end_time, start_time, units = "hours"), "hours\n")

# Stop cluster
stopCluster(cl)

# Save simulation results
write_csv(results, file.path(output_dir, "simulated_results_cv_tuned.csv"))