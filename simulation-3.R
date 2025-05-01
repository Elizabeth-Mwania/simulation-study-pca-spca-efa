# Monte Carlo Simulation: Comparing PCA, Sparse PCA, and EFA
# Based on provided guidelines

# Required packages
if (!require("parallel")) install.packages("parallel")
if (!require("doParallel")) install.packages("doParallel")
if (!require("foreach")) install.packages("foreach")
if (!require("psych")) install.packages("psych")
if (!require("sparsepca")) install.packages("sparsepca")
if (!require("GPArotation")) install.packages("GPArotation")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")
if (!require("corrplot")) install.packages("corrplot")
if (!require("vegan")) install.packages("vegan")  # Package containing procrustes function

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
library(vegan)  # Load vegan package for procrustes function

# Set up output directory
output_dir <- "simulation_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Setup parallel processing
n_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Function to create loading matrix with mixed structure
create_loading_matrix <- function(p, m, alpha, gamma) {
  L <- matrix(0, nrow = p, ncol = m)
  
  # Set the sparsity structure
  alpha_p <- floor(alpha * p)
  
  # Variables 1 to αp: load only on factor 1 (sparse)
  L[1:alpha_p, 1] <- runif(alpha_p, 0.7, 0.9)
  
  # Variables αp+1 to 2αp: load only on factor 2
  if (2*alpha_p <= p) {
    L[(alpha_p+1):(2*alpha_p), 2] <- runif(alpha_p, 0.7, 0.9)
  }
  
  # Remaining variables: small cross-loadings on all factors (dense)
  if (2*alpha_p < p) {
    dense_indices <- (2*alpha_p+1):p
    L[dense_indices, ] <- matrix(runif(length(dense_indices) * m, 0.1, 0.4), 
                                nrow = length(dense_indices), ncol = m)
  }
  
  # Apply eigenvalue decay
  if (gamma > 0) {
    for (j in 1:m) {
      decay_factor <- j^(-gamma)
      L[, j] <- L[, j] * sqrt(decay_factor)
    }
  }
  
  return(L)
}

# Generate data based on factor model
generate_data <- function(n, p, m, alpha, gamma, sigma_sq) {
  # Generate loading matrix
  L <- create_loading_matrix(p, m, alpha, gamma)
  
  # Generate factors
  F <- matrix(rnorm(n * m), nrow = n, ncol = m)
  
  # Generate error
  E <- matrix(rnorm(n * p, 0, sqrt(sigma_sq)), nrow = n, ncol = p)
  
  # Generate data
  X <- F %*% t(L) + E
  
  return(list(X = X, L = L, F = F, E = E))
}

# Apply PCA and reconstruct data
apply_pca <- function(X, m) {
  # Center data 
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  
  # Apply PCA with error handling
  tryCatch({
    pca_result <- prcomp(X_centered, center = FALSE)
    
    # Check if we have enough components
    if(ncol(pca_result$rotation) < m) {
      warning("PCA could not extract enough components. Using available components.")
      m_actual <- ncol(pca_result$rotation)
    } else {
      m_actual <- m
    }
    
    loadings <- pca_result$rotation[, 1:m_actual, drop = FALSE]
    scores <- pca_result$x[, 1:m_actual, drop = FALSE]
    
    # Apply varimax rotation if we have more than one component
    if(m_actual > 1) {
      rot_result <- varimax(loadings)
      loadings_rotated <- loadings %*% rot_result$rotmat
      scores_rotated <- scores %*% solve(rot_result$rotmat)
    } else {
      loadings_rotated <- loadings
      scores_rotated <- scores
    }
    
    # Pad with NA columns if needed to maintain dimension
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
    
    # Create NA matrices with proper dimensions
    loadings <- matrix(NA, nrow = p, ncol = m)
    scores <- matrix(NA, nrow = n, ncol = m)
    X_reconstructed <- matrix(NA, nrow = n, ncol = p)
    
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  })
}

# Apply Sparse PCA
apply_sparse_pca <- function(X, m, sparsity_level) {
  # Center data
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  means <- colMeans(X)
  
  # Apply sparse PCA with error handling
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
    
    # Create NA matrices with proper dimensions
    loadings <- matrix(NA, nrow = p, ncol = m)
    scores <- matrix(NA, nrow = n, ncol = m)
    X_reconstructed <- matrix(NA, nrow = n, ncol = p)
    
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  })
}

# Apply EFA
apply_efa <- function(X, m) {
  # Apply EFA using maximum likelihood with error handling
  tryCatch({
    efa_result <- fa(X, nfactors = m, rotate = "varimax", fm = "ml")
    
    loadings <- efa_result$loadings
    scores <- efa_result$scores
    
    # Reconstruct data
    X_reconstructed <- scores %*% t(loadings) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
    
    return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
  }, error = function(e) {
    # Fallback to minimum residual if ML fails
    tryCatch({
      efa_result <- fa(X, nfactors = m, rotate = "varimax", fm = "minres")
      
      loadings <- efa_result$loadings
      scores <- efa_result$scores
      
      # Reconstruct data
      X_reconstructed <- scores %*% t(loadings) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
      
      return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
    }, error = function(e) {
      # If both methods fail, return NA values with proper dimensions
      p <- ncol(X)
      n <- nrow(X)
      
      # Create NA matrices with proper dimensions
      loadings <- matrix(NA, nrow = p, ncol = m)
      scores <- matrix(NA, nrow = n, ncol = m)
      X_reconstructed <- matrix(NA, nrow = n, ncol = p)
      
      warning("EFA failed to converge. Returning NA matrices.")
      return(list(loadings = loadings, scores = scores, X_reconstructed = X_reconstructed))
    })
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
  
  # Loading accuracy using Procrustes analysis
  # Align estimated loadings with true loadings using vegan's procrustes function
  tryCatch({
    procrustes_result <- vegan::procrustes(true_L, estimated_L)
    loading_error <- procrustes_result$ss  # Sum of squared differences after Procrustes rotation
  }, error = function(e) {
    warning("Procrustes analysis failed. Using NA for loading error.")
    loading_error <- NA
  })
  
  # Sparsity measures
  estimated_L_is_zero <- abs(estimated_L) < 0.05
  
  # Check if the denominators are zero
  if(sum(true_L_is_zero == TRUE) == 0) {
    fp <- NA
  } else {
    # False positive rate: predicted non-zero when truly zero
    fp <- sum(estimated_L_is_zero == FALSE & true_L_is_zero == TRUE) / sum(true_L_is_zero == TRUE)
  }
  
  if(sum(true_L_is_zero == FALSE) == 0) {
    fn <- NA
  } else {
    # False negative rate: predicted zero when truly non-zero
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

# Find optimal sparsity level for sparse PCA using BIC
find_optimal_sparsity <- function(X, m, true_sparsity_level) {
  # Use a grid of sparsity levels around true value
  alpha_grid <- seq(max(0.1, true_sparsity_level - 0.2), 
                    min(0.9, true_sparsity_level + 0.2), 
                    by = 0.1)
  
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  
  bic_values <- numeric(length(alpha_grid))
  
  for (i in seq_along(alpha_grid)) {
    tryCatch({
      alpha <- alpha_grid[i]
      sparse_pca_result <- sparsepca::spca(X_centered, k = m, alpha = alpha)
      loadings <- sparse_pca_result$loadings
      scores <- X_centered %*% loadings
      X_reconstructed <- scores %*% t(loadings)
      
      # Calculate BIC
      rss <- sum((X_centered - X_reconstructed)^2)
      num_nonzero <- sum(abs(loadings) > 0.05)
      bic <- n * log(rss/n) + log(n) * num_nonzero
      
      bic_values[i] <- bic
    }, error = function(e) {
      bic_values[i] <- Inf
    })
  }
  
  # If all BIC values are Inf, return default value
  if(all(is.infinite(bic_values))) {
    warning("All sparsity levels failed. Using default sparsity level.")
    return(0.5)  # Default sparsity level
  }
  
  # Return alpha with minimum BIC
  return(alpha_grid[which.min(bic_values)])
}

# Main simulation function
run_simulation <- function(params) {
  tryCatch({
    set.seed(params$seed)
    
    # Unpack parameters
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
    
    # Calculate true sparsity for tuning sparse PCA
    true_sparsity <- sum(true_L_is_zero) / (p * m)
    
    # Find optimal sparsity for sparse PCA
    optimal_sparsity <- find_optimal_sparsity(X, m, true_sparsity)
    
    # Apply methods
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
    
    # Create result dataframe with error handling
    results <- data.frame(
      method = c("PCA", "Sparse PCA", "EFA"),
      mse = c(pca_eval$mse, sparse_pca_eval$mse, efa_eval$mse),
      loading_error = c(pca_eval$loading_error, sparse_pca_eval$loading_error, efa_eval$loading_error),
      false_positive_rate = c(pca_eval$false_positive_rate, sparse_pca_eval$false_positive_rate, efa_eval$false_positive_rate),
      false_negative_rate = c(pca_eval$false_negative_rate, sparse_pca_eval$false_negative_rate, efa_eval$false_negative_rate),
      mean_factor_score_cor = c(pca_eval$mean_factor_score_cor, sparse_pca_eval$mean_factor_score_cor, efa_eval$mean_factor_score_cor),
      n = n,
      p = p,
      m = m,
      alpha = alpha,
      gamma = gamma,
      sigma_sq = sigma_sq,
      rep = rep
    )
    
    # Save loadings for deeper diagnostics
    loadings_file <- sprintf("%s/loadings_n%d_p%d_m%d_alpha%.1f_gamma%d_sigma%.1f_rep%d.rds", 
                            output_dir, n, p, m, alpha, gamma, sigma_sq, rep)
    
    loadings_data <- list(
      true_L = true_L,
      pca_L = pca_result$loadings,
      sparse_pca_L = sparse_pca_result$loadings,
      efa_L = efa_result$loadings
    )
    
    saveRDS(loadings_data, loadings_file)
    
    return(results)
  }, error = function(e) {
    # Create a row with error information
    warning("Error in simulation run: ", e$message)
    
    # Return a data frame with NAs
    data.frame(
      method = c("PCA", "Sparse PCA", "EFA"),
      mse = rep(NA, 3),
      loading_error = rep(NA, 3),
      false_positive_rate = rep(NA, 3),
      false_negative_rate = rep(NA, 3),
      mean_factor_score_cor = rep(NA, 3),
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

# Create parameter grid
params <- expand.grid(
  n = c(100, 500),
  p = c(50, 100, 200),
  m = c(3, 5),
  alpha = c(0.2, 0.5, 0.8),
  gamma = c(0, 1, 2),
  sigma_sq = c(0.1, 1, 5),
  rep = 1:100  # 100 repetitions per scenario
)

# Add seed for reproducibility
params$seed <- 1000 + (1:nrow(params))

# Run simulation in parallel
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

# Save full results
write_csv(results, file.path(output_dir, "full_results.csv"))

# Aggregate results
aggregated_results <- results %>%
  group_by(method, n, p, m, alpha, gamma, sigma_sq) %>%
  summarize(
    mean_mse = mean(mse),
    se_mse = sd(mse) / sqrt(n()),
    mean_loading_error = mean(loading_error),
    se_loading_error = sd(loading_error) / sqrt(n()),
    mean_false_positive_rate = mean(false_positive_rate),
    se_false_positive_rate = sd(false_positive_rate) / sqrt(n()),
    mean_false_negative_rate = mean(false_negative_rate),
    se_false_negative_rate = sd(false_negative_rate) / sqrt(n()),
    mean_factor_score_cor = mean(mean_factor_score_cor),
    se_factor_score_cor = sd(mean_factor_score_cor) / sqrt(n())
  )

write_csv(aggregated_results, file.path(output_dir, "aggregated_results.csv"))

# Create plots
plot_results <- function(data, x_var, y_var, y_se_var, title, y_lab) {
  ggplot(data, aes_string(x = x_var, y = y_var, color = "method", group = "method")) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes_string(ymin = paste0(y_var, "-", y_se_var), 
                            ymax = paste0(y_var, "+", y_se_var)),
                 width = 0.1) +
    facet_grid(n ~ p, labeller = label_both, scales = "free_y") +
    theme_bw() +
    labs(title = title, y = y_lab, x = x_var, color = "Method") +
    theme(legend.position = "bottom")
}

# Plot results by varying parameters
# Plot by alpha
plot_alpha_mse <- plot_results(
  filter(aggregated_results, m == 3, gamma == 1, sigma_sq == 1),
  "alpha", "mean_mse", "se_mse", 
  "Reconstruction Error by Sparsity", "Mean Squared Error"
)
ggsave(file.path(output_dir, "alpha_mse.png"), plot_alpha_mse, width = 10, height = 8)

# Plot by gamma
plot_gamma_loading <- plot_results(
  filter(aggregated_results, m == 3, alpha == 0.5, sigma_sq == 1),
  "gamma", "mean_loading_error", "se_loading_error", 
  "Loading Error by Eigenvalue Decay", "Loading Error"
)
ggsave(file.path(output_dir, "gamma_loading.png"), plot_gamma_loading, width = 10, height = 8)

# Plot by sigma_sq
plot_sigma_fp <- plot_results(
  filter(aggregated_results, m == 3, alpha == 0.5, gamma == 1),
  "sigma_sq", "mean_false_positive_rate", "se_false_positive_rate", 
  "False Positive Rate by Noise Level", "False Positive Rate"
)
ggsave(file.path(output_dir, "sigma_fp.png"), plot_sigma_fp, width = 10, height = 8)

# Plot by number of factors
plot_m_factor_cor <- plot_results(
  filter(aggregated_results, alpha == 0.5, gamma == 1, sigma_sq == 1),
  "m", "mean_factor_score_cor", "se_factor_score_cor", 
  "Factor Score Correlation by Number of Factors", "Mean Factor Score Correlation"
)
ggsave(file.path(output_dir, "m_factor_cor.png"), plot_m_factor_cor, width = 10, height = 8)

# Create comprehensive report
cat("## Simulation Study: Comparing PCA, Sparse PCA, and EFA\n\n", file = file.path(output_dir, "report.md"))
cat("### Simulation Design\n\n", file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Number of variables (p): %s\n", paste(unique(params$p), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Sample sizes (n): %s\n", paste(unique(params$n), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Number of true factors (m): %s\n", paste(unique(params$m), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Sparsity fraction (alpha): %s\n", paste(unique(params$alpha), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Eigenvalue decay (gamma): %s\n", paste(unique(params$gamma), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Noise level (sigma_sq): %s\n", paste(unique(params$sigma_sq), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Repetitions per scenario: %d\n\n", max(params$rep)), file = file.path(output_dir, "report.md"), append = TRUE)

cat("### Summary of Results\n\n", file = file.path(output_dir, "report.md"), append = TRUE)
cat("See generated plots for visual comparisons.\n\n", file = file.path(output_dir, "report.md"), append = TRUE)

cat("Simulation completed. Results saved in", output_dir, "directory.\n")