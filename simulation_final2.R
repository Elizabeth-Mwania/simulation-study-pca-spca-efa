# PCA,SPCA WITH FIXED ALPHA=0.5, 4 EFA VARIANTS, AND REGULARIZED EFA

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
library(glmnet) # For regularized estimation

# Set up output directory
output_dir <- "simulation_final2_100"
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

# --- Modified Sparse PCA with Fixed Alpha = 0.5 ---
apply_sparse_pca_fixed <- function(X, m, fixed_alpha = 0.00000000000000001) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Center the data
  full_means <- colMeans(X)
  X_centered <- scale(X, center = full_means, scale = FALSE)
  
  tryCatch({
    # Fit sparse PCA with fixed alpha
    fit <- sparsepca::spca(X_centered, k = m, alpha = fixed_alpha, center = FALSE)
    
    loadings <- fit$loadings
    scores <- X_centered %*% loadings
    X_reconstructed <- scores %*% t(loadings) + matrix(full_means, n, p, byrow = TRUE)
    
    return(list(
      loadings = loadings,
      scores = scores,
      X_reconstructed = X_reconstructed,
      used_alpha = fixed_alpha  # renamed from optimal_alpha to used_alpha
    ))
  }, error = function(e) {
    warning("Sparse PCA fit failed with alpha = ", fixed_alpha, ". Returning NA matrices. Error: ", e$message)
    return(list(
      loadings = matrix(NA, p, m),
      scores = matrix(NA, n, m),
      X_reconstructed = matrix(NA, n, p),
      used_alpha = fixed_alpha
    ))
  })
}

# NEW FUNCTION: Apply Regularized EFA using glmnet
apply_regularized_efa <- function(X, m, alpha = 0.5, lambda = 0.1) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize loadings matrix and scores
  loadings <- matrix(NA, nrow = p, ncol = m)
  scores <- matrix(NA, nrow = n, ncol = m)
  
  # Center the data
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  
  tryCatch({
    # First get initial factor scores using standard FA method
    initial_fa <- fa(X_centered, nfactors = m, rotate = "none", fm = "minres")
    initial_scores <- initial_fa$scores
    
    # Iterate through each factor and fit a regularized regression
    for (j in 1:m) {
      # For each factor j, regress all variables on the factor score
      fit <- glmnet(y = X_centered, x = as.matrix(initial_scores[, j, drop = FALSE]), 
                   alpha = alpha, lambda = lambda)
      
      # Extract the regularized loadings for factor j
      loadings[, j] <- as.vector(coef(fit)[-1])  # Skip intercept
    }
    
    # Apply varimax rotation to improve interpretability
    rot_result <- varimax(loadings)
    loadings_rotated <- loadings %*% rot_result$rotmat
    
    # Compute factor scores based on rotated loadings
    # Using regression method to get scores
    S <- cov(X_centered)
    scores_rotated <- X_centered %*% solve(S) %*% loadings_rotated
    
    # Reconstruct the data
    X_reconstructed <- scores_rotated %*% t(loadings_rotated) + matrix(rep(colMeans(X), each = n), nrow = n)
    
    return(list(
      loadings = loadings_rotated,
      scores = scores_rotated,
      X_reconstructed = X_reconstructed,
      lambda = lambda,
      alpha = alpha
    ))
  }, error = function(e) {
    warning("Regularized EFA failed. Returning NA matrices. Error: ", e$message)
    return(list(
      loadings = matrix(NA, p, m),
      scores = matrix(NA, n, m),
      X_reconstructed = matrix(NA, n, p),
      lambda = lambda,
      alpha = alpha
    ))
  })
}

# Apply EFA with multiple rotation methods
apply_efa <- function(X, m) {
  # Define rotation methods to compare
  rotation_methods <- c("varimax", "oblimin", "quartimin", "quartimax")
  
  # Initialize results list
  efa_results <- list()
  
  for (rot in rotation_methods) {
    # Apply EFA using maximum likelihood with error handling
    tryCatch({
      efa_result <- fa(X, nfactors = m, rotate = rot, fm = "ml")
      
      loadings <- efa_result$loadings
      scores <- efa_result$scores
      
      # Reconstruct data
      X_reconstructed <- scores %*% t(loadings) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
      
      efa_results[[rot]] <- list(
        loadings = loadings,
        scores = scores,
        X_reconstructed = X_reconstructed,
        rotation = rot,
        convergence = TRUE
      )
      
    }, error = function(e) {
      # Fallback to minimum residual if ML fails
      tryCatch({
        efa_result <- fa(X, nfactors = m, rotate = rot, fm = "minres")
        
        loadings <- efa_result$loadings
        scores <- efa_result$scores
        
        # Reconstruct data
        X_reconstructed <- scores %*% t(loadings) + matrix(rep(colMeans(X), each = nrow(X)), nrow = nrow(X))
        
        efa_results[[rot]] <- list(
          loadings = loadings,
          scores = scores,
          X_reconstructed = X_reconstructed,
          rotation = rot,
          convergence = TRUE
        )
        
      }, error = function(e) {
        # If both methods fail, return NA values with proper dimensions
        p <- ncol(X)
        n <- nrow(X)
        
        # Create NA matrices with proper dimensions
        loadings <- matrix(NA, nrow = p, ncol = m)
        scores <- matrix(NA, nrow = n, ncol = m)
        X_reconstructed <- matrix(NA, nrow = n, ncol = p)
        
        warning(paste("EFA with rotation", rot, "failed to converge. Returning NA matrices."))
        efa_results[[rot]] <- list(
          loadings = loadings,
          scores = scores,
          X_reconstructed = X_reconstructed,
          rotation = rot,
          convergence = FALSE
        )
      })
    })
  }
  
  return(efa_results)
}

# Evaluation metrics with rotation support
evaluate_results <- function(true_X, reconstructed_X, true_L, estimated_L, true_F, estimated_scores, true_L_is_zero) {
  # Check for NA values
  if(any(is.na(estimated_L)) || any(is.na(estimated_scores)) || any(is.na(reconstructed_X))) {
    return(list(
      mse = NA,
      loading_error = NA,
      loading_alignment = NA,
      false_positive_rate = NA,
      false_negative_rate = NA,
      mean_factor_score_cor = NA,
      factor_correlations = NA
    ))
  }
  
  # 1. Reconstruction error (MSE)
  mse <- mean((true_X - reconstructed_X)^2)
  
  # 2. Loading matrix evaluation
  loading_metrics <- tryCatch({
    # Procrustes rotation to align estimated and true loadings
    proc <- vegan::procrustes(true_L, estimated_L)
    aligned_loadings <- proc$Yrot
    
    # Sum of squared differences after alignment
    loading_error <- proc$ss
    
    # Tucker's congruence coefficient for factor alignment
    congruence <- diag(cor(true_L, aligned_loadings))
    mean_congruence <- mean(abs(congruence))
    
    list(
      loading_error = loading_error,
      loading_alignment = mean_congruence,
      factor_correlations = congruence
    )
  }, error = function(e) {
    warning("Loading matrix evaluation failed: ", e$message)
    list(
      loading_error = NA,
      loading_alignment = NA,
      factor_correlations = rep(NA, ncol(true_L))
    )
  })
  
  # 3. Sparsity pattern evaluation
  sparsity_metrics <- tryCatch({
    estimated_L_is_zero <- abs(estimated_L) < 0.05  # Threshold for "zero" loadings
    
    # False positive rate (type I error)
    fp <- ifelse(sum(true_L_is_zero) > 0,
                sum(estimated_L_is_zero == FALSE & true_L_is_zero == TRUE) / sum(true_L_is_zero),
                NA)
    
    # False negative rate (type II error)
    fn <- ifelse(sum(!true_L_is_zero) > 0,
                sum(estimated_L_is_zero == TRUE & true_L_is_zero == FALSE) / sum(!true_L_is_zero),
                NA)
    
    list(
      false_positive_rate = fp,
      false_negative_rate = fn
    )
  }, error = function(e) {
    warning("Sparsity pattern evaluation failed: ", e$message)
    list(
      false_positive_rate = NA,
      false_negative_rate = NA
    )
  })
  
  # 4. Factor score evaluation
  score_metrics <- tryCatch({
    # Correlation between true and estimated factor scores
    factor_cors <- diag(cor(true_F, estimated_scores))
    mean_factor_cor <- mean(factor_cors)
    
    list(
      mean_factor_score_cor = mean_factor_cor,
      factor_correlations = factor_cors
    )
  }, error = function(e) {
    warning("Factor score evaluation failed: ", e$message)
    list(
      mean_factor_score_cor = NA,
      factor_correlations = rep(NA, ncol(true_F))
    )
  })
  
  # Return all metrics
  return(list(
    mse = mse,
    loading_error = loading_metrics$loading_error,
    loading_alignment = loading_metrics$loading_alignment,
    false_positive_rate = sparsity_metrics$false_positive_rate,
    false_negative_rate = sparsity_metrics$false_negative_rate,
    mean_factor_score_cor = score_metrics$mean_factor_score_cor,
    factor_correlations = score_metrics$factor_correlations,
    loading_correlations = loading_metrics$factor_correlations
  ))
}

# Function to determine true zero loadings (threshold)
get_true_zero_loadings <- function(L) {
  return(abs(L) < 1e-8)  # More stringent threshold for true zeros
}

# Modified run_simulation function with consistent output structure and fixed alpha
run_simulation <- function(params) {
  # Define all possible columns upfront
  base_cols <- c("method", "mse", "loading_error", "false_positive_rate", 
                "false_negative_rate", "mean_factor_score_cor", "used_alpha",
                "rotation_method", "n", "p", "m", "alpha", "gamma", "sigma_sq", 
                "rep", "error_message", "lambda")
  
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
    fixed_alpha <- 0.00000000000000001  # Using fixed alpha = 0.5 for Sparse PCA
    reg_alpha <- 0.00000000000000001   # Using alpha = 0.5 for regularized EFA (elasticnet balance)
    reg_lambda <- 0.1   # Using lambda = 0.1 for regularized EFA
    
    # Generate data
    data <- generate_data(n, p, m, alpha, gamma, sigma_sq)
    X <- data$X
    true_L <- data$L
    true_F <- data$F
    
    # Identify true zero loadings
    true_L_is_zero <- get_true_zero_loadings(true_L)
    
    # Apply methods
    pca_result <- apply_pca(X, m)
    sparse_pca_result <- apply_sparse_pca_fixed(X, m, fixed_alpha)
    efa_results <- apply_efa(X, m)
    
    # Apply regularized EFA
    reg_efa_result <- apply_regularized_efa(X, m, reg_alpha, reg_lambda)
    
    # Evaluate results with consistent structure
    evaluate_with_structure <- function(true_X, reconstructed_X, true_L, estimated_L, 
                                      true_F, estimated_scores, true_L_is_zero) {
      res <- evaluate_results(true_X, reconstructed_X, true_L, estimated_L, 
                            true_F, estimated_scores, true_L_is_zero)
      # Ensure all expected metrics exist
      for (metric in c("mse", "loading_error", "false_positive_rate", 
                      "false_negative_rate", "mean_factor_score_cor")) {
        if (is.null(res[[metric]])) res[[metric]] <- NA
      }
      return(res)
    }
    
    pca_eval <- evaluate_with_structure(X, pca_result$X_reconstructed, true_L, 
                                      pca_result$loadings, true_F, pca_result$scores, 
                                      true_L_is_zero)
    
    sparse_pca_eval <- evaluate_with_structure(X, sparse_pca_result$X_reconstructed, true_L, 
                                             sparse_pca_result$loadings, true_F, 
                                             sparse_pca_result$scores, true_L_is_zero)
    
    # Evaluate all EFA rotation methods
    efa_evals <- lapply(efa_results, function(efa_res) {
      evaluate_with_structure(X, efa_res$X_reconstructed, true_L, 
                            efa_res$loadings, true_F, efa_res$scores, 
                            true_L_is_zero)
    })
    
    # Evaluate regularized EFA
    reg_efa_eval <- evaluate_with_structure(X, reg_efa_result$X_reconstructed, true_L,
                                          reg_efa_result$loadings, true_F,
                                          reg_efa_result$scores, true_L_is_zero)
    
    # Create consistent result rows for each method
    create_result_row <- function(method_name, eval_res, used_alpha = NA, rot_method = NA, lambda = NA) {
      data.frame(
        method = method_name,
        mse = eval_res$mse,
        loading_error = eval_res$loading_error,
        false_positive_rate = eval_res$false_positive_rate,
        false_negative_rate = eval_res$false_negative_rate,
        mean_factor_score_cor = eval_res$mean_factor_score_cor,
        used_alpha = used_alpha,  # renamed from optimal_alpha to used_alpha
        rotation_method = rot_method,
        n = n,
        p = p,
        m = m,
        alpha = alpha,
        gamma = gamma,
        sigma_sq = sigma_sq,
        rep = rep,
        lambda = lambda,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }
    
    # Combine all results
    results <- rbind(
      create_result_row("PCA", pca_eval),
      create_result_row("Sparse PCA (Fixed)", sparse_pca_eval, sparse_pca_result$used_alpha),
      create_result_row("EFA (varimax)", efa_evals$varimax, rot_method = "varimax"),
      create_result_row("EFA (oblimin)", efa_evals$oblimin, rot_method = "oblimin"),
      create_result_row("EFA (quartimin)", efa_evals$quartimin, rot_method = "quartimin"),
      create_result_row("EFA (quartimax)", efa_evals$quartimax, rot_method = "quartimax"),
      create_result_row("Regularized EFA", reg_efa_eval, used_alpha = reg_alpha, rot_method = "varimax", lambda = reg_lambda)
    )
    
    # Save loadings
    file_prefix <- sprintf("%s/results_n%d_p%d_m%d_alpha%.1f_gamma%d_sigma%.1f_rep%d", 
                         output_dir, n, p, m, alpha, gamma, sigma_sq, rep)
    
    loadings_data <- list(
      true_L = true_L,
      pca_L = pca_result$loadings,
      sparse_pca_L = sparse_pca_result$loadings,
      efa_varimax_L = efa_results$varimax$loadings,
      efa_oblimin_L = efa_results$oblimin$loadings,
      efa_quartimin_L = efa_results$quartimin$loadings,
      efa_quartimax_L = efa_results$quartimax$loadings,
      reg_efa_L = reg_efa_result$loadings
    )
    
    saveRDS(loadings_data, paste0(file_prefix, "_loadings.rds"))
    
    return(results)
    
  }, error = function(e) {
    warning("Error in simulation run: ", e$message)
    # Return a consistent data frame structure even on error
    methods <- c("PCA", "Sparse PCA (Fixed)", "EFA (varimax)", "EFA (oblimin)", 
                "EFA (quartimin)", "EFA (quartimax)", "Regularized EFA")
    rotation_methods <- c(NA, NA, "varimax", "oblimin", "quartimin", "quartimax", "varimax")
    lambdas <- c(NA, NA, NA, NA, NA, NA, 0.1)
    
    data.frame(
      method = methods,
      mse = NA,
      loading_error = NA,
      false_positive_rate = NA,
      false_negative_rate = NA,
      mean_factor_score_cor = NA,
      used_alpha = NA,  # renamed from optimal_alpha to used_alpha
      rotation_method = rotation_methods,
      n = params$n,
      p = params$p,
      m = params$m,
      alpha = params$alpha,
      gamma = params$gamma,
      sigma_sq = params$sigma_sq,
      rep = params$rep,
      lambda = lambdas,
      error_message = e$message,
      stringsAsFactors = FALSE
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
  rep = 1:50  # 1 repetition per scenario (change to more if needed)
)

# Add seed for reproducibility
params$seed <- 1000 + (1:nrow(params))

# Run simulation in parallel with proper result combining
start_time <- Sys.time()
cat("Starting simulation at", as.character(start_time), "\n")

results_list <- foreach(i = 1:nrow(params), .packages = c("psych", "sparsepca", "GPArotation", "dplyr", "vegan", "glmnet"),
                       .errorhandling = "pass") %dopar% {
  cat("Running scenario", i, "of", nrow(params), "\n")
  run_simulation(params[i, ])
}

# Combine results safely
results <- tryCatch({
  do.call(rbind, results_list)
}, error = function(e) {
  warning("Error combining results: ", e$message)
  # Fallback - manually combine results
  final_results <- data.frame()
  for (res in results_list) {
    if (is.data.frame(res)) {
      final_results <- rbind(final_results, res)
    }
  }
  return(final_results)
})

end_time <- Sys.time()
cat("Simulation completed at", as.character(end_time), "\n")
cat("Total time:", difftime(end_time, start_time, units = "hours"), "hours\n")

# Stop cluster
stopCluster(cl)

# Ensure results is a proper data frame before saving
if (!is.data.frame(results)) {
  warning("Results is not a data frame - attempting to convert")
  results <- as.data.frame(results)
}

# Save full results
tryCatch({
  write.csv(results, file.path(output_dir, "full_results.csv"), row.names = FALSE)
}, error = function(e) {
  warning("Error writing CSV: ", e$message)
  saveRDS(results, file.path(output_dir, "full_results.rds"))
})

# Aggregate results
aggregated_results <- results %>%
  group_by(method, n, p, m, alpha, gamma, sigma_sq, rotation_method, lambda) %>%
  summarize(
    mean_mse = mean(mse, na.rm = TRUE),
    se_mse = sd(mse, na.rm = TRUE) / sqrt(sum(!is.na(mse))),
    mean_loading_error = mean(loading_error, na.rm = TRUE),
    se_loading_error = sd(loading_error, na.rm = TRUE) / sqrt(sum(!is.na(loading_error))),
    mean_false_positive_rate = mean(false_positive_rate, na.rm = TRUE),
    se_false_positive_rate = sd(false_positive_rate, na.rm = TRUE) / sqrt(sum(!is.na(false_positive_rate))),
    mean_false_negative_rate = mean(false_negative_rate, na.rm = TRUE),
    se_false_negative_rate = sd(false_negative_rate, na.rm = TRUE) / sqrt(sum(!is.na(false_negative_rate))),
    mean_factor_score_cor = mean(mean_factor_score_cor, na.rm = TRUE),
    se_factor_score_cor = sd(mean_factor_score_cor, na.rm = TRUE) / sqrt(sum(!is.na(mean_factor_score_cor)))
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
cat("## Simulation Study: Comparing PCA, Sparse PCA (Fixed Alpha=0.5), and EFA\n\n", file = file.path(output_dir, "report.md"))
cat("### Simulation Design\n\n", file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Number of variables (p): %s\n", paste(unique(params$p), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Sample sizes (n): %s\n", paste(unique(params$n), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Number of true factors (m): %s\n", paste(unique(params$m), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Sparsity fraction (alpha): %s\n", paste(unique(params$alpha), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Eigenvalue decay (gamma): %s\n", paste(unique(params$gamma), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Noise level (sigma_sq): %s\n", paste(unique(params$sigma_sq), collapse = ", ")), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Repetitions per scenario: %d\n\n", max(params$rep)), file = file.path(output_dir, "report.md"), append = TRUE)
cat(sprintf("- Fixed Alpha for Sparse PCA: 0.5\n\n"), file = file.path(output_dir, "report.md"), append = TRUE)

cat("### Summary of Results\n\n", file = file.path(output_dir, "report.md"), append = TRUE)
cat("See generated plots for visual comparisons.\n\n", file = file.path(output_dir, "report.md"), append = TRUE)

cat("Simulation completed. Results saved in", output_dir, "directory.\n")