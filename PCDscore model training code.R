# =============================================================================
# R Script:  Development and Evaluation of the PCDscore
# Author: [Xianwen Lin]
# Date: 2026-02-19
# R Version: 4.5.1
# =============================================================================


# 1. Environment Setup####


# Load all necessary libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(ggsci)
library(survival)
library(randomForestSRC)
library(snowfall)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(devtools)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(ggbreak)
library(tidyr)
library(edgeR)
library(limma)
library(survminer)
library(stringi)
library(tidyverse)
library(ggpubr)
library(beepr)
library(pheatmap)
library(ggsignif)
library(future.apply)
library(gplots)
library(DESeq2)
library(ggrepel)
library(Rcpp)
library(rms)
library(pec)
library(ggDCA)
library(foreign)
library(regplot)
library(timeROC)
library(caret)
library(obliqueRSF)
library(remotes)
library(aorsf)
library(xgboost)
library(party)
library(partykit)
library(openxlsx)
library(truncnorm)




# 2. Global Constants and Helper Functions####

# Set global options
options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
# Define file paths (modify as needed)
BASE_DIR <- "D:/MR-CRC"
COHORT_DIR <- file.path(BASE_DIR, "cohort")
CLI_DIR <- "D:/MR-CRC/cli"

# Set random seed for reproducibility
MODEL_SEED <- 123
SET_SEED <- 123





# 3. Data Loading and Preprocessing####


#' Normalize expression matrix using various methods
#'
#' @param exp_mat Expression matrix (genes x samples)
#' @param method Normalization method: "zscore", "minmax", "quantile"
#' @param ref_mean Reference mean vector (for test set)
#' @param ref_sd Reference sd vector (for test set)
#' @return Normalized expression matrix with attributes
normalize_expression_data <- function(exp_mat,
                                      method = "zscore",
                                      ref_mean = NULL,
                                      ref_sd = NULL) {
  # Input validation and conversion
  if (is.list(exp_mat) && !is.data.frame(exp_mat)) {
    exp_mat <- tryCatch(
      as.matrix(exp_mat),
      error = function(e) stop("Cannot convert list to matrix: ", e$message)
    )
  }
  if (!is.matrix(exp_mat)) exp_mat <- as.matrix(exp_mat)
  
  # Replace infinite/NaN with NA
  exp_mat[is.infinite(exp_mat) | is.nan(exp_mat)] <- NA
  
  # Remove constant or all-NA genes
  valid_genes <- apply(exp_mat, 1, function(x) {
    !all(is.na(x)) && sd(x, na.rm = TRUE) > 1e-10
  })
  if (sum(valid_genes) < nrow(exp_mat)) {
    warning(sprintf("Removed %d invalid genes", nrow(exp_mat) - sum(valid_genes)))
    exp_mat <- exp_mat[valid_genes, , drop = FALSE]
  }
  
  # Perform normalization
  normalized_mat <- switch(
    method,
    "zscore" = {
      if (is.null(ref_mean) || is.null(ref_sd)) {
        # Training mode
        gene_means <- rowMeans(exp_mat, na.rm = TRUE)
        gene_sds   <- apply(exp_mat, 1, sd, na.rm = TRUE)
        gene_sds[gene_sds == 0] <- 1
        norm_mat <- (exp_mat - gene_means) / gene_sds
        attr(norm_mat, "ref_params") <- list(mean = gene_means, sd = gene_sds,
                                             method = "zscore")
      } else {
        # Test mode
        norm_mat <- (exp_mat - ref_mean) / ref_sd
      }
      norm_mat
    },
    "minmax" = {
      gene_mins <- apply(exp_mat, 1, min, na.rm = TRUE)
      gene_maxs <- apply(exp_mat, 1, max, na.rm = TRUE)
      range_zero <- gene_maxs - gene_mins == 0
      if (any(range_zero)) {
        warning(sprintf("%d genes have constant values", sum(range_zero)))
        gene_maxs[range_zero] <- gene_mins[range_zero] + 1
      }
      norm_mat <- (exp_mat - gene_mins) / (gene_maxs - gene_mins)
      attr(norm_mat, "ref_params") <- list(min = gene_mins, max = gene_maxs,
                                           method = "minmax")
      norm_mat
    },
    "quantile" = {
      if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        stop("Package 'preprocessCore' is required for quantile normalization")
      }
      norm_mat <- preprocessCore::normalize.quantiles(exp_mat)
      rownames(norm_mat) <- rownames(exp_mat)
      colnames(norm_mat) <- colnames(exp_mat)
      attr(norm_mat, "ref_params") <- list(method = "quantile")
      norm_mat
    },
    stop("Unsupported normalization method: ", method)
  )
  
  # Stabilize output
  normalized_mat[is.na(normalized_mat) | is.infinite(normalized_mat)] <- 0
  return(normalized_mat)
}

#' Build a normalized dataset from expression and clinical data
#'
#' @param exp_mat Expression matrix (genes x samples) or combined data.frame
#' @param cli_df Clinical data.frame (ignored if data_type = "combined")
#' @param valid_genes Vector of gene symbols to retain
#' @param normalize_method Normalization method
#' @param ref_params Reference parameters from training set (optional)
#' @param data_type "separate" (exp_mat and cli_df separate) or "combined"
#' @return Normalized dataset (samples x genes) with OS.time and OS columns
build_normalized_dataset <- function(exp_mat,
                                     cli_df = NULL,
                                     valid_genes,
                                     normalize_method = "zscore",
                                     ref_params = NULL,
                                     data_type = "separate") {
  if (data_type == "combined") {
    # Combined data: first two columns are OS.time and OS
    combined <- exp_mat
    # Identify time and status columns
    time_col <- grep("^OS\\.time$|^time$|^Time$", colnames(combined), value = TRUE)[1]
    status_col <- grep("^OS$|^status$|^Status$|^OS\\.status$", colnames(combined), value = TRUE)[1]
    if (is.na(time_col) || is.na(status_col)) {
      warning("Standard column names not found, using first two columns as time and status")
      time_col <- colnames(combined)[1]
      status_col <- colnames(combined)[2]
    }
    cli_data <- data.frame(
      ID = rownames(combined),
      OS.time = combined[[time_col]],
      OS = combined[[status_col]],
      stringsAsFactors = FALSE
    )
    exp_cols <- setdiff(colnames(combined), c(time_col, status_col))
    exp_mat <- t(combined[, exp_cols, drop = FALSE])
    cli_df <- cli_data
  }
  
  # Validate inputs
  if (is.null(exp_mat) || (is.matrix(exp_mat) && nrow(exp_mat) == 0)) {
    stop("Expression matrix is empty")
  }
  if (is.null(cli_df) || nrow(cli_df) == 0) {
    stop("Clinical data is empty")
  }
  if (!is.matrix(exp_mat)) exp_mat <- as.matrix(exp_mat)
  
  # Subset to valid genes
  common_genes <- intersect(valid_genes, rownames(exp_mat))
  if (length(common_genes) == 0) {
    stop("No valid genes found")
  }
  exp_mat <- exp_mat[common_genes, , drop = FALSE]
  
  # Normalize expression
  if (!is.null(ref_params)) {
    cat(sprintf("Normalizing %s dataset using reference parameters (method: %s)\n",
                ifelse(data_type == "combined", "validation", "test"),
                ref_params$method))
    normalized_exp <- normalize_expression_data(exp_mat,
                                                method = ref_params$method,
                                                ref_mean = ref_params$mean,
                                                ref_sd = ref_params$sd)
  } else {
    cat(sprintf("Training set normalization (method: %s)\n", normalize_method))
    normalized_exp <- normalize_expression_data(exp_mat, method = normalize_method)
    ref_params <- attr(normalized_exp, "ref_params")
  }
  
  # Transpose to sample x gene
  exp_df <- as.data.frame(t(normalized_exp))
  exp_df$ID <- rownames(exp_df)
  
  # Merge with clinical data
  # Identify time and status columns in cli_df
  avail_cols <- colnames(cli_df)
  time_col <- if ("OS.time" %in% avail_cols) {
    "OS.time"
  } else if ("time" %in% avail_cols) {
    "time"
  } else if ("Time" %in% avail_cols) {
    "Time"
  } else {
    grep("time", avail_cols, ignore.case = TRUE, value = TRUE)[1]
  }
  status_col <- if ("OS" %in% avail_cols) {
    "OS"
  } else if ("OS.status" %in% avail_cols) {
    "OS.status"
  } else if ("status" %in% avail_cols) {
    "status"
  } else if ("Status" %in% avail_cols) {
    "Status"
  } else {
    grep("status", avail_cols, ignore.case = TRUE, value = TRUE)[1]
  }
  if (is.na(time_col) || !time_col %in% avail_cols) {
    stop("Survival time column not found")
  }
  if (is.na(status_col) || !status_col %in% avail_cols) {
    stop("Survival status column not found")
  }
  
  cli_sub <- cli_df[, c("ID", time_col, status_col)]
  colnames(cli_sub) <- c("ID", "OS.time", "OS")
  cli_sub$OS.time <- as.numeric(as.character(cli_sub$OS.time))
  cli_sub$OS      <- as.numeric(as.character(cli_sub$OS))
  
  merged_df <- merge(cli_sub, exp_df, by = "ID", all.x = TRUE)
  rownames(merged_df) <- merged_df$ID
  merged_df$ID <- NULL
  
  # Attach normalization parameters if this is training set
  if (is.null(ref_params)) {
    attr(merged_df, "normalization_params") <- ref_params
  }
  
  return(merged_df)
}


cat("\n===== 2. Loading and preparing data =====\n")

# 3.1 Read key genes
key_gene_file <- file.path(BASE_DIR, "key_gene.xlsx")
crc_keygene <- openxlsx::read.xlsx(key_gene_file,
                                   sheet = 1,
                                   startRow = 1,
                                   colNames = TRUE,
                                   rowNames = FALSE)
key_genes <- unique(crc_keygene$Gene_Symbol)

# 3.2 Load TCGA data
load(file.path(COHORT_DIR, "TCGA_CRC_EXP.Rdata"))
load(file.path(COHORT_DIR, "TCGA_CRC_CLI.Rdata"))
colnames(TCGA_CRC_CLI)[35] <- "OS.time"
colnames(TCGA_CRC_CLI)[36] <- "OS"

# 3.3 Load NFYY and GEO test datasets
load(file.path(COHORT_DIR, "NFYY_CRC_test_dataset.Rdata"))
load(file.path(COHORT_DIR, "GSE17536_CRC_test_dataset.Rdata"))
load(file.path(COHORT_DIR, "GSE17537_CRC_test_dataset.Rdata"))
load(file.path(COHORT_DIR, "GSE29621_CRC_test_dataset.Rdata"))

# 3.4 Split TCGA into training and test (75/25)
set.seed(SET_SEED)
train_indices <- sample(ncol(TCGA_CRC_EXP),
                        size = floor(0.75 * ncol(TCGA_CRC_EXP)))
train_exp <- TCGA_CRC_EXP[, train_indices, drop = FALSE]
train_cli <- TCGA_CRC_CLI[train_indices, , drop = FALSE]
test_exp  <- TCGA_CRC_EXP[, -train_indices, drop = FALSE]
test_cli  <- TCGA_CRC_CLI[-train_indices, , drop = FALSE]

# 3.5 Intersect key genes with expression matrix
valid_genes <- intersect(key_genes, rownames(train_exp))
cat("Number of valid key genes:", length(valid_genes), "\n")

# 3.6 Build normalized datasets
TCGA_CRC_train_dataset <- build_normalized_dataset(
  exp_mat = as.matrix(train_exp),
  cli_df  = train_cli,
  valid_genes = valid_genes,
  normalize_method = "zscore",
  data_type = "separate"
)
norm_params <- attr(TCGA_CRC_train_dataset, "normalization_params")

TCGA_CRC_test_dataset <- build_normalized_dataset(
  exp_mat = as.matrix(test_exp),
  cli_df  = test_cli,
  valid_genes = valid_genes,
  normalize_method = "zscore",
  ref_params = norm_params,
  data_type = "separate"
)

# 3.7 Normalize external validation sets (combined format)
NFYY_CRC_test_dataset_norm <- build_normalized_dataset(
  exp_mat = NFYY_CRC_test_dataset,
  cli_df  = NULL,
  valid_genes = valid_genes,
  normalize_method = "zscore",
  ref_params = norm_params,
  data_type = "combined"
)

GSE17536_CRC_test_dataset_norm <- build_normalized_dataset(
  exp_mat = GSE17536_CRC_test_dataset,
  cli_df  = NULL,
  valid_genes = valid_genes,
  normalize_method = "zscore",
  ref_params = norm_params,
  data_type = "combined"
)

GSE17537_CRC_test_dataset_norm <- build_normalized_dataset(
  exp_mat = GSE17537_CRC_test_dataset,
  cli_df  = NULL,
  valid_genes = valid_genes,
  normalize_method = "zscore",
  ref_params = norm_params,
  data_type = "combined"
)

GSE29621_CRC_test_dataset_norm <- build_normalized_dataset(
  exp_mat = GSE29621_CRC_test_dataset,
  cli_df  = NULL,
  valid_genes = valid_genes,
  normalize_method = "zscore",
  ref_params = norm_params,
  data_type = "combined"
)

# 3.8 Create list of all datasets and extract training set as data.frame
trainlist <- list(
  Train           = TCGA_CRC_train_dataset,
  Test_TCGA       = TCGA_CRC_test_dataset,
  Test_NFYY       = NFYY_CRC_test_dataset_norm,
  Test_GSE17536   = GSE17536_CRC_test_dataset_norm,
  Test_GSE17537   = GSE17537_CRC_test_dataset_norm,
  Test_GSE29621   = GSE29621_CRC_test_dataset_norm
)
train <- as.data.frame(TCGA_CRC_train_dataset)




# 4. Utility Functions for Model Evaluation####



#' Remove rows with infinite values in the RS column
remove_inf <- function(df) {
  if ("RS" %in% names(df)) {
    if (is.numeric(df$RS)) {
      df <- df[is.finite(df$RS), , drop = FALSE]
    } else {
      rs_char <- as.character(df$RS)
      df <- df[!rs_char %in% c("Inf", "-Inf"), , drop = FALSE]
    }
  }
  df
}

#' Safely calculate C-index using survcomp::concordance.index
calculate_cindex <- function(data) {
  required <- c("OS.time", "OS", "RS")
  if (!all(required %in% colnames(data))) {
    warning("Missing required columns: OS.time, OS, or RS")
    return(NA_real_)
  }
  if (nrow(data) == 0) {
    warning("Empty dataset")
    return(NA_real_)
  }
  valid_data <- data[is.finite(data$RS) & !is.na(data$RS), , drop = FALSE]
  if (nrow(valid_data) == 0) {
    warning("All RS values are invalid")
    return(NA_real_)
  }
  if (sum(valid_data$OS) < 2) {
    warning("Insufficient events (<2)")
    return(NA_real_)
  }
  if (any(valid_data$OS.time <= 0)) {
    warning("Non-positive survival times exist; removing those rows")
    valid_data <- valid_data[valid_data$OS.time > 0, , drop = FALSE]
    if (nrow(valid_data) == 0) return(NA_real_)
  }
  tryCatch({
    cindex_obj <- survcomp::concordance.index(
      x = valid_data$RS,
      surv.time = valid_data$OS.time,
      surv.event = valid_data$OS,
      method = "noether"
    )
    round(cindex_obj$c.index, 3)
  }, error = function(e) {
    warning("C-index calculation failed: ", e$message)
    NA_real_
  })
}

#' Safe subsetting of data frame
safe_subset <- function(data, required_vars) {
  missing <- setdiff(required_vars, colnames(data))
  if (length(missing) > 0) {
    warning("Missing variables: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  data[, required_vars, drop = FALSE]
}

#' Create a template data.frame of NA results
create_na_results <- function(model_name, trainlist) {
  na_cc <- rep(NA_real_, length(trainlist))
  names(na_cc) <- names(trainlist)
  data.frame(ID = names(na_cc),
             Cindex = na_cc,
             Model = model_name,
             stringsAsFactors = FALSE)
}


# 5. Parameter Configuration####

# Set seed for reproducibility
set.seed(MODEL_SEED)

# Define a small grid of hyperparameter combinations.
param_grid <- expand.grid(
  subsample       = c(0.1, 0.5, 1.0),
  max_depth       = c(1, 5, 10),
  min_child_weight = c(1, 5, 10),
  gamma           = c(0, 0.5, 1),
  eta             = c(0.01, 0.2, 0.5),
  lambda          = c(0, 3, 10),
  alpha           = c(0, 0.1, 1)
)

# Remove duplicate rows if any
param_grid <- unique(param_grid)

# Number of boosting rounds (fixed during tuning)
n_rounds <- 100

# Create 10-fold indices (Simple random folds)
folds <- createFolds(1:nrow(train), k = 10, list = TRUE, returnTrain = FALSE)

# Print header for tuning process with descriptive column names
cat("\n", rep("=", 90), "\n", sep="")
cat("   XGBoost Hyperparameter Tuning via 10-Fold Cross-Validation\n")
cat(rep("=", 90), "\n\n", sep="")
cat(sprintf("Total parameter combinations: %d\n", nrow(param_grid)))
cat(sprintf("Cross-validation folds: %d\n", length(folds)))
cat(sprintf("Total models to train: %d\n\n", nrow(param_grid) * length(folds)))

# Column headers
header <- sprintf("%-4s | %-14s | %-14s | %-16s | %-10s | %-10s | %-12s | %-11s | %-12s",
                  "No.", "Subsample_Ratio", "Max_Tree_Depth", "Min_Child_Weight",
                  "Gamma", "Eta", "Lambda_L2", "Alpha_L1", "Avg_Cindex_CV")
cat(header, "\n")
cat(rep("-", nchar(header)), "\n", sep="")

# Initialize counters for progress reporting
total_combos <- nrow(param_grid)
combos_done <- 0

# Vector to store average C-index for each parameter set
avg_cindex <- numeric(nrow(param_grid))

# Loop over each parameter combination
for (i in 1:nrow(param_grid)) {
  # Current parameters (the two fixed objectives are added later)
  current_params <- as.list(param_grid[i, ])
  current_params$objective   <- "survival:cox"
  current_params$eval_metric <- "cox-nloglik"
  
  fold_cindex <- numeric(length(folds))
  
  # Cross-validation inner loop
  for (j in 1:length(folds)) {
    # Split data into training and validation for this fold
    val_idx  <- folds[[j]]
    train_idx <- setdiff(1:nrow(train), val_idx)
    
    train_data <- train[train_idx, ]
    val_data   <- train[val_idx, ]
    
    # Prepare labels: positive time for events, negative time for censored
    train_labels <- ifelse(train_data$OS == 1,
                           train_data$OS.time,
                           -train_data$OS.time)
    
    # Create DMatrix for training
    dtrain <- xgb.DMatrix(data = as.matrix(train_data[, -c(1, 2)]),
                          label = train_labels)
    
    # Train model on this fold
    model <- xgb.train(params = current_params,
                       data   = dtrain,
                       nrounds = n_rounds,
                       verbose = 0)
    
    # Predict on validation set
    pred_rs <- predict(model, newdata = as.matrix(val_data[, -c(1, 2)]))
    
    # Build data frame for C-index calculation
    val_df <- data.frame(OS      = val_data$OS,
                         OS.time = val_data$OS.time,
                         RS      = pred_rs)
    
    # Remove any infinite values (should not occur, but safe)
    val_df <- remove_inf(val_df)
    
    # Compute C-index for this fold
    fold_cindex[j] <- calculate_cindex(val_df)
  }
  
  # Average C-index across all 10 folds
  avg_cindex[i] <- mean(fold_cindex, na.rm = TRUE)
  
  # Print a formatted row for this combination
  row_str <- sprintf("%-4d | %-14.2f | %-14d | %-16.1f | %-10.2f | %-10.2f | %-12.2f | %-11.2f | %-12.4f",
                     i,
                     current_params$subsample,
                     current_params$max_depth,
                     current_params$min_child_weight,
                     current_params$gamma,
                     current_params$eta,
                     current_params$lambda,
                     current_params$alpha,
                     avg_cindex[i])
  cat(row_str, "\n")
  
  # Update progress counter
  combos_done <- combos_done + 1
  
  # Print a progress summary every 10 combinations (and at the end)
  # This avoids cluttering the output with per‑fold progress.
  if (combos_done %% 10 == 0 || combos_done == total_combos) {
    cat(sprintf("Progress: %d/%d combinations (%.1f%%) completed\n",
                combos_done, total_combos, 100 * combos_done / total_combos))
  }
  
  # Repeat column headers every 10 rows to keep them visible when scrolling
  if ((i - 1) %% 10 == 0 && i > 1 && i < nrow(param_grid)) {
    cat("\n")  # blank line for separation
    cat(header, "\n")
    cat(rep("-", nchar(header)), "\n", sep="")
  }
}

# Select the parameter combination with the highest average C-index
best_idx <- which.max(avg_cindex)
best_params <- as.list(param_grid[best_idx, ])
best_params$objective   <- "survival:cox"
best_params$eval_metric <- "cox-nloglik"

# Print the best parameters with a star marker and re-print column headers
cat("\n", rep("=", 90), "\n", sep="")
cat("   Best Parameter Combination (★ indicates chosen set)\n")
cat(rep("=", 90), "\n\n", sep="")

# Re-print the column headers for clarity
cat(header, "\n")
cat(rep("-", nchar(header)), "\n", sep="")

# Print the best row with a star
row_best <- sprintf(" ★ %-2d | %-14.2f | %-14d | %-16.1f | %-10.2f | %-10.2f | %-12.2f | %-11.2f | %-12.4f",
                    best_idx,
                    best_params$subsample,
                    best_params$max_depth,
                    best_params$min_child_weight,
                    best_params$gamma,
                    best_params$eta,
                    best_params$lambda,
                    best_params$alpha,
                    avg_cindex[best_idx])
cat(row_best, "\n\n")

# Also create a summary data frame for user reference
tuning_results <- cbind(
  Combination_No = 1:nrow(param_grid),
  param_grid,
  Avg_Cindex_CV = avg_cindex
)

cat("\nFull tuning results (data frame with column names):\n")
print(tuning_results, row.names = FALSE)

# The best parameters become the final configuration for the full model
xgb_params_config <- best_params


xgb_params_config <- list(
  # Core objective: Cox proportional hazards for survival data
  objective = "survival:cox",
  # Evaluation metric: negative log-likelihood for Cox model
  eval_metric = "cox-nloglik",
  # Sample subsampling ratio per tree (prevents overfitting)
  subsample = 0.5,
  # Maximum tree depth (shallow trees for small datasets)
  max_depth = 5,
  # Minimum sum of instance weight (hessian) needed in a child node
  min_child_weight = 5,
  # Minimum loss reduction required to make a further partition
  gamma = 0.5,
  # Learning rate (shrinkage)
  eta = 0.2,
  # L2 regularization term on weights
  lambda = 3,
  # L1 regularization term on weights (encourages sparsity)
  alpha = 0.1
)

# Number of boosting rounds
n_rounds <- 100

rf_nodesize <- 5   # Node size for Random Survival Forest



# 6. Initialize Results Data Frame####


result <- data.frame()


# 7. Model Training####
                     



# 7.1 RSF (Random Survival Forest) family####


#### 1-1. RSF (Random Survival Forest) ####
cat("\n===== Training RSF =====\n")
set.seed(MODEL_SEED)
fit <- rfsrc(
  Surv(OS.time, OS) ~ .,
  data = trainlist$Train,
  ntree = 1000,
  nodesize = rf_nodesize,
  splitrule = 'logrank',
  importance = TRUE,
  proximity = TRUE,
  forest = TRUE,
  seed = MODEL_SEED
)

# Select optimal number of trees based on error rate
best <- which.min(fit$err.rate)
set.seed(MODEL_SEED)
fit <- rfsrc(
  Surv(OS.time, OS) ~ .,
  data = trainlist$Train,
  ntree = best,
  nodesize = rf_nodesize,
  splitrule = 'logrank',
  importance = TRUE,
  proximity = TRUE,
  forest = TRUE,
  seed = MODEL_SEED
)

# Predict on all cohorts
rs <- lapply(trainlist, function(x) {
  cbind(x[, 1:2], RS = predict(fit, newdata = x)$predicted)
})
rs <- lapply(rs, remove_inf)
cc <- sapply(rs, calculate_cindex)

cc_df <- data.frame(Cindex = cc) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF'
result <- rbind(result, cc_df)

#### 1-2. RSF + Elastic Net (Enet) ####
cat("\n===== RSF + Enet =====\n")
set.seed(MODEL_SEED)

# Extract variable importance from the RSF model and normalize
vi <- data.frame(imp = vimp.rfsrc(fit)$importance)
vi$imp <- (vi$imp - min(vi$imp)) / (max(vi$imp) - min(vi$imp))
vi$ID <- rownames(vi)

# Select genes with importance > 0.01
rid <- rownames(vi)[vi$imp > 0.01]
train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
trainlist2 <- lapply(trainlist, function(x) {
  x[, c('OS.time', 'OS', rid), drop = FALSE]
})

x1 <- as.matrix(train2[, rid, drop = FALSE])
x2 <- as.matrix(Surv(train2$OS.time, train2$OS))

# Loop over alpha values
for (alpha in seq(0, 1, 0.1)) {
  set.seed(MODEL_SEED)
  fit_cv <- cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
  rs_alpha <- lapply(trainlist2, function(x) {
    pred <- predict(fit_cv, type = 'link',
                    newx = as.matrix(x[, -c(1, 2), drop = FALSE]),
                    s = fit_cv$lambda.min)
    cbind(x[, 1:2], RS = as.numeric(pred))
  })
  rs_alpha <- lapply(rs_alpha, remove_inf)
  cc_alpha <- sapply(rs_alpha, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_alpha) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('RSF + Enet[α=', alpha, ']')
  result <- rbind(result, cc_df)
}

# Use cross-validation to select optimal lambda
set.seed(MODEL_SEED)
modelexp <- as.matrix(trainlist2$Train[, -c(1, 2), drop = FALSE])
modelstat <- Surv(trainlist2$Train$OS.time, trainlist2$Train$OS)
Enetmodel_cv <- cv.glmnet(modelexp, modelstat, family = 'cox', nfolds = 10)
model_fit <- glmnet(modelexp, modelstat, family = 'cox',
                    lambda = Enetmodel_cv$lambda.min)

rs_best <- lapply(trainlist2, function(x) {
  pred <- predict(model_fit, type = 'link',
                  newx = as.matrix(x[, -c(1, 2), drop = FALSE]),
                  s = model_fit$lambda)
  cbind(x[, 1:2], RS = as.numeric(pred))
})
rs_best <- lapply(rs_best, remove_inf)
cc_best <- sapply(rs_best, calculate_cindex)

cc_df <- data.frame(Cindex = cc_best) %>%
  rownames_to_column('ID')
cc_df$Model <- paste0('RSF + Enet[lambda=', round(Enetmodel_cv$lambda.min, 3), ']')
result <- rbind(result, cc_df)

#### 1-3. RSF + Stepwise Cox ####
cat("\n===== RSF + StepCox =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train2),
                   direction = direction)
  rs_step <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2], RS = predict(fit_step, type = 'risk', newdata = x))
  })
  rs_step <- lapply(rs_step, remove_inf)
  cc_step <- sapply(rs_step, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_step) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('RSF + StepCox[', direction, ']')
  result <- rbind(result, cc_df)
}

#### 1-4. RSF + CoxBoost ####
cat("\n===== RSF + CoxBoost =====\n")
set.seed(MODEL_SEED)

# Optimize penalty parameter
pen <- optimCoxBoostPenalty(
  time = train2$OS.time,
  status = train2$OS,
  x = as.matrix(train2[, -c(1, 2)]),
  trace = TRUE,
  start.penalty = 500,
  parallel = TRUE
)

# Cross-validation to select step number
cv_res <- cv.CoxBoost(
  time = train2$OS.time,
  status = train2$OS,
  x = as.matrix(train2[, -c(1, 2)]),
  maxstepno = 500,
  K = 10,
  type = "verweij",
  penalty = pen$penalty
)

# Fit final model
fit_cb <- CoxBoost(
  time = train2$OS.time,
  status = train2$OS,
  x = as.matrix(train2[, -c(1, 2)]),
  stepno = cv_res$optimal.step,
  penalty = pen$penalty
)

rs_cb <- lapply(trainlist2, function(x) {
  pred <- predict(fit_cb,
                  newdata = x[, -c(1, 2), drop = FALSE],
                  newtime = x[, 1],
                  newstatus = x[, 2],
                  type = "lp")
  cbind(x[, 1:2], RS = as.numeric(pred))
})
rs_cb <- lapply(rs_cb, remove_inf)
cc_cb <- sapply(rs_cb, calculate_cindex)

cc_df <- data.frame(Cindex = cc_cb) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + CoxBoost'
result <- rbind(result, cc_df)

#### 1-5. RSF + plsRcox ####
cat("\n===== RSF + plsRcox =====\n")
set.seed(MODEL_SEED)

# Safe wrapper for plsRcox training
safe_plsRcox <- function(train_data) {
  model_exp <- train_data[, -c(1, 2), drop = FALSE]
  model_time <- train_data$OS.time
  model_stat <- train_data$OS
  n_samples <- nrow(model_exp)
  n_events <- sum(model_stat)
  
  if (n_samples < 5 || n_events < 2) {
    warning("Insufficient samples or events for plsRcox")
    return(NULL)
  }
  
  # Cross-validation for number of components
  cv.model <- tryCatch({
    nt_max <- min(10, floor(n_samples / 2), n_events - 1)
    nfold <- min(10, max(3, floor(n_samples / 3)))
    cv.plsRcox(
      list(x = model_exp, time = model_time, status = model_stat),
      nt = nt_max,
      nfold = nfold,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("plsRcox CV failed: ", e$message)
    NULL
  })
  
  if (is.null(cv.model)) return(NULL)
  
  n_components <- tryCatch({
    comp <- cv.model$lambda.min5
    if (is.na(comp)) min(5, ncol(model_exp), n_samples - 1) else comp
  }, error = function(e) min(5, ncol(model_exp), n_samples - 1))
  
  n_components <- max(1, min(n_components, n_samples - 1, n_events - 1))
  
  tryCatch({
    plsRcox(
      Xplan = model_exp,
      time = model_time,
      event = model_stat,
      nt = n_components,
      alpha.pvals.expli = 0.05,
      sparse = TRUE,
      pvals.expli = TRUE
    )
  }, error = function(e) {
    warning("plsRcox fitting failed: ", e$message)
    NULL
  })
}

# Safe prediction for plsRcox
safe_predict_pls <- function(model, newdata) {
  if (is.null(model)) return(rep(NA, nrow(newdata)))
  tryCatch({
    as.numeric(predict(model, type = "lp", newdata = newdata[, -c(1, 2)]))
  }, error = function(e) {
    warning("plsRcox prediction failed: ", e$message)
    rep(NA, nrow(newdata))
  })
}

# Train model
model_pls <- safe_plsRcox(train2)
rs_pls <- lapply(trainlist2, function(x) {
  predictions <- safe_predict_pls(model_pls, x)
  cbind(x[, 1:2], RS = predictions)
})
rs_pls <- lapply(rs_pls, remove_inf)
cc_pls <- sapply(rs_pls, calculate_cindex)

cc_df <- data.frame(Cindex = cc_pls) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + plsRcox'
result <- rbind(result, cc_df)

#### 1-6. RSF + SuperPC ####
cat("\n===== RSF + SuperPC =====\n")
set.seed(MODEL_SEED)

# Safe wrapper for SuperPC training
safe_superpc <- function(train_data) {
  n_features <- ncol(train_data) - 2
  n_samples <- nrow(train_data)
  n_events <- sum(train_data$OS)
  
  if (n_samples < 5 || n_events < 2 || n_features < 5) {
    warning("Insufficient data for SuperPC")
    return(NULL)
  }
  
  data_list <- list(
    x = t(train_data[, -c(1, 2), drop = FALSE]),
    y = train_data$OS.time,
    censoring.status = train_data$OS,
    featurenames = colnames(train_data)[-c(1, 2)]
  )
  
  fit <- tryCatch({
    superpc.train(data = data_list, type = 'survival', s0.perc = 0.5)
  }, error = function(e) {
    warning("SuperPC training failed: ", e$message)
    NULL
  })
  if (is.null(fit)) return(NULL)
  
  cv.fit <- tryCatch({
    n_threshold <- min(20, n_features)
    n_fold <- min(10, floor(n_samples / 2))
    min_feat <- max(1, min(5, n_features))
    max_feat <- max(min_feat, min(n_features, n_samples))
    superpc.cv(
      fit, data_list,
      n.threshold = n_threshold,
      n.fold = n_fold,
      n.components = min(3, n_features),
      min.features = min_feat,
      max.features = max_feat,
      compute.fullcv = TRUE,
      compute.preval = TRUE
    )
  }, error = function(e) {
    warning("SuperPC CV failed: ", e$message)
    NULL
  })
  
  if (is.null(cv.fit)) return(NULL)
  
  list(fit = fit, cv.fit = cv.fit)
}

# Safe prediction for SuperPC
safe_predict_superpc <- function(superpc_result, newdata) {
  if (is.null(superpc_result)) return(rep(NA, nrow(newdata)))
  
  tryCatch({
    test_data <- list(
      x = t(newdata[, -c(1, 2), drop = FALSE]),
      y = newdata$OS.time,
      censoring.status = newdata$OS,
      featurenames = colnames(newdata)[-c(1, 2)]
    )
    
    threshold <- if (!is.null(superpc_result$cv.fit)) {
      scor <- superpc_result$cv.fit$scor
      if (length(scor) > 0 && nrow(scor) > 0) {
        superpc_result$cv.fit$thresholds[which.max(scor[1, ])]
      } else {
        median(superpc_result$cv.fit$thresholds, na.rm = TRUE)
      }
    } else {
      0.5
    }
    
    pred <- superpc.predict(
      superpc_result$fit,
      data = NULL,
      test = test_data,
      threshold = threshold,
      n.components = 1
    )
    as.numeric(pred$v.pred)
  }, error = function(e) {
    warning("SuperPC prediction failed: ", e$message)
    rep(NA, nrow(newdata))
  })
}

superpc_result <- safe_superpc(train2)
rs_sp <- lapply(trainlist2, function(x) {
  predictions <- safe_predict_superpc(superpc_result, x)
  cbind(x[, 1:2], RS = predictions)
})
rs_sp <- lapply(rs_sp, remove_inf)
cc_sp <- sapply(rs_sp, calculate_cindex)

cc_df <- data.frame(Cindex = cc_sp) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + SuperPC'
result <- rbind(result, cc_df)

#### 1-7. RSF + GBM ####
cat("\n===== RSF + GBM =====\n")
set.seed(MODEL_SEED)
fit_gbm <- gbm(
  Surv(OS.time, OS) ~ .,
  data = train2,
  distribution = 'coxph',
  n.minobsinnode = 10,
  n.cores = 1,
  n.trees = 1000,
  shrinkage = 0.005,
  interaction.depth = 2,
  cv.folds = 5
)
best_gbm <- which.min(fit_gbm$cv.error)

set.seed(MODEL_SEED)
fit_gbm <- gbm(
  Surv(OS.time, OS) ~ .,
  data = train2,
  distribution = 'coxph',
  n.trees = best_gbm,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.001,
  n.cores = 1
)

rs_gbm <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(fit_gbm, x, n.trees = best_gbm, type = 'link')))
})
rs_gbm <- lapply(rs_gbm, remove_inf)
cc_gbm <- sapply(rs_gbm, calculate_cindex)

cc_df <- data.frame(Cindex = cc_gbm) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + GBM'
result <- rbind(result, cc_df)

#### 1-8. RSF + survival-SVM ####
cat("\n===== RSF + survival-SVM =====\n")
set.seed(MODEL_SEED)
fit_svm <- survivalsvm(
  Surv(OS.time, OS) ~ .,
  data = train2,
  gamma.mu = 2
)
rs_svm <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(fit_svm, x)$predicted))
})
rs_svm <- lapply(rs_svm, remove_inf)
cc_svm <- sapply(rs_svm, calculate_cindex)

cc_df <- data.frame(Cindex = cc_svm) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + survival-SVM'
result <- rbind(result, cc_df)

#### 1-9. RSF + Ridge (logistic regression as proxy) ####
cat("\n===== RSF + Ridge =====\n")
set.seed(MODEL_SEED)
modelexp <- as.matrix(train2[, -c(1, 2), drop = FALSE])

for (alpha in seq(0, 1, 0.1)) {
  set.seed(MODEL_SEED)
  model_cv <- cv.glmnet(
    modelexp,
    train2$OS,
    family = 'binomial',
    alpha = alpha,
    nfolds = 10
  )
  fit_ridge <- glmnet(
    modelexp,
    train2$OS,
    family = 'binomial',
    alpha = alpha,
    lambda = model_cv$lambda.min
  )
  
  rs_ridge <- lapply(trainlist2, function(x) {
    pred <- predict(fit_ridge,
                    newx = as.matrix(x[, -c(1, 2), drop = FALSE]),
                    type = "response")
    cbind(x[, 1:2], RS = as.numeric(pred))
  })
  rs_ridge <- lapply(rs_ridge, remove_inf)
  cc_ridge <- sapply(rs_ridge, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_ridge) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('RSF + Ridge[α=', alpha, ']')
  result <- rbind(result, cc_df)
}

#### 1-10. RSF + obliqueRSF ####
cat("\n===== RSF + obliqueRSF =====\n")
set.seed(MODEL_SEED)
model_orsf <- orsf(
  data = train2,
  n_tree = 100,
  formula = Surv(OS.time, OS) ~ .
)
rs_orsf <- lapply(trainlist2, function(x) {
  pred <- predict(model_orsf, new_data = x, pred_type = "risk")[, 1]
  cbind(x[, 1:2], RS = as.numeric(pred))
})
rs_orsf <- lapply(rs_orsf, remove_inf)
cc_orsf <- sapply(rs_orsf, calculate_cindex)

cc_df <- data.frame(Cindex = cc_orsf) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + obliqueRSF'
result <- rbind(result, cc_df)

#### 1-11. RSF + xgboost ####
cat("\n===== RSF + xgboost =====\n")
set.seed(MODEL_SEED)

# Safe wrapper for XGBoost with feature alignment
safe_xgboost_analysis <- function(train_data, test_data_list) {
  if (nrow(train_data) < 10) {
    warning("Insufficient training samples for XGBoost")
    return(list(success = FALSE, results = NULL))
  }
  n_features <- ncol(train_data) - 2
  if (n_features < 1) {
    warning("No features for XGBoost")
    return(list(success = FALSE, results = NULL))
  }
  
  feature_names <- colnames(train_data)[-c(1, 2)]
  
  # Labels: positive time for events, negative time for censored
  labels_vector <- ifelse(train_data$OS == 1,
                          train_data$OS.time,
                          -train_data$OS.time)
  
  train_matrix <- tryCatch({
    xgb.DMatrix(
      data = as.matrix(train_data[, -c(1, 2)]),
      label = labels_vector,
      feature_names = feature_names
    )
  }, error = function(e) {
    warning("Failed to create xgb.DMatrix: ", e$message)
    NULL
  })
  if (is.null(train_matrix)) return(list(success = FALSE, results = NULL))
  
  model <- tryCatch({
    xgb.train(
      params = xgb_params_config,
      data = train_matrix,
      nrounds = n_rounds,
      verbose = 0
    )
  }, error = function(e) {
    warning("XGBoost training failed: ", e$message)
    NULL
  })
  if (is.null(model)) return(list(success = FALSE, results = NULL))
  
  results <- lapply(test_data_list, function(test_data) {
    # Align features
    missing_features <- setdiff(feature_names, colnames(test_data))
    if (length(missing_features) > 0) {
      warning("Test data missing features: ",
              paste(missing_features, collapse = ", "), "; filling with 0")
      for (f in missing_features) test_data[[f]] <- 0
    }
    test_features <- test_data[, feature_names, drop = FALSE]
    test_matrix <- tryCatch(as.matrix(test_features), error = function(e) NULL)
    if (is.null(test_matrix)) {
      return(cbind(test_data[, 1:2], RS = NA))
    }
    pred <- tryCatch(predict(model, newdata = test_matrix), error = function(e) {
      warning("XGBoost prediction failed: ", e$message)
      rep(NA, nrow(test_data))
    })
    cbind(test_data[, 1:2], RS = as.numeric(pred))
  })
  
  list(success = TRUE, results = results)
}

xgb_result <- safe_xgboost_analysis(train_data = train2,
                                    test_data_list = trainlist2)

if (xgb_result$success) {
  rs_xgb <- xgb_result$results
} else {
  rs_xgb <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2], RS = rep(NA, nrow(x)))
  })
}
rs_xgb <- lapply(rs_xgb, remove_inf)
cc_xgb <- sapply(rs_xgb, calculate_cindex)

cc_df <- data.frame(Cindex = cc_xgb) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + xgboost'
result <- rbind(result, cc_df)

#### 1-12. RSF + Conditional Forest (cforest) ####
cat("\n===== RSF + CForest =====\n")
set.seed(MODEL_SEED)
model_cforest <- party::cforest(
  Surv(OS.time, OS) ~ .,
  data = train2,
  controls = cforest_unbiased(ntree = 50)
)
rs_cforest <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(model_cforest, newdata = x, type = "response")))
})
rs_cforest <- lapply(rs_cforest, remove_inf)
cc_cforest <- sapply(rs_cforest, calculate_cindex)

cc_df <- data.frame(Cindex = cc_cforest) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + CForest'
result <- rbind(result, cc_df)

#### 1-13. RSF + Conditional Tree (ctree) ####
cat("\n===== RSF + CTree =====\n")
set.seed(MODEL_SEED)
model_ctree <- ctree(
  Surv(OS.time, OS) ~ .,
  data = train2
)
rs_ctree <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(model_ctree, newdata = x, type = "response")))
})
rs_ctree <- lapply(rs_ctree, remove_inf)
cc_ctree <- sapply(rs_ctree, calculate_cindex)

cc_df <- data.frame(Cindex = cc_ctree) %>%
  rownames_to_column('ID')
cc_df$Model <- 'RSF + CTree'
result <- rbind(result, cc_df)




# 7.2 Lasso family (using Lasso-selected features)####


#### 2-1. Lasso + Elastic Net ####
cat("\n===== Lasso + Enet =====\n")
set.seed(MODEL_SEED)
modelexp <- as.matrix(train[, -c(1, 2), drop = FALSE])
modelstat <- Surv(train$OS.time, train$OS)

for (alpha in seq(0, 1, 0.1)) {
  set.seed(MODEL_SEED)
  model_cv <- cv.glmnet(modelexp, modelstat, family = 'cox',
                        alpha = alpha, nfolds = 10)
  fit <- glmnet(modelexp, modelstat, family = 'cox',
                alpha = alpha, lambda = model_cv$lambda.min)
  
  rs_enet <- lapply(trainlist, function(x) {
    pred <- predict(fit, type = 'link',
                    newx = as.matrix(x[, -c(1, 2), drop = FALSE]),
                    s = fit$lambda)
    cbind(x[, 1:2], RS = as.numeric(pred))
  })
  rs_enet <- lapply(rs_enet, remove_inf)
  cc_enet <- sapply(rs_enet, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_enet) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('Lasso + Enet[α=', alpha, ']')
  result <- rbind(result, cc_df)
}

#### 2-2. Lasso + RSF ####
cat("\n===== Lasso + RSF =====\n")
set.seed(MODEL_SEED)

# Feature selection via Lasso
fit_lasso <- cv.glmnet(modelexp, modelstat, family = "cox")
coef.min <- coef(fit_lasso, s = "lambda.min")
rid <- coef.min@Dimnames[[1]][coef.min@i + 1]
if (length(rid) == 0) {
  warning("Lasso selected no features; using all features as fallback")
  rid <- colnames(train)[-c(1, 2)]
} else {
  rid <- rid[rid != "(Intercept)"]
}
if (length(rid) == 0) stop("No features available for RSF")

train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
trainlist2 <- lapply(trainlist, function(x) {
  x[, c('OS.time', 'OS', rid), drop = FALSE]
})

# RSF on selected features
set.seed(MODEL_SEED)
fit_rsf <- rfsrc(
  Surv(OS.time, OS) ~ .,
  data = train2,
  ntree = 1000,
  nodesize = rf_nodesize,
  splitrule = 'logrank',
  importance = TRUE,
  proximity = TRUE,
  forest = TRUE,
  MODEL_SEED = MODEL_SEED
)

best <- which.min(fit_rsf$err.rate)
set.seed(MODEL_SEED)
fit_rsf <- rfsrc(
  Surv(OS.time, OS) ~ .,
  data = train2,
  ntree = best,
  nodesize = rf_nodesize,
  splitrule = 'logrank',
  importance = TRUE,
  proximity = TRUE,
  forest = TRUE,
  MODEL_SEED = MODEL_SEED
)

rs_lasso_rsf <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2], RS = predict(fit_rsf, newdata = x)$predicted)
})
rs_lasso_rsf <- lapply(rs_lasso_rsf, remove_inf)
cc_lasso_rsf <- sapply(rs_lasso_rsf, calculate_cindex)

cc_df <- data.frame(Cindex = cc_lasso_rsf) %>%
  rownames_to_column('ID')
cc_df$Model <- 'Lasso + RSF'
result <- rbind(result, cc_df)

#### 2-3. Lasso + StepCox ####
cat("\n===== Lasso + StepCox =====\n")
for (direction in c("both", "backward", "forward")) {
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train2),
                   direction = direction)
  rs_step <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2], RS = predict(fit_step, type = 'risk', newdata = x))
  })
  rs_step <- lapply(rs_step, remove_inf)
  cc_step <- sapply(rs_step, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_step) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('Lasso + StepCox[', direction, ']')
  result <- rbind(result, cc_df)
}

#### 2-4. Lasso + CoxBoost ####
cat("\n===== Lasso + CoxBoost =====\n")
set.seed(MODEL_SEED)

if (nrow(train2) < 10 || ncol(train2) <= 2) {
  warning("Insufficient data for CoxBoost")
  cc_df <- create_na_results('Lasso + CoxBoost', names(trainlist))
  result <- rbind(result, cc_df)
} else {
  train_features <- as.matrix(train2[, -c(1, 2), drop = FALSE])
  
  # Optimize penalty
  pen <- optimCoxBoostPenalty(
    time = train2$OS.time,
    status = train2$OS,
    x = train_features,
    trace = TRUE,
    parallel = TRUE
  )
  
  # Cross‑validation
  cv_res <- cv.CoxBoost(
    time = train2$OS.time,
    status = train2$OS,
    x = train_features,
    maxstepno = 100,
    K = 3,
    type = "verweij",
    penalty = pen$penalty
  )
  
  # Final model
  fit_cb <- CoxBoost(
    time = train2$OS.time,
    status = train2$OS,
    x = train_features,
    stepno = cv_res$optimal.step,
    penalty = pen$penalty
  )
  
  safe_predict_cb <- function(x, fit, train_features) {
    if (is.null(x) || nrow(x) == 0) return(rep(NA, nrow(x)))
    if (ncol(x) <= 2) return(rep(NA, nrow(x)))
    tryCatch({
      pred_data <- as.matrix(x[, -c(1, 2), drop = FALSE])
      # Align features (if names differ)
      if (!is.null(colnames(pred_data))) {
        common <- intersect(colnames(pred_data), colnames(train_features))
        if (length(common) == 0) {
          warning("No common features for prediction")
          return(rep(NA, nrow(x)))
        }
        pred_data <- pred_data[, common, drop = FALSE]
      }
      predict(fit, newdata = pred_data, type = "lp")
    }, error = function(e) {
      warning("CoxBoost prediction failed: ", e$message)
      rep(NA, nrow(x))
    })
  }
  
  rs_cb <- lapply(trainlist2, function(x) {
    pred <- safe_predict_cb(x, fit_cb, train_features)
    cbind(x[, 1:2], RS = as.numeric(pred))
  })
  rs_cb <- lapply(rs_cb, remove_inf)
  cc_cb <- sapply(rs_cb, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_cb) %>%
    rownames_to_column('ID')
  cc_df$Model <- 'Lasso + CoxBoost'
  result <- rbind(result, cc_df)
}

#### 2-5. Lasso + plsRcox ####
cat("\n===== Lasso + plsRcox =====\n")
set.seed(MODEL_SEED)

if (nrow(train2) < 10 || ncol(train2) <= 2) {
  warning("Insufficient data for plsRcox")
  cc_df <- create_na_results('Lasso + plsRcox', names(trainlist))
  result <- rbind(result, cc_df)
} else {
  model_exp <- train2[, -c(1, 2), drop = FALSE]
  model_time <- train2$OS.time
  model_stat <- train2$OS
  
  # Remove columns with >50% NA
  na_ratio <- colMeans(is.na(model_exp))
  keep_cols <- na_ratio <= 0.5
  model_exp <- model_exp[, keep_cols, drop = FALSE]
  
  if (ncol(model_exp) == 0) {
    warning("No valid features after NA removal")
    cc_df <- create_na_results('Lasso + plsRcox', names(trainlist))
    result <- rbind(result, cc_df)
  } else {
    # Cross‑validation for number of components
    cv_pls <- tryCatch({
      cv.plsRcox(
        list(x = model_exp, time = model_time, status = model_stat),
        nt = min(5, ncol(model_exp)),
        verbose = FALSE
      )
    }, error = function(e) {
      warning("plsRcox CV failed: ", e$message)
      NULL
    })
    
    optimal_nt <- if (!is.null(cv_pls) && !is.null(cv_pls$lambda.min5)) {
      min(cv_pls$lambda.min5, ncol(model_exp))
    } else {
      min(3, ncol(model_exp))
    }
    
    # Fit final model
    fit_pls <- tryCatch({
      plsRcox(
        Xplan = model_exp,
        time = model_time,
        event = model_stat,
        nt = optimal_nt,
        alpha.pvals.expli = 0.05,
        sparse = TRUE
      )
    }, error = function(e) {
      warning("plsRcox fitting failed, trying non‑sparse: ", e$message)
      tryCatch({
        plsRcox(
          Xplan = model_exp,
          time = model_time,
          event = model_stat,
          nt = optimal_nt,
          sparse = FALSE
        )
      }, error = function(e2) {
        warning("Non‑sparse also failed: ", e2$message)
        NULL
      })
    })
    
    if (is.null(fit_pls)) {
      cc_df <- create_na_results('Lasso + plsRcox', names(trainlist))
      result <- rbind(result, cc_df)
    } else {
      safe_predict_pls2 <- function(x, fit, ref_exp) {
        if (is.null(x) || nrow(x) == 0) return(rep(NA, nrow(x)))
        tryCatch({
          pred_data <- x[, -c(1, 2), drop = FALSE]
          # Align features
          common <- intersect(colnames(pred_data), colnames(ref_exp))
          if (length(common) == 0) {
            warning("No common features")
            return(rep(NA, nrow(x)))
          }
          pred_data <- pred_data[, common, drop = FALSE]
          # Add missing features with mean values
          missing <- setdiff(colnames(ref_exp), common)
          for (f in missing) {
            pred_data[[f]] <- mean(ref_exp[[f]], na.rm = TRUE)
          }
          pred_data <- pred_data[, colnames(ref_exp), drop = FALSE]
          # Impute NAs
          for (j in seq_len(ncol(pred_data))) {
            if (any(is.na(pred_data[, j]))) {
              col_mean <- mean(ref_exp[, j], na.rm = TRUE)
              pred_data[is.na(pred_data[, j]), j] <- col_mean
            }
          }
          predict(fit, type = "lp", newdata = pred_data)
        }, error = function(e) {
          warning("plsRcox prediction failed: ", e$message)
          rep(NA, nrow(x))
        })
      }
      
      rs_pls <- lapply(trainlist2, function(x) {
        pred <- safe_predict_pls2(x, fit_pls, model_exp)
        cbind(x[, 1:2], RS = as.numeric(pred))
      })
      rs_pls <- lapply(rs_pls, remove_inf)
      cc_pls <- sapply(rs_pls, calculate_cindex)
      
      cc_df <- data.frame(Cindex = cc_pls) %>%
        rownames_to_column('ID')
      cc_df$Model <- 'Lasso + plsRcox'
      result <- rbind(result, cc_df)
    }
  }
}

#### 2-6. Lasso + SuperPC ####
cat("\n===== Lasso + SuperPC =====\n")
set.seed(MODEL_SEED)

if (nrow(train2) < 10 || ncol(train2) <= 2) {
  warning("Insufficient data for SuperPC")
  cc_df <- create_na_results('Lasso + SuperPC', names(trainlist))
  result <- rbind(result, cc_df)
} else {
  safe_superpc_analysis <- function(train_data, test_data_list) {
    n_features <- ncol(train_data) - 2
    if (n_features < 1) return(list(success = FALSE))
    
    data_train <- list(
      x = t(train_data[, -c(1, 2), drop = FALSE]),
      y = train_data$OS.time,
      censoring.status = train_data$OS,
      featurenames = colnames(train_data)[-c(1, 2)]
    )
    
    fit <- tryCatch(
      superpc.train(data = data_train, type = 'survival', s0.perc = 0.5),
      error = function(e) NULL
    )
    if (is.null(fit)) return(list(success = FALSE))
    
    cv_fit <- tryCatch(
      superpc.cv(
        fit, data_train,
        n.threshold = min(20, n_features),
        n.fold = min(10, nrow(train_data)),
        n.components = min(3, n_features),
        min.features = 1,
        max.features = n_features,
        compute.fullcv = TRUE,
        compute.preval = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) return(list(success = FALSE))
    
    best_threshold <- tryCatch({
      cv_fit$thresholds[which.max(cv_fit$scor[1, ])]
    }, error = function(e) median(data_train$x, na.rm = TRUE))
    
    results <- lapply(test_data_list, function(w) {
      data_test <- list(
        x = t(w[, -c(1, 2), drop = FALSE]),
        y = w$OS.time,
        censoring.status = w$OS,
        featurenames = colnames(w)[-c(1, 2)]
      )
      pred <- tryCatch(
        superpc.predict(fit, data_train, data_test,
                        threshold = best_threshold, n.components = 1)$v.pred,
        error = function(e) rep(NA, nrow(w))
      )
      cbind(w[, 1:2], RS = as.numeric(pred))
    })
    
    list(success = TRUE, results = results)
  }
  
  sp_res <- safe_superpc_analysis(train2, trainlist2)
  if (sp_res$success) {
    rs_sp <- sp_res$results
  } else {
    rs_sp <- lapply(trainlist2, function(x) cbind(x[, 1:2], RS = NA))
  }
  rs_sp <- lapply(rs_sp, remove_inf)
  cc_sp <- sapply(rs_sp, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_sp) %>%
    rownames_to_column('ID')
  cc_df$Model <- 'Lasso + SuperPC'
  result <- rbind(result, cc_df)
}

#### 2-7. Lasso + GBM ####
cat("\n===== Lasso + GBM =====\n")
set.seed(MODEL_SEED)
fit_gbm <- gbm(
  Surv(OS.time, OS) ~ .,
  data = train2,
  distribution = 'coxph',
  n.trees = 1000,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.001,
  cv.folds = 5,
  n.cores = 1
)
best_gbm <- which.min(fit_gbm$cv.error)
set.seed(MODEL_SEED)
fit_gbm <- gbm(
  Surv(OS.time, OS) ~ .,
  data = train2,
  distribution = 'coxph',
  n.trees = best_gbm,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.001,
  n.cores = 1
)

rs_gbm <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(fit_gbm, x, n.trees = best_gbm, type = 'link')))
})
rs_gbm <- lapply(rs_gbm, remove_inf)
cc_gbm <- sapply(rs_gbm, calculate_cindex)

cc_df <- data.frame(Cindex = cc_gbm) %>%
  rownames_to_column('ID')
cc_df$Model <- 'Lasso + GBM'
result <- rbind(result, cc_df)

#### 2-8. Lasso + survival-SVM ####
cat("\n===== Lasso + survival-SVM =====\n")
set.seed(MODEL_SEED)

if (nrow(train2) < 10 || ncol(train2) <= 2) {
  warning("Insufficient data for survival-SVM")
  cc_df <- create_na_results('Lasso + survival-SVM', names(trainlist))
  result <- rbind(result, cc_df)
} else {
  train2_clean <- train2[complete.cases(train2), ]
  if (nrow(train2_clean) < 10) {
    warning("Too many NAs in training data")
    cc_df <- create_na_results('Lasso + survival-SVM', names(trainlist))
    result <- rbind(result, cc_df)
  } else {
    fit_svm <- tryCatch(
      survivalsvm(Surv(OS.time, OS) ~ ., data = train2_clean, gamma.mu = 2),
      error = function(e) {
        tryCatch(
          survivalsvm(Surv(OS.time, OS) ~ ., data = train2_clean, gamma.mu = 1),
          error = function(e2) NULL
        )
      }
    )
    
    if (is.null(fit_svm)) {
      cc_df <- create_na_results('Lasso + survival-SVM', names(trainlist))
      result <- rbind(result, cc_df)
    } else {
      safe_predict_svm <- function(x, fit, train_data) {
        if (is.null(x) || nrow(x) == 0) return(rep(NA, nrow(x)))
        tryCatch({
          # Align factors etc.
          pred <- predict(fit, x)$predicted
          as.numeric(pred)
        }, error = function(e) {
          warning("SVM prediction failed: ", e$message)
          rep(NA, nrow(x))
        })
      }
      
      rs_svm <- lapply(trainlist2, function(x) {
        pred <- safe_predict_svm(x, fit_svm, train2_clean)
        cbind(x[, 1:2], RS = pred)
      })
      rs_svm <- lapply(rs_svm, remove_inf)
      cc_svm <- sapply(rs_svm, calculate_cindex)
      
      cc_df <- data.frame(Cindex = cc_svm) %>%
        rownames_to_column('ID')
      cc_df$Model <- 'Lasso + survival-SVM'
      result <- rbind(result, cc_df)
    }
  }
}

#### 2-9. Lasso + Ridge (logistic regression)####
cat("\n===== Lasso + Ridge =====\n")
set.seed(MODEL_SEED)

#' Safely train Elastic Net model for a given alpha and compute C-index
#'
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = lasso)
#' @param modelexp Training feature matrix
#' @param train2 Training data with time, status, and features
#' @param trainlist2 List of test datasets (each with time, status, features)
#' @param seed Random seed for reproducibility
#' @return List with C-index vector, alpha, and success flag
safe_elastic_net <- function(alpha, modelexp, train2, trainlist2, seed) {
  cat("  Testing alpha =", alpha, "\n")
  
  tryCatch({
    set.seed(MODEL_SEED)
    
    # Cross-validation to select lambda.min
    model_cv <- cv.glmnet(
      modelexp,
      train2$OS,
      family = "binomial",
      alpha = alpha,
      nfolds = 10,
      parallel = FALSE
    )
    
    # Final model with lambda.min
    fit <- glmnet(
      modelexp,
      train2$OS,
      family = "binomial",
      alpha = alpha,
      lambda = model_cv$lambda.min
    )
    
    # Safe prediction function for a single dataset
    safe_predict_glmnet <- function(x, fit, ref_exp) {
      if (is.null(x) || nrow(x) == 0) {
        warning("Empty dataset; returning NA")
        return(rep(NA, nrow(x)))
      }
      if (ncol(x) <= 2) {
        warning("Insufficient columns; returning NA")
        return(rep(NA, nrow(x)))
      }
      
      tryCatch({
        # Extract features (exclude time and status)
        pred_data <- as.matrix(x[, -c(1, 2), drop = FALSE])
        
        # Align features with training data
        if (ncol(pred_data) != ncol(ref_exp)) {
          warning("Feature dimension mismatch; attempting alignment")
          common_features <- intersect(colnames(pred_data), colnames(ref_exp))
          if (length(common_features) == 0) {
            warning("No common features; returning NA")
            return(rep(NA, nrow(x)))
          }
          pred_data <- pred_data[, common_features, drop = FALSE]
          
          # Add missing features with mean values from training
          if (ncol(pred_data) < ncol(ref_exp)) {
            warning("Adding missing features with training means")
            missing_features <- setdiff(colnames(ref_exp), common_features)
            for (feat in missing_features) {
              pred_data <- cbind(
                pred_data,
                rep(mean(ref_exp[, feat], na.rm = TRUE), nrow(pred_data))
              )
            }
            colnames(pred_data) <- c(common_features, missing_features)
            pred_data <- pred_data[, colnames(ref_exp), drop = FALSE]
          }
        }
        
        # Impute remaining NAs with column means
        if (any(is.na(pred_data))) {
          warning("Imputing NA values with training column means")
          for (j in seq_len(ncol(pred_data))) {
            na_idx <- is.na(pred_data[, j])
            if (any(na_idx)) {
              col_mean <- mean(ref_exp[, j], na.rm = TRUE)
              pred_data[na_idx, j] <- col_mean
            }
          }
        }
        
        # Predict risk scores (type = "response" gives probabilities)
        as.numeric(predict(fit, type = "response", newx = pred_data))
        
      }, error = function(e) {
        warning("Prediction failed: ", e$message)
        rep(NA, nrow(x))
      })
    }
    
    # Generate risk scores for all cohorts
    rs <- lapply(names(trainlist2), function(cohort_name) {
      x <- trainlist2[[cohort_name]]
      pred <- safe_predict_glmnet(x, fit, modelexp)
      data.frame(
        OS.time = x$OS.time,
        OS = x$OS,
        RS = pred,
        row.names = rownames(x)
      )
    })
    names(rs) <- names(trainlist2)
    
    # Remove infinite values (helper defined globally)
    rs <- lapply(rs, remove_inf)
    
    # Calculate C-index for each cohort
    cc <- sapply(rs, calculate_cindex)
    
    list(cc = cc, alpha = alpha, success = TRUE)
    
  }, error = function(e) {
    warning("Elastic Net failed for alpha = ", alpha, ": ", e$message)
    cc_na <- rep(NA, length(trainlist2))
    names(cc_na) <- names(trainlist2)
    list(cc = cc_na, alpha = alpha, success = FALSE)
  })
}

# ===== Main execution for Lasso + Ridge
tryCatch({
  # Data validation
  if (nrow(train2) < 10) stop("Insufficient samples (<10)")
  if (ncol(train2) <= 3) stop("Need at least 1 feature (time, status, and ≥1 feature)")
  
  # Extract feature matrix
  modelexp <- as.matrix(train2[, -c(1, 2), drop = FALSE])
  
  cat("Data validation:\n")
  cat("  Samples:", nrow(modelexp), "\n")
  cat("  Features:", ncol(modelexp), "\n")
  cat("  Event rate:", round(mean(train2$OS), 3), "\n")
  
  if (ncol(modelexp) < 2) stop("Need at least 2 features for Elastic Net")
  
  # Impute missing values in training features
  if (any(is.na(modelexp))) {
    warning("NA values detected; imputing with column means")
    modelexp <- apply(modelexp, 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    })
  }
  
  # Define alpha values to test
  alpha_values <- seq(0, 1, 0.1)
  results_list <- list()
  
  # Loop over alpha values
  for (alpha in alpha_values) {
    res <- safe_elastic_net(alpha, modelexp, train2, trainlist2, MODEL_SEED)
    results_list[[as.character(alpha)]] <- res
  }
  
  # Compile results into data frame and merge with global 'result'
  for (alpha in alpha_values) {
    alpha_str <- as.character(alpha)
    res_alpha <- results_list[[alpha_str]]
    
    cc_df <- data.frame(
      ID = names(res_alpha$cc),
      Cindex = res_alpha$cc,
      Model = paste0("Lasso + Ridge[α=", alpha, "]"),
      stringsAsFactors = FALSE
    )
    
    # Merge (global 'result' should exist)
    if (exists("result")) {
      result <- tryCatch(
        rbind(result, cc_df),
        error = function(e) {
          warning("Failed to merge for alpha = ", alpha, ": ", e$message)
          result
        }
      )
    } else {
      result <- cc_df
    }
  }
  
  # Summary
  cat("Elastic Net analysis completed.\n")
  successful <- sapply(results_list, function(x) x$success)
  cat("Successful alpha values:", alpha_values[successful], "\n")
  cat("Failed alpha values:    ", alpha_values[!successful], "\n")
  
  if (any(successful)) {
    all_cc <- unlist(lapply(results_list[successful], function(x) x$cc))
    all_cc <- all_cc[!is.na(all_cc)]
    if (length(all_cc) > 0) {
      cat("C-index summary (successful models):\n")
      print(summary(all_cc))
      cat("Number of valid C-index values:", length(all_cc), "\n")
    } else {
      cat("No valid C-index values obtained.\n")
    }
  } else {
    cat("All alpha values failed.\n")
  }
  
}, error = function(e) {
  # Complete failure: assign NA for all alphas
  warning("Lasso + Ridge analysis completely failed: ", e$message)
  cat("Assigning NA for all cohorts and alphas.\n")
  
  alpha_values <- seq(0, 1, 0.1)
  for (alpha in alpha_values) {
    cc_na <- rep(NA, length(trainlist2))
    names(cc_na) <- names(trainlist2)
    cc_df <- data.frame(
      ID = names(cc_na),
      Cindex = cc_na,
      Model = paste0("Lasso + Ridge[α=", alpha, "]"),
      stringsAsFactors = FALSE
    )
    if (exists("result")) {
      result <- tryCatch(
        rbind(result, cc_df),
        error = function(e2) {
          warning("Failed to merge NA results for alpha = ", alpha, ": ", e2$message)
          result
        }
      )
    } else {
      result <- cc_df
    }
  }
})

#### 2-10. Lasso + obliqueRSF ####
cat("\n===== Lasso + obliqueRSF =====\n")
set.seed(MODEL_SEED)
model_orsf <- orsf(
  data = train2,
  n_tree = 100,
  formula = Surv(OS.time, OS) ~ .
)
rs_orsf <- lapply(trainlist2, function(x) {
  pred <- predict(model_orsf, new_data = x, pred_type = "risk")[, 1]
  cbind(x[, 1:2], RS = as.numeric(pred))
})
rs_orsf <- lapply(rs_orsf, remove_inf)
cc_orsf <- sapply(rs_orsf, calculate_cindex)

cc_df <- data.frame(Cindex = cc_orsf) %>%
  rownames_to_column('ID')
cc_df$Model <- 'Lasso + obliqueRSF'
result <- rbind(result, cc_df)

#### 2-11. Lasso + xgboost ####
cat("\n===== Lasso + xgboost =====\n")
set.seed(MODEL_SEED)
labels_vector <- ifelse(train$OS == 1, train$OS.time, -train$OS.time)
model_mat <- xgb.DMatrix(
  data = as.matrix(train[, -c(1, 2)]),
  label = labels_vector
)
model_xgb <- xgb.train(
  params = xgb_params_config,
  data = model_mat,
  nrounds = n_rounds,
  verbose = 0
)
rs_xgb <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(model_xgb, newdata = as.matrix(x[, -c(1, 2)]))))
})
rs_xgb <- lapply(rs_xgb, remove_inf)
cc_xgb <- sapply(rs_xgb, calculate_cindex)

cc_df <- data.frame(Cindex = cc_xgb) %>%
  rownames_to_column('ID')
cc_df$Model <- 'Lasso + xgboost'
result <- rbind(result, cc_df)

#### 2-12. Lasso + CForest ####
cat("\n===== Lasso + CForest =====\n")
set.seed(MODEL_SEED)
model_cforest <- party::cforest(
  Surv(OS.time, OS) ~ .,
  data = train2,
  controls = cforest_unbiased(ntree = 50)
)
rs_cforest <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(model_cforest, newdata = x, type = "response")))
})
rs_cforest <- lapply(rs_cforest, remove_inf)
cc_cforest <- sapply(rs_cforest, calculate_cindex)

cc_df <- data.frame(Cindex = cc_cforest) %>%
  rownames_to_column('ID')
cc_df$Model <- 'Lasso + CForest'
result <- rbind(result, cc_df)

#### 2-13. Lasso + CTree ####
cat("\n===== Lasso + CTree =====\n")
set.seed(MODEL_SEED)
model_ctree <- ctree(
  Surv(OS.time, OS) ~ .,
  data = train2
)
rs_ctree <- lapply(trainlist2, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(model_ctree, newdata = x, type = "response")))
})
rs_ctree <- lapply(rs_ctree, remove_inf)
cc_ctree <- sapply(rs_ctree, calculate_cindex)

cc_df <- data.frame(Cindex = cc_ctree) %>%
  rownames_to_column('ID')
cc_df$Model <- 'Lasso + CTree'
result <- rbind(result, cc_df)


# 8. StepCox Family (Stepwise Cox Regression)####


#### 3-1. StepCox alone ####
cat("\n===== StepCox alone =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rs_step <- lapply(trainlist, function(x) {
    cbind(x[, 1:2], RS = predict(fit_step, type = 'risk', newdata = x))
  })
  rs_step <- lapply(rs_step, remove_inf)
  cc_step <- sapply(rs_step, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_step) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, ']')
  result <- rbind(result, cc_df)
}

#### 3-2. StepCox + RSF ####
cat("\n===== StepCox + RSF =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 1) {
    warning(sprintf("Skipped direction '%s': no features selected", direction))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  set.seed(MODEL_SEED)
  fit_rsf <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = train2,
    ntree = 1000,
    nodesize = rf_nodesize,
    splitrule = 'logrank',
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE,
    MODEL_SEED = MODEL_SEED
  )
  best <- which.min(fit_rsf$err.rate)
  set.seed(MODEL_SEED)
  fit_rsf <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = train2,
    ntree = best,
    nodesize = rf_nodesize,
    splitrule = 'logrank',
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE,
    MODEL_SEED = MODEL_SEED
  )
  
  rs <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2], RS = predict(fit_rsf, newdata = x)$predicted)
  })
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, '] + RSF')
  result <- rbind(result, cc_df)
}

#### 3-3. StepCox + Elastic Net ####
cat("\n===== StepCox + Enet =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 2) {
    warning(sprintf("Skipped direction '%s': %d features selected (<2)", direction, length(rid)))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  x1 <- as.matrix(train2[, rid])
  x2 <- Surv(train2$OS.time, train2$OS)
  
  for (alpha in seq(0, 1, 0.1)) {
    if (ncol(x1) < 2) next
    set.seed(MODEL_SEED)
    fit_cv <- tryCatch(
      cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10),
      error = function(e) NULL
    )
    if (is.null(fit_cv)) next
    
    rs <- lapply(trainlist2, function(x) {
      newx <- as.matrix(x[, rid])
      pred <- predict(fit_cv, type = 'link', newx = newx, s = fit_cv$lambda.min)
      cbind(x[, 1:2], RS = as.numeric(pred))
    })
    rs <- lapply(rs, remove_inf)
    cc <- sapply(rs, calculate_cindex)
    
    cc_df <- data.frame(Cindex = cc) %>%
      rownames_to_column('ID')
    cc_df$Model <- sprintf('StepCox[%s] + Enet[α=%.1f]', direction, alpha)
    result <- rbind(result, cc_df)
  }
}

#### 3-4. StepCox + CoxBoost ####
cat("\n===== StepCox + CoxBoost =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 2) {
    warning(sprintf("Skipped direction '%s': %d features selected (<2)", direction, length(rid)))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  x1 <- as.matrix(train2[, rid])
  x2 <- Surv(train2$OS.time, train2$OS)
  
  for (alpha in seq(0, 1, 0.1)) {
    set.seed(MODEL_SEED)
    fit_cv <- tryCatch(
      cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10),
      error = function(e) NULL
    )
    if (is.null(fit_cv)) next
    
    rs <- lapply(trainlist2, function(x) {
      newx <- as.matrix(x[, rid])
      pred <- predict(fit_cv, type = 'link', newx = newx, s = "lambda.min")
      cbind(x[, 1:2], RS = as.numeric(pred))
    })
    rs <- lapply(rs, remove_inf)
    cc <- sapply(rs, calculate_cindex)
    
    cc_df <- data.frame(Cindex = cc) %>%
      rownames_to_column('ID')
    cc_df$Model <- sprintf('StepCox[%s] + CoxBoost(α=%.1f)', direction, alpha)
    result <- rbind(result, cc_df)
  }
}

#### 3-5. StepCox + plsRcox ####
cat("\n===== StepCox + plsRcox =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- tryCatch(
    step(coxph(Surv(OS.time, OS) ~ ., data = train), direction = direction),
    error = function(e) NULL
  )
  if (is.null(fit_step)) next
  rid <- names(coef(fit_step))
  if (length(rid) < 2) {
    warning(sprintf("Skipped direction '%s': %d features selected (<2)", direction, length(rid)))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  pls_data <- list(
    x = as.matrix(train2[, rid]),
    time = train2$OS.time,
    status = train2$OS
  )
  
  set.seed(MODEL_SEED)
  cv_pls <- tryCatch(
    cv.plsRcox(pls_data, nt = 10, nfold = 10, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(cv_pls)) next
  
  best_nt <- as.numeric(cv_pls[5])
  fit_pls <- tryCatch(
    plsRcox(Xplan = pls_data$x, time = pls_data$time,
            event = pls_data$status, nt = best_nt),
    error = function(e) NULL
  )
  if (is.null(fit_pls)) next
  
  rs <- lapply(trainlist2, function(x) {
    newdata <- as.matrix(x[, rid])
    pred <- predict(fit_pls, type = "lp", newdata = newdata)
    cbind(x[, 1:2], RS = as.numeric(pred))
  })
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, '] + plsRcox')
  result <- rbind(result, cc_df)
}

#### 3-6. StepCox + SuperPC ####
cat("\n===== StepCox + SuperPC =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 1) {
    warning(sprintf("Skipped direction '%s': no features selected", direction))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  prepare_sp_data <- function(df, ref_features = NULL) {
    features <- df[, -c(1, 2), drop = FALSE]
    if (!is.null(ref_features)) {
      missing <- setdiff(ref_features, colnames(features))
      for (f in missing) features[[f]] <- 0
      features <- features[, ref_features, drop = FALSE]
    }
    x_mat <- t(as.matrix(features))
    if (is.null(rownames(x_mat))) rownames(x_mat) <- colnames(features)
    list(x = x_mat, y = df$OS.time, censoring.status = df$OS,
         featurenames = rownames(x_mat))
  }
  
  data_train <- tryCatch(
    prepare_sp_data(train2),
    error = function(e) NULL
  )
  if (is.null(data_train)) next
  train_features <- rownames(data_train$x)
  
  set.seed(MODEL_SEED)
  fit_sp <- tryCatch(
    superpc.train(data = data_train, type = 'survival', s0.perc = 0.5),
    error = function(e) NULL
  )
  if (is.null(fit_sp)) next
  
  cv_sp <- tryCatch(
    superpc.cv(
      fit_sp, data_train,
      n.threshold = min(20, length(rid)),
      n.fold = min(5, nrow(train2)),
      n.components = min(3, length(rid)),
      min.features = 1,
      max.features = length(rid),
      compute.fullcv = TRUE,
      compute.preval = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(cv_sp)) next
  
  best_threshold <- tryCatch(
    cv_sp$thresholds[which.max(cv_sp$scor[1, ])],
    error = function(e) 0
  )
  
  rs <- lapply(trainlist2, function(w) {
    data_test <- tryCatch(
      prepare_sp_data(w, ref_features = train_features),
      error = function(e) NULL
    )
    if (is.null(data_test)) return(cbind(w[, 1:2], RS = NA))
    pred <- tryCatch(
      superpc.predict(fit_sp, data_train, data_test,
                      threshold = best_threshold, n.components = 1)$v.pred,
      error = function(e) rep(NA, nrow(w))
    )
    cbind(w[, 1:2], RS = as.numeric(pred))
  })
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, '] + SuperPC')
  result <- rbind(result, cc_df)
}

#### 3-7. StepCox + GBM ####
cat("\n===== StepCox + GBM =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 1) {
    warning(sprintf("Skipped direction '%s': no features selected", direction))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  set.seed(MODEL_SEED)
  fit_gbm <- gbm(
    Surv(OS.time, OS) ~ .,
    data = train2,
    distribution = 'coxph',
    n.trees = 1000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 5,
    n.cores = 1
  )
  best_gbm <- which.min(fit_gbm$cv.error)
  set.seed(MODEL_SEED)
  fit_gbm <- gbm(
    Surv(OS.time, OS) ~ .,
    data = train2,
    distribution = 'coxph',
    n.trees = best_gbm,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    n.cores = 1
  )
  
  rs <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2],
          RS = as.numeric(predict(fit_gbm, x, n.trees = best_gbm, type = 'link')))
  })
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, '] + GBM')
  result <- rbind(result, cc_df)
}

#### 3-8. StepCox + survival-SVM ####
cat("\n===== StepCox + survival-SVM =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- tryCatch(
    step(coxph(Surv(OS.time, OS) ~ ., data = train), direction = direction, trace = 0),
    error = function(e) NULL
  )
  if (is.null(fit_step)) next
  rid <- setdiff(names(coef(fit_step)), "(Intercept)")
  if (length(rid) == 0) next
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  train2_clean <- train2[complete.cases(train2), ]
  if (nrow(train2_clean) < 10) next
  trainlist2 <- lapply(trainlist, function(x) {
    safe_subset(x, c('OS.time', 'OS', rid))
  })
  if (any(sapply(trainlist2, is.null))) next
  
  # Try several gamma values
  for (gm in c(1, 2, 0.5)) {
    fit_svm <- tryCatch(
      survivalsvm(Surv(OS.time, OS) ~ ., data = train2_clean, gamma.mu = gm),
      error = function(e) NULL
    )
    if (is.null(fit_svm)) next
    
    safe_predict_svm2 <- function(fit, newdata) {
      if (is.null(fit)) return(rep(NA, nrow(newdata)))
      tryCatch(
        as.numeric(predict(fit, newdata)$predicted),
        error = function(e) rep(NA, nrow(newdata))
      )
    }
    
    rs <- lapply(trainlist2, function(x) {
      if (is.null(x)) return(data.frame(OS.time = NA, OS = NA, RS = NA))
      pred <- safe_predict_svm2(fit_svm, x)
      data.frame(OS.time = x$OS.time, OS = x$OS, RS = pred)
    })
    rs <- lapply(rs, remove_inf)
    cc <- sapply(rs, calculate_cindex)
    
    cc_df <- data.frame(Cindex = cc) %>%
      rownames_to_column('ID')
    cc_df$Model <- sprintf('StepCox[%s] + survSVM[γ=%.1f]', direction, gm)
    result <- rbind(result, cc_df)
  }
}

#### 3-9. StepCox + Ridge (logistic) ####
cat("\n===== StepCox + Ridge =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 2) {
    warning(sprintf("Skipped direction '%s': %d features selected (<2)", direction, length(rid)))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  modelexp <- as.matrix(train2[, -c(1, 2), drop = FALSE])
  for (alpha in seq(0, 1, 0.1)) {
    set.seed(MODEL_SEED)
    model_cv <- cv.glmnet(modelexp, train2$OS, family = 'binomial',
                          alpha = alpha, nfolds = 10)
    fit <- glmnet(modelexp, train2$OS, family = 'binomial',
                  alpha = alpha, lambda = model_cv$lambda.min)
    
    rs <- lapply(trainlist2, function(x) {
      pred <- predict(fit, newx = as.matrix(x[, -c(1, 2), drop = FALSE]),
                      type = "response")
      cbind(x[, 1:2], RS = as.numeric(pred))
    })
    rs <- lapply(rs, remove_inf)
    cc <- sapply(rs, calculate_cindex)
    
    cc_df <- data.frame(Cindex = cc) %>%
      rownames_to_column('ID')
    cc_df$Model <- sprintf('StepCox[%s] + Ridge[α=%.1f]', direction, alpha)
    result <- rbind(result, cc_df)
  }
}

#### 3-10. StepCox + obliqueRSF ####
cat("\n===== StepCox + obliqueRSF =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 1) {
    warning(sprintf("Skipped direction '%s': no features selected", direction))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  model_orsf <- orsf(
    data = train2,
    n_tree = 100,
    formula = Surv(OS.time, OS) ~ .
  )
  rs <- lapply(trainlist2, function(x) {
    pred <- predict(model_orsf, new_data = x, pred_type = "risk")[, 1]
    cbind(x[, 1:2], RS = as.numeric(pred))
  })
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, '] + obliqueRSF')
  result <- rbind(result, cc_df)
}

#### 3-11. StepCox + xgboost ####
cat("\n===== StepCox + xgboost =====\n")
set.seed(MODEL_SEED)

# Precompute xgboost DMatrix for the full training set (independent of stepwise direction)
# Correct label construction: positive time for events, negative time for censored
labels_vector <- ifelse(train$OS == 1, train$OS.time, -train$OS.time)
xgb_train_matrix <- xgb.DMatrix(
  data = as.matrix(train[, -c(1, 2)]),
  label = labels_vector
)

for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  
  # Stepwise Cox regression to select features
  step_fit <- step(
    coxph(Surv(OS.time, OS) ~ ., data = train),
    direction = direction
  )
  selected_features <- names(coef(step_fit))
  
  # Subset data to selected features
  train_sub <- train[, c("OS.time", "OS", selected_features), drop = FALSE]
  test_list_sub <- lapply(trainlist, function(x) {
    x[, c("OS.time", "OS", selected_features), drop = FALSE]
  })
  
  # Train xgboost model (using the full training set, not the subset)
  # Note: xgboost uses all genes; stepwise only selects features for interpretation
  xgb_model <- xgb.train(
    params = xgb_params_config,
    data = xgb_train_matrix,
    nrounds = n_rounds,
    verbose = 0
  )
  
  # Predict risk scores on all cohorts
  rs <- lapply(test_list_sub, function(x) {
    pred <- predict(xgb_model, newdata = as.matrix(x[, -c(1, 2)]))
    data.frame(
      OS.time = x$OS.time,
      OS = x$OS,
      RS = as.numeric(pred),
      row.names = rownames(x)
    )
  })
  
  # Remove infinite values and calculate C-index
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  # Compile results
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column("ID")
  cc_df$Model <- paste0("StepCox[", direction, "] + xgboost")
  
  # Append to global results
  result <- rbind(result, cc_df)
}


#### 3-12. StepCox + CForest ####
cat("\n===== StepCox + CForest =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 1) {
    warning(sprintf("Skipped direction '%s': no features selected", direction))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  model_cforest <- party::cforest(
    Surv(OS.time, OS) ~ .,
    data = train2,
    controls = cforest_unbiased(ntree = 50)
  )
  rs <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2],
          RS = as.numeric(predict(model_cforest, newdata = x, type = "response")))
  })
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, '] + CForest')
  result <- rbind(result, cc_df)
}

#### 3-13. StepCox + CTree ####
cat("\n===== StepCox + CTree =====\n")
for (direction in c("both", "backward", "forward")) {
  set.seed(MODEL_SEED)
  fit_step <- step(coxph(Surv(OS.time, OS) ~ ., data = train),
                   direction = direction)
  rid <- names(coef(fit_step))
  if (length(rid) < 1) {
    warning(sprintf("Skipped direction '%s': no features selected", direction))
    next
  }
  
  train2 <- train[, c('OS.time', 'OS', rid), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid), drop = FALSE]
  })
  
  model_ctree <- ctree(Surv(OS.time, OS) ~ ., data = train2)
  rs <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2],
          RS = as.numeric(predict(model_ctree, newdata = x, type = "response")))
  })
  rs <- lapply(rs, remove_inf)
  cc <- sapply(rs, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc) %>%
    rownames_to_column('ID')
  cc_df$Model <- paste0('StepCox[', direction, '] + CTree')
  result <- rbind(result, cc_df)
}


#### 4-1. CoxBoost alone ####
cat("\n===== CoxBoost alone =====\n")
set.seed(MODEL_SEED)

# Optimize penalty
pen_cb <- optimCoxBoostPenalty(
  time = train$OS.time,
  status = train$OS,
  x = as.matrix(train[, -c(1, 2)]),
  trace = TRUE,
  start.penalty = 500,
  parallel = TRUE
)

# Cross-validation
cv_cb <- cv.CoxBoost(
  time = train$OS.time,
  status = train$OS,
  x = as.matrix(train[, -c(1, 2)]),
  maxstepno = 500,
  K = 10,
  type = "verweij",
  penalty = pen_cb$penalty
)

# Final model
fit_cb <- CoxBoost(
  time = train$OS.time,
  status = train$OS,
  x = as.matrix(train[, -c(1, 2)]),
  stepno = cv_cb$optimal.step,
  penalty = pen_cb$penalty
)

rs_cb_alone <- lapply(trainlist, function(x) {
  pred <- predict(
    fit_cb,
    newdata = x[, -c(1, 2)],
    newtime = x[, 1],
    newstatus = x[, 2],
    type = "lp"
  )
  cbind(x[, 1:2], RS = as.numeric(pred))
})
rs_cb_alone <- lapply(rs_cb_alone, remove_inf)
cc_cb_alone <- sapply(rs_cb_alone, calculate_cindex)

cc_df <- data.frame(Cindex = cc_cb_alone) %>%
  rownames_to_column('ID')
cc_df$Model <- 'CoxBoost'
result <- rbind(result, cc_df)

#### 4-2. CoxBoost + Elastic Net ####
cat("\n===== CoxBoost + Enet =====\n")
set.seed(MODEL_SEED)

# Features with non-zero coefficients from CoxBoost
cb_coef <- coef(fit_cb)
rid_cb <- names(cb_coef)[which(cb_coef != 0)]

if (length(rid_cb) == 0) {
  warning("CoxBoost selected no features; skipping CoxBoost+Enet")
  for (alpha in seq(0, 1, 0.1)) {
    cc_df <- create_na_results(paste0('CoxBoost + Enet[α=', alpha, ']'), names(trainlist))
    result <- rbind(result, cc_df)
  }
} else {
  train2 <- train[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  })
  
  x1 <- as.matrix(train2[, rid_cb])
  x2 <- as.matrix(Surv(train2$OS.time, train2$OS))
  
  for (alpha in seq(0, 1, 0.1)) {
    if (ncol(x1) < 2) {
      warning(sprintf("Skipped alpha=%.1f: insufficient features", alpha))
      cc_df <- create_na_results(paste0('CoxBoost + Enet[α=', alpha, ']'), names(trainlist))
      result <- rbind(result, cc_df)
      next
    }
    set.seed(MODEL_SEED)
    fit_enet_cb <- tryCatch(
      cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10),
      error = function(e) NULL
    )
    if (is.null(fit_enet_cb)) {
      cc_df <- create_na_results(paste0('CoxBoost + Enet[α=', alpha, ']'), names(trainlist))
      result <- rbind(result, cc_df)
      next
    }
    
    rs_enet_cb <- lapply(trainlist2, function(x) {
      newx <- as.matrix(x[, rid_cb])
      pred <- predict(fit_enet_cb, type = 'link', newx = newx, s = fit_enet_cb$lambda.min)
      cbind(x[, 1:2], RS = as.numeric(pred))
    })
    rs_enet_cb <- lapply(rs_enet_cb, remove_inf)
    cc_enet_cb <- sapply(rs_enet_cb, calculate_cindex)
    
    cc_df <- data.frame(Cindex = cc_enet_cb) %>%
      rownames_to_column('ID')
    cc_df$Model <- paste0('CoxBoost + Enet[α=', alpha, ']')
    result <- rbind(result, cc_df)
  }
}

#### 4-3. CoxBoost + StepCox ####
cat("\n===== CoxBoost + StepCox =====\n")
if (length(rid_cb) > 0) {
  train2 <- train[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  })
  
  for (direction in c("both", "backward", "forward")) {
    fit_step_cb <- step(coxph(Surv(OS.time, OS) ~ ., data = train2),
                        direction = direction)
    rs_step_cb <- lapply(trainlist2, function(x) {
      cbind(x[, 1:2], RS = predict(fit_step_cb, type = 'risk', newdata = x))
    })
    rs_step_cb <- lapply(rs_step_cb, remove_inf)
    cc_step_cb <- sapply(rs_step_cb, calculate_cindex)
    
    cc_df <- data.frame(Cindex = cc_step_cb) %>%
      rownames_to_column('ID')
    cc_df$Model <- paste0('CoxBoost + StepCox[', direction, ']')
    result <- rbind(result, cc_df)
  }
} else {
  for (direction in c("both", "backward", "forward")) {
    cc_df <- create_na_results(paste0('CoxBoost + StepCox[', direction, ']'), names(trainlist))
    result <- rbind(result, cc_df)
  }
}

#### 4-4. CoxBoost + RSF ####
cat("\n===== CoxBoost + RSF =====\n")
if (length(rid_cb) > 0) {
  train2 <- train[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  })
  
  set.seed(MODEL_SEED)
  fit_rsf_cb <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = train2,
    ntree = 1000,
    nodesize = rf_nodesize,
    splitrule = 'logrank',
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE,
    MODEL_SEED = MODEL_SEED
  )
  best <- which.min(fit_rsf_cb$err.rate)
  set.seed(MODEL_SEED)
  fit_rsf_cb <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = train2,
    ntree = best,
    nodesize = rf_nodesize,
    splitrule = 'logrank',
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE,
    MODEL_SEED = MODEL_SEED
  )
  
  rs_rsf_cb <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2], RS = predict(fit_rsf_cb, newdata = x)$predicted)
  })
  rs_rsf_cb <- lapply(rs_rsf_cb, remove_inf)
  cc_rsf_cb <- sapply(rs_rsf_cb, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_rsf_cb) %>%
    rownames_to_column('ID')
  cc_df$Model <- 'CoxBoost + RSF'
  result <- rbind(result, cc_df)
} else {
  cc_df <- create_na_results('CoxBoost + RSF', names(trainlist))
  result <- rbind(result, cc_df)
}

#### 4-5. CoxBoost + plsRcox ####
cat("\n===== CoxBoost + plsRcox =====\n")
if (length(rid_cb) >= 2) {
  train2 <- train[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  })
  
  pls_data <- list(
    x = as.matrix(train2[, rid_cb]),
    time = train2$OS.time,
    status = train2$OS
  )
  
  set.seed(MODEL_SEED)
  cv_pls_cb <- tryCatch(
    cv.plsRcox(pls_data, nt = 10, nfold = 10, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(cv_pls_cb)) {
    cc_df <- create_na_results('CoxBoost + plsRcox', names(trainlist))
    result <- rbind(result, cc_df)
  } else {
    best_nt <- as.numeric(cv_pls_cb[5])
    fit_pls_cb <- tryCatch(
      plsRcox(Xplan = pls_data$x, time = pls_data$time,
              event = pls_data$status, nt = best_nt),
      error = function(e) NULL
    )
    if (is.null(fit_pls_cb)) {
      cc_df <- create_na_results('CoxBoost + plsRcox', names(trainlist))
      result <- rbind(result, cc_df)
    } else {
      rs_pls_cb <- lapply(trainlist2, function(x) {
        newdata <- as.matrix(x[, rid_cb])
        pred <- predict(fit_pls_cb, type = "lp", newdata = newdata)
        cbind(x[, 1:2], RS = as.numeric(pred))
      })
      rs_pls_cb <- lapply(rs_pls_cb, remove_inf)
      cc_pls_cb <- sapply(rs_pls_cb, calculate_cindex)
      
      cc_df <- data.frame(Cindex = cc_pls_cb) %>%
        rownames_to_column('ID')
      cc_df$Model <- 'CoxBoost + plsRcox'
      result <- rbind(result, cc_df)
    }
  }
} else {
  cc_df <- create_na_results('CoxBoost + plsRcox', names(trainlist))
  result <- rbind(result, cc_df)
}

#### 4-6. CoxBoost + SuperPC ####
cat("\n===== CoxBoost + SuperPC =====\n")
if (length(rid_cb) >= 1) {
  train2 <- train[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  })
  
  prepare_sp_cb <- function(df) {
    x_mat <- t(as.matrix(df[, rid_cb, drop = FALSE]))
    if (is.null(rownames(x_mat))) rownames(x_mat) <- rid_cb
    list(x = x_mat, y = df$OS.time, censoring.status = df$OS,
         featurenames = rownames(x_mat))
  }
  
  data_train <- tryCatch(
    prepare_sp_cb(train2),
    error = function(e) NULL
  )
  if (is.null(data_train)) {
    cc_df <- create_na_results('CoxBoost + SuperPC', names(trainlist))
    result <- rbind(result, cc_df)
  } else {
    set.seed(MODEL_SEED)
    fit_sp_cb <- tryCatch(
      superpc.train(data = data_train, type = 'survival', s0.perc = 0.5),
      error = function(e) NULL
    )
    if (is.null(fit_sp_cb)) {
      cc_df <- create_na_results('CoxBoost + SuperPC', names(trainlist))
      result <- rbind(result, cc_df)
    } else {
      cv_sp_cb <- tryCatch(
        superpc.cv(
          fit_sp_cb, data_train,
          n.threshold = 20,
          n.fold = 5,
          n.components = 3,
          min.features = 1,
          max.features = length(rid_cb),
          compute.fullcv = TRUE,
          compute.preval = TRUE
        ),
        error = function(e) NULL
      )
      if (is.null(cv_sp_cb)) {
        cc_df <- create_na_results('CoxBoost + SuperPC', names(trainlist))
        result <- rbind(result, cc_df)
      } else {
        best_threshold <- tryCatch(
          cv_sp_cb$thresholds[which.max(cv_sp_cb$scor[1, ])],
          error = function(e) 0
        )
        
        rs_sp_cb <- lapply(trainlist2, function(w) {
          data_test <- tryCatch(
            prepare_sp_cb(w),
            error = function(e) NULL
          )
          if (is.null(data_test)) return(cbind(w[, 1:2], RS = NA))
          pred <- tryCatch(
            superpc.predict(fit_sp_cb, data_train, data_test,
                            threshold = best_threshold, n.components = 1)$v.pred,
            error = function(e) rep(NA, nrow(w))
          )
          cbind(w[, 1:2], RS = as.numeric(pred))
        })
        rs_sp_cb <- lapply(rs_sp_cb, remove_inf)
        cc_sp_cb <- sapply(rs_sp_cb, calculate_cindex)
        
        cc_df <- data.frame(Cindex = cc_sp_cb) %>%
          rownames_to_column('ID')
        cc_df$Model <- 'CoxBoost + SuperPC'
        result <- rbind(result, cc_df)
      }
    }
  }
} else {
  cc_df <- create_na_results('CoxBoost + SuperPC', names(trainlist))
  result <- rbind(result, cc_df)
}

#### 4-7. CoxBoost + GBM ####
cat("\n===== CoxBoost + GBM =====\n")
if (length(rid_cb) >= 1) {
  train2 <- train[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  trainlist2 <- lapply(trainlist, function(x) {
    x[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  })
  
  set.seed(MODEL_SEED)
  fit_gbm_cb <- gbm(
    Surv(OS.time, OS) ~ .,
    data = train2,
    distribution = 'coxph',
    n.trees = 1000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 3,
    n.cores = 1
  )
  best_gbm <- which.min(fit_gbm_cb$cv.error)
  set.seed(MODEL_SEED)
  fit_gbm_cb <- gbm(
    Surv(OS.time, OS) ~ .,
    data = train2,
    distribution = 'coxph',
    n.trees = best_gbm,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    n.cores = 1
  )
  
  rs_gbm_cb <- lapply(trainlist2, function(x) {
    cbind(x[, 1:2],
          RS = as.numeric(predict(fit_gbm_cb, x, n.trees = best_gbm, type = 'link')))
  })
  rs_gbm_cb <- lapply(rs_gbm_cb, remove_inf)
  cc_gbm_cb <- sapply(rs_gbm_cb, calculate_cindex)
  
  cc_df <- data.frame(Cindex = cc_gbm_cb) %>%
    rownames_to_column('ID')
  cc_df$Model <- 'CoxBoost + GBM'
  result <- rbind(result, cc_df)
} else {
  cc_df <- create_na_results('CoxBoost + GBM', names(trainlist))
  result <- rbind(result, cc_df)
}

#### 4-8. CoxBoost + survival-SVM ####
cat("\n===== CoxBoost + survival-SVM =====\n")
if (length(rid_cb) >= 1 && nrow(train2) >= 10) {
  train2 <- train[, c('OS.time', 'OS', rid_cb), drop = FALSE]
  train2_clean <- train2[complete.cases(train2), ]
  if (nrow(train2_clean) >= 10) {
    trainlist2 <- lapply(trainlist, function(x) {
      x[, c('OS.time', 'OS', rid_cb), drop = FALSE]
    })
    
    # Try multiple gamma values
    for (gm in c(0.1, 1, 2)) {
      fit_svm_cb <- tryCatch(
        survivalsvm(Surv(OS.time, OS) ~ ., data = train2_clean, gamma.mu = gm),
        error = function(e) NULL
      )
      if (is.null(fit_svm_cb)) next
      
      safe_predict_svm_cb <- function(fit, newdata) {
        tryCatch(
          as.numeric(predict(fit, newdata)$predicted),
          error = function(e) rep(NA, nrow(newdata))
        )
      }
      
      rs_svm_cb <- lapply(trainlist2, function(x) {
        pred <- safe_predict_svm_cb(fit_svm_cb, x)
        cbind(x[, 1:2], RS = pred)
      })
      rs_svm_cb <- lapply(rs_svm_cb, remove_inf)
      cc_svm_cb <- sapply(rs_svm_cb, calculate_cindex)
      
      cc_df <- data.frame(Cindex = cc_svm_cb) %>%
        rownames_to_column('ID')
      cc_df$Model <- sprintf('CoxBoost + survSVM[γ=%.1f]', gm)
      result <- rbind(result, cc_df)
    }
  } else {
    cc_df <- create_na_results('CoxBoost + survival-SVM', names(trainlist))
    result <- rbind(result, cc_df)
  }
} else {
  cc_df <- create_na_results('CoxBoost + survival-SVM', names(trainlist))
  result <- rbind(result, cc_df)
}

# 9. Other Single Models####


#### 5. plsRcox alone ####
cat("\n===== plsRcox =====\n")
set.seed(MODEL_SEED)

if (sum(train$OS) < 2) {
  warning("Insufficient events for plsRcox")
  cc_df <- create_na_results('plsRcox', names(trainlist))
  result <- rbind(result, cc_df)
} else {
  cv_pls <- tryCatch(
    cv.plsRcox(
      list(x = train[, -c(1, 2)], time = train$OS.time, status = train$OS),
      nt = 10,
      nfold = 10,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(cv_pls)) {
    cc_df <- create_na_results('plsRcox', names(trainlist))
    result <- rbind(result, cc_df)
  } else {
    n_comp <- as.numeric(cv_pls[5])
    if (is.na(n_comp)) n_comp <- 1
    fit_pls <- tryCatch(
      plsRcox(
        X = train[, -c(1, 2)],
        time = train$OS.time,
        event = train$OS,
        nt = n_comp
      ),
      error = function(e) NULL
    )
    if (is.null(fit_pls)) {
      cc_df <- create_na_results('plsRcox', names(trainlist))
      result <- rbind(result, cc_df)
    } else {
      rs_pls <- lapply(trainlist, function(x) {
        pred <- tryCatch(
          as.numeric(predict(fit_pls, type = "lp", newdata = x[, -c(1, 2)])),
          error = function(e) rep(NA, nrow(x))
        )
        cbind(x[, 1:2], RS = pred)
      })
      rs_pls <- lapply(rs_pls, remove_inf)
      cc_pls <- sapply(rs_pls, calculate_cindex)
      
      cc_df <- data.frame(Cindex = cc_pls) %>%
        rownames_to_column('ID')
      cc_df$Model <- 'plsRcox'
      result <- rbind(result, cc_df)
    }
  }
}

#### 6. SuperPC alone ####
cat("\n===== SuperPC =====\n")
set.seed(MODEL_SEED)

if (ncol(train) <= 2) {
  warning("Insufficient features for SuperPC")
  cc_df <- create_na_results('SuperPC', names(trainlist))
  result <- rbind(result, cc_df)
} else {
  prepare_sp <- function(df) {
    features <- df[, -c(1, 2), drop = FALSE]
    if (ncol(features) == 0) return(NULL)
    x_mat <- t(as.matrix(features))
    if (is.null(rownames(x_mat))) rownames(x_mat) <- colnames(features)
    list(x = x_mat, y = df$OS.time, censoring.status = df$OS,
         featurenames = rownames(x_mat))
  }
  
  data_train <- prepare_sp(train)
  if (is.null(data_train)) {
    cc_df <- create_na_results('SuperPC', names(trainlist))
    result <- rbind(result, cc_df)
  } else {
    train_features <- rownames(data_train$x)
    fit_sp <- tryCatch(
      superpc.train(data = data_train, type = 'survival', s0.perc = 0.5),
      error = function(e) NULL
    )
    if (is.null(fit_sp)) {
      cc_df <- create_na_results('SuperPC', names(trainlist))
      result <- rbind(result, cc_df)
    } else {
      cv_sp <- tryCatch(
        superpc.cv(
          fit_sp, data_train,
          n.threshold = 20,
          n.fold = 5,
          n.components = 3,
          min.features = 1,
          max.features = nrow(data_train$x),
          compute.fullcv = TRUE,
          compute.preval = TRUE
        ),
        error = function(e) NULL
      )
      if (is.null(cv_sp)) {
        cc_df <- create_na_results('SuperPC', names(trainlist))
        result <- rbind(result, cc_df)
      } else {
        best_threshold <- tryCatch(
          cv_sp$thresholds[which.max(cv_sp$scor[1, ])],
          error = function(e) median(data_train$x, na.rm = TRUE)
        )
        
        rs_sp <- lapply(trainlist, function(w) {
          data_test <- prepare_sp(w)
          if (is.null(data_test)) return(cbind(w[, 1:2], RS = NA))
          pred <- tryCatch(
            superpc.predict(fit_sp, data_train, data_test,
                            threshold = best_threshold, n.components = 1)$v.pred,
            error = function(e) rep(NA, nrow(w))
          )
          cbind(w[, 1:2], RS = as.numeric(pred))
        })
        rs_sp <- lapply(rs_sp, remove_inf)
        cc_sp <- sapply(rs_sp, calculate_cindex)
        
        cc_df <- data.frame(Cindex = cc_sp) %>%
          rownames_to_column('ID')
        cc_df$Model <- 'SuperPC'
        result <- rbind(result, cc_df)
      }
    }
  }
}

#### 7. GBM alone ####
cat("\n===== GBM =====\n")
set.seed(MODEL_SEED)
fit_gbm <- gbm(
  Surv(OS.time, OS) ~ .,
  data = train,
  distribution = 'coxph',
  n.trees = 1000,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.001,
  cv.folds = 5,
  n.cores = 1
)
best_gbm <- which.min(fit_gbm$cv.error)
set.seed(MODEL_SEED)
fit_gbm <- gbm(
  Surv(OS.time, OS) ~ .,
  data = train,
  distribution = 'coxph',
  n.trees = best_gbm,
  interaction.depth = 3,
  n.minobsinnode = 10,
  shrinkage = 0.001,
  n.cores = 1
)

rs_gbm <- lapply(trainlist, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(fit_gbm, x, n.trees = best_gbm, type = 'link')))
})
rs_gbm <- lapply(rs_gbm, remove_inf)
cc_gbm <- sapply(rs_gbm, calculate_cindex)

cc_df <- data.frame(Cindex = cc_gbm) %>%
  rownames_to_column('ID')
cc_df$Model <- 'GBM'
result <- rbind(result, cc_df)

#### 8. survival-SVM alone ####
cat("\n===== survival-SVM =====\n")
set.seed(MODEL_SEED)
train_clean <- train[complete.cases(train), ]
if (nrow(train_clean) < 10) {
  warning("Insufficient complete cases for survival-SVM")
  cc_df <- create_na_results('survival-SVM', names(trainlist))
  result <- rbind(result, cc_df)
} else {
  fit_svm <- tryCatch(
    survivalsvm(Surv(OS.time, OS) ~ ., data = train_clean, gamma.mu = 2),
    error = function(e) NULL
  )
  if (is.null(fit_svm)) {
    cc_df <- create_na_results('survival-SVM', names(trainlist))
    result <- rbind(result, cc_df)
  } else {
    rs_svm <- lapply(trainlist, function(x) {
      pred <- tryCatch(
        as.numeric(predict(fit_svm, x)$predicted),
        error = function(e) rep(NA, nrow(x))
      )
      cbind(x[, 1:2], RS = pred)
    })
    rs_svm <- lapply(rs_svm, remove_inf)
    cc_svm <- sapply(rs_svm, calculate_cindex)
    
    cc_df <- data.frame(Cindex = cc_svm) %>%
      rownames_to_column('ID')
    cc_df$Model <- 'survival-SVM'
    result <- rbind(result, cc_df)
  }
}

#### 9. xgboost alone ####
cat("\n===== xgboost =====\n")
set.seed(MODEL_SEED)
labels_vec <- ifelse(train$OS == 1, train$OS.time, -train$OS.time)
model_mat <- xgb.DMatrix(
  data = as.matrix(train[, -c(1, 2)]),
  label = labels_vec
)
model_xgb <- xgb.train(
  params = xgb_params_config,
  data = model_mat,
  nrounds = n_rounds,
  verbose = 0
)
rs_xgb <- lapply(trainlist, function(x) {
  cbind(x[, 1:2],
        RS = as.numeric(predict(model_xgb, newdata = as.matrix(x[, -c(1, 2)]))))
})
rs_xgb <- lapply(rs_xgb, remove_inf)
cc_xgb <- sapply(rs_xgb, calculate_cindex)

cc_df <- data.frame(Cindex = cc_xgb) %>%
  rownames_to_column('ID')
cc_df$Model <- 'xgboost'
result <- rbind(result, cc_df)

# 10. Combine and Process Results####


cat("\n===== All models completed =====\n")
cat("Total models evaluated:", nrow(result), "\n")


# 10.1 Prepare data for summarization####


#' Combine results from multiple seeds (here only one seed)
#'
#' @param seed Random seed used for the run
#' @param current_results Data frame with current results
#' @param all_results Accumulated results (optional)
#' @return Combined data frame
combine_seed_results <- function(seed, current_results, all_results = NULL) {
  current_results$seed <- seed
  if (is.null(all_results)) return(current_results)
  
  cols_current <- names(current_results)
  cols_all    <- names(all_results)
  
  if (!identical(cols_current, cols_all)) {
    missing_in_all   <- setdiff(cols_current, cols_all)
    missing_in_current <- setdiff(cols_all, cols_current)
    
    if (length(missing_in_all) > 0) {
      for (col in missing_in_all) all_results[[col]] <- NA
    }
    if (length(missing_in_current) > 0) {
      for (col in missing_in_current) current_results[[col]] <- NA
    }
    all_results <- all_results[, cols_current, drop = FALSE]
  }
  
  rbind(all_results, current_results)
}

# Rename cohorts and round C-index
results_renamed <- result %>%
  dplyr::select("ID", "Cindex", "Model") %>%
  dplyr::rename(Cohort = "ID") %>%
  mutate(
    Cohort = case_when(
      Cohort == "Train"          ~ "TCGA-CRC-Trainset",
      Cohort == "Test_TCGA"      ~ "TCGA-CRC-Testset",
      Cohort == "Test_NFYY"      ~ "NFYY-CRC-Testset",
      Cohort == "Test_GSE17536"  ~ "GSE17536-CRC-Testset",
      Cohort == "Test_GSE17537"  ~ "GSE17537-CRC-Testset",
      Cohort == "Test_GSE29621"  ~ "GSE29621-CRC-Testset",
      TRUE ~ as.character(Cohort)
    ),
    Cindex = round(Cindex, 3)
  )

# Split into training and test sets
train_results <- results_renamed %>%
  filter(grepl("Trainset", Cohort, ignore.case = TRUE)) %>%
  as.data.frame()

test_results <- results_renamed %>%
  filter(grepl("Testset", Cohort, ignore.case = TRUE)) %>%
  as.data.frame()

# Pivot to wide format
train_wide <- train_results %>%
  distinct(Model, Cohort, .keep_all = TRUE) %>%
  mutate(column_name = paste0(Cohort, "_Cindex")) %>%
  dplyr::select(Model, column_name, Cindex) %>%
  pivot_wider(names_from = column_name, values_from = Cindex) %>%
  as.data.frame() %>%
  mutate(across(-Model, as.numeric))

test_wide <- test_results %>%
  distinct(Model, Cohort, .keep_all = TRUE) %>%
  mutate(column_name = paste0(Cohort, "_Cindex")) %>%
  dplyr::select(Model, column_name, Cindex) %>%
  pivot_wider(names_from = column_name, values_from = Cindex) %>%
  as.data.frame() %>%
  mutate(across(-Model, as.numeric))

# Combine and compute average C-index values
cindex_combined <- train_wide %>%
  full_join(test_wide, by = "Model") %>%
  mutate(
    na_count   = rowSums(is.na(across(-Model))),
    has_na     = na_count > 0,
    allset_mean = round(rowMeans(across(-c(Model, na_count, has_na)), na.rm = TRUE), 3),
    testset_mean = round(rowMeans(across(matches("Testset")), na.rm = TRUE), 3)
  ) %>%
  dplyr::select(-c(na_count, has_na))

# Merge with seed info (only one seed here)
all_results <- combine_seed_results(MODEL_SEED, cindex_combined)

# Remove rows with all NA (if any)
all_results <- na.omit(all_results)

# Sort by test set average, then overall average
sorted_results <- all_results %>%
  arrange(desc(testset_mean), desc(allset_mean))


# 10.2 Identify identical models (same numeric fingerprint)####


identify_identical_models <- function(sorted_results) {
  # Exclude non‑numeric columns
  numeric_cols <- names(sorted_results)[sapply(sorted_results, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, "seed")
  
  # Create a fingerprint by concatenating rounded numeric values
  fingerprints <- apply(
    sorted_results[, numeric_cols, drop = FALSE],
    1,
    function(x) paste(round(as.numeric(x), 6), collapse = "|")
  )
  
  # Group models by fingerprint
  groups <- split(sorted_results$Model, fingerprints)
  duplicate_groups <- groups[sapply(groups, length) > 1]
  
  if (length(duplicate_groups) == 0) {
    cat("No identical models found.\n")
    return(invisible(NULL))
  }
  
  # Extract performance for each group
  group_metrics <- lapply(names(duplicate_groups), function(fp) {
    idx <- which(fingerprints == fp)[1]
    c(
      testset_mean = sorted_results$testset_mean[idx],
      allset_mean  = sorted_results$allset_mean[idx]
    )
  })
  
  group_info <- data.frame(
    fingerprint   = names(duplicate_groups),
    testset_mean  = sapply(group_metrics, `[`, "testset_mean"),
    allset_mean   = sapply(group_metrics, `[`, "allset_mean"),
    model_count   = sapply(duplicate_groups, length),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(testset_mean), desc(allset_mean))
  
  # Print report
  cat("\n", rep("═", 50), "\n", sep = "")
  cat("Identical model groups detected:\n")
  for (i in seq_len(nrow(group_info))) {
    cat(sprintf("\nGroup %d (testset_mean = %.3f, allset_mean = %.3f, %d models):\n",
                i, group_info$testset_mean[i], group_info$allset_mean[i],
                group_info$model_count[i]))
    cat("  ", paste(duplicate_groups[[group_info$fingerprint[i]]], collapse = "\n  "), "\n")
  }
  cat("\n", rep("═", 50), "\n", sep = "")
  
  invisible(list(
    groups = duplicate_groups,
    info   = group_info
  ))
}

identical_models <- identify_identical_models(sorted_results)


# 10.3 Manual adjustment of model order (for identical models)####


#' Adjust order of two models that are functionally identical
#'
#' @param sorted_results Data frame sorted by performance
#' @return Adjusted data frame
adjust_model_order <- function(sorted_results) {
  adjusted <- sorted_results
  
  # Define pairs to swap: first should come before second
  # Here: xgboost and StepCox[forward] + xgboost are identical
  manual_pairs <- list(
    c("xgboost", "StepCox[forward] + xgboost",
      "StepCox[forward] selected no features, identical to xgboost")
  )
  
  cat("\n", rep("─", 60), "\n", sep = "")
  cat("Manual model order adjustment\n")
  
  success_count <- 0
  skip_count    <- 0
  
  for (pair in manual_pairs) {
    model_first  <- pair[1]
    model_second <- pair[2]
    reason       <- pair[3]
    
    cat(sprintf("\nPair: '%s' → '%s'\n  Reason: %s\n", model_first, model_second, reason))
    
    idx_first  <- which(adjusted$Model == model_first)
    idx_second <- which(adjusted$Model == model_second)
    
    if (length(idx_first) == 0 || length(idx_second) == 0) {
      missing <- c()
      if (length(idx_first) == 0)  missing <- c(missing, model_first)
      if (length(idx_second) == 0) missing <- c(missing, model_second)
      cat("  Skipped: model(s) not found:", paste(missing, collapse = ", "), "\n")
      skip_count <- skip_count + 1
      next
    }
    
    # Check that they are indeed identical
    row_first  <- adjusted[idx_first, ]
    row_second <- adjusted[idx_second, ]
    
    numeric_cols <- names(adjusted)[sapply(adjusted, is.numeric)]
    numeric_cols <- setdiff(numeric_cols, "seed")
    identical_perf <- all(sapply(numeric_cols, function(col) {
      isTRUE(all.equal(row_first[[col]], row_second[[col]], tolerance = 1e-6))
    }))
    
    if (!identical_perf) {
      cat("  Skipped: models are not numerically identical\n")
      skip_count <- skip_count + 1
      next
    }
    
    if (idx_first < idx_second) {
      cat("  Already in correct order\n")
      next
    }
    
    # Perform swap
    cat("  Swapping...\n")
    # Remove the two rows
    temp_df <- adjusted[-c(idx_first, idx_second), , drop = FALSE]
    
    # Insert model_first before model_second
    # Determine insertion point: maintain overall sorting
    insert_pos <- find_insert_position(temp_df, row_first, row_second)
    
    if (insert_pos == 1) {
      adjusted <- rbind(row_first, row_second, temp_df)
    } else if (insert_pos > nrow(temp_df)) {
      adjusted <- rbind(temp_df, row_first, row_second)
    } else {
      adjusted <- rbind(
        temp_df[1:(insert_pos - 1), , drop = FALSE],
        row_first,
        row_second,
        temp_df[insert_pos:nrow(temp_df), , drop = FALSE]
      )
    }
    rownames(adjusted) <- NULL
    success_count <- success_count + 1
    cat("  Done\n")
  }
  
  cat("\n", rep("─", 60), "\n", sep = "")
  cat("Adjustment summary:\n")
  cat("  Total pairs :", length(manual_pairs), "\n")
  cat("  Successful  :", success_count, "\n")
  cat("  Skipped     :", skip_count, "\n")
  cat(rep("─", 60), "\n", sep = "")
  
  adjusted
}

#' Helper: find where to insert a pair to maintain sorting
find_insert_position <- function(temp_df, row_first, row_second) {
  if (nrow(temp_df) == 0) return(1)
  
  test_mean_first <- row_first$testset_mean
  all_mean_first  <- row_first$allset_mean
  
  for (i in seq_len(nrow(temp_df))) {
    if (temp_df$testset_mean[i] < test_mean_first) {
      return(i)
    } else if (temp_df$testset_mean[i] == test_mean_first &&
               temp_df$allset_mean[i] < all_mean_first) {
      return(i)
    }
  }
  nrow(temp_df) + 1
}

# Apply manual adjustment
sorted_results <- adjust_model_order(sorted_results)


# 11. Visualization of C-index Results####



# 11.1 Color palette exploration (optional, can be removed)####

# The following lines demonstrate various colour palettes.
# They are kept for reference but can be commented out in production.
if (FALSE) {
  library(paletteer)
  paletteer_c("scico::berlin", n = 10)
  paletteer_d("RColorBrewer::Paired")
  paletteer_d("ggsci::nrc_npg")
  paletteer_d("ggsci::default_jco")
  paletteer_d("ggsci::default_nejm")
  paletteer_d("ggsci::lanonc_lancet")
  paletteer_d("ggsci::default_jama")
  paletteer_dynamic("cartography::green.pal", 5)
  paletteer_dynamic("cartography::green.pal", 10)
  
  library(cols4all)
  mycol1 <- c4a('pastel', 6)
  c4a_plot(mycol1)
  mycol2 <- c4a('classic_blue_red12', 6)
  c4a_plot(mycol2)
  mycol3 <- c4a('dark2', 6)
  c4a_plot(mycol3)
  mycol4 <- c4a('accent', 6)
  c4a_plot(mycol4)
  
  library(ggsci)
  library(scales)
  npg_colors <- pal_npg("nrc", alpha = 0.7)(9)
  show_col(npg_colors)
  aaas_colors <- pal_aaas("default", alpha = 1)(5)
  show_col(aaas_colors)
  d3_colors <- pal_d3("category10")(10)
  show_col(d3_colors)
}

# Colour palette for heatmap columns
heatmap_colors <- c(
  "#BD3C29", "#AE7000", "#925E9FFF", "#00468BFF",
  "#008280FF", "#B09C85FF"
)


# 11.2 Full heatmap (all models)####


# Prepare matrix: exclude seed and average columns
numeric_cols <- names(sorted_results)[sapply(sorted_results, is.numeric)]
numeric_cols <- setdiff(numeric_cols, c("seed", "allset_mean", "testset_mean"))
heatmap_mat <- as.matrix(sorted_results[, numeric_cols, drop = FALSE])
rownames(heatmap_mat) <- sorted_results$Model

# Desired column order
desired_order <- c(
  "TCGA-CRC-Trainset_Cindex",
  "TCGA-CRC-Testset_Cindex",
  "NFYY-CRC-Testset_Cindex",
  "GSE17536-CRC-Testset_Cindex",
  "GSE17537-CRC-Testset_Cindex",
  "GSE29621-CRC-Testset_Cindex"
)
desired_order <- intersect(desired_order, colnames(heatmap_mat))
heatmap_mat <- heatmap_mat[, desired_order, drop = FALSE]

# Column annotation colours
col_anno <- heatmap_colors[seq_along(desired_order)]
names(col_anno) <- desired_order

# Average vectors for right annotation
avg_test <- sorted_results$testset_mean
avg_all  <- sorted_results$allset_mean

hm_full <- ComplexHeatmap::Heatmap(
  heatmap_mat,
  name = "C-index",
  col = circlize::colorRamp2(
    c(0.50, 0.65, 0.80),
    c("#20ACBD", "#FFFFFF", "#F69896")
  ),
  row_gap = unit(1, "mm"),
  column_gap = unit(3, "mm"),
  rect_gp = gpar(col = "grey", lwd = 1),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  width = unit(ncol(heatmap_mat) + 3, "cm"),
  height = unit(nrow(heatmap_mat) / 2, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    legend_height = unit(4, "cm"),
    grid_width = unit(0.6, "cm")
  ),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  top_annotation = HeatmapAnnotation(
    Dataset = colnames(heatmap_mat),
    col = list(Dataset = col_anno),
    show_annotation_name = FALSE,
    simple_anno_size = unit(0.5, "cm")
  ),
  right_annotation = rowAnnotation(
    "Testset\nAvg" = anno_barplot(
      avg_test,
      bar_width = 0.8,
      border = TRUE,
      gp = gpar(
        fill = sapply(avg_test, function(x) {
          rng <- range(avg_test, na.rm = TRUE)
          if (diff(rng) == 0) intensity <- 0.5
          else intensity <- (x - rng[1]) / diff(rng)
          rgb(colorRamp(c("#FFCCCC", "#B84D64"))(intensity), maxColorValue = 255)
        }),
        col = "grey"
      ),
      add_numbers = TRUE,
      numbers_offset = unit(-8, "mm"),
      numbers_gp = gpar(fontsize = 8, col = "white"),
      width = unit(2, "cm")
    ),
    "Allset\nAvg" = anno_barplot(
      avg_all,
      bar_width = 0.8,
      border = TRUE,
      gp = gpar(
        fill = sapply(avg_all, function(x) {
          rng <- range(avg_all, na.rm = TRUE)
          if (diff(rng) == 0) intensity <- 0.5
          else intensity <- (x - rng[1]) / diff(rng)
          rgb(colorRamp(c("#E0E0E0", "#666666"))(intensity), maxColorValue = 255)
        }),
        col = "grey"
      ),
      add_numbers = TRUE,
      numbers_offset = unit(-8, "mm"),
      numbers_gp = gpar(fontsize = 8, col = "white"),
      width = unit(2, "cm")
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 10),
    gap = unit(2, "mm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", heatmap_mat[i, j]), x, y, gp = gpar(fontsize = 7))
  }
)

draw(hm_full)


# 11.3 Top N models heatmap####


n_show <- 15
if (nrow(heatmap_mat) > n_show) {
  top_mat <- heatmap_mat[1:n_show, , drop = FALSE]
  top_avg_test <- avg_test[1:n_show]
  top_avg_all  <- avg_all[1:n_show]
} else {
  top_mat <- heatmap_mat
  top_avg_test <- avg_test
  top_avg_all  <- avg_all
}

hm_top <- ComplexHeatmap::Heatmap(
  top_mat,
  name = "C-index",
  col = circlize::colorRamp2(
    c(0.50, 0.65, 0.80),
    c("#20ACBD", "#FFFFFF", "#F69896")
  ),
  row_gap = unit(1, "mm"),
  column_gap = unit(3, "mm"),
  rect_gp = gpar(col = "grey", lwd = 1),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  width = unit(ncol(top_mat) + 3, "cm"),
  height = unit(nrow(top_mat) / 1.8, "cm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 12),
    legend_height = unit(4, "cm"),
    grid_width = unit(0.6, "cm")
  ),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  top_annotation = HeatmapAnnotation(
    Dataset = colnames(top_mat),
    col = list(Dataset = col_anno),
    show_annotation_name = FALSE,
    simple_anno_size = unit(0.5, "cm")
  ),
  right_annotation = rowAnnotation(
    "Testset\nAvg" = anno_barplot(
      top_avg_test,
      bar_width = 0.8,
      border = TRUE,
      gp = gpar(
        fill = sapply(top_avg_test, function(x) {
          rng <- range(top_avg_test, na.rm = TRUE)
          if (diff(rng) == 0) intensity <- 0.5
          else intensity <- (x - rng[1]) / diff(rng)
          rgb(colorRamp(c("#FFCCCC", "#B84D64"))(intensity), maxColorValue = 255)
        }),
        col = "grey"
      ),
      add_numbers = TRUE,
      numbers_offset = unit(-8, "mm"),
      numbers_gp = gpar(fontsize = 8, col = "white"),
      width = unit(2, "cm")
    ),
    "Allset\nAvg" = anno_barplot(
      top_avg_all,
      bar_width = 0.8,
      border = TRUE,
      gp = gpar(
        fill = sapply(top_avg_all, function(x) {
          rng <- range(top_avg_all, na.rm = TRUE)
          if (diff(rng) == 0) intensity <- 0.5
          else intensity <- (x - rng[1]) / diff(rng)
          rgb(colorRamp(c("#E0E0E0", "#666666"))(intensity), maxColorValue = 255)
        }),
        col = "grey"
      ),
      add_numbers = TRUE,
      numbers_offset = unit(-8, "mm"),
      numbers_gp = gpar(fontsize = 8, col = "white"),
      width = unit(2, "cm")
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 10),
    gap = unit(2, "mm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", top_mat[i, j]), x, y, gp = gpar(fontsize = 7))
  }
)

draw(hm_top)






# 12. Session Information####

cat("\n===== Session Information =====\n")
sessionInfo()
# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26100)
# 
# Matrix products: default
# LAPACK version 3.12.1
# 
# locale:
#   [1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8
# [4] LC_NUMERIC=C                                LC_TIME=Chinese (Simplified)_China.utf8    
# 
# time zone: Asia/Shanghai
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] risksetROC_1.0.4.1          openxlsx_4.2.8              partykit_1.2-24             libcoin_1.0-10              party_1.3-18                strucchange_1.5-4          
# [7] sandwich_3.1-1              zoo_1.8-14                  modeltools_0.2-24           mvtnorm_1.3-3               xgboost_1.7.11.1            aorsf_0.1.5                
# [13] remotes_2.5.0               obliqueRSF_0.1.2            caret_7.0-1                 lattice_0.22-7              timeROC_0.4                 regplot_1.1                
# [19] foreign_0.8-90              ggDCA_1.2                   pec_2023.04.12              Hmisc_5.2-3                 Rcpp_1.0.14                 ggrepel_0.9.6              
# [25] DESeq2_1.48.1               SummarizedExperiment_1.38.1 Biobase_2.68.0              MatrixGenerics_1.20.0       matrixStats_1.5.0           GenomicRanges_1.60.0       
# [31] GenomeInfoDb_1.44.0         IRanges_2.42.0              S4Vectors_0.46.0            BiocGenerics_0.54.0         generics_0.1.4              gplots_3.2.0               
# [37] future.apply_1.20.0         future_1.58.0               ggsignif_0.6.4              pheatmap_1.0.13             beepr_2.0                   lubridate_1.9.4            
# [43] forcats_1.0.0               stringr_1.5.1               purrr_1.0.4                 readr_2.1.5                 tidyverse_2.0.0             stringi_1.8.7              
# [49] survminer_0.5.0             ggpubr_0.6.1                edgeR_4.6.2                 limma_3.64.1                tidyr_1.3.1                 ggbreak_0.1.5              
# [55] BART_2.9.9                  nlme_3.1-168                tibble_3.3.0                dplyr_1.1.4                 survivalsvm_0.0.6           CoxBoost_1.5               
# [61] prodlim_2025.04.28          devtools_2.4.5              usethis_3.1.0               gbm_2.2.2                   superpc_1.12                plsRcox_1.7.7              
# [67] glmnet_4.1-9                Matrix_1.7-3                snowfall_1.84-6.3           snow_0.4-4                  randomForestSRC_3.4.1       survival_3.8-3             
# [73] ggsci_4.1.0                 ggplot2_3.5.2               data.table_1.17.6           RColorBrewer_1.1-3          circlize_0.4.16             ComplexHeatmap_2.24.1      
# 
# loaded via a namespace (and not attached):
#   [1] coin_1.4-3              dichromat_2.0-0.1       urlchecker_1.0.1        rARPACK_0.11-0          nnet_7.3-20             TH.data_1.1-3           vctrs_0.6.5            
# [8] digest_0.6.37           png_0.1-8               corpcor_1.6.10          shape_1.4.6.1           parallelly_1.45.0       permute_0.9-8           magick_2.8.7           
# [15] MASS_7.3-65             plsRglm_1.5.1           reshape2_1.4.4          httpuv_1.6.16           foreach_1.5.2           withr_3.0.2             xfun_0.52              
# [22] ggfun_0.1.9             doRNG_1.8.6.2           ellipsis_0.3.2          memoise_2.0.1           MatrixModels_0.5-4      profvis_0.4.0           GlobalOptions_0.1.2    
# [29] gtools_3.9.5            Formula_1.2-5           ellipse_0.5.0           promises_1.3.3          httr_1.4.7              rstatix_0.7.2           globals_0.18.0         
# [36] rstudioapi_0.17.1       UCSC.utils_1.4.0        miniUI_0.1.2            inum_1.0-5              missForest_1.5          base64enc_0.1-3         bootstrap_2019.6       
# [43] SuppDists_1.1-9.9       fields_16.3.1           randomForest_4.7-1.2    quadprog_1.5-8          GenomeInfoDbData_1.2.14 SparseArray_1.8.0       xtable_1.8-4           
# [50] pracma_2.4.4            doParallel_1.0.17       rms_8.0-0               evaluate_1.0.4          S4Arrays_1.8.1          hms_1.1.3               colorspace_2.1-1       
# [57] visNetwork_2.1.2        magrittr_2.0.3          later_1.4.2             SparseM_1.84-2          class_7.3-23            pillar_1.11.0           iterators_1.0.14       
# [64] sna_2.8                 caTools_1.18.3          compiler_4.5.1          RSpectra_0.16-2         rmeta_3.0               gower_1.0.2             minqa_1.2.8            
# [71] plyr_1.8.9              crayon_1.5.3            abind_1.4-8             mixOmics_6.32.0         gridGraphics_0.5-1      sm_2.2-6.0              locfit_1.5-9.12        
# [78] codetools_0.2-20        multcomp_1.4-28         recipes_1.3.1           GetoptLong_1.0.5        mime_0.13               splines_4.5.1           survivalROC_1.0.3.1    
# [85] quantreg_6.1            lars_1.3                bipartite_2.21          knitr_1.50              clue_0.3-66             lme4_1.1-37             itertools_0.1-3        
# [92] fs_1.6.6                listenv_0.9.1           checkmate_2.3.2         Rdpack_2.6.4            pkgbuild_1.4.8          ggplotify_0.1.2         statmod_1.5.0          
# [99] tzdb_0.5.0              pkgconfig_2.0.3         network_1.19.0          tools_4.5.1             cachem_1.1.0            rbibutils_2.3           viridisLite_0.4.2      
# [106] numDeriv_2016.8-1.1     fastmap_1.2.0           rmarkdown_2.29          scales_1.4.0            audio_0.1-11            broom_1.0.8             patchwork_1.3.1        
# [113] coda_0.19-4.1           dotCall64_1.2           carData_3.0-5           rpart_4.1.24            farver_2.1.2            reformulas_0.4.1        mgcv_1.9-3             
# [120] DiagrammeR_1.0.11       ggthemes_5.1.0          cli_3.6.5               lifecycle_1.0.4         lava_1.8.1              kernlab_0.9-33          sessioninfo_1.2.3      
# [127] backports_1.5.0         BiocParallel_1.42.1     timechange_0.3.0        gtable_0.3.6            rjson_0.2.23            pROC_1.18.5             parallel_4.5.1         
# [134] jsonlite_2.0.0          bitops_1.0-9            yulab.utils_0.2.0       vegan_2.7-1             zip_2.3.3               polspline_1.1.25        survMisc_0.5.6         
# [141] timeDate_4041.110       set_1.2                 shiny_1.11.1            collapse_2.1.2          htmltools_0.5.8.1       KMsurv_0.1-6            data.tree_1.1.0        
# [148] glue_1.8.0              timereg_2.0.6           spam_2.11-1             XVector_0.48.0          gridExtra_2.3           boot_1.3-31             survcomp_1.58.0        
# [155] igraph_2.1.4            R6_2.6.1                km.ci_0.5-6             rngtools_1.5.2          cluster_2.1.8.1         pkgload_1.4.0           survAUC_1.3-0          
# [162] beanplot_1.3.1          aplot_0.2.8             ipred_0.9-15            nloptr_2.2.1            statnet.common_4.12.0   vioplot_0.5.1           DelayedArray_0.34.1    
# [169] tidyselect_1.2.1        htmlTable_2.4.3         maps_3.4.3              car_3.1-3               ModelMetrics_1.2.2.2    KernSmooth_2.23-26      htmlwidgets_1.6.4      
# [176] pls_2.8-5               rlang_1.1.6             Cairo_1.6-2             hardhat_1.4.1     



# cat("\nKey package versions:\n")
# cat("survival:", packageVersion("survival"), "\n")
# cat("randomForestSRC:", packageVersion("randomForestSRC"), "\n")
# cat("glmnet:", packageVersion("glmnet"), "\n")
# cat("xgboost:", packageVersion("xgboost"), "\n")
# cat("ComplexHeatmap:", packageVersion("ComplexHeatmap"), "\n")

# End of script










# 13. Save the trained xgb model and training dataset locally####

# Save the trained xgb model and training dataset locally
cat("\n===== Saving xgb model as PCDscore and training dataset as TCGA-CRC-Trainset =====\n")
setwd(BASE_DIR)
getwd()
# Save the trained xgb model to a local file named PCDscore.model
xgb.save(model_xgb, "PCDscore.model")
cat("Model saved successfully as PCDscore.model\n")

# Save the training dataset for future standardization purposes
saveRDS(train_exp, "TCGA-CRC-Trainset.rds")
cat("Training dataset saved successfully as TCGA-CRC-Trainset.rds\n")

