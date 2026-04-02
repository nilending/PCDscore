# =============================================================================
# R Script: Load training dataset for standardization, generate test data and evaluate PCDscore model
# Author: [Xianwen Lin]
# Date: 2026-02-19
# R Version: 4.5.1
# =============================================================================




# Load training dataset for standardization, generate test data and evaluate PCDscore model####
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

cat("\n===== Loading training dataset for standardization, generating test data and evaluating PCDscore model =====\n")
setwd("D:/MR-CRC")
getwd
# Load the training dataset for standardization parameters
TCGA_CRC_Trainset <- readRDS("TCGA-CRC-Trainset.rds")
cat("Training dataset loaded successfully for standardization parameters\n")

# Extract gene expression data from training set to calculate standardization parameters
TCGA_CRC_Trainset_t <- t.data.frame(TCGA_CRC_Trainset)
gene_names <- c("OPRL1", "CDKN1A", "CD274", "SPI1", "ELANE", "TGFB1")
train_gene_data <- TCGA_CRC_Trainset_t[, gene_names, drop = FALSE]
gene_means <- apply(train_gene_data, 2, mean, na.rm = TRUE)
gene_sds <- apply(train_gene_data, 2, sd, na.rm = TRUE)

# Generate a new standardized test dataset following the format of GSE29621_CRC_test_dataset
set.seed(123)  # For reproducible results

# Define gene names that should match those in the training set
available_genes <- colnames(train_gene_data)
example_genes <- head(available_genes, 6)  # Take first 6 genes as example

# Create a demonstration test dataset with 50 samples
num_samples <- 50
new_test_dataset <- data.frame(
  OS.time = runif(num_samples, min = 5, max = 120),  # Survival time in months
  OS = rbinom(num_samples, 1, 0.3),                  # Event indicator (0=censored, 1=event)
  stringsAsFactors = FALSE
)

# Add gene expression values (randomly generated but following realistic patterns)
for (gene in example_genes) {
  new_test_dataset[[gene]] <- rtruncnorm(num_samples, a = 0, mean = 30, sd = 15)
}


# Assign row names to match expected format 
rownames(new_test_dataset) <- paste0("Sample", sprintf("%06d", 100000:(100000+num_samples-1)))

# Display structure of new test dataset
cat("New test dataset structure:\n")
print(head(new_test_dataset))

# Load the PCDscore model from local storage
PCDscore_model <- xgb.load("PCDscore.model")
cat("PCDscore model loaded successfully\n")

# Standardize the new test dataset using parameters from the training set
new_test_gene_data <- new_test_dataset[, example_genes, drop = FALSE]

# Apply z-score transformation using training set's mean and SD
new_test_standardized <- as.data.frame(matrix(NA, nrow = nrow(new_test_gene_data), ncol = ncol(new_test_gene_data)))
colnames(new_test_standardized) <- colnames(new_test_gene_data)
rownames(new_test_standardized) <- rownames(new_test_gene_data)

for (gene in example_genes) {
  if (gene %in% names(gene_means) && gene %in% names(gene_sds)) {
    if (gene_sds[gene] != 0) {
      new_test_standardized[[gene]] <- (new_test_gene_data[[gene]] - gene_means[gene]) / gene_sds[gene]
    } else {
      new_test_standardized[[gene]] <- new_test_gene_data[[gene]] - gene_means[gene]
    }
  } else {
    warning(paste("Gene", gene, "not found in training set parameters"))
    new_test_standardized[[gene]] <- new_test_gene_data[[gene]]  # Keep original if no parameters
  }
}

# Make predictions using the loaded model
new_test_matrix <- as.matrix(new_test_standardized)
predictions <- predict(PCDscore_model, newdata = new_test_matrix)

# Add predictions to the dataset
new_test_dataset_with_predictions <- cbind(new_test_dataset[, c("OS.time", "OS")], 
                                           RS = as.numeric(predictions))

# Remove any rows with invalid risk scores
valid_data <- new_test_dataset_with_predictions[
  is.finite(new_test_dataset_with_predictions$RS) & 
    !is.na(new_test_dataset_with_predictions$RS), ]

# Calculate C-index using survcomp package
if (nrow(valid_data) > 0 && sum(valid_data$OS) >= 2) {
  c_index_result <- survcomp::concordance.index(
    x = valid_data$RS,
    surv.time = valid_data$OS.time,
    surv.event = valid_data$OS,
    method = "noether"
  )
  
  cat("C-index for new test dataset using PCDscore model:", round(c_index_result$c.index, 3), "\n")
} else {
  cat("Not enough valid samples to calculate C-index\n")
  if (nrow(valid_data) == 0) {
    cat("All predicted risk scores were invalid\n")
  } else if (sum(valid_data$OS) < 2) {
    cat("Insufficient events (", sum(valid_data$OS), ") to calculate C-index\n", sep = "")
  }
}

cat("Evaluation completed successfully\n")



# Session Information
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
