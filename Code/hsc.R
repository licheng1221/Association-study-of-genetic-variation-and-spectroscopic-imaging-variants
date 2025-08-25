set.seed(66)

library(kernlab)
library(factoextra)
library(paran)
library(tidyverse)

# --------------------------
# User paths (EDIT THESE)
# --------------------------
spec_rds <- "/path/to/Pheno_corr.rds"
out_dir  <- "path/to/your/output/directory"


# Read spectra data
Spec = readRDS(spec_rds)
# Remove the "Taxa" column and spectral of 350-399 nm before calculating distances
Spec_no_taxa <- Spec[,-c(1:51)]

# Transpose it because dist() calculates distances between rows, not columns
Spec_t <- t(Spec_no_taxa)

# Initialize a global variable to keep track of the total number of segments
total_segments <- 0

perform_pa_and_pca <- function(spec, level) {
  # Increment the global variable
  total_segments <<- total_segments + 1
  segment <- total_segments
  
  # Perform parallel analysis
  pa <- paran(spec, iterations = 100, centile = 0)
  N_PCs = pa$Retained
  
  # Perform PCA
  PCA <- prcomp(spec, scale. = TRUE)
  
  # Get variances
  variances <- PCA$sdev^2
  # Get proportion of variance explained
  prop_variance <- variances / sum(variances)
  # Get cumulative proportion of variance explained
  cumulative_prop_variance <- cumsum(prop_variance)
  
  # Write the pcs to a file
  write.table(PCA$x, paste0("Seg", segment, "_pcs.txt"), quote = FALSE, row.names = FALSE)
  
  # Create a data frame with the results
  results_df <- data.frame(segment = segment, 
                           level = level, 
                           N_PCs = N_PCs, 
                           Per_expl = cumulative_prop_variance[N_PCs])

  # Create a data frame with the segment information
  seg_df <- data.frame(wavelength = colnames(spec), 
                       level = level, 
                       segment = segment)

  return(list(results_df = results_df, seg_df = seg_df, N_PCs = N_PCs))
}

recursive_spectral_clustering_PA <- function(data, original_data, current_level = 1, cluster_number = 1, results_df, seg_df) {
  # Perform spectral clustering
  specc_result <- specc(as.matrix(data), centers = 2)
  clusters <- specc_result@.Data
  
  # Separate data into two matrices
  data1 <- data[clusters == 1, ]
  data2 <- data[clusters == 2, ]
  
  # Corresponding subsets of the original data
  original_data1 <- original_data[, rownames(data1)]
  original_data2 <- original_data[, rownames(data2)]
  
  
  # Perform parallel analysis and PCA on each subset of the original data
  list1 <- perform_pa_and_pca(original_data1, current_level)
  list2 <- perform_pa_and_pca(original_data2, current_level)
  
  # Append results to results dataframe
  results_df <- rbind(results_df, list1$results_df, list2$results_df)

  # Append segment information to seg_df
  seg_df <- rbind(seg_df, list1$seg_df, list2$seg_df)

  # If N_PCs is larger than 1, perform spectral clustering recursively on each subset
  if (list1$N_PCs > 1) {
    list_results <- recursive_spectral_clustering_PA(data1, original_data1, current_level + 1, 2 * cluster_number - 1, results_df, seg_df)
    results_df <- list_results$results_df
    seg_df <- list_results$seg_df
  }
  
  if (list2$N_PCs > 1) {
    list_results <- recursive_spectral_clustering_PA(data2, original_data2, current_level + 1, 2 * cluster_number, results_df, seg_df)
    results_df <- list_results$results_df
    seg_df <- list_results$seg_df
  }
  
  return(list(results_df = results_df, seg_df = seg_df))
}

# Initialize results_df
results_df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(results_df) <- c("segment", "level", "N_PCs", "Per_expl")

# Initialize seg_df
seg_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(seg_df) <- c("wavelength", "level", "segment")

# Run the combined function
results <- recursive_spectral_clustering_PA(Spec_t, Spec_no_taxa, results_df = results_df, seg_df = seg_df)

# Save results_df and seg_df
results_df <- results$results_df
seg_df <- results$seg_df

# Write results_df
write.csv(results_df,
          file = file.path(out_dir, "results.csv"),
          row.names = FALSE, quote = FALSE)
saveRDS(results_df,
        file = file.path(out_dir, "results.rds"))

# Write seg_df
write.csv(seg_df,
          file = file.path(out_dir, "segments.csv"),
          row.names = FALSE, quote = FALSE)
saveRDS(seg_df,
        file = file.path(out_dir, "segments.rds"))
