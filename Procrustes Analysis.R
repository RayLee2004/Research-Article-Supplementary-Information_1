# Display Procrustes analysis and plot for only one group here as an example.
# Other groups can be analyzed likely after modifying the imported dataset.

# Load required packages
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)

########## PART 1 PROCRUSTES ANALYZING #########

# Read in merged dataset.
load("datasets/merged_datasets.Rdata")

# Set fish as object to analyze.
# Create a function to transpose tables to meet computational format.
transpose <- function(y) {
  current_df <- as.data.frame(y)
  col_names <- current_df[, 1]
  current_df <- t(current_df[, -1])
  colnames(current_df) <- col_names
  return(as.data.frame(current_df))
}

# Execute the function
fish <- transpose(fish)
invertebrate <- transpose(invertebrate)

# Perform Hellinger transformation for each group.
fish_hel <- decostand(fish, method = 'hellinger')
invertebrate_hel <- decostand(invertebrate, method = 'hellinger')

# Perform PCA for each group.
fish_pca <- rda(fish_hel, scale = TRUE)
invertebrate_pca <- rda(invertebrate_hel, scale = FALSE)

# Align the two dimensionality reduction results and perform a Procrustes analysis.
proc <- procrustes(X = fish_pca, Y = invertebrate_pca, symmetric = TRUE)

# Summary Procrustes analyzing results.
summary(proc)

# Conduct significance tests.
set.seed(123) # Set random seeds to ensure the reproducibility of the results
prot <- protest(fish_pca, invertebrate_pca, permutations = 999) # Fish as X variable.

# Show the significance tests results.
print(prot$ss)
print(prot$signif)

########## PART 2 ANALYSIS RESULTS PLOTING #########
# Figure 1: Procrustes Errors

# The format and size of plot output can be modified here.
png("Procrustes Errors.png", width = 2000, height = 1000)

# Plot margin (bottom→left→top→right)
par(mar = c(6.5, 6.5, 4.3, 5)) 
# Distance between axis.lab/axis.content and axis.
par(mgp = c(4, 1.4, 0)) 

plot(
  proc, 
  kind = 2, 
  cex.axis = 1.75,
  cex.lab = 2.6,
  main = "",  # Leave the title blank for separate setting.
  xlab = "Sample Index", 
  ylab = "Procrustes residual"
)

# Add title.
title(
  main = "Procrustes Errors", 
  outer = TRUE, 
  cex.main = 3, 
  font.main = 2, # Bold face.
  line = -2.9 # Adjust title's vertical position.
) 

# Add subtitle.
title(
  main = "Fish - Invertebrate", 
  outer = FALSE, 
  cex.main = 2, 
  font.main = 1,
  line = -2.9
) 

dev.off() # Save the plot.
