# Display Procrustes analysis and plot for only one group here as an example.
# Other groups can be analyzed likely after modifying the imported dataset.

# Load required packages
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)

########## PART 1 PROCRUSTES ANALYZING #########

# Read in merged dataset.
load('datasets/merged datasets.Rdata')

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
prot <- protest(fish_pca, invertebrate_pca, permutations = 999) # Fish as group1.

# Print the M2 statistic.
print(prot$ss)

# Print the significance (p value).
print(prot$signif)

########## PART 2 ANALYSIS RESULTS PLOTING #########

##### Figure 1: Procrustes Errors #####

# Remove redundant objects in R environment.
rm(list = setdiff(ls(), list('proc', 'prot')))

# The format and size of plot output can be modified here.
png('Procrustes Errors.png', width = 2000, height = 1000)

# Plot margin (bottom→left→top→right)
par(mar = c(5, 6, 5, 6))
# Distance between axis.lab/axis.text and axis.
par(mgp = c(4, 1.4, 0))

plot(
  proc,
  kind = 2,
  cex.axis = 1.75,
  cex.lab = 2.25,
  font.lab = 2,
  main = '',  # Leave the title blank for separate setting.
  xlab = 'Sample Index',
  ylab = 'Procrustes Residual'
)

# Add title.
title(
  main = 'Procrustes Errors',
  outer = TRUE,
  cex.main = 3,
  font.main = 2, # Bold face.
  line = -3 # Adjust title's vertical position. (Positive number means up)
)

# Add subtitle.
mtext(
  'Fish - Invertebrate',
  outer = FALSE,  # Place the subtitle within the figure.
  cex = 2,
  font = 1,
  col = 'red',
  line = -2.5,
  adj = 0.5 # Place the subtitle on the mediate side.
)

dev.off() # Save the plot.

##### Figure 2 Procrustes Correlation #####

# Create a data.frame used for plotting including Begin and End Coordinates from analysis result.
plotdf <- data.frame(
  sample = rownames(proc$X), # Sample name
  X_begin = proc$X[, 1], # group 1 x
  X_end = proc$Yrot[, 1], # group 2 x
  Y_begin = proc$X[, 2], # group 1 y
  Y_end = proc$Yrot[, 2] # group 2 y
)

# Create a data.frame to store group names for adding legend in figure.
plotdf_group <- data.frame(
  group = c('Fish', 'Invertebrate')
)

# The format and size of plot output can be modified here.
png('Procrustes Analysis.png', width = 1500, height = 1200)

### 1 ggplot() part (to import plotdf)
ggplot(plotdf) +
  
### 2 geom_segment() part (to specify the ploting object as a line segment)
  geom_segment(
    aes(
      x = X_begin,
      y = Y_begin,
      xend = X_end,
      yend = Y_end
    ), # Map the coordinates to the figure
    linewidth = 0.5, # Line width
    color = 'black' # Line color
  ) +
  
### 3 lab() part (to define axis.labs and figure title)
  labs(
    x = 'Dimension 1',
    y = 'Dimension 2',
    title = 'Procrustes Analysis',
    color = 'Group'
  ) + # Set the legend name 'Group'.
  
### 4 Specify the corresponding colors for each group's name.
  scale_color_manual(
    values = c('red', 'blue'),  # (red for 'Fish', blue for 'Invertebrate')
    labels = plotdf_group$group
  ) +
  
### 5 geom_point() part (to add the points at both ends of the line segment)
  geom_point(aes(x = X_begin, y = Y_begin, color = 'Fish'), size = 5) + # points of group 1
  # Unify the color of the points to the color corresponding to the string "Fish".
  geom_point(aes(x = X_end, y = Y_end, color = 'Invertebrate'), size = 5) + # points of group 2
  # Unify the color of the points to the color corresponding to the string "Invertebrate".
  
### 6 geom_text() part (to add point sample name)
  geom_text(
    aes(x = X_begin, y = Y_begin, label = sample), # show the point name of group 1
    color = 'black',
    size = 6.3,
    vjust = -1
  ) + # Positive number means up
  
### 7 Add annotation of procrutes M2 statistic and p value.
  annotate(
    'text',
    x = 0.33, # Place annotation on a specific coordinate. 
    y = 0.38, # Here you should be checked and modified this.
    label = paste0('M2 = ', round(proc$ss, 3), 'P value = ', prot$signif),  # Read M2 and p value and paste.
    size = 8,
    color = 'black'
  ) +
  
### 8 theme() part (to further polish)
  theme(
    ## Figure margin section
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = 'pt'),
    
    ## Figure background section
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    ## Legend section
    legend.title = element_text(
      size = 27,
      vjust = 2
    ), # Legend title
    legend.text = element_text(
      size = 24
    ), # Legend text
    legend.key.size = unit(1, 'cm'), # circle's size
    legend.spacing = unit(0.6, 'cm'), # distance between adjacent legends
    
    ## Title section
    plot.title = element_text(
      hjust = 0.6, # Adjust title's horizontal position.
      vjust = 0.1, # Adjust title's vertical position. (positive number means up)
      size = 48,
      face = 'bold', # bold face
      margin = margin(b = 15)
    ),
    
    ## Axis lab section
    axis.title = element_text(
      size = 30, # axis.lab size
      face = 'bold', # axis.lab bold face
      vjust = 2
    ),
    
    ## Axis line section
    axis.line = element_line(
      color = 'black',
      size = 0.9
    ),
    
    ## Axis text section
    axis.text = element_text(
      face = 'bold',
      size = 20,
      vjust = -2.5
    )
  )

dev.off() # Save the plot.
