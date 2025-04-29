# Read fish abundance dataset.
load("datasets/merged datasets.Rdata")

# Remove redundant objects in R environment.
rm(list = setdiff(ls(), "fish"))

# Load required packages.
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

# Split the dataset.
fish2022 <- cbind(fish$species, fish[, (2:19)])  # 2022 dataset
fish2023 <- cbind(fish$species, fish[, (20:37)])  # 2023 dataset
fish2024 <- cbind(fish$species, fish[, (38:55)])  # 2024 dataset

# Delete Sampling time number.
colnames(fish2022) <- sub("\\_[1]", "", colnames(fish2022))
colnames(fish2023) <- sub("\\_[2]", "", colnames(fish2023))
colnames(fish2024) <- sub("\\_[3]", "", colnames(fish2024))

# Only one of the dominant fish species is plotted here as an example. 
# Identify top20 fishes.
fish2024$total_abundance <- rowSums(fish2024[, -1])
top20 <- fish2024[order(-fish2024$total_abundance), ][1:20, ]

# Merge the remaining species.
others <- fish2024[!(fish2024$species %in% top20$species), ]
others_sum <- colSums(others[, -c(1, ncol(others) - 1)])
others_row <- c("others", as.numeric(others_sum), sum(others_sum))

# Merge two tables.
final_table <- rbind(top20[, -ncol(top20) - 1], others_row)[, (1:19)]
colnames(final_table)[1] <- 'species'

# Convert table to long data to meet plotting demands.
longdata <- final_table %>%
  pivot_longer(
    cols = -species,
    names_to = "sampling_site",
    values_to = "abundance"
  )

# The format and size of plot output can be modified here.
png('Top 20 Fish Species in 2024.png', width = 2000, height = 1000)

### 1 ggplot() module
ggplot(longdata) +
  
### 2 geom_bar() module (to specify the plotting type)
  geom_bar(
    aes(
      x = sampling_site,
      y = as.numeric(abundance),
      fill = species
    ),
    stat = "identity",
    position = "fill",  # Use relative value as column height.
    width = 0.55,       # Column width.
    show.legend = TRUE  # Show legend.
  ) +
  
### 3 labs() module
  labs(
    x = "Sampling Site",
    y = "Relative Abundance",
    fill = "Species",  # Specify colors based on species.
    title = "Top 20 Fish Species in 2024"
  ) +
  
### 4 Specify the colors of species.
  scale_fill_manual(
    values = c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
      "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
      "peachpuff"
    )
  ) +
  
### 5 theme() module
  theme(
    ## Figure margin section
    plot.margin = margin(t = 30, r = 30, b = 20, l = 50),
    
    ## Figure background section
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank(),      # Hide axis line. 
    axis.ticks = element_blank(),     # Hide axis tick.
    
    ## Title section 
    plot.title = element_text(
      size = 40,
      face = "bold",
      hjust = 0.5,
      vjust = 0.8  # Positive number means up.
    ),
    
    ## Set the position of the title globally.
    plot.title.position = "plot",
    
    ## Legend section
    legend.position = "right",
    legend.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0,
      vjust = 3.5
    ),
    legend.text = element_text(
      size = 16.5,
      margin = margin(l = 10)
    ),  # Text shifts to the right.
    legend.key.width = unit(0.85, "cm"),
    legend.key.height = unit(1.27, "cm"),
    legend.spacing = unit(0.2, "cm"),  # Vertical distance of the legend.
    
    ## Axis text section
    axis.text.x = element_text(
      size = 19,
      vjust = 8,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 19,
      margin = margin(r = 0),
      face = "bold"
    ),  # Use 'axis.text' to set X and Y axis simultaneously.
    
    ## Axis lab section
    axis.title.x = element_text(
      size = 29,
      vjust = 3,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 29,
      vjust = 4.5,
      face = "bold"
    )   # Use 'axis.text' to set X and Y axis simultaneously.
  ) +
  
### 6 Set legend parameters.  
  guides(
    fill = guide_legend(
      override.aes = list(shape = 16),  # Legend shape.
      ncol = 1 # Keep the legend in one column.
    )
  )

dev.off()  # Save the plot.
