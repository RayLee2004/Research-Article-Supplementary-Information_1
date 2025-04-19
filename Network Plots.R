# Only one network is plot here as an example. 
# Other networks can be drawn by modifying the imported 'edges and nodes' file.

# Load required packages.
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)

# Read edges(links) file and nodes file, 
# in which species with a relative abundance less than 0.001 were excluded.
edges <- read.csv("datasets/Network Edges & Nodes/edges2024.csv")
nodes <- read.csv("datasets/Network Edges & Nodes/nodes2024.csv")

# Delete the serial number column.
edges <- edges[, -1]
nodes <- nodes[, -1]

# Convert table to graph.
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE) 

# The format and size of plot output can be modified here.
png("2024.png", width = 2000, height = 1800) 

# Specify the color of the nodes of each group.
unique_groups <- unique(nodes$group) # get groups.
unique_groups <- sort(unique_groups) # Sort by the first letter.

group_colors_manual <- c("darkseagreen", "tomato", "gold", "thistle", "lightskyblue") 

node_colors <- sapply(V(graph)$group, function(x) {
  group_colors_manual[which(unique_groups == x)]
})

V(graph)$color <- node_colors # Node color
V(graph)$size <- 4.2 # Node size

E(graph)$color <- adjustcolor("black", alpha.f = 0.3)  # Edge color and alpha 
E(graph)$width <- 0.26  

par(mar = c(1.6, 1.6, 6, 1.6), xaxs = 'i', yaxs = 'i') # Plot margin

plot(
  graph, 
  vertex.label = NA, # Hide node labels
  edge.arrow.size = 0  # Hide arrows
)

title(
  main = "Trophic Interaction Network (2024) ",
  mar = c(1, 1, 1, 1), 
  cex.main = 5.5
) 

legend(
  "topright",          
  legend = unique_groups, 
  col = group_colors_manual,
  cex = 2.7,           
  bty = "none",        
  pch = 16,            
  pt.cex = 2.7        
) # Add legend on topright

dev.off() # Save the plot.
