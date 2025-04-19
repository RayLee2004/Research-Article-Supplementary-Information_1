# R version: R-4.5.0
########## PART 1 COMMUNITY INDICES CALCULATION #########

# Read in merged - dataset.
load("datasets/merged_datasets.Rdata")

# Load required packages.
library(vegan)
library(tidyr)
library(dplyr)
library(codyn)
library(openxlsx)
library(tibble)

# Create a 'collection' containing all taxonomic groups' datasets.
all_groups <- list(
  algae = merged_dataset_algae[, (7:61)],
  bacteria = merged_dataset_bacteria[, (7:61)],
  invertebrate = merged_dataset_benthic_invertebrate[, (7:61)],
  fish = merged_dataset_fish[, (8:62)],
  fungi = merged_dataset_fungi[, (7:61)],
  insect = merged_dataset_insect[, (7:61)],
  zooplankton = merged_dataset_zooplankton[, (7:61)]
) # Remove the needless columns.

# Unify the name of first column of all tables.
all_groups <- lapply(all_groups, function(df) {
  df <- as.data.frame(df)
  colnames(df)[1] <- "species"
  return(df)
})

# Remove redundant objects in R environment.
rm(list = setdiff(ls(), "all_groups"))

# Create a function to convert absolute abundance to relative abundance (data standardization).
standardize <- function(x) {
  species <- as.data.frame(x[, 1])
  colnames(species)[1] <- "species"
  x <- x[, -1]
  col_sums <- colSums(x)
  x <- sweep(x, 2, col_sums, FUN = "/")
  x <- cbind(species, x)
  return(x)
}

# Execute the function.
all_groups_standardized <- lapply(all_groups, standardize)  # ('an' means 'all normalized')

# Create a function to transpose tables to meet computational format.
transpose <- function(x) {
  current_df <- as.data.frame(x)
  col_names <- current_df[, 1]
  current_df <- t(current_df[, -1])
  current_df<-as.data.frame(current_df)
  colnames(current_df) <- col_names
  return(current_df)
}

# Execute the function.
all_groups_transposed <- lapply(all_groups, transpose)  # ('at' means 'all transposed')

###### Community Indices For Single Taxonomic Groups ######

# Create a function to calculate community richness. The output is a data.frame.
richness <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]])) # Read sampling sites and sampling time (or order).
  for (y in names(x)) {
    richness_values <- vegan::specnumber(x[[y]]) # If the package is not marked, the function will not run normally because the defined function's name has existed in R packages.
    temp_df <- data.frame(
      richness = as.numeric(richness_values)
    ) # Calculate by traversing each table in the list.
    colnames(temp_df)[1] <- substitute(y)
    result_df <- cbind(result_df, temp_df)
  } # Merge the sample name list with the richness value column.
  return(result_df)
}

# Execute the function.
richness_single_groups <- richness(all_groups_transposed)

# Create a function to calculate community diversity.
diversity <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  for (y in names(x)) {
    diversity_values <- vegan::diversity(x[[y]])
    temp_df <- data.frame(
      diversity = as.numeric(diversity_values)
    )
    result_df <- cbind(result_df, temp_df)
  }
  return(result_df)
}

# Execute the function.
diversity_single_groups <- diversity(all_groups_transposed)

# Create a function to convert all tables to long data to meet computational format.
convert_to_longdata <- function(df) {
  df %>%
    pivot_longer(
      cols = -1, # The first column remains unchanged.
      names_to = c("sampling_site", "sampling_order"),  # Next step requires these two columns.
      names_sep = "_", # Read sampling site from 'sampling site + order".
      values_to = "abundance" 
    ) %>%
    rename(species = 1) %>% # Unify the name of first column of all long data tables.
    mutate(sampling_order = as.integer(sampling_order))
}

# Execute the function.
all_groups_longdata <- lapply(all_groups, convert_to_longdata)

# Create a function to calculate community synchrony.
synchrony <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    synchrony_values <- codyn::synchrony(
      x[[y]],
      species.var = "species",
      time.var = "sampling_order",
      abundance.var = "abundance",
      replicate.var = "sampling_site"
    ) # Correspond the parameters to columns of the long data table.
    colnames(synchrony_values)[2] <- substitute(y)
    result_df <- cbind(result_df, synchrony_values[2])
  }
  return(result_df)
}

# Execute the function.
synchrony_single_groups <- synchrony(all_groups_longdata)

# Create a function to calculate community stability.
stability <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    stability_values <- codyn::community_stability(
      x[[y]],
      time.var = "sampling_order",
      abundance.var = "abundance",
      replicate.var = "sampling_site"
    )
    colnames(stability_values)[2] <- substitute(y)
    result_df <- cbind(result_df, stability_values[2])
  }
  return(result_df)
}

# Execute the function.
stability_single_groups <- stability(all_groups_longdata)

###### Community Indices For Multiple Taxonomic Groups ######

# Merge datasets of each groups.
multi_groups <- rbind(
  all_groups$fish,
  all_groups$zooplankton,
  all_groups$invertebrate,
  all_groups$insect,
  all_groups$algae
)

# Standardize and transpose the dataset of multiple groups.
multi_groups_transposed <- transpose(multi_groups)
multi_groups_longdata <- convert_to_longdata(multi_groups)

# Calculate the four community indices like before.
richness_multi_groups <- data.frame(multi_groups = vegan::specnumber(multi_groups_transposed))
diversity_multi_groups <- data.frame(multi_groups = vegan::diversity(multi_groups_transposed))
synchrony_multi_groups <- synchrony(list(multi_groups = multi_groups_longdata))
stability_multi_groups <- stability(list(multi_groups = multi_groups_longdata))

# Now, we should have 8 result tables, respectively from the single groups and multi groups.
# Next, merge all the results into one table and output/save it.
collection_richness <- cbind(richness_single_groups, richness_multi_groups)
collection_diversity <- cbind(diversity_single_groups, diversity_multi_groups)
collection_synchrony <- cbind(synchrony_single_groups, multi_groups = synchrony_multi_groups[, 2])
collection_stability <- cbind(stability_single_groups, multi_groups = stability_multi_groups[, 2])
community_indices_collection <- list(
  richness = collection_richness,
  diversity = collection_diversity,
  synchrony = collection_synchrony,
  stability = collection_stability
)

# Each sheet stores the result set of one index.
write.xlsx(community_indices_collection, "Community Indices.xlsx")

########## PART 2 NETWORK INDICES CALCULATION ##########
###### WARNING: WILL CLEAR ALL OBJECTS IN R ENVIRONMENT
rm(list = ls())

# Load required packages
library(openxlsx)
library(igraph)
library(bipartite)
library(tibble)
library(purrr)

# Read the species occurrence table for the nodes of each subnet.
load("datasets/species_occurrences.Rdata")

# Store them in a list for subsequent network processing.
all_objects <- ls()
occurrences <- list()
for (x in all_objects) {
  if (is.data.frame(get(x))) {
    occurrences[[x]] <- get(x)
  }
}
occurrences <- lapply(occurrences, as.data.frame)

# Delete redundant tables in the list.
occurrences <- occurrences[!names(occurrences) %in% c("sample2022", "sample2023", "sample2024")]

# Remove redundant objects in R environment.
rm(list = setdiff(ls(), "occurrences"))

# Load total network which includes all species in 3 years.
total_network <- read.xlsx("datasets/trophic interaction 0-1 adjacent matrix.xlsx")

original_network <- total_network[,-(1:2)]
original_network <- as.matrix(original_network)
diag(original_network) <- 0 # Ignore intraspecific self-feeding.
original_network <- as.data.frame(original_network)
rownames(original_network) <- colnames(original_network)

# Create three lists to store dfs(data.frames), matrices and graphs separately
dfs <- list()
matrices <- list()
graphs <- list()

# Define a function to get subnets from original network by occurrence lists.
get_df <- function(x) {
  occurrence_list <- x[,ncol(x)]
  current_df <- original_network[occurrence_list, occurrence_list] 
  rownames(current_df) <- colnames(current_df)
  dfs <<- c(dfs, list(current_df))
}

# Execute the function.
lapply(occurrences, get_df)

# Define a function to convert df to matrix and graph format.
get_matrix_graph <- function(x) {
  matrix <- as.matrix(x) # Transform the data frame into a matrix.
  graph <- graph_from_adjacency_matrix(
    matrix,
    mode = "directed",
    weighted = TRUE,
    diag = FALSE
  ) # Transform the matrix into a graph.
  isolated_nodes <- V(graph)[degree(x) == 0] # Identify isolated nodes.
  graph <- delete_vertices(graph, isolated_nodes) # Delete isolated nodes.
  matrices <<- c(matrices, list(matrix)) # Store the matrix into list.
  graphs <<- c(graphs, list(graph)) # Store the graph into list.
}

# Execute the function.
lapply(dfs, get_matrix_graph)

# Name each subnet.
names(dfs) <- names(occurrences)
names(matrices) <- names(occurrences)
names(graphs) <- names(occurrences)

# Define a function to calculate net-chain indices simultaneously.
chain_indices <- function(x) {
  result_df <- data.frame(
    sampling_site_order = character(0),
    connectance = numeric(0),
    max_path_length = numeric(0),
    mean_path_length = numeric(0),
    no_nodes = numeric(0)
  ) # Specify the format of the result tables.
  
  for (y in names(x)) {
    current_graph <- x[[y]]
    metrics <- data.frame(
      sampling_site_order = y,
      connectance = edge_density(current_graph),
      max_path_length = diameter(current_graph, directed = TRUE),
      mean_path_length = mean_distance(current_graph, directed = TRUE),
      no_nodes = vcount(current_graph)
    ) # Traverse each subnet, calculate and output the result values.
    
    result_df <- rbind(result_df, metrics)
  }
  return(result_df)
}

# Execute the function.
chain_indices <- chain_indices(graphs)

# Define a function to calculate network-structure indices simultaneously.
# a matrix and a graph need to be input.
structure_indices <- function(x, y) { 
  # Section 1 Modularity
  result_df1 <- data.frame(
    sampling_site_order = character(0),
    modularity = numeric(0)
  ) # table to store modularity values.
  for (z in names(y)) {
    current_community <- cluster_walktrap(y[[z]])
    mod_df <- data.frame(
      sampling_site_order = z,
      modularity = modularity(current_community)
    )
    result_df1 <- rbind(result_df1, mod_df)
  }
  
  # Section 2 Vulnerability,nestedness,robustness.
  result_df2 <- data.frame(
    nestedness = numeric(0),
    vulnerability = numeric(0),
    robustness = numeric(0)
  ) # the second table to store other indices values.
  vulnerability <- function(n) {
    gen_vul <- networklevel(n, index = "vulnerability")
    gen_vul[!grepl("generality", names(gen_vul))]
  }
  for (w in names(x)) {
    current_matrix <- x[[w]]
    net_metrics <- data.frame(
      nestedness = nested(current_matrix, "NODF"),
      vulnerability = vulnerability(current_matrix),
      robustness = robustness(
        second.extinct(
          current_matrix,
          participant = "higher",
          method = "random",
          nrep = 10, # nrep > 50 is better.
          details = FALSE
        )
      )
    )
    result_df2 <- rbind(result_df2, net_metrics)
  }
  return(cbind(result_df1, result_df2))
}

# Execute the function.
network_indices <- structure_indices(matrices, graphs)

# Merge result tables into one table and output.
collection_network_indices <- cbind(chain_indices, network_indices[, -1])
write.xlsx(collection_network_indices, "Network Indices.xlsx")
