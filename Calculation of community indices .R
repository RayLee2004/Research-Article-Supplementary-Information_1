# Read in dataset (absolute abundance table of each taxa)
load("datasets/merged_datasets.Rdata")

# Load required packages
library(vegan)
library(tidyr)
library(dplyr)
library(codyn)
library(openxlsx)
library(tibble)

### Part 1: Calculate indices for each taxon ####
# Create taxa collection (number indicates sampling order)
all_groups <- list(
  algae = merged_dataset_algae[, (7:61)],
  bacteria = merged_dataset_bacteria[, (7:61)],
  invertebrate = merged_dataset_benthic_invertebrate[, (7:61)],
  fish = merged_dataset_fish[, (8:62)],
  fungi = merged_dataset_fungi[, (7:61)],
  insect = merged_dataset_insect[, (7:61)],
  zooplankton = merged_dataset_zooplankton[, (7:61)]
)

all_groups <- lapply(all_groups, function(df) {
  df<-as.data.frame(df)
  colnames(df)[1] <- "species"
  return(df)
})

rm(list = setdiff(ls(), "all_groups"))

# Function to transpose abundance tables
normalize <- function(x) {
  species <- as.data.frame(x[, 1])
  colnames(species)[1] <- "species"
  x <- x[, -1]
  col_sums <- colSums(x)
  x <- sweep(x, 2, col_sums, FUN = "/")
  x <- cbind(species, x)
}

all_groups_normalized <- lapply(all_groups, normalize)  # ('an' means all normalized)

# Apply transposition
transpose <- function(y) {
  current_df <- as.data.frame(y)
  col_names <- current_df[, 1]
  current_df <- t(current_df[, -1])
  colnames(current_df) <- col_names
  as.data.frame(current_df)
}

all_groups_transposed <- lapply(all_groups, transpose)  # ('at' means all transposed)

# Calculate richness
Calc_richness <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  for (y in names(x)) {
    richness_values <- specnumber(x[[y]])
    temp_df <- data.frame(
      richness = as.numeric(richness_values)
    )
    colnames(temp_df)[1] <- substitute(y)
    result_df <- cbind(result_df, temp_df)
  }
  return(result_df)
}

richness_single_groups <- Calc_richness(all_groups_transposed)

# Calculate diversity
Calc_diversity <- function(x) {
  result_df <- data.frame(sampling_site_order = rownames(x[[1]]))
  for (y in names(x)) {
    diversity_values <- diversity(x[[y]])
    temp_df <- data.frame(
      diversity = as.numeric(diversity_values)
    )
    colnames(temp_df)[1] <- substitute(y)
    result_df <- cbind(result_df, temp_df)
  }
  return(result_df)
}

diversity_single_groups <- Calc_diversity(all_groups_transposed)

convert_to_longdata <- function(df) {
  df %>%
    pivot_longer(
      cols = -1,
      names_to = c("sampling_site", "sampling_order"),
      names_sep = "_",
      values_to = "abundance"
    ) %>%
    rename(species = 1) %>%
    mutate(sampling_order = as.integer(sampling_order))
}

all_groups_longdata <- lapply(all_groups, convert_to_longdata)

# Calculate synchrony
Calc_synchrony <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    synchrony_values <- synchrony(
      x[[y]],
      species.var = "species",
      time.var = "sampling_order",
      abundance.var = "abundance",
      replicate.var = "sampling_site"
    )
    colnames(synchrony_values)[2] <- substitute(y)
    result_df <- cbind(result_df, synchrony_values[2])
  }
  return(result_df)
}

synchrony_single_groups <- Calc_synchrony(all_groups_longdata)

# Calculate stability
Calc_stability <- function(x) {
  result_df <- data.frame(sampling_site = (x[[1]][(1:18), 2]))
  for (y in names(x)) {
    stability_values <- community_stability(
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

stability_single_groups <- Calc_stability(all_groups_longdata)

### Part 2: Calculate indices for multitaxon ####
# Merge normalized data
multi_groups <- rbind(
  all_groups$fish,
  all_groups$zooplankton,
  all_groups$invertebrate,
  all_groups$insect,
  all_groups$algae
)

multi_groups_transposed <- transpose(multi_groups)
multi_groups_longdata <- convert_to_longdata(multi_groups)

# Calculate indices
richness_multi_groups <- data.frame(multi_groups = specnumber(multi_groups_transposed))
diversity_multi_groups <- data.frame(multi_groups = diversity(multi_groups_transposed))
synchrony_multi_groups <- Calc_synchrony(list(multi_groups = multi_groups_longdata))
stability_multi_groups <- Calc_stability(list(multi_groups = multi_groups_longdata))

# Consolidate results
richness_collection <- cbind(richness_single_groups, richness_multi_groups)
diversity_collection <- cbind(diversity_single_groups, diversity_multi_groups)
synchrony_collection <- cbind(synchrony_single_groups, multi_groups = synchrony_multi_groups[, 2])
stability_collection <- cbind(stability_single_groups, multi_groups = stability_multi_groups[, 2])

community_indices_collection <- list(
  richness = richness_collection,
  diversity = diversity_collection,
  synchrony = synchrony_collection,
  stability = stability_collection
)

write.xlsx(community_indices_collection, "resultsets/community indices.xlsx")