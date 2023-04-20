#
# This script creates stacked bar charts of clonevol (sub)clone proportions.
# I'm using them as a starting point in drawing the Jellyfish plots.
#
# Author: Kari Lavikka
#

library(tidyverse)
library(svglite)

# Takes an array of indices indicating the parent node of each element
# Returns the elements (their indices) in depth-first order.
dfs <- function(parents) {
  order <- integer()
  
  do_dfs <- function(parents, current) {
    order <<- c(order, current)
    for (child in which(parents == current)) {
      do_dfs(parents, child)
    }
  }
  
  do_dfs(parents, 1)
  
  order
}

selected_trees <- read_tsv("~/XXX/mutTree_selected_models_20210311.csv")

for (i in seq_len(nrow(selected_trees))) {
  patient <- selected_trees$patient[i]
  model <- selected_trees$model[i]
  
  print(paste(i, patient))
  
  dir <- str_glue("~/YYY/mutTrees/{patient}_v2/")
  y_dir_name <- list.files(dir, "vaf_*")
  
  y <- readRDS(list.files(file.path(dir, y_dir_name), "*_y.rds", full.names = T))
  
  # Tree contains all clusters and their clonal ordering
  tree <- y$matched$merged.trees[[model]]
  
  expanded_parents <- rep(NA, max(as.integer(tree$lab)))
  expanded_parents[as.integer(tree$lab)] <- as.integer(tree$parent)
  
  # Calculate depth-first order for the clusters
  dfs_order <- data.frame(lab = as.character(seq_along(expanded_parents)),
                          dfs.order = NA)
  d <- dfs(expanded_parents)
  dfs_order$dfs.order[d] <- seq_along(d)
  
  # Add a column indicates the depth-first order
  tree <- tree %>%
    right_join(dfs_order)
  
  # Split the cluster table into separate sample-specific fractions
  all <- matrix(ncol = 4)
  colnames(all) <- c("cluster", "sample", "lower", "upper")
  for (cluster in tree$lab) {
    fracs <- tree %>%
      filter(lab == cluster) %>%
      pull(sample.with.cell.frac.ci)
  
    samples <- str_split(fracs, ",")[[1]]
    print(samples)
    matched <- cbind(cluster, 
                     matrix(str_match(samples, "_([A-Za-z0-9_]+).* : (-?[0-9.]+)-([0-9.]+)")[,(2:4)], ncol=3))
    
    all <- rbind(all, matched)
  }
  
  # Join the sample-specific clusters to the tree so that we have the depth-first
  # order for all samples.
  # The order is always exactly the same, but not all clusters all present.
  all_df <- as.data.frame(all) %>%
    filter(!is.na(cluster)) %>%
    mutate(sample = str_replace(sample, "_DNA[0-9]", ""),
           lower = as.numeric(lower) / 100,
           upper = as.numeric(upper) / 100,
           frac = (lower + upper) / 2) %>%
    # Filter out clusters that are too small
    filter(frac > 0.02) %>%
    inner_join(tree %>% transmute(cluster = lab, dfs.order)) 
  
  
  colors <- tree$color
  names(colors) <- tree$lab
  
  p <- ggplot(all_df, aes(x = sample, y = frac, fill = cluster, group = dfs.order)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = colors) +
    theme_classic()
  ggsave(p, file = paste0(patient, ".svg"))
}
