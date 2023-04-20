#!/usr/bin/env Rscript

suppressMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(tibble)
  library(magrittr)
  library(stringr)
  library(optparse)
})
  
# Relative location from component directories.
# This is ugly. TODO: better way to refer ascat-algorithm.R
source("../../ascat-algorithm.R")

# #################

args <- OptionParser(
  usage = "%prog [options] segment_model_file",
  option_list = list(
    make_option("--snp-purity", type="character",
                help="Optional truncal SNP (TP53, for example) for for additional purity evidence. Format: chrom,pos,altCount,refCount"),
    make_option("--sample", type="character", default="sample", 
                help="The sample identifier"),
    make_option("--segment-output", type="character",
                help="Output tsv file for the segments"),
    make_option("--estimate-output", type="character",
                help="Output tsv file for purity, ploidy, and goodness of fit"),
    make_option("--sunrise-output", type="character",
                help="Output file for the sunrise plot"))
  ) %>% parse_args(positional_arguments = 1)


# TODO: Add some command-line option validation

segmentsFile <- args$args
options <- args$options
sampleId <- options$sample

if (!is.null(options$`snp-purity`)) {
  parts <- str_split(options$`snp-purity`, ",", simplify = TRUE)
  SNP <- list(chrom = parts[, 1],
              pos = as.integer(parts[, 2]),
              altCount = as.integer(parts[, 3]),
              refCount = as.integer(parts[, 4]))
} else {
  SNP <- NULL
}

###### 1. Read the GATK segment model #####

# Use only chr1-chrX
chroms <- paste0("chr", c(1:22, "X"))

modeledSegments <- read_tsv(segmentsFile, comment = "@", na = c("NA", "NaN")) %>%
  filter(chr %in% chroms)

###### 2. Run the ASCAT algorithm #####

result <- ASCAT(modeledSegments, SNP)

###### 3. Write the results

if (!is.null(options$`segment-output`)) {
  result$segments %>%
    mutate(sample = sampleId) %>%
    relocate(sample) %>%
    rename(baf = BAF) %>%
    mutate_if(is.numeric, round, digits = 4) %>%
    write_tsv(options$`segment-output`)
}

if (!is.null(options$`estimate-output`)) {
  candidates <- result$candidates
  
  # Add a single row with NAs.
  if (nrow(candidates) < 1) {
    candidates <- candidates %>%
      add_row(rank = 1)
  }
  
  # Note: We include all candidates here. The top candidate can be picked afterwards
  candidates %>%
    mutate(sample = sampleId,
           aberrant = result$aberrant) %>%
    relocate(sample) %>%
    mutate_if(is.numeric, round, digits = 4) %>%
    write_tsv(options$`estimate-output`)
}


if (!is.null(options$`sunrise-output`)) {
  sunrise <- result$sunrise
  
  if (sampleId != "sample") {
    sunrise <- sunrise + ggtitle(sampleId)
  }
  
  ggsave(options$`sunrise-output`,
         plot = sunrise,
         width = 7, height = 7, dpi = 100)
}

