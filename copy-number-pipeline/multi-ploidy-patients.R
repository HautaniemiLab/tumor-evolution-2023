#
# Plot ploidy statistics of patients that have samples with heterogeneous ploidies
#

library(magrittr)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(readxl)

conclusionLevels <- c("default", "adjust", "discard")
conclusionColors <- c("black", "#00e050", "#F04030")

args <- commandArgs(trailingOnly = TRUE)
allCandidatesFilename <- args[1]
ascatParameterOverridesFilename <- args[2]
summaryOutputFilename <- args[3]
patientDir <- args[4]

allCandidates <- read_tsv(allCandidatesFilename) %>%
  filter(aberrant & !is.na(ploidy)) %>%
  mutate(patient = str_extract(sample, "^[^_]+"))

aberrantSamples <- allCandidates %>%
  filter(rank == 1) %>%
  select(sample, patient, ploidy, penalizedGoodnessOfFit) %>%
  left_join(read_xlsx(ascatParameterOverridesFilename) %>% rename(adjustedPloidy = ploidy)) %>%
  mutate(conclusion = factor(
      ifelse(!is.na(discardModels), "discard", ifelse(!is.na(adjustedPloidy), "adjust", "default")),
      levels = conclusionLevels))

find_peaks <- function(d) d$x[c(F, diff(diff(d$y) >= 0) < 0, F) & d$y > 0.01]
find_ploidy_peaks <- function(ploidy) find_peaks(density(ploidy, from = 0, to = 8, bw = 0.2))
find_ploidy_peak_count <- function(ploidy) length(find_ploidy_peaks(ploidy))

patientsWithMultiplePloidies <- aberrantSamples %>%
  group_by(patient) %>%
  summarise(ploidyPeaks = find_ploidy_peak_count(ploidy),
            n = n()) %>%
  ungroup() %>%
  filter(ploidyPeaks > 1)

p <- aberrantSamples %>%
  mutate(penalGoodOfFit = cut(penalizedGoodnessOfFit, seq(60, 100, 10), include.lowest = TRUE)) %>%
  filter(patient %in% patientsWithMultiplePloidies$patient) %>%
  ggplot(aes(x = ploidy)) +
  geom_dotplot(aes(fill = penalGoodOfFit, colour = conclusion),
               binwidth = 0.4, stackgroups = TRUE, binpositions = "all") +
  xlim(c(0, 8)) +
  facet_wrap(vars(patient)) +
  scale_fill_brewer() +
  scale_colour_manual(values = conclusionColors) +
  labs(fill = "Penalized\ngoodness\nof fit") +
  ggtitle("Patients with more than one ploidy solution in their aberrant samples (uncurated ploidies shown)")

ggsave(summaryOutputFilename, plot = p, width = 40, height = 25, units = "cm")


for (patientId in unique(patientsWithMultiplePloidies$patient)) {
  p <- allCandidates %>% 
    filter(patient == patientId) %>%
    group_by(sample) %>%
    slice_max(penalizedGoodnessOfFit, n = 3) %>%
    ungroup() %>%
    mutate(penalGoodOfFit = cut(penalizedGoodnessOfFit, seq(60, 100, 10), include.lowest = TRUE)) %>%
    ggplot() +
    geom_point(aes(x = ploidy,
                   y = sample,
                   fill = penalGoodOfFit,
                   size = rank == 1),
               shape = 21) +
    geom_point(aes(x = adjustedPloidy, y = sample),
               data = aberrantSamples %>%
                  filter(patient == patientId & !is.na(adjustedPloidy)),
               pch = 23,
               size = 3,
               fill = "#00e040",
               alpha = 0.7) +
    xlim(c(1, 8)) +
    scale_size_manual(values = c(4, 7)) +
    scale_fill_brewer(drop = F) +
    labs(fill = "Penalized\ngoodness\nof fit") +
    ggtitle(patientId)
  
  ggsave(file.path(patientDir, paste0(patientId, ".png")), width = 9, height = 6, dpi = 100)
}
