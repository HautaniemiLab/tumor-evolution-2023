# Finds copy number statuses for all genes
# Written by Kari, Yilin, Amjad

segments <- table1
estimates <- table2
genes <- table3

ploidy <- estimates$ploidy

suppressMessages({
  library(GenomicRanges)
  library(dplyr)
  library(magrittr)
})

if (!(TRUE %in% estimates$aberrant)) {
  # Don't output anything
  genes <- genes %>%
    slice_head(n = 0)
}

# Ignore chromosomes Y and MT
genes <- genes %>%
  filter(chr %in% c(1:22, "X")) %>%
  mutate(chr = paste0("chr", chr))

segs_gr = with(segments, GRanges(chr, IRanges(startpos, endpos)))
genes_gr = with(genes, GRanges(chr, IRanges(start, end)))

overlaps = findOverlaps(genes_gr, segs_gr)

genes_and_segments <- data.frame(
  genes[from(overlaps), ],
  segments[to(overlaps), ] %>% select(-chr, -startpos, -endpos),
  width = width(pintersect(genes_gr[from(overlaps)], segs_gr[to(overlaps)])))

# Compute the number of breaks and the max intra-gene copy-number differences
breaks <- genes_and_segments %>%
  group_by(ID) %>%
  summarise(minPurifiedLogR = min(purifiedLogR),
	    maxPurifiedLogR = max(purifiedLogR),
            breaksInGene = n() - 1) %>%
  ungroup()

table.out <- genes_and_segments %>%
  group_by(ID) %>%
  # Find the largest overlap for each gene and compute the final copy-number status based on that
  top_n(1, width) %>%
  ungroup() %>%
  select(-width) %>%
  mutate(
    # Determine if amplification, deletion or normal copy number
    CNstatus = ifelse(
      nAraw + nBraw > ploidy * 2,
      "AMP",
      ifelse(nAraw + nBraw < ploidy * 0.5, "DEL", "Normal")
    ),
    # Determine if Loss of Heterozygosity or Heterozygous. TODO: This is trivial, consider removal
    LOHstatus = ifelse(nMinor == 0, "LOH", "HET")) %>%
  inner_join(breaks) %>%
  arrange(as.numeric(gsub("chr", "", gsub("X", "23", chr))), start)

