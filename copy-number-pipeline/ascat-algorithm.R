#
# The ASCAT algorithm rewritten and optimized for data that has been segmented using external tools such as GATK.
# Also uses TP53 VAF as additional evidence for purity in model selection.
#
# TODO: Make thresholds, factors, etc. configurable.
#
# Author: Kari Lavikka
#

# Partially based on:
# ASCAT 2.5.2
# author: Peter Van Loo
# license: GPLv3

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggfx)
  library(binom)
})

#' Eq. S7
calculateNA <- function(logR, BAF, psi, rho, gamma = 1) {
  (rho - 1 - (BAF - 1) * 2 ^ (logR / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
}

#' Eq. S8
calculateNB <- function(logR, BAF, psi, rho, gamma = 1) {
  (rho - 1 + (BAF    ) * 2 ^ (logR / gamma) * ((1 - rho) * 2 + rho * psi)) / rho
}

#' Eq. S7 + S8
calculateN <- function(logR, BAF, psi, rho, gamma = 1) {
  calculateNA(logR, BAF, psi, rho, gamma) + calculateNB(logR, BAF, psi, rho, gamma) 
}

#' Purify R, based on: http://www.nature.com/articles/ng.2760
#' But there's an error, which is fixed here: https://github.com/lima1/PureCN/issues/40
#'
#' @param purity Purity
#' @param ploidy Ploidy of the cancer cells
#' @param R Copy ratio (not logR but 2^logR)
purifyR <- function(purity, ploidy, R) (purity*ploidy*R + 2*(1-purity)*(R-1)) / (purity*ploidy)

#' Purify logR
#' Because logarithm is undefined for negative numbers and not all tools understand
#' -Inf, we clamp it to a very small value (-10).
#'
#' @param purity Purity
#' @param ploidy Ploidy of the cancer cells
#' @param logR Log2 copy ratio
purifyLogR <- function(purity, ploidy, logR) log2(pmax(2^-10, purifyR(purity, ploidy, 2^logR)))

#' Purify baf, based on equations S2, S7, and S8 of
#'
#' @param purity Purity
#' @param ploidy Ploidy of the cancer cells
#' @param baf B-allele frequency, [0, 0.5]
#' @param R Copy ratio (not logR but 2^logR)
purifyBaf <- function(purity, ploidy, baf, R) {
  f <- function(af) purity - 1 + R*af*(2*(1 - purity) + purity*ploidy)
  f(baf) / (f(1 - baf) + f(baf))
}

#' Compute distances for the given psi and rho candidates using the Eq. 3.
computeDistances <- function(segments, psiCandidates, rhoCandidates) {
  crossing(psi = psiCandidates, rho = rhoCandidates) %>%
    rowwise() %>%
    mutate(distance = {
      nA <- calculateNA(segments$logR, segments$BAF, psi, rho)
      nB <- calculateNB(segments$logR, segments$BAF, psi, rho)
      
      # choose the minor allele (why?)
      nMinor = NULL
      if (sum(nA, na.rm=T) < sum(nB, na.rm=T)) {
        nMinor = nA
      } else {
        nMinor = nB
      }
      
      # Note: Eq. 3. in the paper takes both alleles into account. For some reason ASCAT
      # actually uses only the minor allele.
      # We skip the BAF-based weights (w) because they are probably unnecessary with our data
      sum(abs(nMinor - pmax(round(nMinor), 0)) ^ 2 * segments$length, na.rm=T)
    })
}

distancesToMatrix <- function(distances) {
  rho <- sort(unique(distances$rho))
  psi <- sort(unique(distances$psi))
  
  m <- distances %>%
    arrange(psi, rho) %>%
    pull(distance) %>%
    matrix(ncol = length(psi),
           nrow = length(rho))
  
  colnames(m) <- psi 
  rownames(m) <- rho
  m
}

#' Find local minima from a matrix
findValleys <- function(distanceMatrix, radius = 3) {
  d <- distanceMatrix
  r <- radius
  
  minima <- data.frame(i = numeric(), j = numeric(), value = numeric())
  
  for (i in seq_len(ncol(d))) { # psi
    for (j in seq_len(nrow(d))) { # rho
      m <- d[j, i]
      horiz <- (i-r):(i+r)
      vert <- (j-r):(j+r)

      seld <- d[
        ifelse(vert > 0 & vert <= nrow(d), vert, NA),
        ifelse(horiz > 0 & horiz <= ncol(d), horiz, NA)
      ]
      seld[1+r, 1+r] = NA
      
      if (min(seld, na.rm = T) > m) {
        minima <- rbind(minima, data.frame(i, j, value = m))
      }
    }
  }
  
  minima
}

#' Find cancidate solutions, i.e., the local minima
findCandidates <- function(distances) {
  rho <- sort(unique(distances$rho))
  psi <- sort(unique(distances$psi))
  
  distancesToMatrix(distances) %>%
    findValleys() %>%
    transmute(psi = psi[i], rho = rho[j], distance = value)
}

#' Annotate the candidates with goodnessOfFit, etc., that can be used for
#' ranking the cancidates.
annotateCandidates <- function(candidates, segments) {
  totalLen <- sum(segments$length)
  theoretMaxDist <- totalLen * 0.25
    
  candidates %>%
    rowwise() %>%
    mutate(
      goodnessOfFit = (1 - distance / theoretMaxDist) * 100,
      {
        nA <- calculateNA(segments$logR, segments$BAF, psi, rho)
        nB <- calculateNB(segments$logR, segments$BAF, psi, rho)
        
        tibble(
          # ploidy is recalculated based on results, to avoid bias (due to differences in normalization of LogR)
          ploidy = sum((nA + nB) * segments$length) / totalLen,
          percentZero = (sum((round(nA) == 0) * segments$length) +
                           sum((round(nB) == 0) * segments$length)) / totalLen,
          percentOddEven = sum((
            round(nA) %% 2 == 0 & round(nB) %% 2 == 1 |
              round(nA) %% 2 == 1 & round(nB) %% 2 == 0) * segments$length) / totalLen
        )
      }
    ) %>%
    ungroup()
}

#' Calculate penalties for candidates.
#' 
#' @param candidates A data frame of candidates
#' @param ploidyPenalty
#' @param vafPenaltyScale
#' @param percentZeroPenalty
#' 
penalizeCandidates <- function(candidates, ploidyPenalty = 0.01, vafPenaltyScale = 0.45, percentZeroPenalty = 0.05) {
  candidates %>%
    mutate(penalty.vafFactor = {
             # Prefer solutions that agree with TP53 VAF-based purity
             sd <- (TP53.purity.upper - TP53.purity.lower) / vafPenaltyScale
             vf <- dnorm(rho, mean = TP53.purity.mean, sd) / dnorm(0, 0, sd = sd)
             ifelse(is.na(vf), 1, vf)
             },
           penalty.ploidyFactor = pmax(0, 1 - abs(ploidy - 2) * ploidyPenalty),
           # The threshold 0.02 comes from ASCAT. No idea what it is based on.
           penalty.percentZeroFactor = ifelse(percentZero <= 0.02, 1 - percentZeroPenalty, 1)) %>%
    mutate(penalizedGoodnessOfFit = goodnessOfFit *
             penalty.vafFactor *
             penalty.ploidyFactor *
             penalty.percentZeroFactor)
}

addAlleleCountsToSegments <- function(segments, psi, rho, gamma = 1) {
  if (!is.na(rho) && rho == 0) {
    rho <- 1
  }
  
  # We are dealing only with XX. Thus, sex chromosomes need no special handling here.
  # If BAF is NA, the nAraw/nBraw etc will be NA too.
  segments <- segments %>%
    mutate(nAraw = calculateNA(logR, BAF, psi, rho, gamma),
           nBraw = calculateNB(logR, BAF, psi, rho, gamma)) %>%
    rowwise() %>%
    mutate({
      if (is.na(nAraw) | is.na(nBraw)) {
        tibble(nAraw = NA, nBraw = NA)
      } else if (nAraw + nBraw < 0) {
        # correct for negative values:
        tibble(nAraw = 0, nBraw = 0)
      } else if (nAraw < 0) {
        tibble(nAraw = 0, nBraw = nAraw + nBraw)
      } else if (nBraw < 0) {
        tibble(nAraw = nAraw + nBraw, nBraw = 0)
      } else {
        tibble(nAraw = nAraw, nBraw = nBraw)
      }
    }) %>%
    ungroup()

  # Handle samples that have no model. PurifiedLogR etc. equal to raw logR.
  if (is.na(rho)) {
    rho <- 1
    psi <- 2
  }

  segments %>%
    mutate(nMajor = round(nAraw),
           nMinor = round(nBraw),
           purifiedLogR = purifyLogR(rho, psi, logR),
           purifiedBaf = purifyBaf(rho, psi, BAF, 2^logR),
           purifiedLoh = ifelse(nProbesAf <= 0, NA, abs(purifiedBaf - 0.5) * 2))
}


#' Classify a sample as aberrant if it has enough loss of heterozygosity and large enough
#' fraction of the genome with R deviating at least 10%.
checkAberrancy <- function(segments) {
  segments %>%
    mutate(length = endpos - startpos) %>%
    summarise(aberrantLogR = sum(ifelse(abs(logR) > log(1.1), length, 0)) / sum(length),
              meanLoh = sum(abs(BAF - 0.5) * 2 * length, na.rm = TRUE) / sum(ifelse(is.na(BAF), 0, length))) %>%
    # I chose the thresholds by looking at the log-transformed histogram and choosing points
    # from the bottoms of the valleys.
    # A more scientific approach would be to use Otsu's method, for example.
    with(aberrantLogR > exp(-5) && meanLoh >= 0.025)
}

calculateSnpPurity <- function(copies, VAF) 2 / ((copies / VAF) - (copies - 2))

#' Fits a model, outputs the results as tables and a sunrise plot
#' 
#' @param segments A data frame with at least: chr, startpos, endpos, logR, BAF
#' @param TP53 A list: chrom, pos, altCount, refCount
#' @return A list with dataframes: candidates, segments, and a ggplot object: sunrise
#' 
ASCAT <- function(segments, TP53 = NULL) {
  MIN_GOODNESS_OF_FIT <- 70
  MIN_PENALIZED_GOODNESS_OF_FIT <- 55
  
  # Use only chr1-chrX
  chroms <- paste0("chr", c(1:22, "X"))
  
  psi <- seq(1, 8, 0.05)
  rho <- seq(0.05, 1, 0.005)
  
  # ##############
  
  psiExtremes <- c(min(psi), max(psi))
  
  aberrant <- checkAberrancy(segments)
  
  estimationSegments <- segments %>%
    filter(!is.na(BAF)) %>%
    # Only include autosomes in model fitting
    filter(chr %in% paste0("chr", 1:22)) %>%
    mutate(length = endpos - startpos)
  
  distances <- estimationSegments %>%
    computeDistances(psi, rho)
  
  candidates <- findCandidates(distances) %>%
    annotateCandidates(estimationSegments) %>%
    filter(!(psi %in% psiExtremes) & goodnessOfFit >= MIN_GOODNESS_OF_FIT)
  
  if (is.list(TP53)) {
    TP53 <- c(TP53,
              with(TP53, binom.confint(altCount, refCount + altCount, conf.level = 0.95, methods = "wilson")) %>%
                select(mean, lower, upper) %>%
                as.list())
    
    # Find the segment that contains the variant TP53 SNP
    TP53_segment <- segments %>%
      filter(chr == TP53$chrom & startpos < TP53$pos & endpos > TP53$pos) %>%
      as.list()
    
    candidates <- candidates %>%
      mutate(TP53.estimated.copies = calculateN(TP53_segment$logR, TP53_segment$BAF, ploidy, rho)) %>%
      mutate(TP53.purity.lower = calculateSnpPurity(TP53.estimated.copies, TP53$lower),
             TP53.purity.upper = calculateSnpPurity(TP53.estimated.copies, TP53$upper),
             TP53.purity.mean = calculateSnpPurity(TP53.estimated.copies, TP53$mean))
    
  } else {
    candidates <- candidates %>%
      mutate(TP53.estimated.copies = NA,
             TP53.purity.lower = NA,
             TP53.purity.upper = NA,
             TP53.purity.mean = NA)
  }
    
  candidates <- candidates %>%
    penalizeCandidates() %>%
#    filter(penalizedGoodnessOfFit >= MIN_PENALIZED_GOODNESS_OF_FIT) %>%
    arrange(-penalizedGoodnessOfFit) %>%
    mutate(rank = row_number())
  
  if (aberrant & nrow(candidates) > 0) {
    rho <- candidates$rho[1]
    psi <- candidates$psi[1]
    
    segments <- segments %>%
      addAlleleCountsToSegments(psi, rho) %>%
      mutate_if(is.numeric, round, digits = 4)
    
  } else {
    segments <- segments %>%
      mutate(nAraw = NA, nBraw = NA, nMajor = NA, nMinor = NA)
  }
  
  p <- ggplot((candidates %>% slice_head(n = 3)), aes(psi, rho)) +
    geom_tile(aes(fill = distance),
              data = distances,
              show.legend = FALSE) +
    scale_fill_gradientn(trans = "log",
                         colours = rev(RColorBrewer::brewer.pal(10, "RdBu"))) +
    scale_x_continuous(breaks = 1:8) +
    xlab("psi (Biased ploidy)") +
    ylab("rho (Aberrant cell fraction)")
  
  if (aberrant & nrow(candidates) > 0) {
    if (is.list(TP53)) {
      p <- p + with_shadow(
        geom_errorbar(aes(x = psi, ymin = TP53.purity.lower, ymax = TP53.purity.upper),
                      alpha = 0.4,
                      width = 0.2),
        x_offset = 0, y_offset = 0, colour = "white", sigma = 2) +
      geom_point(aes(x = psi, y = TP53.purity.mean),
                 alpha = 0.4) +
      geom_segment(aes(xend = psi, yend = TP53.purity.mean),
                   alpha = 0.2,
                   linetype = "dashed")
    }
      
    p <- p + with_shadow(
      geom_point(aes(psi, rho),
                 colour = "white",
                 shape = 4,
                 size = 3),
      x_offset = 0, y_offset = 0, sigma = 2.5) +
    with_shadow(
      geom_text(aes(psi,
                    rho,
                    label = paste0("(", rank, ".)\n",
                                   round(goodnessOfFit, 3), "%\n",
                                   round(penalizedGoodnessOfFit, 3), "%\n") ),
                colour = "white",
                size = 3,
                hjust = 0,
                nudge_x = 0.1),
      x_offset = 0, y_offset = 0, sigma = 2.5) +
    geom_point(data = slice_head(candidates),
               colour = "white",
               shape = 1, size = 5)
  }
    
  list(
    aberrant = aberrant,
    candidates = candidates %>% select(-distance),
    segments = segments,
    sunrise = p
  )
}

