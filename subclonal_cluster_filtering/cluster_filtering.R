#Script for filtering certain subclonal clusters in evolution project.
#Generates a table of retained clusters for signature analysis.
#Also generates a table with decisions and metrics/scores used in it.

# --------------------------------------------------
# Packages
# --------------------------------------------------

library(dplyr)

# --------------------------------------------------
# Preprocess data
# --------------------------------------------------

#Read data
dir_prefix = "signatures/results/attributions/filt_set_validation"
dir_suffix = "/subclones/reference/contributions_sbs_subclone_proportions.tsv"
set_ids = c(".no_filt", ".no_filt.arts", "", ".arts")

exps_raw = lapply(
    set_ids,
    function(set_id) {
        read.table(
            paste0(dir_prefix, set_id, dir_suffix),
            header=T, sep="\t", as.is=T
        )
    }
)

#Summarise the exposure tables (sum of artefact proportions)
art_sigs = paste0("SBS", c("7c", 34, 41, 47, 57, 58, 90))
cols_old = c("sample_id", "mut_count", "propArt", "cos_sim")
cols_new = c("cluster", "size", "propArt", "cosSim")

exps_processed = lapply(
    exps_raw,
    function(tbl) {
        tbl %>%
            mutate( propArt = rowSums(select(., colnames(.)[colnames(.) %in% art_sigs])) ) %>%
            select(all_of(cols_old)) %>%
            rename(setNames(cols_old, cols_new))
    }
)

#Merge the tables
exps_processed_merged = exps_processed[[1]] %>%
    merge(exps_processed[[2]] %>% select(cluster, artFit.propArt=propArt, artFit.cosSim=cosSim), sort=F) %>%
    merge(exps_processed[[3]] %>% select(cluster, filt.propArt=propArt, filt.cosSim=cosSim), sort=F) %>%
    merge(exps_processed[[4]] %>% select(cluster, filt.size=size, filt.artFit.propArt=propArt, filt.artFit.cosSim=cosSim), sort=F) %>%
    select(
        cluster, size, filt.size,
        propArt, artFit.propArt, filt.propArt, filt.artFit.propArt,
        cosSim, artFit.cosSim, filt.cosSim, filt.artFit.cosSim,
    ) %>%
    mutate( propFilt = 1 - filt.size / size )

# --------------------------------------------------
# Quick visualisation
# --------------------------------------------------

titles = c("", "artFit", "filt", "filt.artFit")

#Cosine similarity vs size
par(mfrow=c(2,2))
    for (i in 1:4) {
        plot(
            exps_processed[[i]]$size, exps_processed[[i]]$cosSim,
            xlim = c(0, 500),
            main = titles[i], xlab = "cosine similarity", ylab = "size"
        )
        lines(smooth.spline(exps_processed[[i]]$size, exps_processed[[i]]$cosSim), col="red")
        abline(v=96, lty=2)
        abline(h=.7, lty=2)
    }
par(mfrow=c(1,1))

#Forced artefacts vs none
par(mfrow=c(2,1))
    for (i in c(1,3)) {
        plot(
            exps_processed[[i]]$cosSim, exps_processed[[i+1]]$cosSim,
            main = paste0("cosine similarities - ", titles[i+1], " vs ", titles[i]),
            xlab = "cosine similarity", ylab = "cosine similarity",
            col = c("black", "blue")[(exps_processed[[i+1]]$propArt > .2 ) + 1]
        )
        abline(0, 1, col="red")
    }
par(mfrow=c(1,1))

#Forced artefacts vs none
par(mfrow=c(2,2))
    #i=1 (no forced), i=2 (forced)
    for (i in c(1,2)) {
        #1. neither (grey), 2. size (red), 3. art (blue), 4. both (magenta)
        point_cols = c("#00000066", "#FF000066", "#0000FF66", "#FF00FF66")[
            matrix(
                c(
                    exps_processed[[i+2]]$size / exps_processed[[i]]$size < .6, #Heavily filtered (decrease by filtering)
                    exps_processed[[2]]$propArt > .2                            #high artefact to begin with (artefact fit, no filt)
                ),
                ncol=2
            ) %*% 1:2 + 1
        ]
    
        plot(
            exps_processed[[i]]$cosSim, exps_processed[[i+2]]$cosSim,
            main = paste0("cosine similarities - ", titles[i+2], " vs ", titles[i]),
            xlab = "cosine similarity", ylab = "cosine similarity",
            col = point_cols, pch=16
        )
        abline(0, 1, col="orange")
    
        plot(
            exps_processed[[i]]$propArt, exps_processed[[i+2]]$propArt,
            main = paste0("%artefacts - ", titles[i+2], " vs ", titles[i]),
            xlab = "%artefacts", ylab = "%artefacts",
            col = c("#00000066", "#0000FF66")[(exps_processed[[i+2]]$size / exps_processed[[i]]$size < .6 ) + 1], pch=16
        )
        abline(0, 1, col="orange")
    }
par(mfrow=c(1,1))

#cosine similarity change vs %artefact change
par(mfrow=c(2,1))
    for (i in c(1,3)) {
        y = exps_processed[[i+1]]$cosSim - exps_processed[[i]]$cosSim
        x = exps_processed[[i+1]]$propArt - exps_processed[[i]]$propArt
        plot(
            x, y,
            main = paste0("cos.sim diff vs %artefact diff - ", titles[i+1], " vs ", titles[i]),
            xlab = "delta %artefact", ylab = "delta cosine similarity"
        )
        lines(smooth.spline(x, y, spar=1.1), col="red")
    }
par(mfrow=c(1,1))

# --------------------------------------------------
# Test filtering
# --------------------------------------------------

#Target clusters to filter

#View data
exps_processed_merged %>%
    mutate(score = propArt + artFit.propArt * pmax(artFit.cosSim / cosSim - 1, 0)) %>%
    mutate(filt.score = filt.propArt + filt.artFit.propArt * pmax(filt.artFit.cosSim / filt.cosSim - 1, 0)) %>%
    mutate(propBad = propFilt + (1 - propFilt) * filt.artFit.propArt) %>%
    select(cluster, size, filt.size, propFilt, propBad, score, filt.score, everything()) %>%
    View()

#Test filtering with some scores computed
test_filtered = exps_processed_merged %>%
    mutate(score = (propArt * cosSim + artFit.propArt * artFit.cosSim) / (cosSim + artFit.cosSim)) %>%
    mutate(
        filt.score = (filt.propArt * filt.cosSim + filt.artFit.propArt * filt.artFit.cosSim) /
            (filt.cosSim + filt.artFit.cosSim)
    ) %>%
    mutate(propBad = propFilt + (1 - propFilt) * filt.artFit.propArt) %>%
    select(cluster, size, filt.size, propFilt, propBad, score, filt.score, everything()) %>%
    filter(propBad >= .50 & score >= .25 | filt.size < 30) %>%
    pull(cluster)

#Extra
sum(test_filtered %in% target_clusters) / length(test_filtered)
test_filtered[! test_filtered %in% target_clusters]

#Missing
sum(target_clusters %in% test_filtered) / length(target_clusters)
target_clusters[! target_clusters %in% test_filtered]

# --------------------------------------------------
# Filter and produce output tables
# --------------------------------------------------

#Perform actual filtering
#Scores are weighted (cosine similarity) sums of artefact proportions
#propBad adds proportion of original mutations filtered with remainder's artefact proportion
filter_tbl = exps_processed_merged %>%
    mutate(score = (propArt * cosSim + artFit.propArt * artFit.cosSim) / (cosSim + artFit.cosSim)) %>%
    mutate(propBad = propFilt + (1 - propFilt) * filt.artFit.propArt) %>%
    select(cluster, size, filt.size, propFilt, propBad, score, everything()) %>%
    mutate( reason = ifelse(filt.size < 30, "small_30", "") ) %>%
    mutate(
        reason = ifelse(
            propBad >= .50 & score >= .25,
            sub(",$", "", paste0("propBad_50_score_25,", reason)),
            reason
        )
    ) %>%
    mutate( reason = ifelse(reason == "", NA, reason) ) %>%
    mutate(filtered = ifelse(is.na(reason), F, T)) %>%
    select(cluster, filtered, reason, everything())

#Table of retained clusters to use in signature analysis
subclones_to_use_tbl = filter_tbl %>%
    filter(!filtered) %>%
    mutate(patient = sub("_\\d+$", "", cluster)) %>%
    mutate(cluster_id = sub("^\\w+_", "", cluster)) %>%
    mutate(clusters = ".") %>%
    mutate(time_points = NA) %>%
    select(patient, cluster_id, clusters, time_points)

# --------------------------------------------------
# Write outputs
# --------------------------------------------------

write.table(filter_tbl, "cluster_filter_status.validation.csv", row.names=F, quote=F, sep="\t")
write.table(subclones_to_use_tbl, "signature.subclones.validation.csv", row.names=F, quote=F, sep="\t")


