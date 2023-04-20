#!/usr/bin/env Rscript

#R script for computing mutational spectra for subclonal clusters from variant
#tables with SAC columns and variant cluster tracks. This is done for super-
#samples (constituting select clusters) and for individual subclonal clusters
#(variants specific to those clusters).
#The outputs are two catalogue tables; supersamples in one and single clusters
#in the other.

# --------------------------------------------------
# Sourcing
# --------------------------------------------------

#Load reference signatures
source("src/features/f0_load_reference_signatures.R")

#Load spectrum building functions
source("src/features/f0_spectrum_construction_functions.R")

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include: input subclone and input variant table directories, subclone
#and variant table suffices, mutation class and output file.
#Possible option for
# - number of cores to use
# - subclone info table to specify retained clusters and combination clusters
    #Default options
    setOption("mc.cores", max(detectCores() - 1, 1))
    in_subclone_info = ""

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--input-subclone-directory", "-s")) {
            #input subclones directory
            i = i + 1
            in_subclones_dir = args[i]
        } else if (args[i] %in% c("--input-subclone-suffix", "-in-s-suf")) {
            #input subclones table suffix
            i = i + 1
            in_subclones_suf = args[i]
        } else if (args[i] %in% c("--subclone-info", "-s-info")) {
            #input supersample subclones info table
            i = i + 1
            in_subclone_info = args[i]
        } else if (args[i] %in% c("--input-variants-directory", "-v")) {
            #input variants directory
            i = i + 1
            in_vars_dir = args[i]
        } else if (args[i] %in% c("--input-variants-suffix", "-in-v-suf")) {
            #input variant file suffix
            i = i + 1
            in_vars_suf = args[i]
        } else if (args[i] %in% c("--mutation-class", "-m")) {
            #type (SBS, DBS or ID)
            i = i + 1
            mut_class = args[i]
        } else if (args[i] %in% c("--output-directory", "-o")) {
            #output directory
            i = i + 1
            out_dir = args[i]
        } else if (args[i] %in% c("--output-suffix", "-out-suf")) {
            #output suffix
            i = i + 1
            out_suf = args[i]
        } else if (args[i] %in% c("--cores", "-c")) {
            #number of cores to use
            i = i + 1
            setOption("mc.cores", max(min(as.integer(args[i]), detectCores() - 1), 1))
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (
        !exists("in_subclones_dir") | !exists("in_subclones_suf") |
        !exists("in_vars_dir") | !exists("in_vars_suf") |
        !exists("mut_class") | !exists("out_dir") | !exists("out_suf")
    ) {
        stop(
            "Please supply input subclone and input variant table directories, input subclone and variant table suffices, mutation class (SBS/DBS/ID), output directory and file suffix",
            call. = F
        )
    }

# --------------------------------------------------
# Helper functions
# --------------------------------------------------

#Function for creating variant join key (default CHROM-POS) for each variant line
    joinKey = function(tbl, columns=c("CHROM", "POS")) {
        data.frame(
            "join_key" = apply(
                tbl,
                1,
                function(x) paste(
                    gsub(" ", "", x[columns]),
                    collapse="-"
                )
            ),
            stringsAsFactors = F
        )
    }

# --------------------------------------------------
# Read in variants (w/o chrM)
# --------------------------------------------------

#Read the mutation type-subtype tables
    in_vars_files = grep(
        paste0(in_vars_suf, "$"),
        paste0(in_vars_dir, dir(in_vars_dir)),
        value = T
    )
    variants_by_patient = readVarFile(in_vars_files)

# --------------------------------------------------
# Prepare variables for capturing clusters of supersample and combinations
# --------------------------------------------------

#Default with no cluster filtering or combination clusters if such table is not given
if (in_subclone_info == "") {
    #All clusters are retained
    supersample_regex = ""

    #No combination clusters
    combination_clusters = list()
    combination_clusters_patient = c()
} else {
    #Read table
    subclone_info = read.table(in_subclone_info, header=T, as.is=T)

    #Create regex to capture relevant clusters for supersample
    supersample_regex = unique(unlist(
        lapply(
            unique(subclone_info$patient),
            function(patient) {
                paste(
                    patient,
                    unlist(
                        strsplit(as.character(subclone_info[subclone_info$patient == patient, "cluster_id"]), ",")
                    ),
                    sep = "_"
                )
            }
        )
    ))
    supersample_regex = paste0("^(", paste(supersample_regex, collapse="|"), ").SAC$")

    #Create list for creating combination clusters
    combination_clusters = lapply(
        grep(",", subclone_info$clusters),
        function(i) {
            paste(
                subclone_info[i, "patient"],
                strsplit(as.character(subclone_info[i, "clusters"]), ",")[[1]],
                sep = "_"
            )
        }
    )
    names(combination_clusters) = apply(
        subclone_info[grep(",", subclone_info$clusters),],
        1,
        function(x) {
            paste(x["patient"], x["cluster_id"], sep="_")
        }
    )

    combination_clusters_patient = data.frame(
        "patient" = subclone_info[grep(",", subclone_info$clusters), "patient"],
        row.names = names(combination_clusters)
    )
}

# --------------------------------------------------
# Join variants with subclonal clusters
# --------------------------------------------------

#Prepare list of subclone data
    in_subclone_files = grep(
        paste0(in_subclones_suf, "$"),
        paste0(in_subclones_dir, dir(in_subclones_dir)),
        value = T
    )
    names(in_subclone_files) = sub(".*\\/", "", sub(in_subclones_suf, "", in_subclone_files))

#Analyse only those with variant data
    in_subclone_files = in_subclone_files[names(in_subclone_files) %in% names(variants_by_patient)]

#Join variants with subclonal clusters converted to SACs
    variants_joined = mclapply(
        names(in_subclone_files),
        function(patient) {
            #Load subclone data
            variants_clusters = read.table(in_subclone_files[patient], header=T, as.is=T)

            #Prepare cluster variants for joining with SAC for presence absence
            vars_clusters_join = cbind(
                joinKey(variants_clusters),
                sapply(
                    sort(unique(variants_clusters$cluster)),
                    function(cluster) {
                        c("0,0,0,0", "1,1,1,1")[(variants_clusters$cluster == cluster) + 1]
                    }
                )
            )
            colnames(vars_clusters_join)[-1] = paste0(patient, "_", sort(unique(variants_clusters$cluster)), ".SAC")

            #Prepare classified variants without SAC
            vars_mut_types_join = cbind(
                joinKey(variants_by_patient[[patient]]),
                variants_by_patient[[patient]][,!endsWith(colnames(variants_by_patient[[patient]]), ".SAC")]
            )

            #Join the tables and output without join key
            merge(vars_mut_types_join, vars_clusters_join, by = "join_key")[-1]
        }
    )
    names(variants_joined) = names(in_subclone_files)

# --------------------------------------------------
# Compute catalogues
# --------------------------------------------------

#Choose correct reference signatures based on mutation class
    if (mut_class == "SBS") ref_sigs = ref_sigs_sbs
    if (mut_class == "DBS") ref_sigs = ref_sigs_dbs
    if (mut_class == "ID") ref_sigs = ref_sigs_id

#Compute supersample catalogues and combine into one table
    supersample_catalogues = buildCatalogue(
        variants_joined, ref_sigs, mut_class,
        supersample_regex
    )

#Compute single subclone catalogues
    cluster_catalogues = buildCatalogueCohortSamples(variants_joined, ref_sigs, mut_class)

    #Add combination clusters
    if (length(combination_clusters) > 0) {
        combination_cluster_catalogues = sapply(
            combination_clusters,
            function(x) {
                rowSums(cluster_catalogues$spectra[rownames(ref_sigs), x, drop=F])
            }
        )

        cluster_catalogues$metadata = rbind(cluster_catalogues$metadata, combination_clusters_patient)
        cluster_catalogues$spectra = cbind(cluster_catalogues$spectra, combination_cluster_catalogues)
    }

    #Add filtering and combination status metadata
    filtered_cluster = (
        !grepl(supersample_regex, paste0(colnames(cluster_catalogues$spectra), ".SAC")) &
        !colnames(cluster_catalogues$spectra) %in% names(combination_clusters)
    )

    cluster_catalogues$metadata = data.frame(
        cluster_catalogues$metadata,
        filtered = filtered_cluster,
        combination = c(
            rep(F, ncol(cluster_catalogues$spectra) - length(combination_clusters)),
            rep(T, length(combination_clusters))
        )
    )

    #Sort cluster catalogues by patient with combination clusters at the back
    cluster_catalogues$spectra = cluster_catalogues$spectra[, order(cluster_catalogues$metadata$patient), drop=F]
    cluster_catalogues$metadata = cluster_catalogues$metadata[order(cluster_catalogues$metadata$patient),,drop=F]

    #Combine cluster metadata and catalogues into a table
    cluster_catalogues_out = data.frame(
        cluster_catalogues$metadata,
        t(cluster_catalogues$spectra),
        check.names = F
    )

# --------------------------------------------------
# Write the catalogues
# --------------------------------------------------

#Write the catalogues
    #Supersamples
    write.table(
        t(supersample_catalogues),
        paste0(out_dir, "supersample_spectra", out_suf),
        quote = F, sep = "\t", col.names = NA
    )
    #Single subclones
    write.table(
        cluster_catalogues_out,
        paste0(out_dir, "subclone_spectra", out_suf),
        quote = F, sep = "\t", col.names = NA
    )

# --------------------------------------------------
# End of file
# --------------------------------------------------


