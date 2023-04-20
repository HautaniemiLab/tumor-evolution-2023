#!/usr/bin/env Rscript

#R script for computing mutational spectra for sample tree branches from variant
#tables with SAC columns and the tree branches themselves. This is done for
#supersamples (constituting variants from primary branches i.e. 1+) and for
#single tree branches (variants specific to those branches).
#The outputs are two catalogue tables; supersamples in one and single branches
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

#Options include: input sample tree and input variant table directories,
#input sample tree and variant table suffices, mutation class and output file.
#Possible option for number of cores to use
    #Default options
    setOption("mc.cores", max(detectCores() - 1, 1))

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--input-sample-tree-directory", "-s")) {
            #input sample tree RData directory
            i = i + 1
            in_sample_trees_dir = args[i]
        } else if (args[i] %in% c("--input-sample-tree-suffix", "-in-st-suf")) {
            #input sample_tree RData suffix
            i = i + 1
            in_sample_trees_suf = args[i]
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
        !exists("in_sample_trees_dir") | !exists("in_sample_trees_suf") |
        !exists("in_vars_dir") | !exists("in_vars_suf") |
        !exists("mut_class") | !exists("out_dir") | !exists("out_suf")
    ) {
        stop(
            "Please supply input sample tree and input variant table directories, input sample tree and variant table suffices, mutation class (SBS/DBS/ID), output directory and file suffix",
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
# Join variants with tree branches
# --------------------------------------------------

#Prepare list of sample tree data
    in_sample_tree_files = grep(
        paste0(in_sample_trees_suf, "$"),
        paste0(in_sample_trees_dir, dir(in_sample_trees_dir)),
        value = T
    )
    names(in_sample_tree_files) = sub(".*\\/", "", sub(in_sample_trees_suf, "", in_sample_tree_files))

#Join variants with sample tree branches converted to SACs
    variants_joined = mclapply(
        names(variants_by_patient),
        function(patient) {
            #Load sample tree data
            load(in_sample_tree_files[patient])

            #Prepare branch variants for joining with SAC for presence absence
            vars_branches_join = cbind(
                joinKey(variants_branches),
                sapply(
                    sort(unique(variants_branches$branch)),
                    function(branch) {
                        c("0,0,0,0", "1,1,1,1")[(variants_branches$branch == branch) + 1]
                    }
                )
            )
            colnames(vars_branches_join)[-1] = paste0(patient, "_", sort(unique(variants_branches$branch)), ".SAC")

            #Prepare classified variants without SAC
            vars_mut_types_join = cbind(
                joinKey(variants_by_patient[[patient]]),
                variants_by_patient[[patient]][,!endsWith(colnames(variants_by_patient[[patient]]), ".SAC")]
            )

            #Join the tables and output without join key
            merge(vars_mut_types_join, vars_branches_join, by = "join_key")[-1]
        }
    )
    names(variants_joined) = names(variants_by_patient)

# --------------------------------------------------
# Compute catalogues
# --------------------------------------------------

#Choose correct reference signatures based on mutation class
    if (mut_class == "SBS") ref_sigs = ref_sigs_sbs
    if (mut_class == "DBS") ref_sigs = ref_sigs_dbs
    if (mut_class == "ID") ref_sigs = ref_sigs_id

#Compute supersample catalogues excluding branches -2, -1, 0, then combine into one table
    supersample_catalogues = buildCatalogue(
        variants_joined, ref_sigs, mut_class,
        "[A-Za-z]+[0-9]+_(-2|-1|0)", invert = T
    )

#Compute single sample tree branch catalogues, incl. branches -2, -1, 0
    branch_catalogues = buildCatalogueCohortSamples(variants_joined, ref_sigs, mut_class)

    #Combine single sample metadata and catalogues into a table
    branch_catalogues_out = data.frame(
        branch_catalogues$metadata,
        t(branch_catalogues$spectra),
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
    #Single sample tree branches
    write.table(
        branch_catalogues_out,
        paste0(out_dir, "sample_tree_branch_spectra", out_suf),
        quote = F, sep = "\t", col.names = NA
    )

# --------------------------------------------------
# End of file
# --------------------------------------------------


