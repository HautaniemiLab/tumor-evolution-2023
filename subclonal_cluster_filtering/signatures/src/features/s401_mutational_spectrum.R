#!/usr/bin/env Rscript

#R script for computing mutational spectra from variant tables with SAC columns.
#This is done for supersamples (variants added up from all samples of a patient)
#and for single samples (only variants with ALT reads in the sample are
#included).
#The outputs are two catalogue tables; supersamples in one and single samples in
#the other.

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

#Options include: input directory, input variant table suffix, mutation class
#and output file.
#Possible options for normal sample regex match-string and number of cores to
#use
    #Default options
    regex_normals = ""
    setOption("mc.cores", max(detectCores() - 1, 1))

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--input-directory", "-i")) {
            #input directory
            i = i + 1
            in_dir = args[i]
        } else if (args[i] %in% c("--input-suffix", "-in-suf")) {
            #input file suffix
            i = i + 1
            in_suf = args[i]
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
        } else if (args[i] %in% c("--regex-normals", "-r")) {
            #match-regex string for capturing normal samples
            i = i + 1
            regex_normals = args[i]
            if (is.na(regex_normals)) regex_normals = ""
        } else if (args[i] %in% c("--cores", "-c")) {
            #number of cores to use
            i = i + 1
            setOption("mc.cores", max(min(as.integer(args[i]), detectCores() - 1), 1))
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (!exists("in_dir") | !exists("in_suf") | !exists("mut_class") | !exists("out_dir") | !exists("out_suf")) {
        stop(
            "Please supply input directory, input variant table suffix, mutation class (SBS/DBS/ID), output directory and file suffix",
            call. = F
        )
    }

# --------------------------------------------------
# Read in variants (w/o chrM)
# --------------------------------------------------

#Read the mutation type-subtype tables
    in_files = grep(
        paste0(in_suf, "$"),
        paste0(in_dir, dir(in_dir)),
        value = T
    )
    variants = readVarFile(in_files)

# --------------------------------------------------
# Compute catalogues
# --------------------------------------------------

#Choose correct reference signatures based on mutation class
    if (mut_class == "SBS") ref_sigs = ref_sigs_sbs
    if (mut_class == "DBS") ref_sigs = ref_sigs_dbs
    if (mut_class == "ID") ref_sigs = ref_sigs_id

#Compute supersample catalogues excluding normals, then combine into one table
    supersample_catalogues = buildCatalogue(
        variants, ref_sigs, mut_class,
        regex_normals, nchar(regex_normals) > 0
    )

#Compute single sample catalogues, incl. normals (for troubleshooting)
    singlesample_catalogues = buildCatalogueCohortSamples(variants, ref_sigs, mut_class)

    #Combine single sample metadata and catalogues into a table
    singlesample_catalogues_out = data.frame(
        singlesample_catalogues$metadata,
        t(singlesample_catalogues$spectra),
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
    #Single samples
    write.table(
        singlesample_catalogues_out,
        paste0(out_dir, "single_sample_spectra", out_suf),
        quote = F, sep = "\t", col.names = NA
    )

# --------------------------------------------------
# End of file
# --------------------------------------------------


