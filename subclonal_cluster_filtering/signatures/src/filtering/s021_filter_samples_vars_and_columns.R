#!/usr/bin/env Rscript

#R script for filtering TSV variant tables' samples and variants.

#Samples are filtered using a meta data table and R/dplyr style
#filters given via arguments. Outputs a table with minimal
#columns for identifying variants (chrom, pos, ref, alt) and
#sample-specific strand allele count (SAC) columns.

#SAC can be computed from F2R1/F1R2 FORMAT columns. AF, if
#missing, can be computed from either column. If both SAC and
#F2R1/F1R2 columns are missing, filtering cannot proceed.

#Variants are filtered after sample filtering based on (ignoring
#regex matched normals):
# - minimum 1 ALT read in retained samples
# - exceed minimum VAF in at least one sample using AF columns
# - exceed minimum forward ALT read in any sample
# - exceed minimum reverse ALT read in any sample
# - R/dplyr style filters given via arguments

#Variants can additionally be filtered based on maximum VAF
#threshold in given samples (e.g. zero purity samples) via
#a maximum VAF filter table and patient ID. This uses AF columns.

# --------------------------------------------------
# Packages
# --------------------------------------------------

#dplyr
suppressMessages(require(dplyr))

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include: variants table with INFO and FORMAT annotations, sample (metadata) table and output file paths
#Possible options for
# - sample column name in metadata table
# - sample dplyr filters as a single string
# - variant dplyr filters as a single string
# - normal sample regex match-string
# - variant minimum VAF in any sample
# - variant minimum forward read count in any sample
# - variant minimum reverse read count in any sample
# - maximum VAF filter table
# - patient ID for VAF filter table
    #Default options
    sample_col_name = "sample"
    sample_filters_concatenated = ""
    variant_filters_concatenated = ""
    regex_normals = ""
    vaf_thres = 0.05
    fwd_thres = 1
    rev_thres = 1
    vaf_filter_table_file = ""
    patient = ""

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--input", "-i")) {
            #input variants table
            i = i + 1
            in_file = args[i]
        } else if (args[i] %in% c("--sample-list", "-s")) {
            #sample list
            i = i + 1
            sample_list_file = args[i]
        } else if (args[i] %in% c("--output", "-o")) {
            #output table path
            i = i + 1
            out_file = args[i]
        } else if (args[i] %in% c("--sample-column-name")) {
            #sample column name
            i = i + 1
            sample_col_name = args[i]
        } else if (args[i] %in% c("--sample-filters-concatenated", "-sf")) {
            #sample dplyr filters concatenated as a single string
            i = i + 1
            sample_filters_concatenated = args[i]
        } else if (args[i] %in% c("--variant-filters-concatenated", "-vf")) {
            #variant dplyr filters concatenated as a single string
            i = i + 1
            variant_filters_concatenated = args[i]
        } else if (args[i] %in% c("--regex-normals", "-r")) {
            #match-regex string for capturing normal samples
            i = i + 1
            regex_normals = args[i]
            if (is.na(regex_normals)) regex_normals = ""
        } else if (args[i] %in% c("--vaf")) {
            #vaf threshold
            i = i + 1
            vaf_thres = as.numeric(args[i])
        } else if (args[i] %in% c("--fwd")) {
            #forward threshold
            i = i + 1
            fwd_thres = as.numeric(args[i])
        } else if (args[i] %in% c("--rev")) {
            #reverse threshold
            i = i + 1
            rev_thres = as.numeric(args[i])
        } else if (args[i] %in% c("--vaf-filter-table")) {
            #maximum vaf filter table file
            i = i + 1
            vaf_filter_table_file = args[i]
        } else if (args[i] %in% c("--patient", "-p")) {
            #patient ID
            i = i + 1
            patient = args[i]
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (!exists("in_file") | !exists("sample_list_file") | !exists("out_file")) {
        stop(
            "Please supply input variants table, metadata table and output table paths",
            call. = F
        )
    }

# --------------------------------------------------
# Read in data
# --------------------------------------------------

#Read variants
    vars_tbl = read.table(in_file, header=T, sep="\t", quote="", fill=F, as.is=T)

#Read sample list
    #Try proper header metadata table, otherwise it is the old format
    sample_list = tryCatch({
        read.table(sample_list_file, header=T, sep="\t", as.is=T)
    }, error = function(e) {
        #Read and keep non-empty lines not starting with '#', '>' or '-'
        vec = readLines(sample_list_file, -1)
        out = data.frame(
            grep("^#|^>|^-|^$", test, perl=T, value=T, invert=T),
            stringsAsFactors = F
        )
        colnames(out) = NULL
        out
    })

#Read VAF filter table if given
    if (vaf_filter_table_file != "") {
        vaf_filter_table = read.table(vaf_filter_table_file, header=T, as.is=T)
    }

# --------------------------------------------------
# Remove NA lines
# --------------------------------------------------

#Remove NA lines from variants table
    vars_tbl = vars_tbl[!is.na(vars_tbl[,1]),,]

# --------------------------------------------------
# Compute strand read counts, and AFs if missing
# --------------------------------------------------

#Compute REF and ALT read counts for each strand from F2R1 and F1R2 or SAC (in order)
    #Sample names
    samples = unique(sub("\\..*$", "", grep("\\.(GT|DP|AD|AF|SAC|F1R2|F2R1)", colnames(vars_tbl), value=T)))

    #F2R1 and F1R2 and SAC columns
    f2r1_cols = grep("\\.F2R1$", colnames(vars_tbl), value=T)
    f1r2_cols = grep("\\.F1R2$", colnames(vars_tbl), value=T)
    sac_cols = grep("\\.SAC$", colnames(vars_tbl), value=T)

    #Attempt F2R1/F1R2, else SAC, else fail
    if (length(f2r1_cols) > 0 & length(f1r2_cols) > 0) {
        #Reference and ALT read counts
        ref_and_alt_counts = lapply(
            samples,
            function(sample) {
                f2r1 = vars_tbl[,grep(paste0("^", sample, "\\."), f2r1_cols, value=T)]
                f1r2 = vars_tbl[,grep(paste0("^", sample, "\\."), f1r2_cols, value=T)]

                cbind(
                    matrix(as.numeric(unlist(strsplit(f2r1, ","))), ncol=2, byrow=T),
                    matrix(as.numeric(unlist(strsplit(f1r2, ","))), ncol=2, byrow=T)
                )[,c(1,3,2,4)]
            }
        )
        names(ref_and_alt_counts) = samples
    } else if (length(sac_cols) > 0) {
        #Reference and ALT read counts
        ref_and_alt_counts = lapply(
            samples,
            function(sample) {
                sac = vars_tbl[,grep(paste0("^", sample, "\\."), sac_cols, value=T)]

                matrix(as.numeric(unlist(strsplit(sac, ","))), ncol=4, byrow=T)
            }
        )
        names(ref_and_alt_counts) = samples
    } else {
        stop(paste("Missing both F2R1/F1R2 and SAC FORMAT columns from, in_file, - either SAC or both F2R1/F1R2 is needed!"))
    }

#Compute AF from read counts if missing
    if (length(grep("\\.AF$", colnames(vars_tbl))) == 0) {
        #Compute AF for each sample
        afs_tbl = lapply(
            samples,
            function(sample) {
                apply(
                    ref_and_alt_counts[[sample]],
                    1,
                    function(x) sum(x[3:4]) / sum(x)
                )
            }
        )
        names(afs_tbl) = paste0(samples, ".AF")

        #Add AF columns to the variants table
        vars_tbl = cbind(vars_tbl, do.call(cbind, afs_tbl))
    }

# --------------------------------------------------
# Filter by max VAF with low purity sample
# --------------------------------------------------

#Filter with sample's VAFs against maximum VAF of samples present in VAF filter table
    if (vaf_filter_table_file != "") {
        #Choose only samples marked for the patient in VAF table
        if (patient != "") {
            vaf_filter_table = vaf_filter_table[vaf_filter_table$patient == patient,, drop=F]
        }
        #Only apply filter if VAF filter table contains a sample in variant table
        if (any(samples %in% vaf_filter_table[,sample_col_name])) {
            #Retrieve AF tables with appropriate sample name
            afs_tbl = vars_tbl[,grep("\\.AF$", colnames(vars_tbl), value=T)]
            colnames(afs_tbl) = sub("\\.AF$", "", colnames(afs_tbl))

            #Apply max VAF filter
            vaf_filter_idx = apply(
                apply(
                    vaf_filter_table[vaf_filter_table[,sample_col_name] %in% samples,,drop=F],
                    1,
                    function(x) {
                        afs_tbl[,x[sample_col_name]] < x["thres"] | x["thres"] >= 1
                    }
                ),
                1,
                all,
                na.rm = T
            )

            vars_tbl = vars_tbl[vaf_filter_idx,, drop=F]
            for (sample in samples) {
                ref_and_alt_counts[[sample]] = ref_and_alt_counts[[sample]][vaf_filter_idx,, drop=F]
            }
        }
    }

# --------------------------------------------------
# Filter samples
# --------------------------------------------------

#Compute the samples to retain using filters
    #Check if sample column name is in sample list to determine type of metadata table
    if (sample_col_name %in% colnames(sample_list)) {
        #Apply Dplyr filter on metadata table to choose samples
        for (filter_string in strsplit(sample_filters_concatenated, " -sf ")[[1]][-1]) {
            sample_list = eval(parse(
                text = paste("sample_list %>%", filter_string)
            ))
        }

        samples_retained = sample_list[,sample_col_name]
    } else {
        #Keep all samples
        samples_retained = sample_list[,1]
    }

#Filter sample FORMAT columns not among samples to retain
    #Get info and format column indices
    info_columns_idx = grep("\\.", colnames(vars_tbl), invert=T)
    format_columns_idx = grep("\\.", colnames(vars_tbl))

    #Filter FORMAT columns by sample name
    format_column_sample = sub("\\..*$", "", colnames(vars_tbl)[format_columns_idx])
    format_columns_idx = format_columns_idx[
        format_column_sample %in% samples_retained | ifelse(
            rep(nchar(regex_normals) > 0, length(format_column_sample)),
            grepl(regex_normals, format_column_sample),
            F
        )
    ]

    #End if no tumour samples are kept (no writing outputs)
    if (! any(grep(regex_normals, format_column_sample, value=T, invert=T) %in% samples_retained) | F) {
        quit("no", 0, F)
    }

    #Apply column filters
    vars_tbl = vars_tbl[,c(info_columns_idx, format_columns_idx)]

# --------------------------------------------------
# Compute missing SAC from read counts
# --------------------------------------------------

#Sample names retained
    samples = unique(sub("\\..*$", "", grep("\\.", colnames(vars_tbl), value=T)))

#Compute SAC from read counts if SAC is missing
    if (length(grep("\\.SAC$", colnames(vars_tbl))) == 0) {
        #Compute SAC for each patient
        sacs_tbl = lapply(
            samples,
            function(sample) {
                apply(
                    ref_and_alt_counts[[sample]],
                    1,
                    paste,
                    collapse = ","
                )
            }
        )
        names(sacs_tbl) = paste0(samples, ".SAC")

        #Add SAC columns to the variants table
        vars_tbl = cbind(vars_tbl, do.call(cbind, sacs_tbl))
    }

# --------------------------------------------------
# Filter variants by VAF and ALT read counts
# --------------------------------------------------

#Apply filters
    vars_tbl = vars_tbl[
        #Any ALT read filter (at least one sample must have ALT read(s))
        apply(
            sapply(
                grep(regex_normals, samples, value=T, invert=T),
                function(sample) {
                    rowSums(ref_and_alt_counts[[sample]][,3:4]) > 0
                }
            ),
            1,
            any
        ) &
        #VAF filter
        apply(
            vars_tbl[,
                grep(regex_normals,
                    grep("\\.AF$", colnames(vars_tbl), value=T),
                    invert=T, value=T
                ),
                drop = F
            ] > vaf_thres,
            1,
            any,
            na.rm=T
        ) &
        #Forward read
        apply(
            sapply(
                grep(regex_normals, samples, value=T, invert=T),
                function(sample) {
                    ref_and_alt_counts[[sample]][,3] >= fwd_thres
                }
            ),
            1,
            any
        ) &
        #Reverse read
        apply(
            sapply(
                grep(regex_normals, samples, value=T, invert=T),
                function(sample) {
                    ref_and_alt_counts[[sample]][,4] >= rev_thres
                }
            ),
            1,
            any
        ),
    ]

# --------------------------------------------------
# Custom variant filters
# --------------------------------------------------

#Apply variant Dplyr filters
    for (filter_string in strsplit(variant_filters_concatenated, " -vf ")[[1]][-1]) {
        vars_tbl = eval(parse(
            text = paste("vars_tbl %>%", filter_string)
        ))
    }

# --------------------------------------------------
# Filter columns and write table
# --------------------------------------------------

#Keep only "CHROM", "POS", "REF", "ALT" and SAC columns
    vars_tbl = vars_tbl[,c("CHROM", "POS", "REF", "ALT", grep("\\.SAC$", colnames(vars_tbl), value=T))]

#Write table
    write.table(vars_tbl, out_file, quote = F, sep = "\t", row.names=F)


