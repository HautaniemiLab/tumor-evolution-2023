#Script for computing frequenly mutated genes for the manuscript

#In summary, dNdScv is run using discovery and validation sets
#separately.

#Variants are filtered to have AF > 0 in selected sample for each
#patient.

#For dNdScv, phased mutations (gap <= 10bp) are combined according to
#the input tables to avoid the same mutation being repeated as
#multiple separate ones.

library(tidyverse)
library(parallel)

options(readr.show_col_types = FALSE)
options(mc.cores = max(detectCores(), 1))

# --------------------------------------------------
# Inputs
# --------------------------------------------------

setwd("/path/to/frequently_mutated_genes")

#Define inputs
input_sets = list(
    "discovery" = list(
        "samples" = "phased_mutations/vcf_files.tsv",
        "vars_dir" = "data/variants_csv/",
        "phased_vars_split" = "phased_mutations/phased_mutations.split.tsv",
        "phased_vars_merged" = "phased_mutations/phased_mutations.merged.tsv"
    ),
    "validation" = list(
        "samples" = "phased_mutations/validation.vcf_files.tsv",
        "vars_dir" = "data/variants_csv.validation/",
        "phased_vars_split" = "phased_mutations/validation.phased_mutations.split.tsv",
        "phased_vars_merged" = "phased_mutations/validation.phased_mutations.merged.tsv"
    )
)

# --------------------------------------------------
# Read and pre-process data
# --------------------------------------------------

#Functions
#---------

#Function for reading in phased and split mutation tables with removal of chr prefix from contig names
readPhasedVars = function(file_path) {
    read_tsv(file_path) %>% mutate(chr = sub("^chr", "", chr))
}

#Function to read and process a single variants table
processVarsInput = function(tbl_path, patient, sample) {
    #Keep only specific annotations and AF column for selected tumour sample. Filter AF > 0
    read_tsv(tbl_path) %>%
        mutate(FILTER = {if ("FILTER" %in% colnames(.)) FILTER else "PASS"}) %>% #Add FILTER column as PASS only if missing
        select(
            CHROM, POS, REF, ALT, FILTER,
            ends_with(".refGene"), starts_with("CLN"), starts_with("CADD_"),
            all_of(paste0(sample, ".AF"))
        ) %>%
        rename_with(
            ~sub(paste0(sample, "\\."), "", .),
            starts_with(sample)
        ) %>%
        filter(AF > 0) %>%
        mutate(
            patient = patient,
            sample = sample
        ) %>%
        select(patient, sample, everything())
}

#Process variants
#----------------

#Read and process mutation tables, extracting chosen sample per each patient
#Also prepares the tables for dNdScv with phased variants merged
mutation_sets = mclapply(
    input_sets,
    function(inputs) {
        #Read sample table
        tbl_samples = read_tsv(inputs$samples)

        #Read and process mutations, including filtering
        tbl_mutations = mclapply(
            tbl_samples %>% group_split(patient),
            function (x) {
                processVarsInput(paste0(inputs$vars_dir, x$patient, ".csv"), x$patient, x$sample)
            },
            mc.cores = round(getOption("mc.cores") / min(getOption("mc.cores"), length(input_sets))) %>%
                max(1)
        ) %>%
            do.call(rbind, .)

        #Read phased variant tables
        tbl_phased_vars_split  = readPhasedVars(inputs$phased_vars_split)
        tbl_phased_vars_merged = readPhasedVars(inputs$phased_vars_merged)

        #Output
        list(
            "mutations" = tbl_mutations,
            "phased_vars_split"  = tbl_phased_vars_split,
            "phased_vars_merged" = tbl_phased_vars_merged
        )
    },
    mc.cores = min(getOption("mc.cores"), length(input_sets))
)

# --------------------------------------------------
# Run dNdScv
# --------------------------------------------------

#Function for converting the mutation table to format desired by dndscv
#Merges nearby phased mutations according to table input
#Adds distance to nearest variant for investigating complex events
mutsToDndscv = function(muts_tbl, phased_vars_split=NULL, phased_vars_merged=NULL, merge_adjacent=T, compute_dist=T) {
    out = muts_tbl %>%
        select(sampleID=sample, chr=CHROM, pos=POS, ref=REF, mut=ALT) %>%
        mutate(chr = sub("^chr", "", chr)) %>%
        filter(chr != "M")

    if (merge_adjacent) {
        phased_muts = merge(
            out %>% mutate(patient = sub("_.*$", "", sampleID)),
            phased_vars_split
        ) %>%
            group_by(patient, index) %>%
            filter(n() > 1) %>%
            ungroup() %>%
            select(sampleID, chr, pos, ref, mut, patient, index)

        out = out %>%
            setdiff(phased_muts %>% select(-patient, -index)) %>%
            bind_rows(
                phased_vars_merged %>%
                    merge(
                        phased_muts %>% select(sampleID, patient, index) %>% distinct(),
                        by = c("patient", "index")
                    ) %>%
                    select(-patient, -index)
            ) %>%
            arrange(sampleID, chr, pos)
    }

    if (compute_dist) {
        out = out %>%
            group_by(sampleID, chr) %>%
            mutate(dist = pmin(c(diff(pos), NA), c(NA, diff(pos)), na.rm=T)) %>%
            ungroup()
    }

    out
}

#Prepare dNdScv
library(dndscv)

#Load hg38 covariates object
load("data/covariates_hg19_hg38_epigenome_pcawg.rda")

#Prepare gene list without chrY genes
load("data/RefCDS_human_GRCh38_GencodeV18_recommended.rda")

chr_genes = sapply(
    RefCDS,
    function(x) c(x$chr, x$gene_name)
) %>%
    t() %>%
    as.data.frame() %>%
    select(chr=V1, gene_name=V2) %>%
    filter(chr != "Y") %>%
    pull(gene_name)

#Run dNdScv for each set
dndscv_sets = mclapply(
    mutation_sets,
    function(mutation_set) {
        #Convert mutation table to desired dNdScv format with phased variants merged
        tbl_dndscv_mutations = mutsToDndscv(
            mutation_set$mutations,
            mutation_set$phased_vars_split,
            mutation_set$phased_vars_merged
        )

        #Run dNdScv
        dndscv(
            tbl_dndscv_mutations, gene_list=chr_genes,
            refdb="data/RefCDS_human_GRCh38_GencodeV18_recommended.rda", cv=covs
        )
    }
)

# --------------------------------------------------
# Outputting
# --------------------------------------------------

#Functions
#---------

#Compute frequency of mutated patients per gene (stop loss is not considered by dndscv)
dndscvMutFreq = function(dndscv_result, sort=F) {
    out = dndscv_result$annotmuts %>%
        filter(! impact %in% c("Stop_loss", "Synonymous")) %>%
        select(sampleID, gene) %>%
        distinct() %>%
        group_by(gene) %>%
        summarise(
            nPatients = n(),
            freqPatients = nPatients / n_distinct(dndscv_result$annotmuts$sampleID)
        ) %>%
        merge(dndscv_result$sel_cv %>% select(gene=gene_name), by="gene", all.y=T) %>%
        mutate_if(is.numeric, ~ifelse(is.na(.), 0, .))

    if (sort) {
        out %>% arrange(match(gene, dndscv_result$sel_cv$gene_name))
    } else {
        out
    }
}

#Write results
#-------------

#Make sure output directory exists
dir.create("output", showWarnings=F)

#Write all sets, including combined and filtered list
tbls_sel_cv = lapply(
    names(mutation_sets),
    function(set_name) {
        #Retrieve repeated variables
        dndscv_set = dndscv_sets[[set_name]]

        #Merge the sel_cv table with mutation frequencies
        tbl_sel_cv = dndscv_set$sel_cv %>%
            merge(dndscvMutFreq(dndscv_set), by.x="gene_name", by.y="gene", sort=F) %>%
            arrange(match(gene_name, dndscv_set$sel_cv$gene_name))

        #Write the sel_cv table as well as the annotated mutations of dndscv as is
        write_tsv(tbl_sel_cv, paste0("output/", set_name, ".dndscv.sel_cv.tsv"))
        write_tsv(dndscv_set$annotmuts, paste0("output/", set_name, ".dndscv.annotmuts.tsv"))

        #Output for merging
        tbl_sel_cv %>%
            select(
                Gene=gene_name, Syn=n_syn, Mis=n_mis, Non=n_non, Spl=n_spl, Ind=n_ind,
                `P-value`=pglobal_cv, `Q-value`=qglobal_cv,
                Patients=nPatients, `Fraction of Patients`=freqPatients
            )
    }
) %>%
    setNames(names(mutation_sets))

#Side-by-side table of discovery and validation runs with NAs as empty fields
lapply(
    c("discovery", "validation"),
    function(set_name) {
        rbind(
            colnames(tbls_sel_cv[[set_name]]),
            tbls_sel_cv[[set_name]] %>%
                filter(`P-value` < 0.05) %>%
                setNames(paste0(set_name, colnames(tbls_sel_cv[[set_name]])))
        ) %>%
        mutate(i = row_number())
    }
) %>%
    reduce(full_join, by="i") %>%
    select(-i) %>%
    setNames(NULL) %>%
    write_tsv("output/S2.merged.dndscv.sel_cv.tsv", na="", col_names=F)

