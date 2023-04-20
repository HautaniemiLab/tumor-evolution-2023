#!/usr/bin/env Rscript

#R script for attributing signatures on multisample spectra, given single sample
#as well as supersample tables. Outputs contributions in both counts as well as
#proportions in a table for supersamples and single samples each.

#Performs signature attribution by the following steps:
# 1. data set wide supersample forward selection to choose common signatures
# 2. Supersample signature attribution with the following steps:
#   a. supersample-wise backward-forward selection from data set common
#      signatures
#   b. single sample level backward-forward selection to get sample signatures
#      to improve sensitivity from combined common and supersample signatures
#   c. supersample-wise backward selection from combined single sample
#      signatures (with rules enforced)
# 3. single sample backward selection (with rules enforced) from supersample
#    signatures to get final signature contributions for individual samples

#Data set and supersample attributions only consider supersamples or samples
#respectively with a minimum number of mutations.

#Supersample-only mode skips all steps involving single samples. This means that
#supersample attribution involves only a single backward-forward model selection
#starting from common signatures.

#Independent samples mode analyses samples as if they were independent samples
#from different patients. In practice runs like supersample-only mode where each
#individual sample becomes a supersample not affecting others outside of finding
#common signatures. The single sample tables are still written for plotting.
#Overrides supersample-only mode.

#Sensitive mode relaxes the final supersample backward selection step as well as
#replaces them with backward-forward selection. Only affects multisample
#patients.

#Poisson model uses Poisson loglikelihood base BIC in stepwise model selection,
#with attributions fitted optimising generalised Kullback-Leibler. It overfits
#with WGS data (thousands of mutations per sample) so should only be considered
#in targeted sequencing, e.g. WES data.

#Data set wide forward selection can be disabled and the set of signatures found
#there or otherwise can be supplemented via options.

#Rules are always enforced in backward-forward selection; it implies forcing
#SBS1 and SBS5 by default as well as connected signatures (mainly APOBEC) in SBS
#data. Similar rules are present in COSMIC signatures. Forced signatures can be
#specified via options.

# --------------------------------------------------
# Sourcing
# --------------------------------------------------

#Load reference signatures
source("src/features/f0_load_reference_signatures.R")

#Load signature attribution rules
source("src/attribution/f0_load_signature_attribution_rules.R")

#Load attribution functions
source("src/attribution/f1_signature_attribution_functions.R")

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include: supersample and single sample input spectra table files,
#mutation class (SBS/DBS/ID) and output directory.
#Possible options for
# - normal sample regex match-string
# - patients to ignore in dataset selection step; regex match-string
# - custom reference signatures, e.g. exome signatures
# - whether data set forward selection to obtain common signature should be performed
# - any known common signatures to add in addition to data set selection
# - custom forced signatures (also to disable default option)
# - whether to analyse supersamples only
# - whether all samples are analysed as if they were different patients
# - whether to use Poisson model for attribution
# - whether to use more lenient selection thresholds
# - sensitive mode's single-to-supersample backward threshold
# - number of cores to use
    #Default options
    regex_normals = ""
    regex_dataset_sel = ""
    dataset_selection = T
    common_sigs = c()
    supersample_only = F
    independent_samples = F
    poisson_model = F
    sensitive_mode = F
    sensitive_mode_bwd_default = c(0.001, 0)    #Normal and Poisson
    sensitive_mode_bwd_user = NA
    setOption("mc.cores", max(detectCores() - 1, 1))

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--supersample", "-super")) {
            #input supersample table
            i = i + 1
            in_file_supersample = args[i]
        } else if (args[i] %in% c("--single-sample", "-single")) {
            #input single sample table
            i = i + 1
            in_file_single_sample = args[i]
        } else if (args[i] %in% c("--mutation-class", "-m")) {
            #type (SBS, DBS or ID)
            i = i + 1
            mut_class = args[i]
        } else if (args[i] %in% c("--output-directory", "-o")) {
            #output directory
            i = i + 1
            out_dir = args[i]
        } else if (args[i] %in% c("--regex-normals", "-r")) {
            #match-regex string for capturing normal samples
            i = i + 1
            regex_normals = args[i]
            if (is.na(regex_normals)) regex_normals = ""
        } else if (args[i] %in% c("--regex-dataset-selection", "-r-sel")) {
            #match-regex string for capturing patients to ignore in dataset selection step
            i = i + 1
            regex_dataset_sel = args[i]
            if (is.na(regex_dataset_sel)) regex_dataset_sel = ""
        } else if (args[i] %in% c("--no-dataset-selection", "-n")) {
            #whether to do dataset forward selection
            dataset_selection = F
        } else if (args[i] %in% c("--reference-signatures", "-ref-sigs")) {
            #custom reference signatures path
            i = i + 1
            ref_sig_path = args[i]
        } else if (args[i] %in% c("--common-signatures", "-c-sigs")) {
            #additional common signatures
            i = i + 1
            common_sigs = strsplit(args[i], ",")[[1]]
            if (length(common_sigs) == 0) common_sigs = NULL
        } else if (args[i] %in% c("--forced-signatures", "-f")) {
            #custom forced signatures
            i = i + 1
            sigs_rules_forced = strsplit(args[i], ",")[[1]]
            if (length(sigs_rules_forced) == 0) sigs_rules_forced = NULL
        } else if (args[i] %in% c("--supersample-only")) {
            #analyse only supersamples, ignore single samples altogether
            supersample_only = T
        } else if (args[i] %in% c("--independent-samples")) {
            #analyse only samples as independent patientd
            independent_samples = T
        } else if (args[i] %in% c("--use-poisson-model", "-p")) {
            #use Poisson model for attribution
            poisson_model = T
        } else if (args[i] %in% c("--sensitive-mode", "-sensitive")) {
            #use more sensitive model selection thresholds (relaxed backward step, allow forward in single samples)
            sensitive_mode = T
        } else if (args[i] %in% c("--sensitive-backward-threshold", "-sens-bwd-thres")) {
            #sensitive mode's single-to-supersample backward threshold
            i = i + 1
            sensitive_mode_bwd_user = as.numeric(args[i])
        } else if (args[i] %in% c("--cores", "-c")) {
            #number of cores to use
            i = i + 1
            setOption("mc.cores", max(min(as.integer(args[i]), detectCores() - 1), 1))
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (!exists("in_file_supersample") | !exists("in_file_single_sample") | !exists("out_dir") | !exists("mut_class")) {
        stop(
            "Please supply supersample and single sample input spectra table files, mutation class (SBS/DBS/ID) and output directory",
            call. = F
        )
    }

# --------------------------------------------------
# Prepare reference signatures
# --------------------------------------------------

#Load custom signatures if specified
if (exists("ref_sig_path")) {
    ref_sigs = loadSignatures(ref_sig_path, mut_class)
} else {
    #Choose correct default reference signatures based on mutation class
    if (mut_class == "SBS") ref_sigs = ref_sigs_sbs
    if (mut_class == "DBS") ref_sigs = ref_sigs_dbs
    if (mut_class == "ID") ref_sigs = ref_sigs_id
    if (mut_class == "COSMIC") ref_sigs = ref_sigs_cosmic
}

#Choose signature attribution rules based on mutation class
    #Remove invalid forced signatures (unavailable in reference signatures)
    if (exists("sigs_rules_forced")) {
        sigs_rules_forced = sigs_rules_forced[sigs_rules_forced %in% colnames(ref_sigs)]
        if (length(sigs_rules_forced) == 0) sigs_rules_forced = NULL
    }

    #Defaults for SBS and COSMIC
    if (mut_class == "SBS") {
        if (! exists("sigs_rules_forced")) sigs_rules_forced = c("SBS1", "SBS5")
        sigs_rules_conn = sigs_sbs_rules_conn
        sigs_rules_tsb = NULL #sigs_sbs_rules_tsb
        sigs_rules_mc = NULL #sigs_sbs_rules_mc_file
    } else if (mut_class == "COSMIC") {
        if (! exists("sigs_rules_forced")) sigs_rules_forced = c("Signature_1", "Signature_5")
        sigs_rules_conn = list("APOBEC" = c("Signature_2", "Signature_13"))
    }

    #Defaulting to otherwise NULL
    if (! exists("sigs_rules_forced")) sigs_rules_forced = NULL
    if (! exists("sigs_rules_conn")) sigs_rules_conn = NULL
    if (! exists("sigs_rules_tsb")) sigs_rules_tsb = NULL
    if (! exists("sigs_rules_mc")) sigs_rules_mc = NULL

# --------------------------------------------------
# Read in the mutational spectra
# --------------------------------------------------

#Read the spectra
    #Supersample
    data_supersample = read.table(in_file_supersample, as.is=T, header=T, check.names=F)
    #Single sample
    data_single_sample = read.table(in_file_single_sample, as.is=T, header=T, check.names=F)

    #Remove blood samples from single sample table
    data_single_sample = data_single_sample[
        grep(
            regex_normals, rownames(data_single_sample),
            invert = nchar(regex_normals) > 0
        ),
        ,
        drop = F
    ]

    #Strip annotations
    catal_supersample = t(data_supersample[,rownames(ref_sigs), drop=F])
    catal_single_sample = t(data_single_sample[,rownames(ref_sigs), drop=F])

#Let single samples behave as supersamples if such option is used
    if (independent_samples) {
        catal_supersample = catal_single_sample
    }

#TODO Exit without error if there are no mutations
    if (sum(catal_supersample) == 0) {
        cat(paste("Skipping", mut_class, "class due to zero mutations!\n"))
        quit("no", 0, F)
    }

# --------------------------------------------------
# Selection model preparation
# --------------------------------------------------

#Default is non-Poisson i.e. least-squares with cosine similarity
    mdl = list()
    mdl$thres_block_fwd = Inf
    if (!poisson_model) {
        #Select signatures based on cosine similarity
        #Optimise with SSE
        mdl$f_sel_stat = selStatAcc
        mdl$f_stop_stat = selStatAcc
        mdl$method = sigSSEOptim
        #Forward and backward thresholds for accuracy i.e. cosine similarity
        mdl$thres_dataset_fwd = 0.01
        mdl$thres_bwd = 0.01
        mdl$thres_fwd = 0.05
        mdl$thres_single_to_super_bwd = mdl$thres_bwd
    } else {
        #Select and keep signatures based on BIC using Poisson log-likelihood
        #Optimise with generalised Kullback-Leibler
        mdl$f_sel_stat = poisLLBIC
        mdl$f_stop_stat = poisLLBIC
        mdl$method = sigKLOptim
        #Forward and backward thresholds for BIC
        mdl$thres_dataset_fwd = 2
        mdl$thres_bwd = 2
        mdl$thres_fwd = 10
        mdl$thres_single_to_super_bwd = mdl$thres_bwd
    }
    #Sensitive mode allows forward selection at final single sample selection
    #It also relaxes single sample to supersample backward selection
    if (sensitive_mode) {
        mdl$thres_block_fwd = mdl$thres_fwd

        #Use user supplied value if available, else default
        mdl$thres_single_to_super_bwd = ifelse(
            is.na(sensitive_mode_bwd_user),
            sensitive_mode_bwd_default[1 + poisson_model],
            sensitive_mode_bwd_user
        )
    }

# --------------------------------------------------
# Data set signatures by forward selection
# --------------------------------------------------

#Perform forward selection to get data set signatures
    if (dataset_selection) {
        #Restrict to supersamples with minimum mutation count of the number of unique mutation types
        #and ignore patients selected with regex
        ignored_supersamples = colSums(catal_supersample) < nrow(ref_sigs)
        if (nchar(regex_dataset_sel) > 0) {
            ignored_supersamples = ignored_supersamples | grepl(regex_dataset_sel, colnames(catal_supersample))
        }

        #Skip if there are no acceptable supersamples
        if (sum(!ignored_supersamples) > 0) {
            #Forward select common signatures
            dataset_exps = fitSigFwd(
                catal_supersample[, !ignored_supersamples, drop = F],
                ref_sigs, stop_thres = mdl$thres_dataset_fwd,
                f_sel_stat = mdl$f_sel_stat, f_stop_stat = mdl$f_stop_stat,
                method = mdl$method
            )

            #Remove signatures with zero attributed mutations
            dataset_sigs = rownames(dataset_exps)[rowSums(dataset_exps) > 0]
        } else {
            dataset_sigs = c()
        }

        #Add input common signatures
        dataset_sigs = unique(c(dataset_sigs, colnames(ref_sigs)[colnames(ref_sigs) %in% common_sigs]))
    } else {
        #Use input common signatures (or empty)
        dataset_sigs = colnames(ref_sigs)[colnames(ref_sigs) %in% common_sigs]
    }

# --------------------------------------------------
# Supersample exposures
# --------------------------------------------------

#Attribution for each supersample via backward-forward,
#followed by its single sample backward-forward and
#finally backward (rules enforced) on supersample using combined signatures of its single samples
    #Get patients
    patients = colnames(catal_supersample)

    #Remove patients with no mutations
    patients = patients[colSums(catal_supersample) > 0]

    #Fit supersamples
    exps_supersample = mclapply(
        patients,
        function(patient) {
            #Supersample Backward-forward from data set signatures
            patient_fit_orig = fitSigBF(
                catal_supersample[, patient, drop=F], ref_sigs, dataset_sigs,
                forced_sigs = sigs_rules_forced, connected_sigs = sigs_rules_conn,
                bwd_thres = mdl$thres_bwd, fwd_thres = mdl$thres_fwd,
                f_sel_stat = mdl$f_sel_stat, f_stop_stat = mdl$f_stop_stat,
                method = mdl$method
            )

            #Return if analysing samples independently or supersamples only
            if (supersample_only | independent_samples) {
                #Remove 0-contribution signatures
                return(patient_fit_orig[patient_fit_orig > 0,, drop=F])
            }

            #Single sample backward-forward from supersample signatures
            #Restrict to samples with minimum mutation count of the number of unique mutation types
            samples_of_patient = rownames(data_single_sample)[data_single_sample$patient == patient]
            samples_of_patient = samples_of_patient[
                colSums(catal_single_sample[, samples_of_patient, drop=F]) >= nrow(ref_sigs)
            ]

            #If there are no samples passing mutation count condition, return the supersample exposure
            if (length(samples_of_patient) == 0) return(patient_fit_orig)

            sample_fit_orig = lapply(
                samples_of_patient,
                function(sample) {
                    fitSigBF(
                        catal_single_sample[, sample, drop=F], ref_sigs,
                        unique(c(
                            rownames(patient_fit_orig)[patient_fit_orig > 0],
                            dataset_sigs
                        )),
                        forced_sigs = sigs_rules_forced, connected_sigs = sigs_rules_conn,
                        bwd_thres = mdl$thres_bwd, fwd_thres = mdl$thres_fwd,
                        f_sel_stat = mdl$f_sel_stat, f_stop_stat = mdl$f_stop_stat,
                        method = mdl$method
                    )
                }
            )

            #Join sample exposures
            exps_joined = joinSigExposures(sample_fit_orig, ref_sigs)

            #Supersample backward with rules from joined sample signatures
            exps_out = fitSigBF(
                catal_supersample[, patient, drop=F], ref_sigs, rownames(exps_joined)[rowSums(exps_joined) > 0],
                forced_sigs = sigs_rules_forced, connected_sigs = sigs_rules_conn,
                bwd_thres = mdl$thres_single_to_super_bwd, fwd_thres=mdl$thres_block_fwd, bwd_stop_n=1,
                f_sel_stat = mdl$f_sel_stat, f_stop_stat = mdl$f_stop_stat,
                method = mdl$method
            )

            #Remove 0-contribution signatures
            exps_out[exps_out > 0,, drop=F]
        }
    )

    #Add patient names
    names(exps_supersample) = patients

#Join supersample exposures
    exps_supersample_joined = joinSigExposures(exps_supersample, ref_sigs)

# --------------------------------------------------
# Single sample backward selection
# --------------------------------------------------

#Attribution for each sample backward from its respective supersample's signatures
#Skipped if analysing samples independently or supersamples only
    if (!supersample_only & !independent_samples) {
        #Fit each patient's samples
        exps_samples_by_patient = mclapply(
            patients,
            function(patient) {
                #Samples of the patient
                samples_of_patient = rownames(data_single_sample)[data_single_sample$patient == patient]
                samples_of_patient = samples_of_patient[
                    colSums(catal_single_sample[, samples_of_patient, drop=F]) > 0
                ]

                #Single sample backward with rules from supersample signatures
                exps_by_sample = lapply(
                    samples_of_patient,
                    function(sample) {
                        #Limit forward selection in sensitive mode to samples with minimum
                        #mutation count of the number of unique mutation types
                        fwd_thres = ifelse(sum(catal_single_sample[,sample]) >= nrow(ref_sigs), mdl$thres_block_fwd, Inf)

                        fitSigBF(
                            catal_single_sample[, sample, drop=F], ref_sigs, rownames(exps_supersample[[patient]]),
                            forced_sigs = sigs_rules_forced, connected_sigs = sigs_rules_conn,
                            bwd_thres = mdl$thres_bwd, fwd_thres=mdl$thres_block_fwd, bwd_stop_n=1,
                            f_sel_stat = mdl$f_sel_stat, f_stop_stat = mdl$f_stop_stat,
                            method = mdl$method
                        )
                    }
                )

                #Join sample exposures
                exps_out = joinSigExposures(exps_by_sample, ref_sigs)

                #Remove 0-contribution signatures
                exps_out[rowSums(exps_out) > 0,, drop=F]
            }
        )
        names(exps_samples_by_patient) = patients

        #Join single sample exposures
        exps_single_sample_joined = joinSigExposures(exps_samples_by_patient, ref_sigs)
    }

#In independent samples, copy supersample results to single sample to keep patient grouping for plotting purposes
    if (independent_samples) {
        exps_single_sample_joined = exps_supersample_joined
    }

# --------------------------------------------------
# Write results
# --------------------------------------------------

#Function for transforming exposure table into written output table
    prepareOutputTable = function(exps, catal, prop=F, catal_annot=NULL) {
        #Convert to proportions if desired
        exps_conv = exps
        if (prop) {
            exps_conv = exps_conv / rep(colSums(exps), each=nrow(exps))
        }

        #Transpose
        exps_transp = t(exps_conv)

        #Prepare extra annotations
        if (!is.null(catal_annot)) {
             annotations = catal_annot[
                rownames(exps_transp),
                !colnames(catal_annot) %in% rownames(ref_sigs),
                drop = F
            ]
        } else {
            annotations = NULL
        }

        #Add sample_id, mut_count, accuracy and other annotation columns
        data.frame(
            cbind(
                "sample_id" = colnames(exps),
                annotations
            ),
            "mut_count" = colSums(exps),
            "cos_sim" = sapply(
                colnames(exps),
                function(sample) {
                    1 - selStatAcc(
                        catal[,sample, drop=F],
                        ref_sigs[,rownames(exps), drop=F],
                        exps[,sample, drop=F]
                    )
                }
            ),
            exps_transp
        )
    }

#Write table function
    writeTable = function(table_to_write, file_name) {
        write.table(
            table_to_write,
            file_name,
            quote = F,
            sep = "\t",
            row.names = F
        )
    }

#Write data set common signatures
    write(dataset_sigs, paste0(out_dir, "common_signatures_", tolower(mut_class), ".txt"))

#Format and write output tables with mutation count and accuracy (1 - cosine similarity)
    #Supersample, original frequencies with '_supersample' suffix in column names
    out_table = prepareOutputTable(exps_supersample_joined, catal_supersample, F)
    out_table[,1] = paste0(out_table[,1], "_supersample")
    writeTable(
        out_table,
        paste0(out_dir, "contributions_", tolower(mut_class), "_supersample.tsv")
    )
    #Supersample, proportions with '_supersample' suffix in column names
    out_table = prepareOutputTable(exps_supersample_joined, catal_supersample, T)
    out_table[,1] = paste0(out_table[,1], "_supersample")
    writeTable(
        out_table,
        paste0(out_dir, "contributions_", tolower(mut_class), "_supersample_proportions.tsv")
    )

    #Single samples, skipped when analysing only supersamples
    if (!supersample_only | independent_samples) {
        #original frequencies
        writeTable(
            prepareOutputTable(exps_single_sample_joined, catal_single_sample, F, data_single_sample),
            paste0(out_dir, "contributions_", tolower(mut_class), "_single_sample.tsv")
        )
        #proportions
        writeTable(
            prepareOutputTable(exps_single_sample_joined, catal_single_sample, T, data_single_sample),
            paste0(out_dir, "contributions_", tolower(mut_class), "_single_sample_proportions.tsv")
        )
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


