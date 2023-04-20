#!/usr/bin/env Rscript

#R script for plotting signature exposures and fit plots.
#Plots both supersample and single samples.

# --------------------------------------------------
# Sourcing
# --------------------------------------------------

#Load reference signatures
source("src/features/f0_load_reference_signatures.R")

#Load attribution functions
source("src/attribution/f1_signature_attribution_functions.R")

#Load plotting functions
source("src/visualisation/f2_signature_plotting_functions.R")

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include: supersample and single sample input spectra table files,
#supersample and single sample input signature contribution table files,
#mutation class (SBS/DBS/ID) and output directory
#Possible options for custom reference signatures, e.g. exome signatures

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--supersample-catalogues", "-super-catal")) {
            #input supersample catalogues table
            i = i + 1
            in_file_catal_supersample = args[i]
        } else if (args[i] %in% c("--single-sample-catalogues", "-single-catal")) {
            #input single sample catalogues table
            i = i + 1
            in_file_catal_single_sample = args[i]
        } else if (args[i] %in% c("--supersample-exposures", "-super-exp")) {
            #input supersample exposures table
            i = i + 1
            in_file_exps_supersample = args[i]
        } else if (args[i] %in% c("--single-sample-exposures", "-single-exp")) {
            #input single sample exposures table
            i = i + 1
            in_file_exps_single_sample = args[i]
        } else if (args[i] %in% c("--mutation-class", "-m")) {
            #type (SBS, DBS or ID)
            i = i + 1
            mut_class = args[i]
        } else if (args[i] %in% c("--output-directory", "-o")) {
            #output directory
            i = i + 1
            out_dir = args[i]
        } else if (args[i] %in% c("--reference-signatures", "-ref-sigs")) {
            #custom reference signatures path
            i = i + 1
            ref_sig_path = args[i]
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (
        !exists("in_file_catal_supersample") | !exists("in_file_catal_single_sample") |
        !exists("in_file_exps_supersample") | !exists("in_file_exps_single_sample") |
        !exists("out_dir") | !exists("mut_class")
    ) {
        stop(
            "Please supply supersample and single sample input spectra table files, supersample and single sample input signature contribution table files, mutation class (SBS/DBS/ID) and output directory",
            call. = F
        )
    }

# --------------------------------------------------
# Read in data
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

#Choose correct plotting parameters based on mutation class
    exp_plot.mar.4.shift = ifelse(mut_class == "COSMIC", 5, 3)

#Read the spectra
    #Supersample
    data_catal_supersample = read.table(in_file_catal_supersample, as.is=T, header=T, check.names=F)
    #Single sample
    data_catal_single_sample = read.table(in_file_catal_single_sample, as.is=T, header=T, check.names=F)

    #Strip annotations from the spectra
    catal_supersample = t(data_catal_supersample[,rownames(ref_sigs), drop=F])
    catal_single_sample = t(data_catal_single_sample[,rownames(ref_sigs), drop=F])

#Read the exposures transformed to correct format with signatures on rows
    #Supersample
    data_exps_supersample = read.table(in_file_exps_supersample, as.is=T, header=T, row.names=1, check.names=F)

    exps_supersample = t(as.matrix(
        data_exps_supersample[,colnames(data_exps_supersample) %in% colnames(ref_sigs), drop=F]
    ))

    #Single sample
    if (file.exists(in_file_exps_single_sample)) {
        data_exps_single_sample = read.table(in_file_exps_single_sample, as.is=T, header=T, row.names=1, check.names=F)

        exps_single_sample = t(as.matrix(
            data_exps_single_sample[,colnames(data_exps_single_sample) %in% colnames(ref_sigs), drop=F]
        ))
    }

#Remove suffix from supersamples
    colnames(exps_supersample) = gsub("_.*", "", colnames(exps_supersample))

#Patients
    patients = unique(colnames(exps_supersample))

#Set y-axis scale to thousands if there are samples with tens of thousands of mutations
    exp_plot.ylab.scale = ifelse(max(colSums(catal_single_sample)) >= 10000, 1000, 1)

# --------------------------------------------------
# Exposure plots
# --------------------------------------------------

#Prepare signature colours
    #Observed signatures from supersample and single sample signatures
    observed_signatures = unique(c(
        rownames(exps_supersample),
        if (exists("exps_single_sample")) rownames(exps_single_sample)
    ))
    observed_signatures = colnames(ref_sigs)[colnames(ref_sigs) %in% observed_signatures]

    #Determine colours
    col_sigs = colorRampPalette(
        brewer.pal(
            max(
                min(length(observed_signatures), 9),
                3
            ),
            "Set1"
        )
    )( length(observed_signatures) )
    names(col_sigs) = observed_signatures

#Draw supersample signature exposures plot
    pdf(paste0(out_dir, "contributions_", tolower(mut_class), "_supersample.pdf"), w=8, h=6)
        #Copy supersample exposures
        exps_supersample_copy = exps_supersample

        #Remove names if there are more than ~100 samples
        if (ncol(exps_supersample_copy) > 100) {
            colnames(exps_supersample_copy) = NULL
        }

        plotSigExposures(
            exps_supersample_copy, col_sigs,
            title = paste("Supersamples", mut_class, "Mutations Contributed", sep = " - "),
            ylab.scale=exp_plot.ylab.scale, lab.cex=min(.7, 35 / ncol(exps_supersample_copy)),
            lab.srt=60, lab.adj=c(.5,0),
            lab.x.shift = -.0015 * ncol(exps_supersample_copy),
            lab.y.shift = -min(.04, .08 - .0008 * ncol(exps_supersample_copy)),
            mar.4.shift = exp_plot.mar.4.shift
        )
    invisible(dev.off())

#Single sample exposure plots
    if (exists("exps_single_sample")) {
        #Draw single sample signature exposure plots for each patient
        pdf(paste0(out_dir, "contributions_", tolower(mut_class), "_samples_by_patient.pdf"), w=8, h=6)
            for (patient in patients) {
                samples_of_patient = rownames(data_exps_single_sample)[data_exps_single_sample$patient == patient]

                #Select only signatures observed in patient
                patient_sample_exps = exps_single_sample[,samples_of_patient, drop=F]
                patient_sample_exps = patient_sample_exps[rowSums(patient_sample_exps) > 0,, drop=F]

                #Keep only site info
                colnames(patient_sample_exps) = sub(
                    paste0("^", patient, "_"), "",
                    sub("_DNA\\d*", "", colnames(patient_sample_exps))
                )
                #Determine suitable font size
                axis_cex = min(
                    1,
                    75 / max(sapply(colnames(patient_sample_exps), nchar)) / (ncol(patient_sample_exps) * 1.25 + .2)
                )

                #Draw exposures
                plotSigExposures(
                    patient_sample_exps, col_sigs,
                    title = paste(paste(patient, "samples"), mut_class, "Mutations Contributed", sep = " - "),
                    ylab.scale=exp_plot.ylab.scale, lab.cex=axis_cex, lab.srt=0, lab.adj=c(.5,.5),
                    lab.x.shift=0, lab.y.shift=0, mar.4.shift=exp_plot.mar.4.shift
                )
            }
        invisible(dev.off())

    #Draw single sample signature exposure plots in one plot
        #Copy single sample exposures
        exps_single_sample_copy = exps_single_sample

        #Remove names if there are more than ~50 samples
        if (ncol(exps_single_sample_copy) > 50) {
            colnames(exps_single_sample_copy) = NULL
        } else if (ncol(exps_supersample) == 1) {
            #Strip the patient prefix if there is only one patient
            colnames(exps_single_sample_copy) = sub(
                paste0("^", patient, "_"), "", colnames(exps_single_sample_copy)
            )
        } else {
            #Remove "DNA\d" suffix from sample names
            colnames(exps_single_sample_copy) = sub("_DNA\\d*", "", colnames(exps_single_sample_copy))
        }

        #Determine suitable font size
        axis_cex = ifelse(
            !is.null(colnames(exps_single_sample_copy)),
            min(.5, 4 / max(sapply(colnames(exps_single_sample_copy), nchar), 1)),
            .5
        )

        #Draw exposures
        pdf(paste0(out_dir, "contributions_", tolower(mut_class), "_all_samples.pdf"), w=8, h=6)
            plotSigExposures(
                exps_single_sample_copy, col_sigs,
                title = paste("Single samples", mut_class, "Mutations Contributed", sep = " - "),
                ylab.scale=exp_plot.ylab.scale, lab.cex=axis_cex, lab.srt=60, lab.adj=c(1,.5),
                lab.x.shift=0.1, lab.y.shift=0.03, mar.4.shift=exp_plot.mar.4.shift
            )
        invisible(dev.off())
    }

# --------------------------------------------------
# Signature fit plots
# --------------------------------------------------

#Maximum number of signatures among patients
exp_rows_super = max(colSums(exps_supersample > 0))

#Draw supersample signature fit barplot
    pdf(paste0(out_dir, "fit_barplot_", tolower(mut_class), "_supersamples.pdf"), w=11, h=9)
        for (patient in patients) {
            #Select only signatures observed in patient
            patient_supersample_exps = exps_supersample[,patient, drop=F]
            patient_supersample_exps = patient_supersample_exps[rowSums(patient_supersample_exps) > 0,, drop=F]

            #Draw signature fit plot
            plotSigFitSampleLarge(
                catal_supersample[,patient, drop=F],
                ref_sigs[,rownames(patient_supersample_exps), drop=F],
                patient_supersample_exps,
                col_sigs,
                main_text = paste(
                    patient, "supersample -", mut_class, "profile -",
                    sum(data_catal_single_sample$patient == patient), "sample(s) -",
                    sum(catal_supersample[,patient, drop=F]), "mutations"
                ),
                type = mut_class,
                min_exp_rows = max(6, exp_rows_super),
                draw_signatures = T
            )
        }
    invisible(dev.off())

#Draw single sample signature fit barplots by patient
    if (exists("exps_single_sample")) {
        #Make sure target directory for the patient files exists
        out_dir_sample_fit = paste0(out_dir, "sample_fit_plots/")
        if (! dir.exists(out_dir_sample_fit)) {
            dir.create(out_dir_sample_fit)
        }

        #Loop over patients
        for (patient in patients) {
            #Get samples of the patient
            samples_of_patient = rownames(data_exps_single_sample)[data_exps_single_sample$patient == patient]

            #Prepare PDF output with all samples of the patient
            pdf(paste0(out_dir_sample_fit, "fit_barplot_", patient, "_", tolower(mut_class), "_samples.pdf"), w=11, h=8)
                for (sample in samples_of_patient) {
                    #Select only signatures observed in patient
                    sample_exps = exps_single_sample[,sample, drop=F]
                    sample_exps = sample_exps[rowSums(sample_exps) > 0,, drop=F]

                    #Draw signature fit plot
                    plotSigFitSampleLarge(
                        catal_single_sample[,sample, drop=F],
                        ref_sigs[,rownames(sample_exps), drop=F],
                        sample_exps,
                        col_sigs,
                        main_text = paste0(
                            sample, " - ", mut_class, " profile - ",
                            sum(catal_single_sample[,sample, drop=F]), " mutations"
                        ),
                        type = mut_class
                    )
                }
            invisible(dev.off())
        }
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


