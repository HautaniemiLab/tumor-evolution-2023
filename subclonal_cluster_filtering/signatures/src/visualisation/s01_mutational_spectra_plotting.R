#!/usr/bin/env Rscript

#R script for plotting mutational profiles from SBS, DBS and ID spectra.
#Normal samples are from outputs single sample specific outputs.

# --------------------------------------------------
# Sourcing
# --------------------------------------------------

#Load reference signatures
source("src/features/f0_load_reference_signatures.R")

#Load plotting functions
source("src/visualisation/f2_signature_plotting_functions.R")

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include:  SBS, DBS and ID input spectra table files, output directory,
#output file suffix.
#Possible option for whether to put patients in separate plots and normal sample
#regex match-string
    #Default options
    plot_separate = F
    regex_normals = ""

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--sbs-input", "-sbs")) {
            #input sbs file
            i = i + 1
            in_file_sbs = args[i]
        } else if (args[i] %in% c("--dbs-input", "-dbs")) {
            #input dbs file
            i = i + 1
            in_file_dbs = args[i]
        } else if (args[i] %in% c("--id-input", "-id")) {
            #input id file
            i = i + 1
            in_file_id = args[i]
        } else if (args[i] %in% c("--output-directory", "-o")) {
            #output directory
            i = i + 1
            out_dir = args[i]
        } else if (args[i] %in% c("--output-suffix", "-out-suf")) {
            #output suffix
            i = i + 1
            out_suf = args[i]
        } else if (args[i] %in% c("--plot-separate-patients", "-p")) {
            #whether to separate patients in output plots
            plot_separate = T
        } else if (args[i] %in% c("--regex-normals", "-r")) {
            #match-regex string for capturing normal samples
            i = i + 1
            regex_normals = args[i]
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (!exists("in_file_sbs") | !exists("in_file_dbs") | !exists("in_file_id") | !exists("out_dir") | !exists("out_suf")) {
        stop(
            "Please supply SBS, DBS and ID input spectra table files, output directory and output file suffix",
            call. = F
        )
    }

# --------------------------------------------------
# Read in the mutational spectra
# --------------------------------------------------

#Function for transforming annotated catalogue table to correctly oriented catalogue table without annotations
    catalogueStripAnnotations = function(catal, mut_class="SBS") {
        if (mut_class == "SBS") mut_types = rownames(ref_sigs_sbs)
        if (mut_class == "DBS") mut_types = rownames(ref_sigs_dbs)
        if (mut_class == "ID") mut_types = rownames(ref_sigs_id)

        t(catal[,mut_types])
    }

#Read the spectra and strip annotations
    #SBS
    data_sbs = read.table(in_file_sbs, as.is=T, header=T, check.names=F)
    catalogues_sbs = catalogueStripAnnotations(data_sbs, "SBS")
    #DBS
    data_dbs = read.table(in_file_dbs, as.is=T, header=T, check.names=F)
    catalogues_dbs = catalogueStripAnnotations(data_dbs, "DBS")
    #ID
    data_id = read.table(in_file_id, as.is=T, header=T, check.names=F)
    catalogues_id = catalogueStripAnnotations(data_id, "ID")

# --------------------------------------------------
# Plotting helper functions
# --------------------------------------------------

#SBS/DBS/ID stackplot drawing function
    plotSignaturesMultitypeStackPlot = function(catal_sbs, catal_dbs, catal_id, name_label="") {
        #Original plot settings
        orig_par = par(c("mfrow", "cex"))

        #SBS
        plotSignaturesSBS(
            catal_sbs,
            paste(name_label, "SBS, total:", sum(catal_sbs)),
            thin_plot = T
        )

        #DBS
        plotSignaturesDBS(
            catal_dbs,
            paste(name_label, "DBS, total:", sum(catal_dbs)),
            thin_plot = T
        )

        #ID
        plotSignaturesID(
            catal_id,
            paste(name_label, "ID, total:", sum(catal_id)),
            thin_plot = T
        )

        #Return original par values
        par(orig_par)
    }

# --------------------------------------------------
# Draw patient profiles
# --------------------------------------------------

#Draw all profiles in one file or for separate patients:
    if (!plot_separate) {
        #Single pdf output file; assumes patient per column
        pdf(paste0(out_dir, out_suf, ".pdf"), h=12, w=14)
            par(mfrow=c(3,1), cex=1)

            #Draw stacked SBS/DBS/ID profiles
            for (i in 1:ncol(catalogues_sbs)) {
                #Extract patient name
                patient_name = colnames(catalogues_sbs)[i]

                #Draw plot
                plotSignaturesMultitypeStackPlot(
                    catalogues_sbs[,i, drop=F], catalogues_dbs[,i, drop=F], catalogues_id[,i, drop=F],
                    patient_name
                )
            }
        invisible(dev.off())
    } else {
        #Find normal samples if normal regex is specified
        if (nchar(regex_normals) > 0 & !is.na(regex_normals)) {
            normals = grepl(regex_normals, colnames(catalogues_sbs))
        } else {
            normals = rep(F, ncol(catalogues_sbs))
        }

        #Plot normal sample profiles in a single file, if specified to
        if (any(normals)) {
            pdf(paste0(out_dir, "normal_", out_suf, ".pdf"), h=12, w=14)
                par(mfrow=c(3,1), cex=1)

                #Draw stacked SBS/DBS/ID profiles
                for (i in which(normals)) {
                    #Extract sample name
                    sample_name = colnames(catalogues_sbs)[i]

                    #Draw plot
                    plotSignaturesMultitypeStackPlot(
                        catalogues_sbs[,i, drop=F], catalogues_dbs[,i, drop=F], catalogues_id[,i, drop=F],
                        sample_name
                    )
                }
            invisible(dev.off())
        }

        #Extract patient names from sample names
        sample_of_patient = data_sbs$patient

        #Group SBS/DBS/ID catalogues by patient, ignoring normals
        catalogues_by_patient = lapply(
            unique(sample_of_patient),
            function(patient) {
                #Find tumour samples belonging to patient
                samples_ind = (sample_of_patient == patient) & !normals

                #Gather SBS/DBS/ID catalogues of patient in a list
                out = list()
                out$sbs = catalogues_sbs[,samples_ind, drop=F]
                out$dbs = catalogues_dbs[,samples_ind, drop=F]
                out$id = catalogues_id[,samples_ind, drop=F]

                out
            }
        )
        names(catalogues_by_patient) = unique(sample_of_patient)

        #One pdf output file per patient
        for (i in 1:length(catalogues_by_patient)) {
            #Extract patient name
            patient_name = names(catalogues_by_patient)[i]

            #Start plotting the patient
            pdf(paste0(out_dir, patient_name, "_", out_suf, ".pdf"), h=12, w=14)
                par(mfrow=c(3,1), cex=1)

                #Do nothing if there are no samples/columns for patient
                if (ncol(catalogues_by_patient[[i]]$sbs) > 0) {
                    #Draw stacked SBS/DBS/ID profiles for each sample
                    for (j in 1:ncol(catalogues_by_patient[[i]]$sbs)) {
                        #Extract sample name
                        sample_name = colnames(catalogues_by_patient[[i]]$sbs)[j]

                        #Draw plot
                        plotSignaturesMultitypeStackPlot(
                            catalogues_by_patient[[i]]$sbs[,j, drop=F],
                            catalogues_by_patient[[i]]$dbs[,j, drop=F],
                            catalogues_by_patient[[i]]$id[,j, drop=F],
                            sample_name
                        )
                    }
                }
            invisible(dev.off())
        }
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


