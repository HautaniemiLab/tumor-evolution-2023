#!/usr/bin/env Rscript

#R script for plotting subclonal cluster signature exposures and fit plots.
#Plots both supersample and single subclones.

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
# Packages
# --------------------------------------------------

#Library for drawing circles
suppressMessages(require(plotrix))

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include: input supersample and subclone input spectra table files,
#supersample and subclonal cluster input signature contribution table files,
#mutation class (SBS/DBS/ID) and output directory
#Possible options for custom reference signatures, e.g. exome signatures,
#and whether cluster indices start from zero
    #Default options
    has_cluster_zero = F

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--supersample-catalogues", "-super-catal")) {
            #input supersample catalogues table
            i = i + 1
            in_file_catal_supersample = args[i]
        } else if (args[i] %in% c("--subclone-catalogues", "-subclone-catal")) {
            #input subclone catalogues table
            i = i + 1
            in_file_catal_subclone = args[i]
        } else if (args[i] %in% c("--supersample-exposures", "-super-exp")) {
            #input supersample exposures table
            i = i + 1
            in_file_exps_supersample = args[i]
        } else if (args[i] %in% c("--subclone-exposures", "-subclone-exp")) {
            #input subclone exposures table
            i = i + 1
            in_file_exps_subclone = args[i]
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
        } else if (args[i] %in% c("--zero-indexed-clusters", "-z")) {
            #whether clusters are zero indexed
            has_cluster_zero = T
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (
        !exists("in_file_catal_supersample") | !exists("in_file_catal_subclone") |
        !exists("in_file_exps_supersample") | !exists("in_file_exps_subclone") |
        !exists("out_dir") | !exists("mut_class")
    ) {
        stop(
            "Please supply supersample and subclone input spectra table files, supersample and subclone input signature contribution table files, mutation class (SBS/DBS/ID) and output directory",
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
    #Single cluster
    data_catal_subclone = read.table(in_file_catal_subclone, as.is=T, header=T, check.names=F)

    #Strip annotations from the spectra
    catal_supersample = t(data_catal_supersample[,rownames(ref_sigs), drop=F])
    catal_subclone = t(data_catal_subclone[,rownames(ref_sigs), drop=F])

#Read the exposures
    #Supersample
    data_exps_supersample = read.table(in_file_exps_supersample, as.is=T, header=T, row.names=1, check.names=F)
    #Subclonal clusters
    data_exps_subclone = read.table(in_file_exps_subclone, as.is=T, header=T, row.names=1, check.names=F)

# --------------------------------------------------
# Prepare exposure plot data
# --------------------------------------------------

#Extract cosine similarities
    exps_cos_sim_all = c(data_exps_supersample$cos_sim, data_exps_subclone$cos_sim)
    names(exps_cos_sim_all) = c(rownames(data_exps_supersample), rownames(data_exps_subclone))

#Transform exposure tables back to proper format with signatures on rows
    exps_supersample = t(as.matrix(
        data_exps_supersample[,colnames(data_exps_supersample) %in% colnames(ref_sigs), drop=F]
    ))
    exps_subclone = t(as.matrix(
        data_exps_subclone[,colnames(data_exps_subclone) %in% colnames(ref_sigs), drop=F]
    ))

#Remove suffix from supersamples
    names(exps_cos_sim_all) = gsub("_supersample", "", names(exps_cos_sim_all))
    colnames(exps_supersample) = gsub("_.*", "", colnames(exps_supersample))

#Set y-axis scale to thousands if there are samples with tens of thousands of mutations
    exp_plot.ylab.scale = ifelse(max(colSums(catal_subclone)) >= 10000, 1000, 1)

#Prepare signature colours
    #Observed signatures from supersample and subclone signatures
    #Those detected only in subclones use a different colour scheme
    supersample_signatures = rownames(exps_supersample)
    supersample_signatures = colnames(ref_sigs)[colnames(ref_sigs) %in% supersample_signatures]

    #Other signatures
    other_signatures = rownames(exps_subclone)[
        !rownames(exps_subclone) %in% supersample_signatures
    ]
    other_signatures = colnames(ref_sigs)[colnames(ref_sigs) %in% other_signatures]

    #Determine colours
    col_sigs = c(
        colorRampPalette(
            brewer.pal(
                max(
                    min(length(supersample_signatures), 9),
                    3
                ),
                "Set1"
            )
        )( length(supersample_signatures) ),
        colorRampPalette(
            brewer.pal(
                max(
                    min(length(other_signatures), 8),
                    3
                ),
                "Pastel2"
            )
        )( length(other_signatures) )
    )
    names(col_sigs) = c(supersample_signatures, other_signatures)

#Prepare cosine similarity colours
    col_cos_sim = colorRampPalette(brewer.pal(9,"PuBu"))(8)

# --------------------------------------------------
# Supersample exposure plot
# --------------------------------------------------

#Fill in patients with zero mutations
    if (any(! colnames(catal_supersample) %in% colnames(exps_supersample))) {
        exps_supersample_plot = cbind(
            exps_supersample,
            sapply(
                colnames(catal_supersample)[! colnames(catal_supersample) %in% colnames(exps_supersample)],
                function(patient) {
                    rep(0, nrow(exps_supersample))
                }
            )
        )
        exps_supersample_plot = exps_supersample_plot[,colnames(catal_supersample), drop=F]
    } else {
        exps_supersample_plot = exps_supersample
    }

#Draw supersample signature exposures plot
    pdf(paste0(out_dir, "contributions_", tolower(mut_class), "_supersample.pdf"), w=8, h=6)
    plotSigExposures(
        exps_supersample_plot, col_sigs,
        title = paste("Supersamples", mut_class, "Mutations Contributed", sep = " - "),
        ylab.scale=exp_plot.ylab.scale, lab.cex=0.7, lab.srt=60, lab.adj=c(1,.5),
        lab.x.shift=0.1, lab.y.shift=0.03, mar.4.shift=exp_plot.mar.4.shift
    )
    invisible(dev.off())

# --------------------------------------------------
# Functions for subclonal cluster exposure plots
# --------------------------------------------------

#Helper function for ...

# --------------------------------------------------
# Subclonal exposure plots
# --------------------------------------------------

#Make sure target directory for the patient files exists
    out_dir_subclone_exp = paste0(out_dir, "subclone_exposures/")
    if (! dir.exists(out_dir_subclone_exp)) {
        dir.create(out_dir_subclone_exp, showWarnings=F)
    }

#Loop over patients
    for (patient in colnames(exps_supersample)) {
        #Prepare exposures table for plotting
            #Select only signatures observed in patient
            patient_supersample_exps = exps_supersample[,patient, drop=F]
            patient_supersample_exps = patient_supersample_exps[rowSums(patient_supersample_exps) > 0,, drop=F]

            #Get subclones of the patient
            subclones_of_patient = rownames(data_exps_subclone)[data_exps_subclone$patient == patient]

            #Select subclones of patient
            subclone_exps = exps_subclone[,colnames(exps_subclone) %in% subclones_of_patient, drop=F]
            subclone_exps = subclone_exps[rowSums(subclone_exps) > 0,, drop=F]

            #Exposures table for plotting with all subclones in spectrum table even if missing (when they have 0 mutations)
            exposures_plot = subclone_exps
            colnames(exposures_plot) = sub(paste0(patient, "_"), "", colnames(exposures_plot))
            for (cluster in sub(paste0(patient, "_"), "", subclones_of_patient)) {
                #Add zero column if missing
                if (! cluster %in% colnames(exposures_plot)) {
                    exposures_plot = cbind(
                        exposures_plot,
                        matrix(0, nrow=nrow(subclone_exps), dimnames=list(NULL, cluster))
                    )
                }
            }

            #Add patient supersample exposures and sort the table
            exposures_plot = exposures_plot[, sub(paste0(patient, "_"), "", subclones_of_patient), drop=F]
            exposures_plot = joinSigExposures(list(patient_supersample_exps, exposures_plot), ref_sigs)

            #Cosine similarities for each case
            cos_sim_plot = exps_cos_sim_all[c(patient, colnames(subclone_exps))]
            names(cos_sim_plot) = sub(paste0(patient, "_"), "", names(cos_sim_plot))
            cos_sim_plot = cos_sim_plot[colnames(exposures_plot)]
            names(cos_sim_plot) = colnames(exposures_plot)

        #Prepare PDF output with all samples of the patient
        w = 0.66 * (ncol(exposures_plot) + 1) + 0.11 + 1.24
        h = nrow(exposures_plot) + 4
        pdf(
            paste0(out_dir_subclone_exp, "exposures_", patient, "_", tolower(mut_class), "_subclones.pdf"),
            width = w, height = h
        )
            #Layout
            par(mfrow=c(2*nrow(exposures_plot) + 7, 1), mar=c(.5, 4, 3, 2) + .1)
            layout(
                matrix(c(
                    1,
                    sapply(
                        2:(nrow(exposures_plot) + 4),
                        function(x) rep(x, 2)
                    ),
                    nrow(exposures_plot) + 5
                ))
            )

            #Title
            par(mar = c(0, 0, 0, 0))
            plot(
                1, type = "n",
                xaxt = "n", yaxt = "n", bty = "n"
            )
            title(
                paste0(patient, " subclone clusters"),
                line = -2.5,
                cex.main = 2
            )

            #Signature plots
                par(mar=c(1, 4, 1.5, 2) + .1)

                #Plot per signature contributions (mutation counts)
                for (i in 1:nrow(exposures_plot)) {
                    #Empty plot with grid lines
                    barplot(
                        rep(0, length(exposures_plot[i,])),
                        ylim = c(0, 1.2 * max(exposures_plot[i,])),
                        border = F, yaxt = "n"
                    )
                    abline(h=axTicks(2), col="#DDDDDD")
                    par(new=T)

                    #Barplot
                    barplot(
                        exposures_plot[i,],
                        xlim = c(0, 1.2 * ncol(exposures_plot) + .2),
                        ylim = c(0, 1.2 * max(exposures_plot[i,])),
                        main = paste0(sub("\\.", " ", rownames(exposures_plot)[i]), " contribution"),
                        ylab = "mutation count",
                        yaxt = "n", mgp = c(3, 0, 0), xaxs = "i",
                        col = col_sigs[rownames(exposures_plot)[i]],
                        border = ifelse(exposures_plot[i,] > 0, par("fg"), NA)
                    )
                    #y axis
                    axis(2, las=1, mgp=c(3, .75, 0))
                    #Horizontal line at 0
                    abline(h=0)
                }

                #Plot absolute signature contribution plot
                    #Empty plot with grid lines
                    x_lab_pos = as.numeric(
                        barplot(
                            rep(0, ncol(exposures_plot)),
                            ylim = c(0, 1.2 * max(colSums(exposures_plot))),
                            border = F, yaxt = "n"
                        )
                    )
                    abline(h=axTicks(2), col="#DDDDDD")
                    par(new=T)

                    #Barplot
                    barplot(
                        exposures_plot[rev(rownames(exposures_plot)),, drop=F],
                        xlim = c(0, 1.2 * ncol(exposures_plot) + .2),
                        ylim = c(0, 1.2 * max(colSums(exposures_plot))),
                        main = paste0("Total signature contributions"),
                        ylab = "mutation count",
                        yaxt = "n", mgp = c(3, 0, 0), xaxs = "i",
                        col = rev(col_sigs[rownames(exposures_plot)])
                    )
                    #y axis
                    axis(2, las=1, mgp=c(3, .75, 0))
                    #Horizontal line at 0
                    abline(h=0)

                #Plot signature exposure fraction plot
                    #Barplot
                    barplot(
                        exposures_plot[rev(rownames(exposures_plot)),, drop=F] / rep(colSums(exposures_plot), each=nrow(exposures_plot)),
                        xlim = c(0, 1.2 * ncol(exposures_plot) + .2),
                        ylim = c(0, 1.2),
                        main = paste0("Fractional signature exposures"),
                        ylab = "fraction",
                        yaxt = "n", mgp = c(3, 0, 0), xaxs = "i",
                        col = rev(col_sigs[rownames(exposures_plot)])
                    )
                    #y axis
                    axis(2, at=.2*0:5, las=1, mgp=c(3, .75, 0))
                    #Horizontal line at 0
                    abline(h=0)

            #Plot cosine similarities
                #Empty plot with grid lines
                x_lab_pos = as.numeric(
                    barplot(
                        rep(0, ncol(exposures_plot)),
                        ylim = c(.7, 1.06),
                        border = F, yaxt = "n"
                    )
                )
                abline(h=seq(.7, 1, .05), col="#DDDDDD")
                par(new=T)

                #Barplot
                barplot(
                    cos_sim_plot,
                    xlim = c(0, 1.2 * ncol(exposures_plot) + .2),
                    ylim = c(.7, 1.06),
                    main = paste0("Cosine Similarities"),
                    ylab = "cosine similarity",
                    yaxt = "n", mgp = c(3, 0, 0), xaxs = "i",
                    col = c(NA, col_cos_sim)[findInterval(cos_sim_plot, c(0, .6 + .05 * 0:8), rightmost.closed=T)],
                    xpd = F
                )
                #Horizontal line at 0
                abline(h=.7)
                #Cosine similarities as text above bars
                par(xpd=T)
                text(
                    x_lab_pos, sapply(cos_sim_plot, max, .7, na.rm=T),
                    labels = round(cos_sim_plot, 3),
                    pos = 3, offset = .25
                )
                par(xpd=F)

                #Colour code on the left side
                    par(mar=c(1, 3.5, 1.5, 2 + par("mar")[2]*par("plt")[2] / par("plt")[1] - par("mar")[2]) + .1, new=T)
                    #Rectangles with barplot
                    barplot(
                        matrix(c(.7, rep(.05, 6))),
                        xlim = c(0, 1) + .2, ylim = c(.7, 1.06),
                        xlab = "", ylab = "",
                        yaxt = "n", yaxs = "i", xaxs = "i",
                        col = c(NA, col_cos_sim[3:8]), border = F, xpd = F
                    )
                    #y axis
                    axis(2, at=seq(0.7, 1, .05), las=1, mgp=c(3, .75, 0))

                    #Box
                    rect(.2, .7, 1.2, 1)

            #Overall x label to show cluster information
                par(mar=c(1.1, 4.1, 0, 2.1))

                #Axis limits
                x_lim = c(0, 1.2 * ncol(exposures_plot) + .2)
                ratio = .5 * (par("pin")[2] - par("mai")[1]) / par("pin")[1]
                y_range = (x_lim[2] - x_lim[1]) * ratio
                y_lim = 1 + .5 * y_range * c(-1,1)

                #Empty plot
                plot(
                    1, type = "n",
                    xlim = x_lim, ylim = y_lim,
                    xaxt = "n", yaxt = "n",
                    xaxs = "i", yaxs = "i",
                    ylab = "",
                    bty = "n",
                    mgp = c(3, 0, 0), asp = 1
                )
                #Axis
                axis(1, x_lab_pos, colnames(exposures_plot), tick=F, line=-1, lwd=0)

                #Draw circles in the colours of the subclones
                r_cluster = .48 * (par("usr")[4] - par("usr")[3])
                par(xpd=T)

                    #Supersample
                    par(xpd=F)
                    draw.circle(
                        x_lab_pos[1], 1,
                        r = 1.25*r_cluster, col = "darkgray"
                    )
                    par(xpd=T)

                    #Singular clusters
                    invisible(sapply(
                        grep("^[0-9]+$", colnames(exposures_plot), value=T),
                        function(cluster) {
                            circle_col = brewer.pal(12, "Paired")[(as.integer(cluster) - 1 + has_cluster_zero) %% 12 + 1]

                            #Draw filtered clusters as circles with line fill on darkgray background and borders
                            if (data_exps_subclone[paste0(patient, "_", cluster), "filtered"]) {
                                draw.circle(
                                    x_lab_pos[colnames(exposures_plot) == cluster], 1,
                                    r = r_cluster, col = "darkgray",
                                    border = "darkgray"
                                )
                                draw.circle(
                                    x_lab_pos[colnames(exposures_plot) == cluster], 1,
                                    r = r_cluster, col = circle_col,
                                    density = 35, border = NA
                                )
                            } else {
                                draw.circle(
                                    x_lab_pos[colnames(exposures_plot) == cluster], 1,
                                    r = r_cluster, col = circle_col
                                )
                            }
                        }
                    ))

                    #Combination clusters
                    invisible(sapply(
                        grep(",", colnames(exposures_plot)),
                        function(i) {
                            #Clusters to combine, angle and radii of drawn circle
                            clusters = unlist(strsplit(colnames(exposures_plot)[i], ","))
                            angle = 2*pi / length(clusters)
                            if (length(clusters) %% 2) {
                                r_c = 8 * r_cluster / (7 + cos(angle*.5))
                                y_pos = par('usr')[4] - r_c
                            } else {
                                r_c = 4 * r_cluster / (3 + cos(angle*.5))
                                y_pos = 1
                            }
                            #Draw each transparent circle 3 times
                            for (k in 1:3) {
                                invisible(sapply(
                                    1:length(clusters),
                                    function(j) {
                                        draw.circle(
                                            x_lab_pos[i] + .25 * r_cluster * sin(pi + angle*(j-.5)),
                                            y_pos + .25 * r_cluster * cos(pi + angle*(j-.5)),
                                            r = .75 * r_cluster,
                                            col = paste0(brewer.pal(12, "Paired"), "99")[(as.integer(clusters[j]) - 1 + has_cluster_zero) %% 12 + 1],
                                            lty = 3,
                                            border = "#00000099"
                                        )
                                    }
                                ))
                            }
                        }
                    ))

                #Add mutation counts to the circles
                text(x_lab_pos, 1, round(colSums(exposures_plot)), cex=1.3)
                par(xpd=F)

        invisible(dev.off())
    }

# --------------------------------------------------
# Signature fit plots
# --------------------------------------------------

#Draw supersample signature fit barplot
    pdf(paste0(out_dir, "fit_barplot_", tolower(mut_class), "_supersamples.pdf"), w=10, h=9)
        for (patient in colnames(exps_supersample)) {
            #Select only signatures observed in patient
            patient_supersample_exps = exps_supersample[,patient, drop=F]
            patient_supersample_exps = patient_supersample_exps[rowSums(patient_supersample_exps) > 0,, drop=F]

            #Draw signature fit plot
            plotSigFitSample(
                catal_supersample[,patient, drop=F],
                ref_sigs[,rownames(patient_supersample_exps), drop=F],
                patient_supersample_exps,
                main_text = paste(
                    patient, "supersample -", mut_class, "profile -",
                    sum(data_exps_subclone$patient == patient & !data_exps_subclone$filtered & !data_exps_subclone$combination),
                    "cluster(s) -",
                    sum(catal_supersample[,patient, drop=F]), "mutations"
                )
            )
        }
    invisible(dev.off())

#Draw subclonal signature fit barplots by patient
    #Make sure target directory for the patient files exists
    out_dir_subclone_fit = paste0(out_dir, "subclone_fit_plots/")
    if (! dir.exists(out_dir_subclone_fit)) {
        dir.create(out_dir_subclone_fit, showWarnings=F)
    }

    #Loop over patients
    for (patient in colnames(exps_supersample)) {
        #Get subclones of the patient
        subclones_of_patient = rownames(data_exps_subclone)[data_exps_subclone$patient == patient]

        #Prepare PDF output with all samples of the patient
        pdf(paste0(out_dir_subclone_fit, "fit_barplot_", patient, "_", tolower(mut_class), "_subclones.pdf"), w=10, h=9)
        for (subclone in subclones_of_patient) {
            #Select only signatures observed in patient
            subclone_exps = exps_subclone[,subclone, drop=F]
            subclone_exps = subclone_exps[rowSums(subclone_exps) > 0,, drop=F]

            #Draw signature fit plot
            plotSigFitSample(
                catal_subclone[,subclone, drop=F],
                ref_sigs[,rownames(subclone_exps), drop=F],
                subclone_exps,
                main_text = paste0(
                    patient, " subclone ",
                    sub(paste0(patient, "_"), "", subclone),
                    " - ", mut_class, " profile - ",
                    sum(catal_subclone[,subclone, drop=F]), " mutations"
                )
            )
        }
        invisible(dev.off())
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


