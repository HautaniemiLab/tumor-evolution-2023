#!/usr/bin/env Rscript

#R script for plotting sample tree branch signature exposures and fit plots.
#Plots both supersample and single tree branches.
#Exposure plots include dendrograms

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

#Options include: input sample tree directory and suffix, supersample and sample
#tree branch input spectra table files, supersample and sample tree branch
#input signature contribution table files, mutation class (SBS/DBS/ID) and
#output directory
#Possible options for custom reference signatures, e.g. exome signatures

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--input-sample-tree-directory", "-in-st-dir")) {
            #input sample tree RData directory
            i = i + 1
            in_sample_trees_dir = args[i]
        } else if (args[i] %in% c("--input-sample-tree-suffix", "-in-st-suf")) {
            #input sample_tree RData suffix
            i = i + 1
            in_sample_trees_suf = args[i]
        } else if (args[i] %in% c("--supersample-catalogues", "-super-catal")) {
            #input supersample catalogues table
            i = i + 1
            in_file_catal_supersample = args[i]
        } else if (args[i] %in% c("--tree-branch-catalogues", "-branch-catal")) {
            #input sample tree branch catalogues table
            i = i + 1
            in_file_catal_single_branch = args[i]
        } else if (args[i] %in% c("--supersample-exposures", "-super-exp")) {
            #input supersample exposures table
            i = i + 1
            in_file_exps_supersample = args[i]
        } else if (args[i] %in% c("--tree-branch-exposures", "-branch-exp")) {
            #input sample tree branch exposures table
            i = i + 1
            in_file_exps_single_branch = args[i]
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
        !exists("in_sample_trees_dir") | !exists("in_sample_trees_suf") |
        !exists("in_file_catal_supersample") | !exists("in_file_catal_single_branch") |
        !exists("in_file_exps_supersample") | !exists("in_file_exps_single_branch") |
        !exists("out_dir") | !exists("mut_class")
    ) {
        stop(
            "Please supply input sample tree directory and suffix, supersample and sample tree branch input spectra table files, supersample and sample tree branch input signature contribution table files, mutation class (SBS/DBS/ID) and output directory",
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
    #Sample tree branches
    data_catal_single_branch = read.table(in_file_catal_single_branch, as.is=T, header=T, check.names=F)

    #Strip annotations from the spectra
    catal_supersample = t(data_catal_supersample[,rownames(ref_sigs), drop=F])
    catal_single_branch = t(data_catal_single_branch[,rownames(ref_sigs), drop=F])

#Read the exposures
    #Supersample
    data_exps_supersample = read.table(in_file_exps_supersample, as.is=T, header=T, row.names=1, check.names=F)
    #Sample tree branches
    data_exps_single_branch = read.table(in_file_exps_single_branch, as.is=T, header=T, row.names=1, check.names=F)

# --------------------------------------------------
# Prepare exposure plot data
# --------------------------------------------------

#Extract cosine similarities
    exps_cos_sim_all = c(data_exps_supersample$cos_sim, data_exps_single_branch$cos_sim)
    names(exps_cos_sim_all) = c(rownames(data_exps_supersample), rownames(data_exps_single_branch))

#Transform exposure tables back to proper format with signatures on rows
    exps_supersample = t(as.matrix(
        data_exps_supersample[,colnames(data_exps_supersample) %in% colnames(ref_sigs), drop=F]
    ))
    exps_single_branch = t(as.matrix(
        data_exps_single_branch[,colnames(data_exps_single_branch) %in% colnames(ref_sigs), drop=F]
    ))

#Remove suffix from supersamples
    names(exps_cos_sim_all) = gsub("_supersample", "", names(exps_cos_sim_all))
    colnames(exps_supersample) = gsub("_.*", "", colnames(exps_supersample))

#Set y-axis scale to thousands if there are samples with tens of thousands of mutations
    exp_plot.ylab.scale = ifelse(max(colSums(catal_single_branch)) >= 10000, 1000, 1)

#Prepare list of sample tree data
    in_sample_tree_files = grep(
        paste0(in_sample_trees_suf, "$"),
        paste0(in_sample_trees_dir, dir(in_sample_trees_dir)),
        value = T
    )
    names(in_sample_tree_files) = sub(".*\\/", "", sub(in_sample_trees_suf, "", in_sample_tree_files))

#Prepare signature colours
    #Observed signatures from supersample and single sample signatures
    #Those detected only in single samples use a different colour scheme
    supersample_signatures = rownames(exps_supersample)
    supersample_signatures = colnames(ref_sigs)[colnames(ref_sigs) %in% supersample_signatures]

    #Other signatures
    other_signatures = rownames(exps_single_branch)[
        !rownames(exps_single_branch) %in% supersample_signatures
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

#Draw supersample signature exposures plot
    pdf(paste0(out_dir, "contributions_", tolower(mut_class), "_supersample.pdf"), w=8, h=6)
    plotSigExposures(
        exps_supersample, col_sigs,
        title = paste("Full trees", mut_class, "Mutations Contributed", sep = " - "),
        ylab.scale=exp_plot.ylab.scale, lab.cex=0.7, lab.srt=60, lab.adj=c(1,.5),
        lab.x.shift=0.1, lab.y.shift=0.03, mar.4.shift=exp_plot.mar.4.shift
    )
    invisible(dev.off())

# --------------------------------------------------
# Functions for sample tree branch exposure plots
# --------------------------------------------------

#Helper function for creating branch/sample binary table
    treeBranchSampleBinTable = function(branches, samples, sample_order) {
        #Create matrix of zeroes
        branch_binary = matrix(
            rep(0, length(samples) * length(branches)),
            nrow = length(samples)
        )
        colnames(branch_binary) = 1:length(branches) - 3
        rownames(branch_binary) = samples

        #Update branch/sample pairs to 1
        for (i in 3:length(branches)) {
            branch_binary[branches[[i]], i] = 1
        }

        #Reorder rows
        if (nrow(branch_binary) >= 2) branch_binary = branch_binary[sample_order,]

        branch_binary
    }

#Helper function for compiling tree information as a list
    treeInfoList = function(branches, variants_branches, sample_order) {
        #Initialise list
        branch_info = list()

        for (i in 1:(length(branches)-3)) {
            branch_info[[i]] = list(
                "samples" = branches[[i+3]],
                "size" = sum(variants_branches$branch == i),
                "children" = NULL
            )
        }

        #Determine child nodes, if available
        if (length(sample_order) >= 2) {
            for (sample_id in names(sample_order)) {
                path = treeGetPath(branches, sample_id) - 3
                for (i in 1:(length(path)-1)) {
                    #Check if child is already included
                    if (!length(which(branch_info[[path[i]]]$children == path[i+1]))) {
                        branch_info[[path[i]]]$children = c(branch_info[[path[i]]]$children, path[i+1])
                    }
                }
            }
        }

        #Return value
        branch_info
    }

#Helper function for finding tallest branch via depth first recursive search
    treeFindTallest = function(branch_info, node, curr_path) {
        #Add new node to current path
        curr_path$nodes = c(curr_path$nodes, node)
        curr_path$height = curr_path$height + branch_info[[node]]$size
        tallest = curr_path

        #Recursive call
        if (length(branch_info[[node]]$children) > 0) {
            #Call function for each child node
            for (child in branch_info[[node]]$children) {
                new_path = treeFindTallest(branch_info, child, curr_path)
                #Update tallest
                if (new_path$height > tallest$height) tallest = new_path
            }
        }

        #Return value
        tallest
    }

#Helper function for drawing dendrogram via depth first recursion
    treeDraw = function(branch_info, branch_binary, node, curr_height) {
        #Add new node to current height
        new_height = curr_height + branch_info[[node]]$size
        x_pos = mean(
            which( rownames(branch_binary) %in% names(branch_info[[node]]$samples) )
        )

        #Draw vertical line
        segments(
            x_pos, curr_height,
            y1 = new_height
        )

        #Recursive call - there are either 0 or 2 children (bottom-up clustering)
        if (length(branch_info[[node]]$children) > 0) {
            #Call function for each child node
            x_children = sapply(
                branch_info[[node]]$children,
                function(child) {
                    #Draw horizontal line if child has same samples
                    if (setequal(branch_info[[node]]$samples, branch_info[[child]]$samples)) {
                        segments(
                            par("usr")[1], new_height, x_pos,
                            col = "#33333333"
                        )
                    }
                    treeDraw(branch_info, branch_binary, child, new_height)
                }
            )
            #Draw horizontal line
            segments(
                x_children[1], new_height,
                x_children[2]
            )
            #Draw branch number
            text(
                x_children, new_height,
                branch_info[[node]]$children,
                pos = c(2,4)[order(x_children)]
            )
        } else if (length(branch_info[[node]]$samples) == 1) {
            #Draw sample name for leafs
            text(
                which(rownames(branch_binary) == names(branch_info[[node]]$samples)),
                new_height,
                sample_names[rownames(branch_binary) == names(branch_info[[node]]$samples)],
                pos = 1, las = 2
            )
        }

        #Return value
        x_pos
    }

# --------------------------------------------------
# Sample tree branch exposure plots
# --------------------------------------------------

#Make sure target directory for the patient files exists
    out_dir_sample_exp = paste0(out_dir, "sample_tree_exposures/")
    if (! dir.exists(out_dir_sample_exp)) {
        dir.create(out_dir_sample_exp, showWarnings=F)
    }

#Loop over patients
    for (patient in colnames(exps_supersample)) {
        #Prepare sample tree data for plotting
            #Load sample tree data
            load(in_sample_tree_files[patient])

            #Create binary sample/branch table
            branch_binary = treeBranchSampleBinTable(branches, samples, sample_order)

            #Create branch info list
            branch_info = treeInfoList(branches, variants_branches, sample_order)

            #Determine height of tallest branch
            curr = list("nodes"=c(), "height"=0)
            tallest = treeFindTallest(branch_info, 1, curr)

            #Short sample names
            sample_names = sub(
                "_DNA\\d?", "",
                sub(
                    ".SAC$", "",
                    sub(
                        paste0(patient, "_"), "",
                        rownames(branch_binary)
                    )
                )
            )

        #Prepare exposures table for plotting
            #Select only signatures observed in patient
            patient_supersample_exps = exps_supersample[,patient, drop=F]
            patient_supersample_exps = patient_supersample_exps[rowSums(patient_supersample_exps) > 0,, drop=F]

            #Get tree branches of the patient
            branches_of_patient = rownames(data_exps_single_branch)[data_exps_single_branch$patient == patient]

            #Select branches of patient
            branch_exps = exps_single_branch[,branches_of_patient, drop=F]
            branch_exps = branch_exps[rowSums(branch_exps) > 0,, drop=F]

            #Exposures table for plotting with all branches incl. -2, -1 and 0 if missing (when they have 0 mutations)
            exposures_plot = branch_exps
            colnames(exposures_plot) = sub(paste0(patient, "_"), "", colnames(exposures_plot))
            for (branch in colnames(branch_binary)) {
                #Add zero column if missing
                if (! branch %in% colnames(exposures_plot)) {
                    exposures_plot = cbind(
                        exposures_plot,
                        matrix(0, nrow=nrow(branch_exps), dimnames=list(NULL, branch))
                    )
                }
            }

            #Add patient supersample exposures and sort the table
            exposures_plot = exposures_plot[,order(as.numeric(colnames(exposures_plot))), drop=F]
            exposures_plot = joinSigExposures(list(patient_supersample_exps, exposures_plot), ref_sigs)
            exposures_plot = exposures_plot[, c(2:4, 1, 5:ncol(exposures_plot)), drop=F]

            #Cosine similarities for each case
            cos_sim_plot = exps_cos_sim_all[c(patient, colnames(branch_exps))]
            names(cos_sim_plot) = sub(paste0(patient, "_"), "", names(cos_sim_plot))
            cos_sim_plot = cos_sim_plot[colnames(exposures_plot)]
            names(cos_sim_plot) = colnames(exposures_plot)

        #Prepare PDF output with all samples of the patient
        pdf(
            paste0(out_dir_sample_exp, "exposures_", patient, "_", tolower(mut_class), "_sample_tree_branches.pdf"),
            width = 0.66 * (ncol(branch_binary) + 1) + 0.11 + 1.24,
            height = nrow(exposures_plot) + 6.5
        )
            #Layout
            par(mfrow=c(2*nrow(exposures_plot) + 2*4 + 5, 1), mar=c(.5, 4, 3, 2) + .1)
            layout(
                matrix(c(
                    rep(1, 5),
                    sapply(
                        2:(nrow(exposures_plot) + 5),
                        function(x) rep(x, 2)
                    )
                ))
            )

            #Plot dendrogram
                #Empty plot
                plot(
                    0, type = "n",
                    xlim = c(1, nrow(branch_binary)) + .2 * c(-1,1),
                    ylim = c(1.1 * tallest$height, 0),
                    main = paste(patient, "sample tree branches -", mut_class, "signatures"),
                    xlab = "", ylab = "mutation count",
                    xaxt = "n", #yaxt = "n",
                    yaxs = "i",
                    las = 1, mgp = c(3, .75, 0),
                    bty = "n"
                )

                #Legend
                legend(
                    "topright",
                    c("-2", "-1", "0", patient),
                    text.width = 1.75 * 1.08 * (nrow(branch_binary) - .6) * 2.1/5 / (0.8 * ncol(branch_binary)),
                    adj = c(1, -1),
                    bty = "n",
                    xpd = T
                )
                legend(
                    "topright",
                    c("low depth", "low confidence", "ambiguous", "full tree"),
                    text.width = 1.7 * 1.08 * (nrow(branch_binary) - .6) * 2.1/5 / (0.8 * ncol(branch_binary)),
                    adj = c(0, -1),
                    bty = "n",
                    xpd = T
                )

                par(xpd=T)

                #Run the recursive function to drwan the dendrogram lines
                invisible(treeDraw(branch_info, branch_binary, 1, 0))

                #Branch text for root
                text(
                    mean(sample_order[names(branch_info[[1]]$samples)]), 0,
                    "1",
                    pos = 2
                )

                par(xpd=F)

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
                    #Mutation counts as text above bars
                    par(xpd=T)
                    text(
                        x_lab_pos, colSums(exposures_plot),
                        labels = round(colSums(exposures_plot)),
                        pos = 3, offset = .25
                    )
                    par(xpd=F)

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
                par(mar=c(2, 4, 0, 2) + .1)
                #Empty plot
                x_lab_pos = as.numeric(barplot(
                    exposures_plot,
                    xlim = c(0, 1.2 * ncol(exposures_plot) + .2),
                    ylim = c(0, nrow(branch_binary)) + .5,
                    ylab = "sample count",
                    yaxt = "n", yaxs = "i", mgp = c(3, 0, 0), xaxs = "i",
                    col = "white", border = F
                ))
                names(x_lab_pos) = colnames(exposures_plot)

                #Coloured rectangles based on number of samples in branch, and sample names
                col_branch = colorRampPalette(brewer.pal(6,"BuGn"))(nrow(branch_binary))
                for (i in 4:ncol(branch_binary)) {
                    #Coloured rectangles
                    rect(
                        x_lab_pos[i+1] - .6,
                        nrow(branch_binary) - min(which(branch_binary[,i] == 1)) + 1.5,
                        x_lab_pos[i+1] + .6,
                        nrow(branch_binary) - max(which(branch_binary[,i] == 1)) + .5,
                        col = col_branch[sum(branch_binary[,i])],
                        border = F
                    )

                    #Write sample names
                    for (j in 1:nrow(branch_binary)) {
                        if (branch_binary[j,i]) {
                            text(
                                x_lab_pos[i+1], nrow(branch_binary) - j + 1,
                                sample_names[j],
                                cex = .8
                            )
                        }
                    }
                }

                #Write "low depth", "low confidence" and "ambiguous" for branches -2, -1 and 0
                #and "full tree" for supersample
                text(
                    x_lab_pos[1:4],
                    mean(c(nrow(branch_binary), 1)),
                    c("low depth", "low confidence", "ambiguous", "full tree"),
                    cex = .8,
                    srt = 90
                )

                #Colour code on the left side
                    par(mar=c(2, 3.5, 0, 2 + par("mar")[2]*par("plt")[2] / par("plt")[1] - par("mar")[2]) + .1, new=T)
                    #Rectangles with barplot
                    barplot(
                        matrix(rep(1, nrow(branch_binary))),
                        xlim = c(0, 1) + .2, ylim = c(0, nrow(branch_binary)),
                        xlab = "", ylab = "",
                        yaxt = "n", yaxs = "i", xaxs = "i",
                        col = col_branch, border = F
                    )
                    #Axis
                    axis_skip = ceiling(nrow(branch_binary)/7)
                    axis_left = axis_skip * 1 : round(nrow(branch_binary)/axis_skip) - axis_skip + 1
                    if (axis_left[length(axis_left)] + axis_skip == nrow(branch_binary)) {
                        axis_left = c(axis_left, nrow(branch_binary))
                    }
                    axis(2, at=1:nrow(branch_binary) - .5, labels=F, las=1, mgp=c(3, .75, 0))
                    axis(2, at=axis_left - .5, axis_left, tick=F, las=1, mgp=c(3, .75, 0))

                    #Box
                    box()

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
                    patient, "full tree -", mut_class, "profile -",
                    sum(data_exps_single_branch$patient == patient), "branch(es) -",
                    sum(catal_supersample[,patient, drop=F]), "mutations"
                )
            )
        }
    invisible(dev.off())

#Draw sample tree branch signature fit barplots by patient
    #Make sure target directory for the patient files exists
    out_dir_sample_fit = paste0(out_dir, "tree_branch_fit_plots/")
    if (! dir.exists(out_dir_sample_fit)) {
        dir.create(out_dir_sample_fit, showWarnings=F)
    }

    #Loop over patients
    for (patient in colnames(exps_supersample)) {
        #Get tree branches of the patient
        branches_of_patient = rownames(data_exps_single_branch)[data_exps_single_branch$patient == patient]

        #Prepare PDF output with all samples of the patient
        pdf(paste0(out_dir_sample_fit, "fit_barplot_", patient, "_", tolower(mut_class), "_sample_tree_branches.pdf"), w=10, h=9)
        for (branch in branches_of_patient) {
            #Select only signatures observed in patient
            branch_exps = exps_single_branch[,branch, drop=F]
            branch_exps = branch_exps[rowSums(branch_exps) > 0,, drop=F]

            #Draw signature fit plot
            plotSigFitSample(
                catal_single_branch[,branch, drop=F],
                ref_sigs[,rownames(branch_exps), drop=F],
                branch_exps,
                main_text = paste0(
                    patient, " branch ",
                    sub(paste0(patient, "_"), "", branch),
                    " - ", mut_class, " profile - ",
                    sum(catal_single_branch[,branch, drop=F]), " mutations"
                )
            )
        }
        invisible(dev.off())
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


