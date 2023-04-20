#!/usr/bin/env Rscript

#R script for drawing sample tree VAF plots with optional variant annotation
#tracks. VAF plots can be further sorted by median VAF within branches.

#Multiple tracks can be specified.

# --------------------------------------------------
# Packages
# --------------------------------------------------

#RColorBrewer
library(RColorBrewer)

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include: sample tree RData and output plot file path.
#Possible options for
# - variant tracks tables (can be specified multiple times)
# - whether to sort heatmap also by VAF
    #Default options
    tracks_files = c()
    sort_by_vaf = F

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--input", "-i")) {
            #input sample tree RData file
            i = i + 1
            in_file_rdata = args[i]
        } else if (args[i] %in% c("--output", "-o")) {
            #output plot path
            i = i + 1
            out_file = args[i]
        } else if (args[i] %in% c("--tracks", "-t")) {
            #input variant tracks tables
            i = i + 1
            tracks_files = c(tracks_files, args[i])
        } else if (args[i] %in% c("--sort-by-vafs", "-s")) {
            #sorting by VAF
            sort_by_vaf = T
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (!exists("in_file_rdata") | !exists("out_file")) {
        stop(
            "Please supply input sample tree RData and output plot path",
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

#Function for breaking down tracks to separate binary tracks
    collapseTracks = function(tbl) {
        #Ignore columns "CHROM", "POS", "REF", "ALT"
        tbl = tbl[,! colnames(tbl) %in% c("CHROM", "POS", "REF", "ALT"), drop=F]

        #Collapse each column
        col_tracks = lapply(
            colnames(tbl),
            function(column) {
                #Get unique values of column
                unique_vals = sort(unique(tbl[,column]))

                #Create binary track for each unique value
                sapply(
                    as.character(unique_vals[!is.na(unique_vals)]),
                    function(val) {
                        col_binary = as.numeric(tbl[,column] == val)
                        col_binary[col_binary == 0] = NA

                        col_binary
                    }
                )
            }
        )

        do.call(cbind, col_tracks)
    }

#Function for reading variant tracks table
    readTracks = function(file_path) {
        tbl = read.table(file_path, as.is=T, header=T, check.names=F)

        cbind(
            joinKey(tbl),
            collapseTracks(tbl)
        )
    }

#Function for joining variant tracks tables by join_key
    joinTracks = function(tbl1, tbl2) {
        merge(tbl1, tbl2, by = "join_key", all = T)
    }

# --------------------------------------------------
# Read in data
# --------------------------------------------------

#Read variants
    load(in_file_rdata)

#Read tracks
    tracks_tables = lapply(
        tracks_files,
        function(x) readTracks(x)
    )

# --------------------------------------------------
# Preparing VAF matrix
# --------------------------------------------------

#Extract patient ID
    patient = sub("(_|-|\\.).*", "", grep("\\.SAC", colnames(variants), value=T)[1])

#Compute VAF matrix
    vafs = matrix(
        apply(
            variants[, endsWith(colnames(variants), ".SAC"), drop=F],
            1,
            function(sacs) {
                #Get read counts
                rcs = sacsToADs(sacs)
                #Set REF read count to minimum of 1 to avoid zero division
                rcs[1,] = rcs[1,] + (colSums(rcs) == 0)
                #Count VAF, and if read count is under threshold of mutation type, returns zero
                any(rcs[2,] >= min_rc_patient) * (rcs[2,] >= min_rc_sample) * rcs[2,] / colSums(rcs)
            }
        ),
        ncol = sum(endsWith(colnames(variants), ".SAC")),
        dimnames = list(NULL, grep("\\.SAC", colnames(variants), value=T)),
        byrow = T
    )

    #Apply sample tree variant order
    vafs = vafs[vars_order, sample_order, drop=F]

# --------------------------------------------------
# Prepare variant track
# --------------------------------------------------

#Perform if tracks tables have been given
    if (length(tracks_files) > 0) {
        #Join any variant tracks tables
        tracks = tracks_tables[[1]]
        i = 2
        while (i <= length(tracks_files)) {
            tracks = joinTracks(tracks, tracks_tables[[i]])
            i = i + 1
        }

        #Join variant tracks with variants
        vars_tracks = merge(joinKey(variants), tracks, by="join_key", all.x=T)

        #Sort tracks to the same order as variants and drop columns
        rownames(vars_tracks) = vars_tracks$join_key
        vars_tracks = vars_tracks[
            unlist(joinKey(variants)),
            colnames(vars_tracks) != "join_key",
            drop = F
        ]
        rownames(vars_tracks) = NULL

        #Apply sample tree variant order
        vars_tracks = vars_tracks[vars_order,, drop=F]
    }

# --------------------------------------------------
# Sorting by VAFs
# --------------------------------------------------

#Branch sizes and starting positions for vertical border positions
    branch_sizes = rep(0, length(branches))
    names(branch_sizes) = 1:length(branches) - 3
    branch_sizes[names(table(variants_branches$branch))] = table(variants_branches$branch)
    branch_borders = cumsum(branch_sizes)

#Reordering branches by VAFs
    if (sort_by_vaf) {
        #Initialise order
        vars_order_by_vaf = 1:nrow(vafs)

        #Loop over main branches
        for (i in 4:length(branch_borders)) {
            #Don't sort if the branch has no variants
            if (branch_borders[i] > branch_borders[i-1]) {
                index_range = (branch_borders[i-1] + 1):branch_borders[i]

                #Reorder branch variants by median VAF
                vars_order_by_vaf[index_range] = vars_order_by_vaf[index_range][
                    order(
                        apply(
                            vafs[index_range, names(branches[[i]]), drop=F],
                            1,
                            median, na.rm=T
                        ),
                        decreasing = T
                    )
                ]
            }
        }

        #Apply new order on vafs and tracks
        vafs = vafs[vars_order_by_vaf,, drop=F]
        if (length(tracks_files) > 0) {
            vars_tracks = vars_tracks[vars_order_by_vaf,, drop=F]
        }
    }

# --------------------------------------------------
# Prepare VAF heatmap data
# --------------------------------------------------

#VAF matrix as x and y vectors for plotting
    vafs_x = rep(1:nrow(vafs), length(samples))
    vafs_y = c(sapply(
        length(samples):1,
        function(x) {
            rep(x, nrow(vafs))
        }
    ))

#VAF heatmap color palette
    vaf_breaks = c(0, 3, 6, 10, 15, 20, 25, 30, 35, 45, 100)/100
    vaf_col_levels = c(
        "#FBFCBE", "#FDD294", "#FDA873", "#F87B5D", "#E95461",
        "#C73D73", "#A22F7D", "#7C2382", "#5A157D", "#331068"
    )

    #Colour values for each data point
    vafs_cols = vaf_col_levels[as.numeric(cut(unlist(vafs), breaks=vaf_breaks))]

#Remove variants with white colour so that they are not drawn
    #vafs_x = vafs_x[-which(vafs_cols == "#FBFCBE")]
    #vafs_y = vafs_y[-which(vafs_cols == "#FBFCBE")]
    #vafs_cols = vafs_cols[-which(vafs_cols == "#FBFCBE")]

# --------------------------------------------------
# Drawing VAF heatmap
# --------------------------------------------------

#PDF file and plot size depending on whether tracks were specified
    if (length(tracks_files) > 0) {
        pdf(out_file, height=.25*max(ncol(vafs),3) + .125*(1 + ncol(vars_tracks)) + 1.84, width=12)
        ylim = c(-.5*(1 + ncol(vars_tracks)), length(samples)) + .5
    } else {
        pdf(out_file, height=.25*max(ncol(vafs),3) + 1.84, width=12)
        ylim = c(0, length(samples)) + .5
    }
        par(new=F, xpd=F, mar=c(5, 4, 4, 12.5) + .1)

        #Main plot
            #Empty plot
            plot(0, type = "n",
                xlim = c(0, nrow(vafs)), ylim = ylim,
                main = paste0(patient, " VAF heatmap - variants ordered by tree branches"),
                xlab = "variant, respective tree branch",
                ylab = "",
                xaxs = "i", yaxs = "i",
                xaxt = "n", yaxt = "n"
            )

            #Background colour
            rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey85")

            #x and y axes; sample names without .SAC suffix
            axis(
                1,
                at = .5 * (branch_borders + c(0, branch_borders[-length(branch_borders)]))[!duplicated(branch_borders)],
                labels = (-2:(length(branches)-3))[!duplicated(branch_borders)]
            )
            axis(
                4,
                at = ncol(vafs):1,
                labels = gsub("\\.SAC$", "", colnames(vafs)),
                las = 1
            )

            #Plot the actual data as rectangles
            rect(
                vafs_x - 1, vafs_y - .5,
                vafs_x, vafs_y + .5,
                col = vafs_cols,
                border = NA
            )

        #Tracks
            if (length(tracks_files) > 0) {
                #Plot each track as lines on their own row
                for (i in 1:ncol(vars_tracks)) {
                    #Background gray color
                    rect(
                        0, -.5 * (i - 1),
                        nrow(vars_tracks), -.5 * i,
                        col = c("#BBBBBB", "#777777")[i%%2 + 1],
                        border = NA
                    )
                    #Rectangles
                    if (sum(!is.na(vars_tracks[,i])) > 0) {
                        rect(
                            which(vars_tracks[,i] == 1) - 1,
                            -.5 * (i - 1),
                            which(vars_tracks[,i] == 1),
                            -.5 * i,
                            col = paste0(brewer.pal(12, "Paired"), "66")[(i-1) %% 12 + 1],
                            border = NA
                        )
                    }
                }

                #Track axis labels
                axis(
                    4,
                    at = -.5 * 1:ncol(vars_tracks) + .25,
                    labels = colnames(vars_tracks),
                    las = 1, cex.axis = .75
                )
            }

        #Vertical branch borders
            abline(
                v = unique(branch_borders[-length(branch_borders)]),
                lty=3, lwd=1, col="#33333399"
            )

        #Dendrogram (not representative of distances used in clustering)
            if (length(samples) >= 2) {
                #Function for getting parent branch index
                getParent = function(curr_index) {
                    for (i in (curr_index-1):3) {
                        if (all(!is.na(branches[[i]][names(branches[[curr_index]])]))) {
                            return(i)
                        }
                    }
                }

                #Initialise verticals (contains y positions of branching points)
                verticals = list()
                verticals[[length(branches)]] = 1

                #Loop down branch indices to draw dendrogram
                for (i in length(branches):5) {
                    #Get y position
                    y_pos = ncol(vafs) + 1 - mean(order(sample_order)[branches[[i]]])
                    #Get parent index
                    parent = getParent(i)

                    #Horizontal line
                    segments(
                        branch_borders[parent], y_pos,
                        branch_borders[i],
                        lwd = 1.5
                    )

                    #Add y position to parent list
                    verticals[[parent]] = c(verticals[[parent]], y_pos)

                    #Draw vertical line if there are 2 y positions
                    if (length(verticals[[parent]]) == 2) {
                        segments(
                            branch_borders[parent],
                            min(verticals[[parent]]),
                            y1 = max(verticals[[parent]]),
                            lwd = 1.5
                        )
                    }
                }
            }

            #Root - Horizontal line
            segments(
                branch_borders[3], mean(branches[[4]]),
                branch_borders[4],
                lwd = 1.5
            )

            #Black borders around the plot
            box()

        #Heatmap colour guide
            par(new=T, xpd=T, mar=c(5, 3, 4, 5*par("pin")[1] + 12.5) + .1)

            #Empty plot
            plot(
                0,
                type = "n",
                xlim = 0:1, ylim = 0:1,
                xlab = "", ylab = "",
                xaxs = "i", yaxs = "i",
                xaxt = "n", las = 1
            )

            #Colour code (rectangles on the side)
            rect(
                0, vaf_breaks[-length(vaf_breaks)],
                1, vaf_breaks[-1],
                col = vaf_col_levels,
                border = F
            )

            #Black borders around the plot
            box()

    invisible(dev.off())

# --------------------------------------------------
# End of file
# --------------------------------------------------


