#Collection of functions related to signature plotting.
#List of contents:
# - Reference/Single signature barplots:
#   - plotSignaturesSBS(values, ...)
#   - plotSignaturesDBS(values, ...)
#   - plotSignaturesID(values, ...)
# - Signature fit barplots:
#   - plotSigFitSample(catalogue, signatures, exposures, ...)
#   - plotSigFitSampleLarge(catalogue, signatures, exposures, ...)
#   - plotSigExposures(exposures, col_exp, ...)
# - Error plots:
#   - plotSigFitSelErrors(x, y1, ...)

# --------------------------------------------------
# Packages
# --------------------------------------------------

#RColorBrewer
suppressMessages(require(RColorBrewer))

# --------------------------------------------------
# Helper functions/definitions
# --------------------------------------------------

#Reverse complement function
    revComp = function(x) {
        chartr(
            "ACGT",
            "TGCA",
            paste0(rev(unlist(strsplit(x, ""))), collapse="")
        )
    }

#DBS preferences
    dbsPrefs = list(
        ref_5 = c("T", "C", "A", "G"),
        alt_5 = list(
            A = c("C", "G", "T"),
            C = c("T", "G", "A"),
            G = c("A", "C", "T"),
            T = c("G", "C", "A")
        )
    )

# --------------------------------------------------
# Mutation type - colour definitions
# --------------------------------------------------

#Bases
    bases = c("A", "C", "G", "T")

#Base substitution 6-type colours
    sbs_types_6 = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    cols_sbs_6 = c("cyan", "black", "red", "gray", "lightgreen", "pink")
    names(cols_sbs_6) = sbs_types_6

#Trinucleotide context substitution with 6 colours
    sbs_types_96 = apply(
        merge(
            bases,
            apply(
                merge(bases, sub(">.", "", sbs_types_6)),
                1,
                paste0, collapse=""
            )
        ),
        1,
        function(x) paste0(x[2], x[1])
    )
    cols_sbs_96 = c(sapply(
        cols_sbs_6,
        function(x) {
            rep(x, 16)
        }
    ))
    names(cols_sbs_96) = sbs_types_96

#DBS signature 78-type colours
    #Create all possible DBS
    dbs_types_78 = c(sapply(
        #All dinucleotides
        apply(merge(bases, bases), 1, paste0, collapse=""),
        function(dinuc) {
            #For each dinucleotide combine with all 3*3 possible ALT dinucleotides
            paste(
                dinuc,
                apply(
                    do.call(
                        merge,
                        lapply(
                            unlist(strsplit(dinuc, "")),
                            function(base) bases[!bases %in% base]
                        )
                    ),
                    1,
                    paste0, collapse=""
                ),
                sep = ">"
            )
        }
    ))
    #Take preferred reverse complements of DBS types
    dbs_types_78 = sort(unique(sapply(
        dbs_types_78,
        function(dbs_type) {
            #Get REF ALT dinucleotides
            dbs = unlist(strsplit(dbs_type, ">"))
            #Get REF ALT +/- 5' bases
            ref_5 = substr(dbs_type, 1, 1)
            alt_5 = substr(dbs_type, 4, 4)
            ref_rev_5 = revComp(substr(dbs_type, 2, 2))
            alt_rev_5 = revComp(substr(dbs_type, 5, 5))
            #Output DBS type based on REF 5' base preference -> ALT 5' base preference
            if (which(dbsPrefs$ref_5 == ref_rev_5) < which(dbsPrefs$ref_5 == ref_5)) {
                paste0(sapply(dbs, revComp), collapse=">")
            } else if (ref_rev_5 != ref_5) {
                dbs_type
            } else if (which(dbsPrefs$alt_5[[ref_5]] == alt_rev_5) < which(dbsPrefs$alt_5[[ref_5]] == alt_5)) {
                paste0(sapply(dbs, revComp), collapse=">")
            } else {
                dbs_type
            }
        }
    )))
    names(dbs_types_78) = NULL
    #Colour
    cols_dbs_78 = rep(
        brewer.pal(10, "Paired"),
        table(sapply(dbs_types_78, function(x) unlist(strsplit(x, ">"))[1]))
    )
    names(cols_dbs_78) = dbs_types_78

#ID signature 83-type colours
    #ID type definitions
    id_types_83 = c(
        sort(c(
            apply(merge(c("DEL_C_1", "DEL_T_1"), c(1:5, "6+")), 1, paste0, collapse="_"),
            apply(merge(c("INS_C_1", "INS_T_1"), c(0:4, "5+")), 1, paste0, collapse="_")
        )),
        sort(c(
            apply(merge(paste0("DEL_repeats_", c(2:4, "5+")), c(1:5, "6+")), 1, paste0, collapse="_"),
            apply(merge(paste0("INS_repeats_", c(2:4, "5+")), c(0:4, "5+")), 1, paste0, collapse="_")
        )),
        paste0(
            "DEL_MH_",
            c(
                unlist(sapply(2:4, function(x) paste(x, 1:(x-1), sep="_")), paste0("5+_", c(1:4, "5+"))),
                paste0("5+_", c(1:4, "5+"))
            )
        )
    )
    #Colour
    cols_id_83 = rep(
        c(
            brewer.pal(8, "Paired")[c(7:8,3:4)],
            colorRampPalette(brewer.pal(9, "Reds")[3:7])(4),
            brewer.pal(11, "RdBu")[7:10],
            brewer.pal(11, "PRGn")[5:2]
        ),
        rle(gsub("_[0-9]\\+?$", "", id_types_83))$length
    )
    names(cols_id_83) = id_types_83
    #ID type labels in signature plots
    id_types_83_plot_labs = t(rbind(
        sapply(
            gsub("(_[CT])?_[0-9]\\+?_[0-9]\\+?$", "", id_types_83),
            function(x) {
                if (x == "DEL") {
                    c(paste("1bp deletion"), NA, "Homopolymer length")
                } else if (x == "INS") {
                    c(paste("1bp insertion"), NA, "Homopolymer length")
                } else if (x == "DEL_repeats") {
                    c(paste(">1bp deletions at repeats"), "(Deletion length)", "Number of repeat units")
                } else if (x == "INS_repeats") {
                    c(paste(">1bp insertions at repeats"), "(Insertion length)", "Number of repeat units")
                } else if (x == "DEL_MH") {
                    c(paste("Deletions with microhomology"), "(Deletion length)", "Microhomology length")
                }
            }
        ),
        gsub("(_1)?_[0-9]\\+?$", "", gsub("^[^_]*_((repeats|MH)_)?", "", id_types_83)),
        gsub("^.*_", "", id_types_83)
    ))
    colnames(id_types_83_plot_labs) = c("toplab", "sublab", "xlabs", "colcode", "barlab")

# --------------------------------------------------
# Reference/Single signature barplots
# --------------------------------------------------

#Function called by each signature plot for gridded barplot with title
    plotSignaturesBase = function(values, y_lim, y_lab, gap, borders, cols, main_text, simple=F, thin_plot=simple, line=2) {
        #Number of bars
        n = length(values)

        #Empty plot with grid lines
        barplot(
            rep(0, n),
            xlim = c(0, (n+1)*gap + n),
            ylim = y_lim,
            border = F,
            las = 2,
            xaxs = "i",
            yaxt = "n"
        )
        if (!simple) {
            axis(
                2,
                axTicks(2),
                if (isTRUE(all.equal(sum(values), 1))) paste0(100*axTicks(2), "%") else axTicks(2),
                las = 2,
                mgp = c(3, .75, 0)
            )
            abline(h=axTicks(2), col="#DDDDDD")
        }
        par(new=T)

        #Barplot
        x_pos = barplot(
            as.numeric(values),
            space = gap,
            xlim = c(0, (n+1)*gap + n),
            ylim = y_lim,
            axisnames = F,
            ylab = y_lab,
            cex.names = 0.5,
            col = cols,
            border = borders,
            mgp = c(3, .5, 0),
            las = 2,
            yaxt = "n",
            xaxs = "i",
            xpd = F
        )
        if (!simple) {
            box()
        } else {
            segments(par("usr")[1], par("usr")[3], x1=par("usr")[2], xpd=T)
            segments(par("usr")[1], par("usr")[3], y1=par("usr")[4], xpd=T)
        }

        #Title/signature name
        if (!thin_plot) {
            title(main_text, line=line, font.main=2, cex.main=1.3)
        } else {
            #Signature name in topleft
            text(
                par('usr')[1] + .5 * strwidth("S"),
                par('usr')[4] - .5 * strheight("S"),
                main_text,
                adj = c(0,1), font=2, cex=1.3
            )
        }

        #Return bar positions
        x_pos
    }

#Function for plotting 96-type SBS signatures
    plotSignaturesSBS = function(values, main_text="", y_lab="", y_range=NULL, labels=sbs_types_96, gap=0.5, borders=NA, simple=F, thin_plot=simple, xlab_new=F, line=2, no_top=simple, no_bottom=simple) {
        #Original plot settings
        orig_par = par(c("mar", "las", "cex.axis", "mgp", "font.axis", "family", "xpd"))
        #Thin plot settings
        if (thin_plot) {
            par(mar=pmax(par("mar") - c(0, 1, 1, 0), rep(0,4)))
            y_lab = ""
        }
        #Margins for no top and no bottom labels
        if (no_top) {
            par(mar=pmax(par("mar") - c(0, 0, 2.5, 0), rep(0,4)))
            line = line - 1
        }
        if (no_bottom) {
            par(mar=pmax(par("mar") - c(2, 0, 0, 0), rep(0,4)))
        }

        #Set ylim
        if (is.null(y_range)) {
            y_max = max(.1, 1.2 * as.numeric(values))
            y_lim = c(0, y_max)
        } else {
            y_lim = y_range
        }

        #Plot bars and y axis
        x_pos = plotSignaturesBase(values, y_lim, y_lab, gap, borders, cols_sbs_96, main_text, simple, thin_plot, line=line)

        #Draw x axis and axis labels
        if (!no_bottom) {
            if (xlab_new) {
                par(las=1, cex.axis=.5, mgp=c(1,0,0))
                #short ticks separating 3' bases
                axis(1, x_pos[-(4 * 1:24)] + .75, F, col=NA, col.ticks=1, lwd.ticks=1, tck=-.025)
                #long ticks separating 5' bases
                axis(1, c(.25, x_pos[4 * 1:24] + .75), F, col=NA, col.ticks=1, lwd.ticks=1, tck=-.06)
                #label for 3' base
                axis(1, x_pos[2*1:48], substr(labels, 3, 3)[2*1:48], tick=F, padj=-1.25)
                axis(1, x_pos[2*1:48-1], substr(labels, 3, 3)[2*1:48-1], tick=F, padj=-1.25)
                #label for 5' base
                axis(1, x_pos[4 * 1:24 - 2] + .75, rep(c("A","C","G","T"), 6), tick=F, padj=.75, cex.axis=1)
            } else {
                #Axis labels with colours based on mutating base
                par(las=2, cex.axis=.75, font.axis=2, family="mono")
                #5' base
                axis(1, x_pos, labels=paste0(substr(labels, 1, 1), "  "), tick=F, mgp=c(1,.5,0), col.axis="#888888")
                #Coloured mutating pyrimidine
                for (i in 1:6) {
                    axis(
                        1,
                        x_pos[1:16 + (i-1)*16],
                        labels = paste0(" ", substr(labels, 2, 2), " ")[1:16 + (i-1)*16],
                        tick = F,
                        mgp = c(1,.4,0),
                        col.axis = cols_sbs_6[i]
                    )
                }
                #3' base
                axis(1, x_pos, labels=paste0("  ", substr(labels, 3, 3)), tick=F, mgp=c(1,.3,0), col.axis="#888888")
            }
            par(las=0, cex.axis=1, mgp=c(3,1,0), font.axis=1, family="")
        }

        #Top guide
        if (!no_top) {
            par(xpd=T)
            #Rectangles
            rect(
                x_pos[16 * 0:5 + 1] - .5,
                y_lim[1] + 1.01 * sum(abs(y_lim)),
                x_pos[16 * 1:6] + .5,
                y_lim[1] + 1.06 * sum(abs(y_lim)),
                border = F,
                col = cols_sbs_6
            )
            #Text
            text(
                .5 * (x_pos[16 * 0:5 + 1] + x_pos[16 * 1:6]),
                y_lim[1] + 1.07 * sum(abs(y_lim)),
                c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                pos = 3, offset = 0
            )
        }

        #Return original par values
        par(orig_par)
    }

#Function for plotting 78-type DBS signatures
    plotSignaturesDBS = function(values, main_text="", y_lab="", y_range=NULL, labels=dbs_types_78, gap=0.5, borders=NA, simple=F, thin_plot=simple, line=2, no_top=simple, no_bottom=simple) {
        #Original plot settings
        orig_par = par(c("mar", "las", "cex.axis", "mgp", "font.axis", "family", "xpd"))
        #Thin plot settings
        if (thin_plot) {
            par(mar=pmax(par("mar") - c(0, 1, 1, 0), rep(0,4)))
            y_lab = ""
        }
        #Margins for no top and no bottom labels
        if (no_top) {
            par(mar=pmax(par("mar") - c(0, 0, 2.5, 0), rep(0,4)))
            line = line - 1
        }
        if (no_bottom) {
            par(mar=pmax(par("mar") - c(1.4, 0, 0, 0), rep(0,4)))
        }

        #Set ylim
        if (is.null(y_range)) {
            y_max = max(.1, 1.2 * as.numeric(values))
            y_lim = c(0, y_max)
        } else {
            y_lim = y_range
        }

        #Plot bars and y axis
        x_pos = plotSignaturesBase(values, y_lim, y_lab, gap, borders, cols_dbs_78, main_text, simple, thin_plot, line=line)

        #Draw x axis and axis labels
        if (!no_bottom) {
            #Axis labels with colours based on mutating base
            par(las=2, cex.axis=.75, font.axis=2, family="mono")
            #Dinucleotide
            axis(1, x_pos, labels=substr(labels, 4, 5), tick=F, mgp=c(1,.3,0), col.axis="#888888")
            par(las=0, cex.axis=1, mgp=c(3,1,0), font.axis=1, family="")
        }

        #Top guide
        if (!no_top) {
            par(xpd=T)
            #Cumsum table
            end_pos = c(0, cumsum(table(substr(labels,1,2))))
            #Rectangles
            rect(
                x_pos[end_pos[-length(end_pos)] + 1] - .5,
                y_lim[1] + 1.01 * sum(abs(y_lim)),
                x_pos[end_pos[-1]] + .5,
                y_lim[1] + 1.06 * sum(abs(y_lim)),
                border = F,
                col = unique(cols_dbs_78)
            )
            #Text
            text(
                .5 * (x_pos[end_pos[-length(end_pos)] + 1] + x_pos[end_pos[-1]]),
                y_lim[1] + 1.07 * sum(abs(y_lim)),
                paste0(unique(substr(labels,1,3)), "NN"),
                pos = 3, offset = 0
            )
        }

        #Return original par values
        par(orig_par)
    }

#Function for plotting 83-type ID signatures
    plotSignaturesID = function(values, main_text="", y_lab="", y_range=NULL, labels=id_types_83_plot_labs, gap=0.75, borders=NA, simple=F, thin_plot=simple, line=3, no_top=simple, no_bottom=simple) {
        #Original plot settings
        orig_par = par(c("mar", "las", "cex.axis", "mgp", "font.axis", "family", "xpd"))
        #Thin plot settings
        if (thin_plot) {
            par(mar=pmax(par("mar") - c(0, 1, 1, 0), rep(0,4)))
            y_lab = ""
        }
        #Margins for no top and no bottom labels
        if (no_top) {
            par(mar=pmax(par("mar") - c(0, 0, 3.8, 0), rep(0,4)))
            line = line - 1
        }
        if (no_bottom) {
            par(mar=pmax(par("mar") - c(4, 0, 0, 0), rep(0,4)))
        }

        #Set ylim
        if (is.null(y_range)) {
            y_max = max(.1, 1.2 * as.numeric(values))
            y_lim = c(0, y_max)
        } else {
            y_lim = y_range
        }

        #Plot bars and y axis
        x_pos = plotSignaturesBase(values, y_lim, y_lab, gap, borders, cols_id_83, main_text, simple, thin_plot, line=line)

        #Color guides
        par(xpd=T)
        #Cumsum table
        end_pos = c(0, cumsum(rle(labels[,"colcode"])$lengths))
        #Rectangles, top and bottom
        colouredRectangels = function(i) {
            rect(
                x_pos[end_pos[-length(end_pos)] + 1] - .5,
                y_lim[i] + c(1,-1)[i] * 1.01 * sum(abs(y_lim)),
                x_pos[end_pos[-1]] + .5,
                y_lim[i] + c(1,-1)[i] * 1.09 * sum(abs(y_lim)),
                border = F,
                col = unique(cols_id_83)
            )
        }
        if (!no_top) colouredRectangels(1)
        if (!no_bottom) colouredRectangels(2)
        #Text
        if (!no_top) {
            text(
                .5 * (x_pos[end_pos[-length(end_pos)] + 1] + x_pos[end_pos[-1]]),
                y_lim[1] + 1.05 * sum(abs(y_lim)),
                rle(labels[,"colcode"])$values,
                col = c("black", "white")[c(rep(0:1,2), rep(c(0,0,0,1),3)) + 1]
            )
        }

        #Top labels + sublabels + xlabs
        #Cumsum table
        end_pos = c(0, cumsum(rle(labels[,"toplab"])$lengths))
        has_sublab = !is.na(labels[end_pos[-length(end_pos)]+1, "sublab"])
        if (!no_top) {
            #Toplab
            text(
                .5 * (x_pos[end_pos[-length(end_pos)] + 1] + x_pos[end_pos[-1]]),
                y_lim[1] + 1.1 * sum(abs(y_lim)) + (2.1 - .75 * !has_sublab) * strheight("(DIlg)"),
                rle(labels[,"toplab"])$values,
                pos = 3, offset = 0,
                cex = 1.1
            )
            #Sublab
            text(
                .5 * (x_pos[end_pos[-length(end_pos)] + 1] + x_pos[end_pos[-1]]),
                y_lim[1] + 1.1 * sum(abs(y_lim)) + 0.3 * strheight("(DI)"),
                labels[end_pos[-length(end_pos)]+1, "sublab"],
                pos = 3, offset = 0,
                cex = 1.1
            )
        }
        if (!no_bottom) {
            #xlab
            text(
                .5 * (x_pos[end_pos[-length(end_pos)] + 1] + x_pos[end_pos[-1]]),
                y_lim[2] - 1.115 * sum(abs(y_lim)) - 2 * strheight("123456+", family="mono", cex=.9),
                labels[end_pos[-length(end_pos)]+1, "xlabs"],
                pos = 1, offset = 0
            )
        }

        #x axis labels
        if (!no_bottom) {
            text(
                x_pos,
                y_lim[2] - 1.115 * sum(abs(y_lim)),
                labels[,"barlab"],
                pos = 1, offset = 0,
                cex = .9,
                family = "mono"
            )
        }

        #Return original par values
        par(orig_par)
    }

# --------------------------------------------------
# Signature fit barplots
# --------------------------------------------------

#Single sample signature fit plot; 3 barplots for original catalogue, fitted catalogue
#and differences. Use type = "SBS", "DBS" or "ID" to select plot type or let the
#function decide based on mutation types
    plotSigFitSample = function(catalogue, signatures, exposures, main_text="", y_lab="fraction", type="", err_text="Cosine similarity", err_fun=cosSim, prob=T) {
        #Determine signature plotting function
        if (type == "SBS" || nrow(catalogue) == 96 && all(gsub("\\[|\\]|>[ACGT]", "", rownames(catalogue)) == sbs_types_96)) {
            plotSigsFun = plotSignaturesSBS
            par_mar_shift = c(0,0,0,0)
            line = 2.5
        } else if (type == "DBS" || nrow(catalogue) == 78 && all(rownames(catalogue) == dbs_types_78)) {
            plotSigsFun = plotSignaturesDBS
            par_mar_shift = c(-.6,0,0,0)
            line = 2.5
        } else if (type == "ID" || nrow(catalogue) == 83 && all(rownames(catalogue) == id_types_83)) {
            plotSigsFun = plotSignaturesID
            par_mar_shift = c(1.2,0,1.3,0)
            line = 4
        }

        #Normalisation of prob=T
        if (prob) {
            catal = prop.table(catalogue)
            exps = exposures / sum(exposures)
        } else {
            catal = catalogue
        }

        #Set mfrow and margin
        orig_par = par(c("mfrow", "mar"))
        par(mfrow=c(3,1), mar=c(2,4,5,2) + par_mar_shift + .1)

        #Get relevant information for plot mains
            #Middle main text (signatures with exposures)
            sigs_to_mention = exps[which(exps > 0),,drop=F]

            mid_text = paste(
                sapply(
                    1:length(sigs_to_mention),
                    function(i) {
                        output = paste(
                            rownames(sigs_to_mention)[i],
                            sprintf("%.3f", sigs_to_mention[i]),
                            sep = " : "
                        )
                        sub("nature", "", output)
                    }
                ),
                collapse = " - "
            )

        #Fitted value
        fit_catal = signatures %*% exps

        #Compute ylim
        y_max = max(catal, fit_catal)

        #Top plot, observed frequencies
        plotSigsFun(catal, main_text, y_lab, 1.2*c(0, y_max), line=line)

        #Mid plot, expected frequencies
        par(mar=c(2.5, 4, 4.5, 2) + par_mar_shift + .1)
        plotSigsFun(fit_catal, mid_text, y_lab, 1.2*c(0, y_max), line=line)

        #Bottom plot, difference
        par(mar=c(3, 4, 4, 2) + par_mar_shift + .1)
        plotSigsFun(
            catal - fit_catal,
            paste0(err_text, " = ", sprintf("%.3f", err_fun(catalogue, sum(catalogue) * fit_catal))),
            "difference",
            1.1 * c(-1, 1) * max(0.03, max(abs(catal - fit_catal))),
            line=line
        )
        abline(h=0)

        #Return original par values
        par(orig_par)
    }

#Large single sample signature fit plot; triple barplots for original catalogue,
#fitted catalogue and differences with signature contributions visualised in the
#left-hand side. Use type = "SBS", "DBS" or "ID" to select plot type or let the
#function decide based on mutation types. Option to also show the signatures
#with contribution to the sample
#Top and bottom visuals of each barplot are shown for only a subset of the plots.
#Use a default value to reserve slots for reference signatures to keep plots consistent in size.
    plotSigFitSampleLarge = function(catalogue, signatures, exposures, col_exp, main_text="", type="", err_text="Cosine similarity", err_fun=cosSim, prob=T, min_exp_rows=0, draw_signatures=F) {
        #Determine signature plotting function
        if (type == "SBS" || nrow(catalogue) == 96 && all(gsub("\\[|\\]|>[ACGT]", "", rownames(catalogue)) == sbs_types_96)) {
            plotSigsFun = plotSignaturesSBS
            par_mar_shift = c(2, 0, 2.5, 0) - c(0, 0, 1, 0) * draw_signatures
            line = 2.5
        } else if (type == "DBS" || nrow(catalogue) == 78 && all(rownames(catalogue) == dbs_types_78)) {
            plotSigsFun = plotSignaturesDBS
            par_mar_shift = c(1.4, 0, 2.5, 0) - c(0, 0, 1, 0) * draw_signatures
            line = 2.5
        } else if (type == "ID" || nrow(catalogue) == 83 && all(rownames(catalogue) == id_types_83)) {
            plotSigsFun = plotSignaturesID
            par_mar_shift = c(4, 0, 3.8, 0) - c(1, 0, 1, 0) * draw_signatures
            line = 4
        }

        #Normalisation of prob=T
        if (prob) {
            catal = prop.table(catalogue)
            exps = exposures / sum(exposures)
        } else {
            catal = catalogue
        }

        #Set mfrow and margin
        orig_par = par(c("mfrow", "mar"))
        if (draw_signatures) {
            exp_rows = max(nrow(exps), min_exp_rows)
            layout(
                rbind(1, cbind(2, 3:5), cbind(6, 7:(6 + exp_rows)), 7 + exp_rows),
                widths = c(.15, .85),
                heights = c((3 + .5 * exp_rows) / 31, rep(1,3), rep(.5, exp_rows), (3 + .5 * exp_rows) / 31)
            )
        } else {
            layout(
                rbind(1, cbind(2, 3:4), cbind(6, 5), 7),
                widths = c(.15, .85),
                heights = c(3 / 31, rep(1,3), 3 / 31)
            )
        }

        #Title
        par(mar=c(0,0,0,0))
        plot.new()
        title(main_text, line=-1.5, font.main=2, cex.main=1.3)

        #Left side
        #Choose signatures with non-zero contribution
        sigs_to_mention = rownames(exposures)[rowSums(exposures) > 0]

        #Legend
        par(mar=c(0, 1, 1, 2) + par_mar_shift * c(1,0,1,0) + .1)
        plot.new()
        legend(
            "bottom",
            sigs_to_mention,
            title = "Signatures",
            fill = col_exp[sigs_to_mention],
            bty = "n",
            cex = 1.3
        )

        #Right side
        #Fitted value
        fit_catal = signatures %*% exps

        #Compute ylim
        y_max = max(catal, fit_catal)

        #Top plot, observed frequencies
        par_mar_shift_max = max(par_mar_shift[c(1,3)])

        par(mar=c(0, 4.1, 1, 2.1) + par_mar_shift + c(max(par_mar_shift_max - par_mar_shift[3], 0), 0, 0, 0))
        plotSigsFun(catal, "Observed Profile", y_range=1.2*c(0, y_max), line=line, thin_plot=T, no_bottom=T)

        #Middle plot, difference
        par(
            mar = c(0, 4.1, 1, 2.1) + par_mar_shift + c(par_mar_shift[1], 0, par_mar_shift[3], 0) + .5 * (
                par_mar_shift_max - sum(par_mar_shift[c(1,3)])
            ) * c(1,0,1,0)
        )
        plotSigsFun(
            catal - fit_catal,
            paste("Difference", "-", err_text, "=", sprintf("%.3f", err_fun(catalogue, sum(catalogue) * fit_catal))),
            y_range = 1.2 * c(-1, 1) * max(0.03, max(abs(catal - fit_catal))),
            line = line, thin_plot = T, no_top = T, no_bottom = T
        )
        abline(h=0)

        #Bottom plot, expected frequencies
        par(mar=c(0, 4.1, 1, 2.1) + par_mar_shift + c(0, 0, max(par_mar_shift_max - par_mar_shift[1], 0), 0))
        plotSigsFun(fit_catal, "Reconstructed Profile", y_range=1.2*c(0, y_max), line=line, thin_plot=T, no_top=T)

        #Bottom

        #Signature exposures
        par(mar=c(0, 4, 0, 4) + .1 + c(par("mar")[1] - .1, 0, -.1 + max(par_mar_shift_max - par_mar_shift[1], 0), 0) * !draw_signatures)
        barplot(
            exps[rev(sigs_to_mention),,drop=F],
            space = 0,
            ylim = c(0,1),
            main = "", xlab = "", ylab = "",
            axisnames = F,
            col = rev(col_exp[sigs_to_mention]),
            xaxs = "i", yaxt = "n"
        )
        axis(
            2,
            axTicks(2),
            paste0(100*axTicks(2), "%"),
            las = 2,
            mgp = c(3, .75, 0)
        )
        box()

        #Scale y coordinates
        plot_y_frac = par("pin")[2] / (par("pin")[2] + sum(par("mai")[c(1,3)]))
        plot_y_min = 0 - par("mai")[1] / par("pin")[2]
        plot_y_max = 1 + par("mai")[3] / par("pin")[2]

        sig_midpoints = filter(cumsum(c(0, rev(exps))), c(.5, .5), sides = 2)[1:nrow(exps)]

        if (draw_signatures) {
            #Scale x coordinates
            plot_x_frac = par("pin")[1] / (par("pin")[1] + par("mai")[4])
            plot_x_out_half = (.5 + .5 * plot_x_frac) / plot_x_frac
            plot_x_out_full = 1 / plot_x_frac

            #Compute y coordinates of signatures
            plot_midpoints = filter(
                seq(plot_y_min, plot_y_max, length.out = 1 + exp_rows),
                c(.5, .5),
                sides = 2
            )[exp_rows - nrow(exps) + 1:nrow(exps)]

            #Draw the segment lines
            for (i in 1:nrow(exps)) {
                segments(1, sig_midpoints[i], plot_x_out_half, xpd=T)
                segments(plot_x_out_half, sig_midpoints[i], plot_x_out_full, plot_midpoints[i], xpd=T)
            }

            #Signatures
            par(mar=c(0, 4, 1, 2) + par_mar_shift + .1)
            for (sig_name in rownames(exps)) {
                #Signature
                plotSigsFun(
                    signatures[,sig_name, drop=F],
                    paste(sig_name, ":", sprintf("%.3f", exps[sig_name,])),
                    "", 1.2*c(0, max(signatures[,sig_name])), simple=T
                )

                #Lines to exposures
                segments(
                    par("usr")[1] - par("mai")[2] / par("pin")[1] * (par("usr") %*% c(-1, 1, 0, 0)),
                    mean(par("usr")[c(3,4)]),
                    par("usr")[1],
                    xpd = T
                )
            }

            par(mar=c(0,0,0,0))
            for (i in nrow(exps):exp_rows) {
                plot.new()
            }
        } else {
            #Draw the contributions
            axis(4, sig_midpoints, sprintf("%.3f", exps[rev(sigs_to_mention),]), las=1, xpd=T)

            par(mar=c(0,0,0,0))
            plot.new()
        }

        #Return original par values
        par(orig_par)
    }

#Multisample signature contribution barplot, 2 tiered
    plotSigExposures = function(exposures, col_exp, title="Mutations contributed", ylab.scale=1000, lab.cex=.8, lab.font=1, lab.srt=0, lab.adj=.5, lab.x.shift=0, lab.y.shift=0, mar.4.shift=2.5) {
        #Set mfrow, margin and las
        orig_par = par(c("mfrow", "mar", "las"))
        par(mfrow=c(2,1), mar=c(0, 4, 4, 3 + mar.4.shift) + .1, las=1, lwd=min(.8, 100 / ncol(exposures)))

        #Empty plot with grid lines
        barplot(
            rep(0,ncol(exposures)),
            ylim = 1.2 * c(0, max(colSums(exposures)) / ylab.scale),
            border = F,
            las = 2
        )
        abline(h=axTicks(2), col="#DDDDDD")
        par(new=T)

        #Scaling ylab
        if (ylab.scale == 1) {
            ylab = "mutation count"
        } else {
            ylab = paste0("mutation count (", ylab.scale, "s)")
        }

        #Absolute signature exposure (mutation counts) per sample
        x_pos = barplot(
            exposures[rev(rownames(exposures)),,drop=F] / ylab.scale,
            xlim = c(0, ncol(exposures)),
            ylim = 1.2 * c(0, max(colSums(exposures)) / ylab.scale),
            axisnames = F,
            main = title,
            ylab = ylab,
            col = rev(col_exp[rownames(exposures)]),
            yaxt = "n",
            xaxs = "i",
            space = 0
        )
        box()
        text(x_pos + lab.x.shift, (-0.05 + lab.y.shift) * par('usr')[4], colnames(exposures), cex=lab.cex, font=lab.font, srt=lab.srt, adj=lab.adj, xpd=NA)

        #Signature exposure (fraction) per sample
        par(mar=c(2, 4, 2, 3 + mar.4.shift) + .1)
        prop_exps = exposures / rep(colSums(exposures), each=nrow(exposures))
        if (length(prop_exps) == 1) prop_exps = matrix(prop_exps, dimnames=dimnames(exposures))
        x_pos = barplot(
            prop_exps[rev(rownames(exposures)),,drop=F],
            xlim = c(0, ncol(exposures)),
            ylim = c(0, 1),
            axisnames = F,
            ylab = "fraction",
            col = rev(col_exp[rownames(exposures)]),
            xaxs = "i",
            space = 0
        )
        box()
        text(x_pos + lab.x.shift, (-0.05 + lab.y.shift) * par('usr')[4], colnames(exposures), cex=lab.cex, font=lab.font, srt=lab.srt, adj=lab.adj, xpd=NA)

        par(mfrow=c(1,1), mar=c(2, 4 + 5*par("pin")[1], 4, .5) + .1, new=T)

        #Right-hand legend
        plot.new()
        legend(
            "left",
            sub("\\.", " ", rownames(exposures)),
            fill = col_exp[rownames(exposures)],
            bty = "n"
        )

        #Return original par values
        par(orig_par)
    }

# --------------------------------------------------
# -- Error plots
# --------------------------------------------------

#Model selection error plot; 2 tracking lines (y axes) + vertical lines
    plotSigFitSelErrors = function(x, y1, y2=NULL, v=NULL, main_text="", y_labs="Cosine Similarity", v_col="black", v_lty=5, v_lwd=1.5, v_legend=NULL) {
        #Set margin
        orig_par = par(c("mar", "las", "cex.axis", "mgp", "font.axis", "family", "xpd"))
        par(mar=c(5, 4, 4, 4) + .1)
        #First plot
            plot(
                rep(c(x,NA), length(y1) / length(x)),
                rbind(as.matrix(y1), NA),
                type = "b",
                col = "#0000FF99",
                xlim = c(min(x), max(x)),
                ylim = c(min(y1), max(y1)),
                main = main_text,
                xlab = "Number of signatures",
                ylab = "",
                xaxt = "n",
                las = 2
            )
            mtext(y_labs[1], 2, line = 3)
            #x axis labels
            axis(side=1, lwd=0, lwd.ticks=1, line=NULL, at=(min(x):max(x))[(1:length(min(x):max(x)) %% 2) == 0])
            axis(side=1, lwd=0, lwd.ticks=1, line=NULL, at=(min(x):max(x))[(1:length(min(x):max(x)) %% 2) == 1])

        #Second plot
        if (!is.null(y2)) {
            par(new=T)
            plot(
                rep(c(x,NA), length(y2) / length(x)),
                rbind(as.matrix(y2), NA),
                type = "b",
                col = "#FF000099",
                xlim = c(min(x), max(x)),
                ylim = c(min(y2), max(y2)),
                xlab = "",
                ylab = "",
                xaxt = "n",
                yaxt = "n",
                las = 2
            )

            #Y axis
            axis(4, las=2)
            mtext(y_labs[2], 4, line = 3)
        }

        #Vertical lines
        abline(v=v, lty=v_lty, col=v_col, lwd=v_lwd)

        #Legend
        legend(
            "right",
            c(y_labs, v_legend),
            lty = c(1, if (!is.null(y2)) 1 else NULL, v_lty),
            col = c("blue", if (!is.null(y2)) "red" else NULL, v_col),
            bg = "white",
            cex = 0.75
        )

        #Return original par values
        par(orig_par)
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


