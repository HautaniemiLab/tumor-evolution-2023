#!/usr/bin/env Rscript

#R script to adjust input reference signatures based on kmer frequencies of
#source and target genomes

#Load plotting functions
source("src/visualisation/f2_signature_plotting_functions.R")

#Command line arguments, contains input signatures file, type of signature (SBS/
#DBS), input and target genome counts file, output and plot file names
    args = commandArgs(trailingOnly=TRUE)
    if (length(args) < 5) {
        stop(
            "Please supply input signatures file, type of signature (SBS/DBS), input and target genome counts file, output and plot file names",
            call. = F
        )
    } else {
        #input signatures file
        in_file = args[1]
        #type (SBS or DBS)
        type = args[2]
        #input original counts
        src_counts_file = args[3]
        #input new counts
        out_counts_file = args[4]
        #output signatures file
        out_file = args[5]
        #output plot file
        plot_file = args[6]
    }

#Read signatures file
    src_signatures = read.table(in_file, header=T, as.is=T)

#Read genome kmer frequency tables
    src_counts = read.table(src_counts_file, header=T, as.is=T)
    out_counts = read.table(out_counts_file, header=T, as.is=T)

#Ratio of kmer frequencies between the builds
    kmer_ratio = out_counts[order(src_counts$kmer), "count"] / src_counts[order(src_counts$kmer), "count"]
    names(kmer_ratio) = src_counts$kmer[order(src_counts$kmer)]

#Get signature types
    col_idx = ifelse(
        type == "SBS",
        which(endsWith(colnames(src_signatures), "Type") & nchar(colnames(src_signatures)) > 4),
        which(endsWith(colnames(src_signatures), "Type"))
    )
    types = unlist(sapply(strsplit(src_signatures[,col_idx], ">"), function(x) x[1]))

#Signature correction
    out_signatures = src_signatures
    out_signatures[,!endsWith(colnames(src_signatures), "Type")] = apply(
        out_signatures[,!endsWith(colnames(src_signatures), "Type")],
        2,
        function(signature) {
            prop.table(kmer_ratio[types] * signature)
        }
    )

#Signature barplots
    pdf(plot_file, width = 12, height = 4)
        par(mar=c(2.1, 4.1, 4.1, 2.1))
        for (sig in colnames(out_signatures)[!endsWith(colnames(out_signatures), "Type")]) {
            if (type == "SBS") plotSignaturesSBS(out_signatures[,sig], sig, "probability")
            if (type == "DBS") plotSignaturesDBS(out_signatures[,sig], sig, "probability")
        }
        par(mar=c(5.1, 4.1, 4.1, 2.1))
    invisible(dev.off())

#Save new signature table
    write.table(
        out_signatures,
        out_file,
        quote = F,
        sep = "\t",
        row.names = F
    )


