#R script for loading in reference signatures in desired format
#Input file names:
IN_DIR = "data/features/reference_signatures/"
COSMIC_SIGS = paste0(IN_DIR, "GRCh38.COSMIC_signatures.tsv")
SIG_PROFILER_SBS = paste0(IN_DIR, "GRCh38.sigProfiler_SBS_signatures.tsv")
SIG_PROFILER_DBS = paste0(IN_DIR, "GRCh38.sigProfiler_DBS_signatures.tsv")
SIG_PROFILER_ID = paste0(IN_DIR, "sigProfiler_ID_signatures.tsv")
#Output variable names:
# - COSMIC:
#   - ref_sigs_cosmic
# - SigProfiler signature:
#   - ref_sigs_sbs
#   - ref_sigs_dbs
#   - ref_sigs_id

# --------------------------------------------------
# Helper function for loading signatures
# --------------------------------------------------

#Function for loading non-default reference signatures
loadSignatures = function(path, mut_class) {
    if (mut_class %in% c("SBS", "COSMIC")) {
        #Read file
        ref_sigs = read.table(path, header = T, as.is = T)

        #Compute the 96 subtypes from the 6 types and 16 trinucleotide contexts
        rownames(ref_sigs) = apply(
            ref_sigs[,c("Type","SubType")],
            1,
            function(x) {
                bases = unlist(strsplit(x[2], ""))
                paste0(bases[1], "[", x[1], "]", bases[3])
            }
        )
        ref_sigs = as.matrix(ref_sigs[,-(1:2)])
    } else {
        #Read file
        ref_sigs = as.matrix(
            read.table(path, header = T, row.names = 1, as.is = T)
        )
    }
    return(ref_sigs)
}

# --------------------------------------------------
# Load COSMIC signatures
# --------------------------------------------------

#SBS only
    ref_sigs_cosmic = loadSignatures(COSMIC_SIGS, "SBS")

# --------------------------------------------------
# Load SigProfiler signatures
# --------------------------------------------------

#SBS signatures
    ref_sigs_sbs = loadSignatures(SIG_PROFILER_SBS, "SBS")

#DBS signatures
    ref_sigs_dbs = loadSignatures(SIG_PROFILER_DBS, "DBS")

#ID signatures
    ref_sigs_id = loadSignatures(SIG_PROFILER_ID, "ID")

# --------------------------------------------------
# End of file
# --------------------------------------------------


