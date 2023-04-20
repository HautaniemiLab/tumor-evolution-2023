#R script for loading in SigProfiler attribution rules in desired format
#Input file names:
IN_DIR = "data/raw/reference_signatures/"
RULES_CONNECTED_SIGS = paste0(IN_DIR, "sigProfiler_SBS_signatures_rules_connected.tsv")
RULES_MUTATION_COUNTS = paste0(IN_DIR, "sigProfiler_SBS_signatures_rules_count.tsv")
RULES_TS_BIAS = paste0(IN_DIR, "sigProfiler_SBS_signatures_rules_tsb.tsv")
#Output variable names:
# - sigs_sbs_rules_conn
# - sigs_sbs_rules_mc
# - sigs_sbs_rules_tsb

# --------------------------------------------------
# Load SigProfiler attribution rules
# --------------------------------------------------

#Connected signature rules
    #Read file
    sigs_sbs_rules_conn_file = read.table(
        RULES_CONNECTED_SIGS,
        header = T, as.is = T
    )

    #Break down the table into list of connected signature setd
    sigs_sbs_rules_conn = apply(
        sigs_sbs_rules_conn_file,
        1,
        function(categ) {
            out = unlist(strsplit(categ[2], ","))
            names(out) = NULL
            out
        }
    )
    names(sigs_sbs_rules_conn) = sigs_sbs_rules_conn_file[,1]

#Signature count rules
    #Read file
    sigs_sbs_rules_mc_file = read.table(
        RULES_MUTATION_COUNTS,
        header = T, as.is = T
    )

    #Break down as a vector of signatures and thresholds
    sigs_sbs_rules_mc = sigs_sbs_rules_mc_file[,3]
    names(sigs_sbs_rules_mc) = sigs_sbs_rules_mc_file[,1]

#Transcription strand bias rules
    #Read file
    sigs_sbs_rules_tsb_file = read.table(
        RULES_TS_BIAS,
        header = T, as.is = T, check.names = F
    )

    #
    sigs_sbs_rules_tsb = sigs_sbs_rules_tsb_file[,-1]
    rownames(sigs_sbs_rules_tsb) = sigs_sbs_rules_tsb_file[,1]

# --------------------------------------------------
# End of file
# --------------------------------------------------


