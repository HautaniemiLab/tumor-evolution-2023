#Adjusts sigProfiler (COSMIC v3) signatures for
# - trinucleotide frequencies in targets vs hs37d5 for SBS signatures
# -  dinucleotide frequencies in targets vs hs37d5 for DBS signatures
#and COSMIC signatures for
# - trinucleotide frequencies in targets vs hs19 for SBS signatures
#and relabels ID deletion at repeat subtypes by incrementing repeat counts by 1.
#Also plots the resulting reference signatures.

#Options
    #Defaults
    TARGETS_FILE="data/raw/reference_genomes/GRCh38.d1.vd1.fa"
    RELABEL_INDELS=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -t | --targets )        shift
                                    TARGETS_FILE=$1 ;;
            --relabel-indels )      RELABEL_INDELS="TRUE"
        esac
        shift
    done

#Name variables
    IN_DIR_SIGS="data/raw/reference_signatures/"
    IN_DIR_KMERS="data/interim/reference_genomes/"
    SUFFIX=".tsv"
    SIG_CORRECTION_SCRIPT="src/features/s111_signature_correction.R"

    HG19="hg19"
    HS37D5="hs37d5"
    TARGET=`echo $TARGETS_FILE | perl -pe 's/.*\///; s/\..*//'`
    DINUC_SUFFIX=".dinuc_freq.tsv"
    TRINUC_SUFFIX=".trinuc_freq.tsv"

    COSMIC="COSMIC_signatures"
    LBA_SBS="sigProfiler_SBS_signatures"
    LBA_DBS="sigProfiler_DBS_signatures"
    LBA_ID="sigProfiler_ID_signatures"

    OUT_DIR_PLOTS="results/plots/reference_signatures/"
    OUT_DIR_SIGS="data/features/reference_signatures/"

#Make sure directories exist
    mkdir -p $OUT_DIR_PLOTS $OUT_DIR_SIGS

#Run R script that does the correction
    #COSMIC signatures
    $SIG_CORRECTION_SCRIPT \
        ${IN_DIR_SIGS}${COSMIC}$SUFFIX "SBS" \
        ${IN_DIR_KMERS}${HG19}$TRINUC_SUFFIX ${IN_DIR_KMERS}${TARGET}$TRINUC_SUFFIX \
        "${OUT_DIR_SIGS}${TARGET}.${COSMIC}$SUFFIX" \
        "${OUT_DIR_PLOTS}${TARGET}.${COSMIC}.pdf"
    #SBS signatures
    $SIG_CORRECTION_SCRIPT \
        ${IN_DIR_SIGS}${LBA_SBS}$SUFFIX "SBS" \
        ${IN_DIR_KMERS}${HS37D5}$TRINUC_SUFFIX ${IN_DIR_KMERS}${TARGET}$TRINUC_SUFFIX \
        "${OUT_DIR_SIGS}${TARGET}.${LBA_SBS}$SUFFIX" \
        "${OUT_DIR_PLOTS}${TARGET}.${LBA_SBS}.pdf"
    #DBS signatures
    $SIG_CORRECTION_SCRIPT \
        ${IN_DIR_SIGS}${LBA_DBS}$SUFFIX "DBS" \
        ${IN_DIR_KMERS}${HS37D5}$DINUC_SUFFIX ${IN_DIR_KMERS}${TARGET}$DINUC_SUFFIX \
        "${OUT_DIR_SIGS}${TARGET}.${LBA_DBS}$SUFFIX" \
        "${OUT_DIR_PLOTS}${TARGET}.${LBA_DBS}.pdf"

#Relabel ID deletion subtypes and plot the signatures
    if [[ $RELABEL_INDELS != "" ]]; then
        perl -ane '
            BEGIN{$header = <>; print $header}
            if ($F[0] =~ m/ ^(DEL_(C|T|repeats)_[0-9]\+?_) ([0-9]) (\+?) /x) {
                $F[0] = "$1" . ($3 + 1) . "$4";
            }
            print(join("\t", @F) . "\n");
        ' ${IN_DIR_SIGS}${LBA_ID}$SUFFIX \
            > ${OUT_DIR_SIGS}${LBA_ID}$SUFFIX
        Rscript -e " \
            source('src/visualisation/f2_signature_plotting_functions.R'); \
            args = commandArgs(trailingOnly=T); \
            sigs_id = read.table(args[1], header=T, as.is=T); \
            pdf(args[2], w=14, h=4); \
                par(mar=c(3.2,4,6,2) + .1); \
                for (i in 2:ncol(sigs_id)) { \
                    plotSignaturesID(sigs_id[,i], colnames(sigs_id)[i], 'probability', thin_plot=F, line=4); \
                }; \
            invisible(dev.off()); \
        " "${OUT_DIR_SIGS}${LBA_ID}$SUFFIX" "${OUT_DIR_PLOTS}${LBA_ID}.pdf"
    fi

