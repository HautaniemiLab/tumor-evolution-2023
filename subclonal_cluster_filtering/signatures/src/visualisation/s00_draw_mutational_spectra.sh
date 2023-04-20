#Script for plotting both the supersample as well as single sample mutational
#spectra.
#Normals are plotted separately from tumour single samples.

#Options
    #Defaults
    NUM_ID="0"
    PLOT_SEPARATE=""
    REGEX_NORMALS=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )         shift
                                    NUM_ID=$1 ;;
            -p | --plot-separate )  PLOT_SEPARATE="-p" ;;
            -r | --regex-normals )  shift
                                    REGEX_NORMALS=$1
        esac
        shift
    done

#Name variables
    PLOT_SCRIPT="src/visualisation/s01_mutational_spectra_plotting.R"

    IN_DIR="data/features/filt_set_$NUM_ID/mutational_spectra/"
    IN_SBS=".tsv"

    OUT_DIR="results/plots/filt_set_$NUM_ID/mutational_profiles/"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*.pdf

#Plot the spectra
    #Supersamples
    $PLOT_SCRIPT -sbs "${IN_DIR}supersample_spectra_sbs.tsv" \
        -dbs "${IN_DIR}supersample_spectra_dbs.tsv" \
        -id "${IN_DIR}supersample_spectra_id.tsv" \
        -o $OUT_DIR -out-suf "supersample_profiles"
    #Single samples
    $PLOT_SCRIPT -sbs "${IN_DIR}single_sample_spectra_sbs.tsv" \
        -dbs "${IN_DIR}single_sample_spectra_dbs.tsv" \
        -id "${IN_DIR}single_sample_spectra_id.tsv" \
        -o $OUT_DIR -out-suf "sample_profiles" $PLOT_SEPARATE -r $REGEX_NORMALS

