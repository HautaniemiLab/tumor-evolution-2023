#Script for plotting both the supersample as well as tree branch mutational
#spectra for sample trees.

#Options
    #Defaults
    NUM_ID="0"
    PLOT_SEPARATE=""

    #Command Line
    while [[ $# -gt 0 ]]; do
        case $1 in
            -n | --num-id )         shift
                                    NUM_ID=$1 ;;
            -p | --plot-separate )  PLOT_SEPARATE="-p"
        esac
        shift
    done

#Name variables
    PLOT_SCRIPT="src/visualisation/s01_mutational_spectra_plotting.R"

    IN_DIR="data/features/filt_set_$NUM_ID/sample_tree/mutational_spectra/"
    IN_SBS=".tsv"

    OUT_DIR="results/plots/filt_set_$NUM_ID/sample_tree/mutational_profiles/"

#Make sure output directory exists
    mkdir -p $OUT_DIR

#Clear previous content if it exists
    rm -f $OUT_DIR/*.pdf

#Plot the spectra
    #Supersamples
    $PLOT_SCRIPT -sbs "${IN_DIR}supersample_spectra_sbs.tsv" \
        -dbs "${IN_DIR}supersample_spectra_dbs.tsv" \
        -id "${IN_DIR}supersample_spectra_id.tsv" \
        -o $OUT_DIR -out-suf "supersample_profiles" &
    #Sample tree branches
    $PLOT_SCRIPT -sbs "${IN_DIR}sample_tree_branch_spectra_sbs.tsv" \
        -dbs "${IN_DIR}sample_tree_branch_spectra_dbs.tsv" \
        -id "${IN_DIR}sample_tree_branch_spectra_id.tsv" \
        -o $OUT_DIR -out-suf "sample_tree_branch_profiles" $PLOT_SEPARATE &

wait

