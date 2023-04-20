#Shared variables
BASE_NAME="discovery"

SAMPLE_LIST="data/discovery_subclones.csv"
CLUSTER_INFO="data/raw/discovery.subclones.csv"

VARS_DIR="data/raw/discovery_vcf"
CLUSTERS_DIR="data/raw/discovery.subclones"
CLUSTERS_SUFFIX="_cluster_variants.csv"

NORMAL_REGEX="^[A-Za-z]+[0-9]+_(BDNA|BBC|ME)"

# No filtering
# ============

# Sample-level prep
src/s100_filter_variants_and_compute_spectra.sh -c 8 \
    -n "$BASE_NAME.no_filt" \
    -i "$VARS_DIR" -s "$SAMPLE_LIST" \
    -r "$NORMAL_REGEX" --vaf 0.001 --fwd 0 --rev 0 \
    -vf 'filter( nchar(REF) == 1 & nchar(ALT) == 1 )'

# Sample-level prep
src/s300_create_sample_trees_and_spectra.sh -c 8 \
    -n "$BASE_NAME.no_filt" \
    -r "$NORMAL_REGEX" -s \
    -min-depth 1 -min-rc-patient 1 -min-rc-sample 1

# No filtering subclones prep
src/s400_create_subclonal_spectra.sh -c 20 \
    -n "$BASE_NAME.no_filt" \
    -i "$CLUSTERS_DIR" -s-suffix "$CLUSTERS_SUFFIX" \
    -s-info "$CLUSTER_INFO"

# No filtering subclones attribution
src/s410_attribute_subclonal_sigs.sh -c 20 \
    -n "$BASE_NAME.no_filt" \
    -s-info "$CLUSTER_INFO" -z \
    --sensitive-mode -c-sigs "SBS1,SBS2,SBS3,SBS5,SBS13,SBS40"



# No filtering artefacts
# ======================

# Use base prep
ln -s "filt_set_$BASE_NAME.no_filt" "data/features/filt_set_$BASE_NAME.no_filt.arts"

# No filtering artefacts subclones attribution
src/s410_attribute_subclonal_sigs.sh -c 20 \
    -n "$BASE_NAME.no_filt.arts" \
    -s-info "$CLUSTER_INFO" -z \
    --sensitive-mode -c-sigs "SBS1,SBS2,SBS3,SBS5,SBS13,SBS40" \
    -f "SBS1,SBS5,SBS34,SBS41,SBS47,SBS57,SBS58,SBS90"



# Base with filtering
# ===================

# Sample-level prep
src/s100_filter_variants_and_compute_spectra.sh -c 8 \
    -n "$BASE_NAME" \
    -i "$VARS_DIR" -s "$SAMPLE_LIST" \
    -r "$NORMAL_REGEX" --vaf 0.001 --fwd 0 --rev 0 \
    --blacklist "data/raw/blacklists/pon.4.1.9.0.sets2-8c.TCGA.vcf.gz" --filter-mbs-as-sbs --dont-keep-nested-sbs \
    --info-query '\t%FILTER\t%NALOD\t%CONTQ' \
    -vf 'filter( NALOD >= 1.0 & CONTQ >= 10 )' \
    -vf 'filter( nchar(REF) == 1 & nchar(ALT) == 1 )'

# Sample trees
src/s300_create_sample_trees_and_spectra.sh -c 8 \
    -n "$BASE_NAME" \
    -r "$NORMAL_REGEX" -s \
    -min-depth 1 -min-rc-patient 1 -min-rc-sample 1

# Base subclones prep
src/s400_create_subclonal_spectra.sh -c 20 \
    -n "$BASE_NAME" \
    -i "$CLUSTERS_DIR" -s-suffix "$CLUSTERS_SUFFIX" \
    -s-info "$CLUSTER_INFO"

# Base subclones attribution
src/s410_attribute_subclonal_sigs.sh -c 20 \
    -n "$BASE_NAME" \
    -s-info "$CLUSTER_INFO" -z \
    --sensitive-mode -c-sigs "SBS1,SBS2,SBS3,SBS5,SBS13,SBS40"



# Artefacts
# =========

# Use base prep
ln -s "filt_set_$BASE_NAME" "data/features/filt_set_$BASE_NAME.arts"

# Artefacts attribution (forced artefacts)
src/s410_attribute_subclonal_sigs.sh -c 20 \
    -n "$BASE_NAME.arts" \
    -s-info "$CLUSTER_INFO" -z \
    --sensitive-mode -c-sigs "SBS1,SBS2,SBS3,SBS5,SBS13,SBS40" \
    -f "SBS1,SBS5,SBS34,SBS41,SBS47,SBS57,SBS58,SBS90"



# Final filtered clusters
# =======================

# Sample-level prep
src/s100_filter_variants_and_compute_spectra.sh -c 20 \
    -n "$BASE_NAME.cluster_filt" \
    -i "$VARS_DIR.mouse_filtered/" -s "$SAMPLE_LIST" \
    -r "$NORMAL_REGEX" --vaf 0.001 --fwd 0 --rev 0 \
    --blacklist "data/raw/blacklists/pon.4.1.9.0.sets2-8c.TCGA.vcf.gz" --filter-mbs-as-sbs --dont-keep-nested-sbs \
    --info-query '\t%FILTER\t%NALOD\t%CONTQ' -vf 'filter( NALOD >= 1.0 & CONTQ >= 10 )'

# Sample-level attribution (for GenomeSpy etc.)
src/s200_attribute_signatures_and_plot.sh -c 20 \
    -n "$BASE_NAME.cluster_filt" \
    -r "$NORMAL_REGEX" \
    -c-sigs "SBS1,SBS2,SBS3,SBS5,SBS13,SBS40,DBS2,DBS4,DBS5,DBS6,DBS9,ID1,ID2,ID4,ID5,ID6,ID8,ID9"



