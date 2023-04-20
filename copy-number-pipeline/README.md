# Copy-Number Pipeline

This [Anduril](https://anduril.org/site/) pipeline (`gatk-ascat.scala`) performs
copy-number segmentation using GATK and estimates the purity and ploidy using a
reimplemented ASCAT algorithm. The pipeline uses _TP53_ VAF as an additional
evidence for purity in the model selection.

The `ascat-parameter-overrides.xlsx` file allows for steering ASCAT in the model
selection.

Note: some files have been excluded from this repository as they contain
sensitive sample information.

## Results

The pipeline outputs the following files and directories into the
result_gatk_ascat/output/ folder:

### ascatProfiles-out/

Plots of the true integer copy numbers of major and minor alleles estimated by ASCAT.

### ascatStats-table.csv

Some statistics computed from the segments. Mainly the number of breakpoints.

### ascatSunrises-out/

Plots of ascat cost surfaces. ASCAT uses these to choose the optimal ploidy/purity combination.
The black points and error bars are _TP53_-based purity estimates.

### cnGenesFolder-out/

Per-gene copy-number statistics for each patient.

### combinedAscatEstimates-out.csv

Purity and ploidy estimates for each sample. Note: these may have some errors in purities
because ASCAT is unable to estimate the aberrant-status for our faked SNP data.

### combinedAscatSegments-out.csv

All copy-number segments for all samples. Each segment also contains ASCAT-estimates
for major and minor allele counts.

### combinedCnGenes-out.csv

All genes for all patients in a single file, which is huge but easy to filter.

### gatkPlotsFolder-out/

Plots of copy-number segmentations produced by GATK.

### multiPloidyPatients-bash-folder1/

Ploidy plots for patients that have samples with discrepant ploidy estimates.

# Author

Kari Lavikka
