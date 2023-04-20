#!/usr/bin/env anduril
//$OPT --wrapper slurm-prefix

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.sequencing._

// Anduril workflow for evolution project filtering of variants caused
// by mouse contamination.
//
// Works by screening variants with ALT reads only in affected samples,
// which will define the target variants to check for contamination.
// Reads of the affected sample aligning around the target variants will
// be extracted as FASTQ and then realigned against a concatenated
// genome of human and mouse reference genomes. Finally, the target
// variants are force-called. Variants that disappear or become un-
// callable will be considered filtered and will make a blacklist for
// that patient.

object FilterMouseContamination {
    //Resource paths
    val RESOURCES_DIR = "/path/to/resources/"
    val PIPELINE_RES_DIR = RESOURCES_DIR + "pipeline_resources/"
    val GATK_RES_DIR     = PIPELINE_RES_DIR + "gatk/"
    val PIPELINE_SCRIPTS = "/path/to/wgs_scripts/alignment/"

    //Reference sequence and PoN
    val reference = INPUT("data/reference_genome/GRCh38.d1.vd1_GRCm39.toplevel.fa.gz")
    val contigs   = INPUT("data/reference_genome/contigs.tsv")
    val pon       = INPUT(PIPELINE_RES_DIR + "pon/pon_discovery.vcf.gz")

    //GATK, BEDtools, Picard and bcftools folders
    val gatk = "/opt/share/gatk-4.1.4.1/"
    val bedtools = "/opt/share/bedtools-2.28.0/bin/"
    val picard   = "/opt/share/picard-2.6/"
    val bcftools = "/opt/share/bcftools-1.11/"

    //Temporary directory base path for sorting
    val TEMP_PATH = "./.tmp/"

    //Input tables
    val listSampleNames     = INPUT("data/sample_names.tsv")
    val listAffectedSamples = INPUT("data/affected_samples.tsv")
    val listBams            = INPUT("data/sample_info_extended.csv")

    //Set 5 read group metadata. Since there is only a single sample per read group, we can get the correct metadata...
    val listReadGroupMetadata = INPUT("/path/to/read_group_metadata.csv")

    //Process input lists
    val listMerged = CSVDplyr(
        csv1 = listSampleNames,
        csv2 = listAffectedSamples,
        csv3 = listBams,
        function1 = """
            mutate(affected = evolName %in% csv2$sample) %>%
                merge(csv3, by.x="currName", by.y="sample")
        """
    )

    //Split BAMs by read groups
    //-------------------------

    //Output variables
    val bamsByReadGroupOut = NamedMap[Any]("bamsByReadGroupOut")

    //Iterate over affected samples
    for ( rowMap <- iterCSV(listMerged).filter(_("affected") == "TRUE") ) {
        val sample = rowMap("currName")
        val bam    = rowMap("bamFile")

        withName(sample) {
            //Split BAMs by read groups and sort by read name
            val bamsByReadGroup = bash"""
                $PIPELINE_SCRIPTS/split_by_readgroups.sh \\
                    --link --do-not-index -F 2304 \\
                    -i $bam --sample-key $sample -o @out1@ -d @folder1@

                function sort_by_name {
                    local FILE=$$1

                    samtools sort -n $$FILE -o $$FILE.temp
                    mv $$FILE.temp $$FILE
                }

                for FILE in @folder1@/*.bam; do
                    sort_by_name $$FILE &
                done

                wait
            """
            bamsByReadGroup._filename("out1", "read_groups.tsv")
            bamsByReadGroup._filename("folder1", "bams")
            bamsByReadGroup._custom("memory") = "14336"
            bamsByReadGroup._execute = "once"

            bamsByReadGroupOut(sample) = bamsByReadGroup.out1
        }
    }

    //Create a table of each read group's read name sorted BAMs with underscore separators
    //Since there is only a single sample per read group, we can get the correct metadata
    //by merging with sample name...
    val readGroupBamsCSV = CSVDplyr(
        csv1 = CSVListJoin(in = bamsByReadGroupOut, fileCol = ""),
        csv2 = listMerged,
        csv3 = listReadGroupMetadata,
        function1 = """
            merge(csv2 %>% select(currName, evolName), by.x="sample", by.y="currName") %>%
                merge(csv3 %>% select(readGroupID, sample, platform, library), by="sample") %>%
                select(readGroupID, sample, platform, library, bam=file)
        """
    )

    //Realign read groups
    //-------------------

    //Output variables
    val sortedOut = NamedMap[Any]("sortedOut")

    //Iterate through each read group
    for ( rowMap <- iterCSV(readGroupBamsCSV) ) {
        val readGroupID = rowMap("readGroupID")
        val sample      = rowMap("sample")
        val platform    = rowMap("platform")
        val library     = rowMap("library")
        val bam         = rowMap("bam")

        withName(readGroupID) {
            //Convert BAM to compressed paired FASTQs
            val bamToFastq = bash"""
                $bedtools/bamToFastq -i $bam \\
                    -fq >(gzip -c > @out1@) -fq2 >(gzip -c > @out2@)
            """
            bamToFastq._filename("out1", readGroupID + "_1.fq.gz")
            bamToFastq._filename("out2", readGroupID + "_2.fq.gz")
            bamToFastq._custom("cpu") = "3"
            bamToFastq._execute = "once"

            //Realign paired-ended reads with BWA MEM against the concatenated reference
            val realigned = bash"""
                bwa mem -M -t 8 -K 10000000 \\
                    -R "@RG\\tID:$readGroupID\\tPL:$platform\\tLB:$library\\tSM:$sample" \\
                    $reference ${bamToFastq.out1} ${bamToFastq.out2} |
                    samtools view -b > @out1@
            """
            realigned._keep = false
            realigned._filename("out1", "aligned.bam")
            realigned._custom("cpu") = "8"
            realigned._custom("memory") = "5632"
            realigned._execute = "once"

            //Sort the aligned reads by coordinate with PICARD
            val sorted = bash"""
                TEMP_DIR="$TEMP_PATH/$readGroupID-sorted-bash"
                mkdir -p $$TEMP_DIR
                trap "rm -rf $$TEMP_DIR" EXIT

                java -Xmx8G -jar $picard/picard.jar SortSam \\
                    VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate \\
                    CREATE_INDEX=true CREATE_MD5_FILE=true \\
                    MAX_RECORDS_IN_RAM=1200000 \\
                    TMP_DIR=$$TEMP_DIR \\
                    INPUT=${realigned.out1} \\
                    OUTPUT=@out1@
            """
            sorted._filename("out1", "sorted.bam")
            sorted._custom("memory") = "3584"
            sorted._execute = "once"

            sortedOut(readGroupID) = sorted.out1
        }
    }

    //Create a table of each sample's read group specific sorted BAMs
    val readGroupBamsRealignedCSV = CSVDplyr(
        csv1 = readGroupBamsCSV,
        csv2 = Array2CSV(sortedOut, _name = "sortedOutCSV"),
        function1 = """
            merge(csv2, by=1, sort=F) %>%
                select(sample, readGroupID, bam=File)
        """
    )

    //Combine read groups to samples and postprocess
    //----------------------------------------------

    //Output variables
    val markedOut = NamedMap[Any]("markedOut")

    //Iterate through each sample
    val samples = iterCSV(readGroupBamsRealignedCSV).map(_("sample")).toList.distinct

    for (sample <- samples) {
        withName(sample) {
            //Generate input parameter of sample's read group BAM files
            val inputMarked = iterCSV(readGroupBamsRealignedCSV)
                .filter(_("sample") == sample)
                .map("INPUT=" + _("bam"))
                .mkString(" ")

            //Merge BAMs and mark duplicates with Picard MarkDuplicates
            val marked = bash"""
                java -Xmx4G -jar $picard/picard.jar MarkDuplicates \\
                    VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \\
                    ASSUME_SORT_ORDER=coordinate \\
                    "$inputMarked" \\
                    OUTPUT=@out1@ METRICS_FILE=@out2@
            """
            marked._filename("out1", "marked.bam")
            marked._filename("out2", "markedMetrics.txt")
            marked._custom("memory") = "4096"
            marked._execute = "once"

            markedOut(sample) = marked.out1
        }
    }

    //Merge the realigned BAMs, replacing the original in the table
    val listBamsRealigned = CSVDplyr(
        csv1 = listMerged,
        csv2 = Array2CSV(markedOut, _name = "markedOutCSV"),
        function1 = """
            merge(csv2, by=1, all.x=T) %>%
                select(
                    patient, evolName, currName, affected, bamFile, realignedBamFile=File, normal,
                    normalBamFile=variantCallingNormalBamFile, normalSample=variantCallingNormalSample
                ) %>%
                mutate(bamFile = ifelse(is.na(realignedBamFile), bamFile, realignedBamFile))
        """
    )

    //Force call variants
    //-------------------

    //Output variables
    val forceCalledFilteredOut = NamedMap[Any]("forceCalledFilteredOut")
    val varsFilteredCSVOut     = NamedMap[Any]("varsFilteredCSVOut")

    //Iterate over patients
    val patients = iterCSV(listBamsRealigned).filter(_("affected") == "TRUE").map(_("patient")).toList.distinct

    for (patient <- patients) {
        withName(patient) {
            //Original patient VCF
            val originalVcf = INPUT("data/variants_vcf/" + patient + ".vcf.gz")

            //Generate tumour and normal inputs for variant calling
            val inputForceCalledTumours = iterCSV(listBamsRealigned)
                .filter(_("patient") == patient)
                .filter(_("normal") == "FALSE")
                .map("-I " + _("bamFile"))
                .mkString(" ")

            val inputForceCalledNormal = iterCSV(listBamsRealigned)
                .filter(_("patient") == patient)
                .filter(_("normal") == "TRUE")
                .map(row => "-I " + row("normalBamFile") + " -normal " + row("normalSample"))
                .toList.distinct(0)

            //Force genotype tumours jointly in multisample setting with normal with GATK Mutect2
            val forceCalled = bash"""
                #// Setup temporary directory
                eval "${Tasks.SET_TEMP_DIR}"

                #// Forced calling
                $gatk/gatk --java-options "-Xmx4G" Mutect2 \\
                    -R $reference -L $originalVcf -ip 250 \\
                    -alleles $originalVcf \\
                    -pon $pon \\
                    --force-call-filtered-alleles \\
                    -emit-lod 1e999 \\
                    "$inputForceCalledTumours" "$inputForceCalledNormal" \\
                    --tmp-dir $$TEMP_DIR -O @out1@
            """
            forceCalled._filename("out1", "out1.vcf.gz")
            forceCalled._custom("cpu") = "4"
            forceCalled._execute = "once"

            //In some rare cases where some variants are not emitted with
            //default settings, retry forced calling for missing alleles using
            //Mutect2 but with '-init-lod -2.15', and add the missing records to
            //the force-called VCF.
            //TODO: Note that -2.15 was chosen as it is the value that happened
            //      to emit all of the missing ones; too low or too high a value
            //      results in missing variants!
            //      In some cases, it is therefore possible that a single run of
            //      salvaging missing variants is insufficient (due to lack of
            //      catch-all '-init-lod' value) - multiples of this component
            //      with different '-init-lod' values would need to be daisy-
            //      chained. In such cases, converting this component to a
            //      function would be advisable
            val forceCalledMissing = bash"""
                BCFTOOLS="$bcftools"
                GATK="$gatk"

                FORCE_CALLED="${forceCalled.out1}"

                #// Setup temporary directory
                eval "${Tasks.SET_TEMP_DIR}"

                #// Check if there are missing records
                $$BCFTOOLS/bcftools isec -c any -n-1 -w1 \\
                    $originalVcf $$FORCE_CALLED \\
                    -Oz -o $$TEMP_DIR/temp.missing_variants.vcf.gz

                tabix -p vcf $$TEMP_DIR/temp.missing_variants.vcf.gz

                #// Force-call missing records and add them if there are any
                if [[ $$($$BCFTOOLS/bcftools view -H $$TEMP_DIR/temp.missing_variants.vcf.gz | wc -l) -gt 0 ]]; then
                    #// Forced genotyping
                    $$GATK/gatk --java-options "-Xmx4G" Mutect2 \\
                        -R $reference -L $$TEMP_DIR/temp.missing_variants.vcf.gz -ip 250 \\
                        -alleles $$TEMP_DIR/temp.missing_variants.vcf.gz \\
                        -pon $pon \\
                        --force-call-filtered-alleles \\
                        -emit-lod 1e999 -init-lod -2.15 \\
                        "$inputForceCalledTumours" "$inputForceCalledNormal" \\
                        --tmp-dir $$TEMP_DIR -O $$TEMP_DIR/temp.missing_variants.force_called.vcf.gz

                    #// Merge the missing records
                    $$BCFTOOLS/bcftools isec -c any -n-1 -w1 \\
                        $$TEMP_DIR/temp.missing_variants.force_called.vcf.gz $$FORCE_CALLED \\
                        -o $$TEMP_DIR/temp.missing_variants.force_called.private.vcf

                    $$GATK/gatk MergeVcfs \\
                        -I $$FORCE_CALLED -I $$TEMP_DIR/temp.missing_variants.force_called.private.vcf \\
                        -O @out1@
                else
                    #// Copy successfully genotyped variants
                    $$GATK/gatk SelectVariants -V $$FORCE_CALLED -O @out1@
                fi
            """
            forceCalledMissing._filename("out1", "out1.vcf.gz")
            forceCalledMissing._custom("cpu") = "4"

            //Affected sample (original name)
            val affectedSampleEvol = iterCSV(listBamsRealigned)
                .filter(_("patient") == patient)
                .filter(_("affected") == "TRUE")
                .map(_("evolName"))
                .toList.distinct(0)

            //Remove sites without desired variants and determine uncalled/mouse contamination variants
            //Only sites originally private to the affected sample are considered
            //Also output as CSV
            val forceCalledFiltered = bash"""
                BCFTOOLS=$bcftools

                #// Setup temporary directory
                eval "${Tasks.SET_TEMP_DIR}"

                #// Determine variants that were originally called from the affected sample
                echo "$affectedSampleEvol" > $$TEMP_DIR/sample.txt

                $$BCFTOOLS/bcftools view \\
                    -i "FORMAT/AD[@$$TEMP_DIR/sample.txt:1] > 0 && FORMAT/AD[@$$TEMP_DIR/sample.txt:1] == SUM(FORMAT/AD[1-:1])" \\
                    $originalVcf \\
                    -Oz -o $$TEMP_DIR/private_sites.vcf.gz
                tabix -p vcf $$TEMP_DIR/private_sites.vcf.gz

                #// Remove non-force-called new variants
                $$BCFTOOLS/bcftools isec -c any -n=2 -w1 \\
                    ${forceCalledMissing.out1} $originalVcf \\
                    -Oz -o $$TEMP_DIR/force_called_variants.vcf.gz
                tabix -p vcf $$TEMP_DIR/force_called_variants.vcf.gz

                #// Determine mouse contamination variants
                $$BCFTOOLS/bcftools annotate \\
                    -a $$TEMP_DIR/private_sites.vcf.gz -m "private_site" \\
                    $$TEMP_DIR/force_called_variants.vcf.gz |
                    $$BCFTOOLS/bcftools filter \\
                        -e 'private_site == 1 && SUM(FORMAT/AD[1-:1]) == 0' \\
                        -s "mouse_contamination" \\
                        -Oz -o @out1@

                tabix -p vcf @out1@

                #// CSV output
                $$BCFTOOLS/bcftools view @out1@ |
                    grep -v "^##" |
                    cut -f-2,4,5,7 |
                    sed 's/^#//' \\
                        > @out2@
            """
            forceCalledFiltered._filename("out1", "mouse_contam_filtered.vcf.gz")
            forceCalledFiltered._filename("out2", "mouse_contam_filtered.csv")

            forceCalledFilteredOut(patient) = forceCalledFiltered.out2

            //Original VCF converted to CSV with contamination variants removed
            val varsFilteredCSV = CSVDplyr(
                csv1 = INPUT("data/variants_csv/" + patient + ".csv"),
                csv2 = forceCalledFiltered.out2,
                function1 = """
                    mutate(index = row_number()) %>%
                        merge(csv2 %>% filter(FILTER == "PASS") %>% select(CHROM, POS))
                """
            )

            varsFilteredCSVOut(patient) = varsFilteredCSV.out
        }
    }

    //Collect filtered variant tables in a folder
    val varsFolderCSV = Array2Folder(forceCalledFilteredOut, fileMode = "@key@.csv")
    val varsFolderFilteredCSV = Array2Folder(varsFilteredCSVOut, fileMode = "@key@.csv")

    //Compute coverages
    //-----------------

    //R script for computing mean coverage as well as total base count over human
    //and mouse contigs (incl. non-primary ones, excl. viruses). Also estimates
    //degree of contamination.
    val estimatedCoveragesRscript = StringInput(
        """
            # Library
            library(tidyverse)

            # Inputs
            genomecov_path = var1
            contigs = table1
            sample_name = param1

            genomecov = read_tsv(
                genomecov_path, col_types="cinid",
                col_names = c("contig", "depth", "count", "contigLength", "fraction")
            ) %>%
                mutate(bases = depth * count)

            # Compute total bases and mean coverages
            selectContigs = function(species, main=F) {
                contigs %>%
                    filter(species == !!species & (primary | !main)) %>%
                    pull(contig)
            }

            genomeLength = function(species, main=F) {
                contigs %>%
                    filter(species == !!species & (primary | !main)) %>%
                    pull(length) %>%
                    sum()
            }

            cov_stats = genomecov %>%
                summarise(
                    human.bases     = bases[contig %in% selectContigs("human")] %>% sum(),
                    humanMain.bases = bases[contig %in% selectContigs("human", T)] %>% sum(),
                    mouse.bases     = bases[contig %in% selectContigs("mouse")] %>% sum(),
                    mouseMain.bases = bases[contig %in% selectContigs("mouse", T)] %>% sum()
                ) %>%
                mutate(
                    human.cov     = human.bases / genomeLength("human"),
                    humanMain.cov = humanMain.bases / genomeLength("human", T),
                    mouse.cov     = mouse.bases / genomeLength("mouse"),
                    mouseMain.cov = mouseMain.bases / genomeLength("mouse", T),
                    mouse.contam     = mouse.bases / (human.bases + mouse.bases),
                    mouseMain.contam = mouseMain.bases / (humanMain.bases + mouseMain.bases)
                )

            # Compute percent coverage
            pctCoverage = function(species, dp, main=F) {
                100 * genomecov %>%
                    filter(
                        depth >= !!dp &
                        contig %in% (contigs %>% filter(species == !!species & (primary | !main)) %>% pull(contig))
                    ) %>%
                    pull(count) %>%
                    sum() /
                    contigs %>% filter(species == !!species & (primary | !main)) %>% pull(length) %>% sum()
            }

            pct_coverages = c()

            for (dp in c(1, 3, 5, 10, 30)) {
                out = c(
                    pctCoverage("human", dp),
                    pctCoverage("human", dp, T),
                    pctCoverage("mouse", dp),
                    pctCoverage("mouse", dp, T)
                )
                names(out) = paste0(c("human", "humanMain", "mouse", "mouseMain"), paste0(".pctCov_x", dp))

                pct_coverages = c(pct_coverages, out)
            }

            # Combine and sort results
            table.out = cbind(
                cov_stats,
                t(pct_coverages)
            ) %>%
                mutate_at(vars(grep("\\.", .)), ~sprintf("%.6f", .)) %>%
                select(starts_with("human"), starts_with("mouse")) %>%
                select(-contains("Main"), contains("Main")) %>%
                select(ends_with("bases"), ends_with("contam"), everything()) %>%
                mutate(sample = sample_name) %>%
                select(sample, everything())
        """
    )

    //Output variables
    val estimatedCoveragesOut = NamedMap[Any]("estimatedCoveragesOut")

    //Iterate through each sample
    for (sample <- samples) {
        withName(sample) {
            //Compute coverage histogram for each contig of genome with bedtools genomecov
            //using non-duplicate reads only, i.e. ignoring reads with 1024 SAM flag
            //Use the CPU annotation for 120%-133% CPU usage (and faster processing) when
            //per job CPU usage is restricted by reservation
            //Note that secondary alignments (256) are not removed. This is because these reads
            //are not always filtered by downstream software, e.g. GATK's GetPileUpSummaries
            //(contamination) as well as CollectAllelicCounts and CollectReadCounts (CNV)
            val genomecov = bash"""
                samtools view -F 1024 -b ${markedOut(sample)} |
                    bedtools genomecov -ibam - \\
                    > @out1@
            """
            genomecov._filename("out1", "coverage.table")
            genomecov._custom("cpu") = "2"
            genomecov._custom("memory") = "2048"
            genomecov._execute = "once"

            //Compute mean coverage over human contigs (incl. non-primary ones, excl. viruses)
            //as well as mouse contigs (same thing) with comparisons to estimate degree of
            //contamination
            val estimatedCoverages = REvaluate(
                var1 = genomecov.out1,
                table1 = contigs,
                param1 = sample,
                script = estimatedCoveragesRscript
            )
            estimatedCoveragesOut(sample) = estimatedCoverages.table
        }
    }

    //Join all samples' coverage stats and remove quotes
    val estimatedCoveragesJoin = CSVListJoin(in = estimatedCoveragesOut, fileCol = "")
    val estimatedCoveragesAll = BashEvaluate(
        var1 = estimatedCoveragesJoin.out,
        script = """sed 's/\"//g' @var1@ > @out1@"""
    )
    estimatedCoveragesAll._filename("out1", "out.csv")
}

//A set of Anduril functions (i.e. components) and string commands to run inside bash components
object Tasks {
    //String of commands for creating temporary directory at /tmp/anduril.<PID>
    //the temporary directory will be marked for deleation on exit
    //The path to temp dir is stored in variable TEMP_DIR
    //eval is a workaround for double quatation marks appearing in the final script '"'
    //Use in bash call with 'eval "${Tasks.SET_TEMP_DIR}"'
    val SET_TEMP_DIR = """
        TEMP_DIR="/tmp/anduril.$$";
        mkdir $TEMP_DIR;
        trap "rm -rf $TEMP_DIR" EXIT
    eval"""
}

