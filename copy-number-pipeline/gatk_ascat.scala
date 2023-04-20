#!/usr/bin/env anduril
//$OPT --wrapper slurm-prefix

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._
import scala.collection.mutable.Map

// Author: Kari Lavikka
// Partially based on work by: Juho and Amjad

// CNV analysis for WGS data after alignment (bam files) and germline variant calling
// Based roughly on GATK4 best practices:
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531092--How-to-part-I-Sensitively-detect-copy-ratio-alterations-and-allelic-segments


// TODO: Check if there are patients having only zero copies on a considerably large fraction of the genome

object GatkAscat {

  val STORAGEBIG_RES = "/XXXX"
  val CSC_RES = "/YYYY"

  val DATAPATH = STORAGEBIG_RES + "/processed_data/HERCULES"
  val PIPELINE_RES = CSC_RES + "/pipeline_resources/ascat"

  val reference = INPUT(CSC_RES + "/references_annotations/reference_genomes/GRCh38d1vd1/GRCh38.d1.vd1.fa")
  val refDict = INPUT(CSC_RES + "/references_annotations/reference_genomes/GRCh38d1vd1/GRCh38.d1.vd1.dict")
  val genes = INPUT(PIPELINE_RES + "/ensembl_genes_v96.csv")

  // List of BAM files that are used for collecting the counts
  val bamsCSV = INPUT(DATAPATH + "/WGSbams/list.csv")

  // Info on sample sequencing platform
  val samples = INPUT(DATAPATH + "/WGSbams/sample_info_extended.csv")

  // Which non-normal samples to include. An empty table means all samples.
  val samplesToInclude = INPUT("samples-to-include.csv")

  // Paths to VCF files that contain the positions of heterozygous SNPs of each patient
  val snpFiles = INPUT(DATAPATH + "/WGS_germline_variants/patient_biallelic_het_snps_vcf_list.csv")

  // TP53 variants for additional purity evidence
  //val tp53File = INPUT("/mnt/storageBig8/web-pub/projects/HERCULES/WGS/variants_allMatched/variants_v4.2/tp53_variants.csv")
  val tp53File = INPUT(DATAPATH + "/tp53_inspection/230123/tp53_filtered.csv")
  // Manually chosen truncal TP53 mutations (from eDuuni)
  val tp53AnnotationFile = INPUT("TP53_annotations_variants4.5.xlsx")

  // Just a list of chromosomes, 1-22,X
  val chromosomes = INPUT(PIPELINE_RES + "/chromosomes.list")

  // Encode blacklist: https://www.nature.com/articles/s41598-019-45839-z
  val blacklist = INPUT(PIPELINE_RES + "/hg38-blacklist.v2.bed")

  // Regions of repetitive DNA or germ-line CNV
  // Based on an earlier segmentation of normal samples
  val cnvBlacklist = INPUT(PIPELINE_RES + "/hercules-cnv-blacklist-v1.bed")

  val excludeFromPon = INPUT("exclude-from-pon.csv")

  // Sets a proper JAVA_HOME and PATH
  val gatk = INPUT("runGatk.sh")

  // Manual curation for select samples
  val ascatParameterOverrides = INPUT("ascat-parameter-overrides.xlsx")

  // Scripts
  // script to get the gene-specific CN st
  val cnGeneScript = INPUT("cn_genes.R")
  val ascatOnGatkScript = INPUT("ascat_on_gatk.R")
  val multiPloidyPatientsScript = INPUT("multi-ploidy-patients.R")
  val filterReadCountsScript = INPUT("filter-read-counts.py")

  // Prepare the intervals for copy-ratio extraction //////////////////////////

  // Generate bins that are used for collecting read counts
  // TODO: Document why OVERLAPPING_ONLY
  val intervalList = bash"""${gatk} PreprocessIntervals \\
      -L ${chromosomes} \\
      -R ${reference} \\
      --sequence-dictionary ${refDict} \\
      -XL ${blacklist} \\
      -XL ${cnvBlacklist} \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      -O @out1@
    """
   intervalList._filename("out1", "preprocessed_intervals.interval_list")

  // Annotate the intervals with GC content, compute it from the reference genome
  val intervalGCAnnot = bash"""${gatk} AnnotateIntervals \\
      -R ${reference} \\
      -L ${intervalList.out1} \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      -O @out1@
      """
  intervalGCAnnot._filename("out1", "intervals_annotated.tsv")

  // Collect read and allelic counts from each sample /////////////////////////

  // Join samples with files that contain patients' germline heterozygous loci
  val samplesWithSNPPos = REvaluate(
    table1 = samples,
    table2 = snpFiles,
    script = """
      library(magrittr)
      library(dplyr)
      
      table.out <- table1 %>%
        inner_join(table2 %>% rename(patient = Key, SNPFile = File))
    """)

  val readCounts = Map[String, TextFile]()
  val allelicCounts = Map[String, TextFile]()

  for (rowMap <- iterCSV(samplesWithSNPPos.table)) {
    val sample = rowMap("sample")
    val bamFile = rowMap("bamFile")

    // Heterozygous loci of the current patient
    val snpFile = rowMap("SNPFile")

    withName(sample) {
      val collectReadCounts = bash"""${gatk} --java-options "-Xmx2560m" CollectReadCounts \\
          -I ${bamFile} \\
          -R ${reference} \\
          -L ${intervalList.out1} \\
          --interval-merging-rule OVERLAPPING_ONLY \\
          --format TSV \\
          -O @out1@
        """
      collectReadCounts._filename("out1", "readCounts_" + sample + ".tsv")
      collectReadCounts._custom("memory") = "2560" 
      readCounts(sample) = collectReadCounts.out1

      // TODO: Why base-quality 20? Seems to be the default - can be removed
      val collectAllelicCounts = bash"""${gatk} --java-options "-Xmx2560m" CollectAllelicCounts \\
          -I ${bamFile} \\
          -R ${reference} \\
          -L ${snpFile} \\
          --minimum-base-quality 20 \\
          -O @out1@
        """
      collectAllelicCounts._filename("out1", "allelicCounts_" + sample + ".tsv")
      collectAllelicCounts._custom("memory") = "2560"
      allelicCounts(sample) = collectAllelicCounts.out1
    }
  }

  // Extends the sample info with sample-specific read and allelic count files.
  // Adds a column that contains the patient-specific NORMAL allelic counts.
  val samplesWithCounts = REvaluate(
    table1 = samples,
    table2 = Array2CSV(readCounts, _name = "readCountsCSV"),
    table3 = Array2CSV(allelicCounts, _name = "allelicCountsCSV"),
    table4 = excludeFromPon,
    script = """
      library(magrittr)
      library(dplyr)
      library(tidyr)

      samples <- table1 %>%
        filter(usable)

      readCounts <- table2 %>%
        select(sample = Key, readCounts = File)

      allelicCounts <- table3 %>%
        select(sample = Key, allelicCounts = File)

      excludeFromPon <- table4

      samplesWithCounts <- samples %>%
        inner_join(readCounts) %>%
        inner_join(allelicCounts)

      # Find any normals for the samples and use them if no same-platform normal can be found
      # Commented out because all samples should now have a same-platform normal
      #fallbackNormals <- samples %>%
      #  inner_join(allelicCounts) %>%
      #  filter(normal) %>%
      #  group_by(patient) %>%
      #  slice_head(n = 1) %>%
      #  select(patient, fallbackNormalAllelicCounts = allelicCounts)

      table.out <- samplesWithCounts %>%
        inner_join(allelicCounts %>%
                    rename(normalSample = sample,
                            normalAllelicCounts = allelicCounts)) %>%
        mutate(pon = usable & normal & !(sample %in% excludeFromPon$sample))
      """)

  

  // Panels of Normals ////////////////////////////////////////////////////////

  // Creates a panel of normals for the given platform
  def createPanelOfNormals(platform: String): BinaryFile = {
      // Concatenate all normal readCount filenames into command line options
      val inputParams = iterCSV(samplesWithCounts.table)
        .filter(_("pon") == "TRUE")
        .filter(_("platform") == platform)
        .map("-I " + _("readCounts"))
        .mkString(" ")

      // TODO: CHECK MINIMUM MEDIAN PERCENTILE VARIABLE --> what is recommended for WGS/our data?
      val ponOut = BashEvaluate(
        var1 = gatk,
        var2 = intervalGCAnnot.out1,
        script = s"""@var1@ --java-options "-Xmx16G" CreateReadCountPanelOfNormals \\
          $inputParams \\
          --annotated-intervals @var2@ \\
          -O @out1@
        """)
    ponOut._filename("out1", "pon.hdf5") // TODO: add platform name to the filename
    ponOut._custom("memory") = "16384"

    return ponOut.out1
  }


  // Panels of normals for each platform
  val panelOfNormals = Map[String, BinaryFile]()

  // Extract all distinct sequencing platforms from the sample list
  val platforms = iterCSV(samples).map(_("platform")).toList.distinct

  // ... and create the PoNs for them
  for (platform <- platforms) {
    panelOfNormals(platform) = withName(platform) {
      createPanelOfNormals(platform)
    }
  }

  //var samplesToProcess = samplesWithCounts

  var samplesToProcess = REvaluate(
    table1 = samplesWithCounts.table,
    table2 = samplesToInclude,
    script = """
      library(magrittr)
      library(dplyr)

      samples <- table1
      samplesToInclude <- table2

      if (nrow(samplesToInclude) > 0) {
        table.out = samples %>%
          inner_join(samplesToInclude)
      } else {
        table.out = samples
      }
      """)


  // Do segmentation for all samples //////////////////////////////////////////

  val modelFilesMap = Map[String, FolderExtractor]()
  val gatkPlotsMap = Map[String, BinaryFile]()

  // Skip normals, TODO: process normals too
  for (rowMap <- iterCSV(samplesToProcess.table)) {

    val sample = rowMap("sample")
    val tumor = rowMap("normal") != "TRUE"

    // Note: We are referring raw file paths here, not ports
    // However, the dependency management works anyway, because they have gone through Array2CSV
    val readCounts = rowMap("readCounts")
    val allelicCounts = rowMap("allelicCounts")
    val normalAllelicCounts = rowMap("normalAllelicCounts")

    // TODO: Figure out which one is proper comparison: null or ""
    val useMatchingNormal = tumor && (normalAllelicCounts != null && normalAllelicCounts != "")

    withName(sample) {
      val denoiseCounts = bash"""${gatk} --java-options "-Xmx8g" DenoiseReadCounts \\
          -I ${readCounts} \\
          --count-panel-of-normals ${panelOfNormals(rowMap("platform"))} \\
          --standardized-copy-ratios @out1@ \\
          --denoised-copy-ratios @out2@
        """
      denoiseCounts._filename("out1", sample + ".standardizedCR.tsv") // Needed for something?
      denoiseCounts._filename("out2", sample + ".denoisedCR.tsv")
      denoiseCounts._custom("memory") = "8192"

      // Because our average coverage is 30, the default minimum-total-allele-count-normal
      // of 30 would drop a half of our snips. Thus, we reduce it to an arbitrary value of 20.
      // TODO: Come up with a better justified value
      // 
      // Using "eval" to get rid of quotes that string interpolation of bash""" adds
      val segmentModel = bash"""eval ${gatk} --java-options "-Xmx8g" ModelSegments \\
        --denoised-copy-ratios ${denoiseCounts.out1} \\
        --allelic-counts ${allelicCounts} \\
        ${if (useMatchingNormal) "--normal-allelic-counts " + normalAllelicCounts else ""} \\
        --output-prefix ${sample} \\
        --minimum-total-allele-count-normal 20 \\
        \\
        --number-of-changepoints-penalty-factor 1 \\
        --kernel-variance-allele-fraction 0 \\
        --kernel-variance-copy-ratio 0.2 \\
        --kernel-scaling-allele-fraction 0.1 \\
        --smoothing-credible-interval-threshold-allele-fraction 2 \\
        --smoothing-credible-interval-threshold-copy-ratio 10 \\
        \\
        --output @folder1@
      """
      segmentModel._custom("memory") = "8192"
      segmentModel._filename("folder1", "model")

      // modelFinal.seg contains the actual segmentation
      // hets.tsv contains the heterozygous loci that were used for computing BAF
      val modelFiles = FolderExtractor(
        folder = segmentModel.folder1,
        filename1 = sample + ".modelFinal.seg",
        filename2 = sample + ".hets.tsv")
      modelFiles._keep = true

      modelFilesMap(sample) = modelFiles


      // Create plots with GATK, gives an overview of the segmentation result

      if (false) {
        // Crashes to some R problem. Disabled for now...
        val segmentPlot = bash"""${gatk} --java-options "-Xmx2g" PlotModeledSegments \\
          --denoised-copy-ratios ${denoiseCounts.out2} \\
          --allelic-counts ${modelFiles.file2} \\
          --segments ${modelFiles.file1} \\
          --sequence-dictionary ${refDict} \\
          --output @folder1@ \\
          --output-prefix ${sample}
        """
        segmentPlot._filename("folder1", "plots")

        val segmentPlotFiles = FolderExtractor(
          folder = segmentPlot.folder1,
          filename1 = sample + ".modeled.png")
        segmentPlotFiles._keep = true
        segmentPlotFiles._filename("file1", sample + ".modeled.png")
        gatkPlotsMap(sample) = segmentPlotFiles.file1
      }
    }
  }

  
  // Run ASCAT for purity/ploidy estimation ///////////////////////////////////

  // Alt and ref read counts for TP53 SNPs. These are used as additional evidence in
  // purity estimation.

  val tp53SNPs = REvaluate(
    table1 = tp53File,
    var1 = tp53AnnotationFile,
    script = """
      library(magrittr)
      library(dplyr)
      library(tidyr)
      library(stringr)
      library(readxl)

      table.out <- read_xlsx(var1) %>%
        filter(truncalDriver == "yes") %>%
        mutate(annotationPos = as.numeric(str_match(Hg38, "^chr17:(\\d+)")[, 2])) %>%
        select(patient, annotationPos) %>%
        # Join to the actual TP53 mutation table
        inner_join(table1) %>%
        mutate(tp53Distance = abs(POS - annotationPos)) %>%
        # The annotation and mutation tables have slightly different positions because
        # of notation differences. Thus, keep only the closest ones
        filter(tp53Distance < 30) %>%
        group_by(patient, annotationPos) %>%
        arrange(tp53Distance) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        select(-annotationPos, -tp53Distance) %>%
        separate_rows(samples, readCounts, sep = ";") %>%
        separate(readCounts, into = c("ref", "alt"), sep = ",", convert = TRUE) %>%
        rename(sample = samples,
              tp53.chrom = CHROM,
              tp53.pos = POS,
              tp53.refCount = ref,
              tp53.altCount = alt) %>%
        group_by(sample) %>%
        # If multiple mutations exist, pick the one with highest alt count
        arrange(-tp53.altCount) %>%
        slice_head(n = 1) %>%
        ungroup()
    """)

  val samplesWithTP53AndStuff = REvaluate(
    table1 = samplesToProcess.table,
    table2 = tp53SNPs.table,
    script = """
      library(magrittr)
      library(dplyr)

      table.out <- table1 %>%
        select(sample) %>%
        left_join(table2)
    """)


  val ascatOnGatkMap = Map[String, BashEvaluate]()
  val topAscatEstimateMap = Map[String, CSV]()
  val segmentsWithAllelesMap = Map[String, CSV]()
  val cnGenesMap = Map[String, CSV]()

  for (rowMap <- iterCSV(samplesWithTP53AndStuff.table)) {
    val sample = rowMap("sample")

    withName(sample) {

      // Feed the GATK segmentation to ASCAT for ploidy and purity estimation
      // If available, use TP53 as additional purity evidence.
      val flags = if (rowMap("tp53.pos") != null && rowMap("tp53.pos") != "NA") {
        "--snp-purity=" + Array("chrom", "pos", "altCount", "refCount").map(x => rowMap("tp53." + x)).mkString(",")
      } else {
        ""
      }

      val segments = bash"""Rscript -e '
        library(magrittr)
        library(dplyr)
        library(readr)

        args <- commandArgs(trailingOnly = TRUE)
        
        read_tsv(args[1], comment = "@", na = c("NA", "NaN")) %>%
          select(chr = CONTIG, startpos = START, endpos = END,
                 nProbesCr = NUM_POINTS_COPY_RATIO, nProbesAf = NUM_POINTS_ALLELE_FRACTION,
                 logR = LOG2_COPY_RATIO_POSTERIOR_50, BAF = MINOR_ALLELE_FRACTION_POSTERIOR_50) %>%
          write_tsv(args[2])
      ' ${modelFilesMap(sample).file1} @out1@"""

      segments._filename("out1", sample + ".gatk.segments.csv")

      val ascatOnGatk = BashEvaluate(
        var1 = ascatOnGatkScript,
        var2 = segments.out1,
        param1 = sample,
        script = s"""@var1@ $flags \\
          --sample=@param1@ \\
          --estimate-output=@out2@ \\
          --sunrise-output=@out3@ \\
          @var2@ \\
        """)

      ascatOnGatk._filename("out2", sample + ".candidates.csv")
      ascatOnGatk._filename("out3", sample + ".sunrise.png")
      ascatOnGatkMap(sample) = ascatOnGatk

      // Picks the best (the highest penalizedGoodnessOfFit) candidate.
      // If manual overrides are available, use them to confine the candidates.
      val topAscatEstimate = bash"""Rscript -e '
        library(magrittr)
        library(dplyr)
        library(readr)
        library(readxl)

        args <- commandArgs(trailingOnly = TRUE)
        sampleId <- args[1]
        ascatParameterOverrides <- read_excel(args[3]) %>%
          transmute(sample = sample,
                    discardModels = !is.na(discardModels),
                    ploidyOverride = ploidy,
                    purityOverride = purity)
        outputFile <- args[5]

        candidates <- read_tsv(args[2]) %>%
          left_join(ascatParameterOverrides)

        candidates <- candidates %>%
          rename(purity = rho)

        # TODO: ascat-algoritm.R could output this. Would simplify the pipeline.
        tp53VAFs <- read_tsv(args[4]) %>%
          mutate(TP53.VAF = tp53.altCount / (tp53.altCount + tp53.refCount)) %>%
          select(sample, TP53.VAF)

        if (!(TRUE %in% candidates$$aberrant)) {
          # Sample is not aberrant
          table.out <- candidates %>%
            slice_head(n = 1) %>%
            mutate(purity = 0,
                   psi = 2,
                   ploidy = 2,
                   TP53.purity.mean = NA,
                   goodnessOfFit = NA,
                   penalizedGoodnessOfFit = NA)
        
        } else if (TRUE %in% candidates$$discardModels) {
          table.out <- candidates %>%
            slice_head(n = 1) %>%
            mutate(purity = NA,
                   psi = NA,
                   ploidy = NA,
                   TP53.purity.mean = NA,
                   goodnessOfFit = NA,
                   penalizedGoodnessOfFit = NA)

        } else if (TRUE %in% is.na(candidates$$ploidy)) {
          # Sample is aberrant but no solution was found. All colums are NA
          table.out <- candidates %>%
            slice_head(n = 1) %>%
            mutate(purity = NA)

        } else {
          table.out <- candidates %>%
            filter(is.na(ploidyOverride) | (ploidy > ploidyOverride - 0.25 & ploidy < ploidyOverride + 0.25)) %>%
            filter(is.na(purityOverride) | (purity > purityOverride - 0.1 & purity < purityOverride + 0.1)) %>%
            slice_max(penalizedGoodnessOfFit, n = 1)

          if (nrow(table.out) == 0) {
            stop("The overridden ploidy did not match any candidate!")
          }
        }

        table.out %>%
          left_join(tp53VAFs) %>%
          select(sample, aberrant, purity, psi, ploidy,
                 TP53.purity.mean, TP53.VAF,
                 goodnessOfFit, penalizedGoodnessOfFit) %>%
          write_tsv(outputFile)
      ' ${sample} ${ascatOnGatk.out2} ${ascatParameterOverrides} ${tp53SNPs.table} @out1@"""

      topAscatEstimateMap(sample) = topAscatEstimate.out1


      // Calculate A and B allele counts using the estimated purity and ploidy
      val segmentsWithAlleles = bash"""Rscript -e '
        library(readr)
        source("../../ascat-algorithm.R")

        args <- commandArgs(trailingOnly = TRUE)
        topEstimate <- read_tsv(args[1])
        segments <- read_tsv(args[2])

        segments %>%
          addAlleleCountsToSegments(topEstimate$$ploidy, topEstimate$$purity) %>%
          mutate(sample = topEstimate$$sample) %>%
          relocate(sample) %>%
          rename(baf = BAF) %>%
          mutate_if(is.numeric, round, digits = 4) %>%
          write_tsv(args[3])
      ' ${topAscatEstimate.out1} ${segments.out1} @out1@"""

      segmentsWithAlleles._filename("out1", sample + ".segments.csv")
  
      segmentsWithAllelesMap(sample) = segmentsWithAlleles.out1


      // Compute per-gene copy number statuses
      val cnGenes = REvaluate(
        table1 = segmentsWithAlleles.out1,
        table2 = topAscatEstimate.out1,
        table3 = genes,
        script = cnGeneScript)
      cnGenes._filename("table", sample + ".cnGenes.csv")
      cnGenesMap(sample) = cnGenes.table

    }
  }

  // All ASCAT solution candidates for all samples
  val combinedAscatCandidates = CSVListJoin(in = ascatOnGatkMap.mapValues(_.out2), fileCol = "")


  // Prepare outputs //////////////////////////////////////////////////////////

  //val gatkPlotsFolder = Array2Folder(in = gatkPlotsMap, fileMode = "@file@")

  val combinedAscatSegments = CSVListJoin(in = segmentsWithAllelesMap, fileCol = "")
  val combinedAscatEstimates = CSVListJoin(in = topAscatEstimateMap, fileCol = "")

  // Calculate some breakpoint statistics
  val ascatStats = REvaluate(
    table1 = combinedAscatSegments,
    script = """
      library(magrittr)
      library(dplyr)
      library(tidyr)

      table.out <- table1 %>%
        group_by(sample, chr) %>%
        mutate(`break` = row_number() != 1) %>%
        mutate(ascatBreak = (nMajor != lag(nMajor) | nMinor != lag(nMinor)) %>%
                 replace_na(TRUE) & `break`) %>%
        ungroup() %>%
        group_by(sample) %>%
        summarise(breaks = sum(`break`),
                  ascatBreaks = sum(ascatBreak)) %>%
        mutate(breakRatio = ascatBreaks / breaks) %>%
        mutate_if(is.numeric, round, digits = 3)
    """)


  def extractFile(folder: BinaryFolder, filename: String): BinaryFile = {
    val extractor = FolderExtractor(folder = folder, filename1 = filename)
    extractor._keep = false
    extractor._filename("file1", filename)
    return extractor.file1
  }

  val ascatSunrises = Array2Folder(in = ascatOnGatkMap.mapValues(_.out3), fileMode = "@file@")

  val cnGenesFolder = Array2Folder(in = cnGenesMap, fileMode = "@file@")

  val combinedCnGenes = CSVListJoin(in = cnGenesMap, fileCol = "")

  // Plot that summarises patients that have multiple different ploidies in their samples
  val multiPloidyPatients = bash"""Rscript ${multiPloidyPatientsScript} ${combinedAscatCandidates} ${ascatParameterOverrides} @out1@ @folder1@"""
  multiPloidyPatients._filename("out1", "patients-with-multiple-ploidies.pdf")
  multiPloidyPatients._filename("folder1", "patients")


  // Do output ////////////////////////////////////////////////////////////////

  // Segmentation results as CSVs
  OUTPUT(combinedAscatSegments)
  OUTPUT(combinedAscatEstimates)
  OUTPUT(ascatStats.table)

  // ASCAT plots
  OUTPUT(ascatSunrises)

  OUTPUT(multiPloidyPatients.out1)
  OUTPUT(multiPloidyPatients.folder1)

  // GATK Segmentation plots
  //OUTPUT(gatkPlotsFolder)

  // Per gene results
  OUTPUT(cnGenesFolder)
  OUTPUT(combinedCnGenes)

  // TODO: combinedCnGenes for "selected genes", BRCA12, CCNE1, etc
  // Note: It's per patient, not per sample.
  // Some pretty plots might be useful: density of absolute copy numbers, etc.
  // https://wiki.helsinki.fi/display/HautaniemiCollaboration/DNA+Sequencing+Selected+Genes

  // TODO: Give proper names for outputs:
  /*
  val result_files=FolderCombiner(
                    files=makeArray( 
                            "imageData.txt" -> perimage,
                            "nucleusData.txt" -> important
                            ),
                    keysAsNames=true)
  */
}
