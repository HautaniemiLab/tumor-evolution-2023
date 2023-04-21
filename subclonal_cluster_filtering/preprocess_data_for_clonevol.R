# Script to prepare data to ClonEvol from Cyclone outputs
# Will filter clusters / mutations in them
  # Small clusters (parameter per patient)
  # Unreliable fitting to cluster 
  # part-of-trunc cluster(s) with CCF not fitting tree: clusters with >50% CCF within every sample
  # 0-3 random cluster until model is found
  # Lowers CCF for clusters above truncal cluster, will skip if too high (parameter overall)



for (package in c('devtools',
                  'ggplot2',
                  'igraph',
                  'gtools',
                  'colorspace',
                  'FactoMineR',
                  'clonevol',
                  'dplyr')) {
    library(package, character.only=T)
}


#############################################################################################################
##                                            SET PARAMETERS                                               ##
#############################################################################################################
######  DATA  ######
auto_parameters <- TRUE 
parametersDir <- "ADD_DIR"
auto.dir <- TRUE

# Input specifications
auto_ignore_clust0 <- FALSE # In old pyclone cluster 0 was unspecific; TRUE if needs to be deleted
mut_probability_threshold = 0.6 # Cluster assignment probability threshold to exclude badly clustered mutations
correctChildCCFmax = 0.05 # How much variation is allowed in CCF values to correct for cluster with CCF higher than trunc?

# Which runs to make
compute_models = TRUE # If FALSE, will only make initial box plots
clonevol_plot = TRUE

# Plotting parameters
clone.colors <- NULL # If set, clones will be these particular colors in the plots.

# Other parameters
driver.genes.list <- c(
  'BRCA1',
  'BRCA2',
  'NF1',
  'TP53'
)

###### OPTIONAL PARAMETERS WHEN auto_parameters=FALSE ######
#patient     <- "H011_test"  # Output name
#patFileName <- "H011"        # As in source files
#founding_cluster = NA         # Original ID!
#outputDir   <- "ADD_DIR"    # output dir (will create subdir with 'patient' name)
#dataDir     <- "/ADD_DIR"     # input
#ignore_samples <- NA # names of samples to ignore
#min_vaf_now        <- 0.05  # min VAF to use
#min_clust_size <- 10
#min_clusters_ignored = 0 
#max_clusters_ignored = 1
#manually_ignored_clusters <- c(1,2) # (options: list or NULL)
#merge_clusters <- NA # List of elements of 2 clusters to merge. Original IDs! Even when doing only 1 merge, it has to be a list (i.e. list(c(1,2))). It's possible to combine 3 or more i.e. list(c(1,2),c(1,3)) <- #combines 1,2,3

###### READ PARAMETERS FOR WHEN auto_parameters=TRUE ######
# exclude marked patients
if (auto_parameters) {
  parameters_list <- read.csv(parametersDir, header = TRUE, sep = '\t')
  parameters_list <- parameters_list[!(parameters_list$include == 'no'),]
  num.patients = length(parameters_list$patient)
} else {num.patients = 1}

# Basic info per patient
patientInfo <- data.frame(patient=character(),
                          ignores=character(),
                          stringsAsFactors=FALSE)
# If redo run, skip manual ignores since no need to check
redo=TRUE
# read coverage data here if wanted, now excluded everywhere

#############################################################################################################
##                                            GET DATA                                                     ##
#############################################################################################################
# run each patient
for (i in 1:num.patients) {
  # get parameters from file per individual
  if (auto_parameters) { 
    # GET parameters
    patient <- as.character(parameters_list$patient[i])  
    patFileName <- as.character(parameters_list$pat.file[i])
    if (is.na(parameters_list$founding[i])) {
    	founding_cluster <- NA
    } else {
      founding_cluster <- as.numeric(parameters_list$founding[i])
    }
    ignore_samples <- unlist(strsplit(as.character(parameters_list$ignore.samples[i]),','))
    min_vaf_now <- parameters_list$min.vaf[i]
    min_clust_size <- parameters_list$min.clust.size[i]
    min_clusters_ignored <- parameters_list$min.ignored[i]
    max_clusters_ignored <- parameters_list$max.ignored[i]
    manually_ignored_clusters <- as.character(parameters_list$manually.ignored[i])
    merge_clusters <- unlist(strsplit(as.character(parameters_list$merge[i]),';'))
    ccf_error <- parameters_list$ccf.error[i]
    if (auto.dir) {
      dataDir <- as.character(parameters_list$input[i])
      outputDir <- as.character(parameters_list$output[i])
    }
    if (ignore_samples[1] == "-") {ignore_samples <- NA}
    if (manually_ignored_clusters[1] == "-") {
        manually_ignored_clusters <- NA
    } else { manually_ignored_clusters <- unlist(strsplit(as.character(manually_ignored_clusters),',')) }
    if (merge_clusters[1] == "-") {merge_clusters            <- NA}
    if (ccf_error == "-") {ccf_error <- 0.05}
    ccf_error <- ccf_error * 100
    correctChildCCFmax <- ccf_error
  }
  
  # DATA PATHS
  # Input data file
  fileDir     <- paste0(dataDir, patFileName, ".csv")
  # Setup output folder
  subDir <- patient
  dir.create(file.path(outputDir, subDir), showWarnings = FALSE)
  setwd(file.path(outputDir, subDir))

  # SAVE PARAMETERS
  HEADING <- c("founding cluster","ignore samples","minimum cluster size", "manually ignored clusters","merged clusters")
  VALUES <- c(founding_cluster, ignore_samples, min_clust_size, paste(manually_ignored_clusters,collapse=","), paste(merge_clusters,collapse = ';') )
  parameters <- data.frame(HEADING, VALUES)

  # READ INPUT
  init_file <- read.csv(fileDir, header = TRUE, sep="\t", stringsAsFactors = FALSE)
  #set missing std values to 0
  if(sum(is.na(init_file$cellular_prevalence_std)) > 0) {
    init_file[is.na(init_file$cellular_prevalence_std),]$cellular_prevalence_std <- 0
  }
  patientInfo[i,]$patient <-  patient
  
  ###########################################################################################################
  ##                                            FILTER AND PLOT DATA                                       ##
  ###########################################################################################################
  
  ####### READ AND SETUP #######
  # Delete ignored samples
  if (!is.na(ignore_samples)) {
    init_file <- init_file[! init_file$sample_id %in% ignore_samples,]
  }
  # ignore mutations less than selected probability threshold, but never remove more than half!
  if (median(init_file$cluster_assignment_prob) < mut_probability_threshold) {
    init_file <- init_file[init_file$cluster_assignment_prob > median(init_file$cluster_assignment_prob),]
  } else {
    init_file <- init_file[init_file$cluster_assignment_prob > mut_probability_threshold,]
  } 
  # list samples
  sample_list <- unique(init_file$sample_id)
  # Identify founding cluster:
  # select highest freq cluster as clonal cluster if not selected by user
  if (is.na(founding_cluster)) {
    bestClusters <- summarise(summarise(group_by(init_file, cluster_id), mean_CCF=mean(cellular_prevalence)),topcluster=which.max(mean_CCF))
    founding_cluster <- unique(init_file$cluster_id)[bestClusters$topcluster[1]]
	}
  # Reformat CCFs: every sample to own column
  x = init_file[init_file$sample_id == sample_list[1],c("mutation_id","cluster_id")]
  colnames(x) <- c("mutation_id","cluster")
  for (j in 1:length(sample_list)) {
    x[sample_list[[j]]] <- init_file$cellular_prevalence[init_file$sample_id == sample_list[j]] * 50
  }
  # save original cellular prevalences if needed
  # add distribution to variants based on SD, mean and error (SD not read separately)
  for (j in unique(x$cluster)) {
    for (samplecol in unique(sample_list)) {
      # if detected, add variation from normal distribution + selected error
      if(mean(x[x$cluster==j,samplecol]) >= min_vaf_now*100 ) { 
        # error can't be larger than frequency:
        used_error <- min(ccf_error,mean(x[x$cluster==j,samplecol]))
        # add variation
        x[x$cluster==j,samplecol] = rnorm(nrow(x[x$cluster==j,]), 
                                          mean = init_file[init_file$sample_id == samplecol & init_file$cluster_id == j,]$cellular_prevalence[1], 
                                          sd = init_file[init_file$sample_id == samplecol & init_file$cluster_id == j,]$cellular_prevalence_std[1] + used_error/100)
        # transform back to VAF as clonevol wants
        x[x$cluster==j,samplecol] = x[x$cluster==j,samplecol] *50
        if (min(x[x$cluster==j,samplecol]) < 0 ) {
          for (k in 1:nrow(x[x$cluster==j,])) {
            if (x[x$cluster==j,samplecol][k] < 0) { x[x$cluster==j,samplecol][k] = 0 }
          }
        }  
      }
      # if too large values, lower because clonevol does not allow them
      if (mean(x[x$cluster==j,samplecol]) > 50) {
        this_difference <- mean(x[x$cluster==j,samplecol]) -50
        x[x$cluster==j,samplecol] = x[x$cluster==j,samplecol] - this_difference -0.1
      }
  } }
  x <- x[order(x$cluster),]

  ####### PLOT ORIGINAL CLUSTERS #######
  # PLOT data with original cluster IDs and without mutation count threshold
  if (clonevol_plot) {
    pdf(paste(patient, 'box.pdf', sep = '_'), useDingbats = FALSE)
    pp <- plot.variant.clusters(x,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = sample_list,
                            vaf.limits = 50,
                            sample.title.size = 6,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.color = clone.colors,
                            jitter.alpha = 1,
                            jitter.center.method = 'median',
                            jitter.center.size = 1,
                            jitter.center.color = 'darkgray',
                            jitter.center.display.value = 'none',
                            highlight = 'is.driver',
                            highlight.shape = 21,
                            highlight.color = 'blue',
                            highlight.fill.color = 'green',
                            highlight.note.col.name = 'gene',
                            highlight.note.size = 1,
                            base_size = 8,
                            plot.margin = 1,
                            vaf.suffix = '',
                            jitter.size = 0.25,
                            jitter.shape = 'circle',
                            order.by.total.vaf = TRUE)
    dev.off()
    rm(pp)
  }

  ####### FILTER CLUSTERS #######
  # Exclude clusters by min number of variants
  filtered_clusters <- vector()
  mutations_filtered = 0
  for (j in unique(x$cluster)) {
    if ( (length(which(x$cluster == j))) <= min_clust_size) {
      filtered_clusters <- c(filtered_clusters, j)
      mutations_filtered <- mutations_filtered + length(which(x$cluster == j))
      x <- x[!(x$cluster == j),]
    } }
  
  # if redo, skip exclusions completely (not within clonevol run)
  if(redo) {
    x <- x[!(x$cluster %in% manually_ignored_clusters),] 
    manually_ignored_clusters <- NA
  }
  
  #exclude 0-cluster if needed
  if (auto_ignore_clust0) { 
    x <- x[!(x$cluster == 0),]
  }
  
  # Rename clusters: numbers from 1 to N where 1 is for founding cluster
  x$cluster_new = NA
  newcode = 2
  for (j in unique(x$cluster)) {
    if(j == founding_cluster) {
      x[x$cluster == j,]$cluster_new = 1
    } else {
      x[x$cluster == j,]$cluster_new = newcode
      newcode = newcode +1
  } }

  # ADD DRIVER gene info
  x$gene <- NA
  for (j in 1:nrow(x)) {
    x[j,]$gene <- unlist(strsplit(as.character(x[j,]$mutation_id),':'))[3]
  }
  x$is.driver <- FALSE
  if(sum(x$gene %in% driver.genes.list) > 0) {
    x[x$gene %in% driver.genes.list,]$is.driver <- TRUE
  }
  
  # get new codes for manual ignores
  clone_codes = unique(x[,c("cluster","cluster_new")])
  manually.ignored.clusters = manually_ignored_clusters
  if(!is.na(manually_ignored_clusters)) {
    manually.ignored.clusters = clone_codes[clone_codes$cluster %in% manually_ignored_clusters,]$cluster_new
  }
  # clone sizes
  clone_codes$size <- 0
  clone_codes$status <- "OK"
  for (thisclone in clone_codes$cluster_new) {
    clone_codes[clone_codes$cluster_new == thisclone,]$size <- nrow(x[x$cluster_new == thisclone,])
    if(thisclone == 1) { clone_codes[clone_codes$cluster_new == thisclone,]$status = "trunc"}
  }
  
  # save cluster freqs (for SubMARine (https://github.com/morrislab/submarine))
  clonefreqs <- clone_codes
  #clonestds <- clone_codes
  for (sample in sample_list) { 
    clonefreqs[,sample] <- 0
    #clonestds[,sample] <- 0
    for (clone in unique(clonefreqs$cluster_new)) {
      clonefreqs[clonefreqs$cluster_new == clone,sample] <- mean(x[x$cluster_new == clone, sample])
  } }
  
  # 1. Check if all freqs 50% or more compared to trunc
  for (clone in clonefreqs[clonefreqs$status == "OK",]$cluster_new) {
    all_over_half = TRUE
    for (sample in sample_list) { 
      if (clonefreqs[clonefreqs$cluster_new == clone,sample] < clonefreqs[clonefreqs$status == "trunc",sample] / 2 ) { all_over_half = FALSE }
    }
    if (all_over_half == TRUE) { clonefreqs[clonefreqs$cluster_new == clone,]$status = "possibly_part_of_trunc" }
  }
  
  #clonefreqs <- clonefreqs[,c("cluster_new",sample_list)]
  colnames(clonefreqs)[1] <- "cluster_orig"
  write.table(clonefreqs, paste0(patient, "_cluster_info.csv"), sep = "\t", row.names = FALSE)

  ## Modify cluster CCFs if within range (based on parameter correctChildCCFmax) and larger than trunc
  for (sample in sample_list) { 
    for (clone in unique(clonefreqs$cluster_new)) {
      this_difference = clonefreqs[clonefreqs$cluster_new == clone,sample] - clonefreqs[clonefreqs$cluster_new == 1,sample]
      if (clone != 1 & this_difference >= 0 & this_difference < correctChildCCFmax) {
        ### if difference is within accepted range, lower child cluster's CCF trunc_CCF - 0.1%
        x[x$cluster_new == clone,sample] <- x[x$cluster_new == clone,sample] - this_difference - 0.1
      } else if (clone != 1 & this_difference > 0 & this_difference > correctChildCCFmax) {
        ### ignore cluster as there is no place in the model
        manually.ignored.clusters = unique(c(clone,manually.ignored.clusters))
      }
    }
  }
  

  ####### PLOT FILTERED CLUSTERS #######
  # Box plot again but with filtered clusters and driver genes
  if (clonevol_plot) {
    pdf(paste(patient, 'box_filtered_drivers.pdf', sep = '_'), useDingbats = FALSE)
    pp <- plot.variant.clusters(x,
                            cluster.col.name = 'cluster_new',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = sample_list,
                            vaf.limits = 50,
                            sample.title.size = 6,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.color = clone.colors,
                            jitter.alpha = 1,
                            jitter.center.method = 'median',
                            jitter.center.size = 1,
                            jitter.center.color = 'darkgray',
                            jitter.center.display.value = 'none',
                            highlight = 'is.driver',
                            highlight.shape = 21,
                            highlight.color = 'blue',
                            highlight.fill.color = 'green',
                            highlight.note.col.name = 'gene',
                            highlight.note.size = 2,
                            base_size = 8,
                            plot.margin = 1,
                            vaf.suffix = '',
                            jitter.size = 0.25,
                            jitter.shape = 'circle',
                            order.by.total.vaf = TRUE)
    dev.off()
    rm(pp)
  }

  x$cluster_old = x$cluster
  x$cluster = x$cluster_new
  x <- x[order(x$cluster),]

  # Save formatted and filtered file
  write.table(x, paste(patient, "formatted_filtered.csv", sep = '_'), sep = "\t", row.names = FALSE)
  patientInfo[i,]$ignores <- manually.ignored.clusters
  }
setwd(file.path("ADD_DIR"))
write.table(patientInfo, "patient_ignores_for_clonevol.csv", sep = "\t", row.names = FALSE)
