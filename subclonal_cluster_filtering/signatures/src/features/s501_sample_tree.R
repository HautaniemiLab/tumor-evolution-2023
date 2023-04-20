#!/usr/bin/env Rscript

#R script for crudely clustering variants by sample to sample trees using
#agglomerative hierarchical clustering (bottom-up) with complete link. Outputs
#as RData for downstream plotting and merging with other forms of variant
#features such as mutation types.

#Variants are first filtered according to read counts:
# 1. Minimum depth in all samples, those failing go to cluster '-2'
# 2. Minimum ALT read count in at least one sample, else in cluster '-1'
# 3. Minimum ALT reads for a sample to have a variant
#There are options to control the thresholds of these filters. Normal samples
#are completely omitted if matching regex is given and found.

#Only variants passing 1. and 2. are used in the agglomerative clustering of
#samples where the similarity is the number of shared mutations and
#dissimilarity is the difference of similarity from the maximum of the
#similarity matrix.

#The clustered samples are used to create a tree with branching order; branches
#are assigned to all variants. Variants not fitting any branch or clusters '-2'
#or '-1' are assigned to cluster '0'.

#There is an option to split the tree trunk into sections by median VAF via
#given number of cutoff points (uniform cutoffs) or specific cutoffs

#Further refactoring possible, i.e. wrapping steps inside functions

# --------------------------------------------------
# Packages
# --------------------------------------------------

#Hash
suppressMessages(library(hash))

# --------------------------------------------------
# Command line arguments
# --------------------------------------------------

#Options include: variants table with VAF and output RData path.
#Possible options for
# - normal sample regex match-string
# - patient minimum ALT read count in any sample
# - sample minimum ALT read count for single sample
# - patient minimum depth in all samples
# - number of cutoffs or cutoff points for splitting the root
    #Default options
    regex_normals = ""
    min_rc_patient = 2
    min_rc_sample = 1
    min_depth = 10
    root_cutoffs = -1

    args = commandArgs(trailingOnly=TRUE)

    i = 1
    while (i <= length(args)) {
        #Look for options
        if (args[i] %in% c("--input", "-i")) {
            #input variants table
            i = i + 1
            in_file = args[i]
        } else if (args[i] %in% c("--output", "-o")) {
            #output table path
            i = i + 1
            out_file = args[i]
        } else if (args[i] %in% c("--regex-normals", "-r")) {
            #match-regex string for capturing normal samples
            i = i + 1
            regex_normals = args[i]
            if (is.na(regex_normals)) regex_normals = ""
        } else if (args[i] %in% c("--minimum-readcount-patient", "-min-rc-patient")) {
            #input single sample table
            i = i + 1
            min_rc_patient = as.numeric(args[i])
        } else if (args[i] %in% c("--minimum-readcount-sample", "-min-rc-sample")) {
            #type (SBS, DBS or ID)
            i = i + 1
            min_rc_sample = as.numeric(args[i])
        } else if (args[i] %in% c("--minimum-depth", "-min-depth")) {
            #output directory
            i = i + 1
            min_depth = as.numeric(args[i])
        } else if (args[i] %in% c("--root-cutoffs", "-cuts")) {
            #output directory
            i = i + 1
            root_cutoffs = args[i]
            #Handle comma separated quantiles or integer number of cuts
            if (grepl(",", root_cutoffs)) {
                root_cutoffs = as.numeric(strsplit(root_cutoffs, ",")[[1]])
                root_cutoff_specific = T
            } else if (grepl(".", root_cutoffs)) {
                root_cutoffs = as.numeric(root_cutoffs)
                root_cutoff_specific = T
            } else {
                root_cutoffs = as.integer(root_cutoffs)
                root_cutoff_specific = F
            }
        }

        i = i + 1
    }

    #Ensure that required arguments are defined
    if (!exists("in_file") | !exists("out_file")) {
        stop(
            "Please supply input variants table and output RData path",
            call. = F
        )
    }

# --------------------------------------------------
# Read in data
# --------------------------------------------------

#Read variants
    variants = read.table(in_file, stringsAsFactors=FALSE, header=T, check.names=F)

#Remove normals
    variants = variants[
        ,
        grep(
            regex_normals, colnames(variants),
            invert = nchar(regex_normals) > 0
        ),
        drop = F
    ]

# --------------------------------------------------
# Filtering variants for clustering
# --------------------------------------------------

#Helper function for transforming SAC vector to sample-wise columns of AD read counts
    sacsToADs = function(sacs) {
        #Set NAs as zeroes
        sacs_fixed = sacs
        sacs_fixed[is.na(sacs)] = "0,0,0,0"

        #Conver to AD table
        matrix(c(1,0,1,0,0,1,0,1), nrow=2) %*% matrix(
            as.numeric(unlist(strsplit(as.character(sacs_fixed), split=','))),
            nrow = 4
        )
    }

#Function for creating variant by sample binary table
    rcBinaryTable = function(vars_tbl, min_rc_patient=1, min_rc_sample=1) {
        matrix(
            apply(
                vars_tbl[, endsWith(colnames(vars_tbl), ".SAC"), drop=F],
                1,
                function(sacs) {
                    #Get read counts
                    rcs = sacsToADs(sacs)
                    #Check that there is sufficient alt read count for at least one
                    #sample (high confidence variant in patient), then for each sample
                    #sufficient alt read count means the sample has the variant
                    any(rcs[2,] >= min_rc_patient) * (rcs[2,] >= min_rc_sample)
                }
            ),
            ncol = sum(endsWith(colnames(variants), ".SAC")),
            dimnames = list(NULL, grep("\\.SAC", colnames(variants), value=T)),
            byrow = T
        )
    }

#Function for
    depthBinaryVector = function(vars_tbl, min_depth=1) {
        apply(
            vars_tbl[, endsWith(colnames(vars_tbl), ".SAC"), drop=F],
            1,
            function(sacs) {
                #Check sufficient coverage of variant in all samples
                all(colSums(sacsToADs(sacs)) >= min_depth)
            }
        )
    }

#Binary table of samples that have variant based on ALT read counts
    vars_tbl_rc_binary = rcBinaryTable(variants, min_rc_patient, min_rc_sample)

#Vector of variants with sufficient depth to be used in clustering
    vars_vec_dp_binary = depthBinaryVector(variants, min_depth)

#Variants with sufficient coverage to be used in clustering
    vars_valid = vars_tbl_rc_binary[vars_vec_dp_binary,, drop=F]

# --------------------------------------------------
# Sample tree clustering
# --------------------------------------------------

#Clustering
    #Similarity matrix; number of variants in common
        #Initialise matrix
        sim_mat = matrix(
            rep(0, (ncol(vars_valid))^2),
            ncol = ncol(vars_valid),
            dimnames = list(colnames(vars_valid), colnames(vars_valid))
        )

        #Compute number of common variants for each pair of samples
        for (i in 1:max(ncol(sim_mat)-1, 1)) {
            for (j in min(i+1, ncol(sim_mat)):ncol(sim_mat)) {
                sim_mat[i,j] = sum(vars_valid[,i] + vars_valid[,j] == 2)
                sim_mat[j,i] = sim_mat[i,j]
            }
        }

    #Check if number of samples is at least 2
    if (nrow(sim_mat) >= 2) {
        #Create distance matrix as difference from maximum
        dissim_mat = max(sim_mat) - sim_mat
        diag(dissim_mat) = 0
        dist_m = as.dist(dissim_mat)

        #Fitting clusters
        fit = hclust(dist_m, method="complete")
    } else {
        #Make fake hclust result with only a single cluster
        fit = list()
        fit$labels = rownames(sim_mat)

        #Replace cutree function to give specific result
        cutree = function(not_used_fit, not_used_k) {
            out = c(1)
            names(out) = fit$labels
            out
        }
    }

#Output table
    variants_branches = cbind(
        variants[, c("CHROM", "POS", "REF", "ALT")],
        "branch" = 0
    )

#Branching - assign variants to branches
    #Variables needed
        curr_i = 4
        already_out = c()   #samples already out
        samples = fit$labels

        #List of branches, each with samples/variants involved; keep last unique branch
        branches = list(-2, -1, 0)
        names(branches[[1]]) = "Low Depth"
        names(branches[[2]]) = "Low Confidence"
        names(branches[[3]]) = "Unknown"

        #Hash for tracking branches to avoid repeats
        h = hash()
        h[["#-2"]] = -2
        h[["#-1"]] = -1
        h[["#0"]] = 0

    #Variants considered low confidence (from insufficient alt read) as -1
    variants_branches[rowSums(vars_tbl_rc_binary) == 0, "branch"] = -1

    #Variants with insufficient depth as -2
    variants_branches[!vars_vec_dp_binary, "branch"] = -2

    #Loop through all tree cuts
    for (k in 1:length(samples)) {
        cutree_indices = cutree(fit, k)
        #Filter out samples already branched out
        if (length(already_out) > 0) {
            samples_remaining = cutree_indices[-already_out]
        } else {
            samples_remaining = cutree_indices
        }

        #Loop over cutree indices of remaining samples
        for (i in max(samples_remaining):1) {
            n_remaining_with_index = length(which(samples_remaining == i))
            #Check if new branch
            if (n_remaining_with_index >= 1) {
                #The samples involved
                samples_involved = which(cutree_indices == i)
                key = paste0("#", samples_involved, collapse="")

                #Do nothing if already in hash
                if (is.null(h[[key]])) {
                    h[[key]] = curr_i
                    #Check if single sample branching out
                    if (n_remaining_with_index == 1) {
                        already_out = c(already_out, samples_involved)
                    }
                    #Assign correct variants (1 rc in involved samples, 0 otherwise) to the branch
                    variants_branches[
                        which(
                            #Correct number of samples (any) have the variant
                            rowSums(vars_tbl_rc_binary) == length(samples_involved) &
                            #Selected samples (in branch) all have the variant
                            rowSums(vars_tbl_rc_binary[,samples_involved, drop=F]) == length(samples_involved) &
                            #Variants that have passed coverage threshold in all samples
                            vars_vec_dp_binary
                        ),
                        "branch"
                    ] = curr_i - 3

                    branches[[curr_i]] = samples_involved
                    curr_i = curr_i + 1
                }
            }
        }
    }

#Helper function that checks for non-contiguous blocks
    treeGetNonContiguous = function(branches, sample_order) {
        #Find nodes with children
        non_leafs = which(sapply(branches, length) >= 2)

        #Find contiguous blocks of samples in non-leaf nodes
        contiguous_blocks = sapply(
            non_leafs,
            function(i) {
                all(sapply(
                    2:length(branches[[i]]),
                    function(j) {
                        sort(sample_order[names(branches[[i]])])[j] - sort(sample_order[names(branches[[i]])])[j-1]
                    }
                ) == 1)
            }
        )

        #Return result
        list("non_leafs" = non_leafs, "contiguous" = contiguous_blocks)
    }

#Helper function that gives order of branches from root to singleton
    treeGetPath = function(branches, sample_id) {
        #Find branch indices the sample is in
        branch_ind = unlist(sapply(
            4:length(branches),
            function(i) {
                if (!is.na(branches[[i]][sample_id])) i
            }
        ))

        #Get the sizes of the branches
        sizes = sapply(
            branch_ind,
            function(i) {
                length(branches[[i]])
            }
        )

        #Return result
        branch_ind[order(sizes, decreasing=T)]
    }

#Reorder branches into a staircase according to already_out order
    if (length(samples) >= 2) {
        #Reorder samples - make sure vertical blocks are contiguous
            curr_sample_order = 1:length(samples)
            names(curr_sample_order) = names(already_out)

            #Fix blocks from furthest (usually smallest) blocks
            blocks = treeGetNonContiguous(branches, curr_sample_order)
            last_block = 0

            while(!all(blocks$contiguous)) {
                block_to_fix = which(blocks$contiguous == F)[length(which(blocks$contiguous == F))]
                samples_in_block = branches[[blocks$non_leafs[block_to_fix]]]
                sib_sorted = sort(curr_sample_order[names(samples_in_block)])
                sib_move_dist = -sib_sorted + c(sib_sorted[-1], sib_sorted[length(sib_sorted)]+1) - 1

                #Moving samples
                cum_move_dist = 0
                for (i in rev(which(sib_move_dist > 0))) {
                    #Update cum_move_dist
                    cum_move_dist = cum_move_dist + sib_move_dist[i]
                    #Reordering sample order
                    curr_sample_order[sib_sorted[i] : (sib_sorted[i]+cum_move_dist)] = c(
                        sib_sorted[i]+cum_move_dist,
                        (sib_sorted[i]+1) : (sib_sorted[i]+cum_move_dist) -1
                    )
                }

                #Update order vector
                curr_sample_order = sort(curr_sample_order)

                #Update blocks
                blocks = treeGetNonContiguous(branches, curr_sample_order)
                blocks
            }

            #Finalise sample order
            sample_order = already_out[names(sort(curr_sample_order))]

        #Reorder branches
            #Get new order of branches
            new_branch_order = c()
            h_branches = hash()

            for (sample_id in names(sample_order)) {
                #Get sample's tree path
                sample_path = treeGetPath(branches, sample_id)
                #Add new unique indices to the hash and the new order
                for (i in sample_path) {
                    if (is.null(h_branches[[as.character(i)]])) {
                        h_branches[[as.character(i)]] = 1
                        new_branch_order = c(new_branch_order, i)
                    }
                }
            }

            #Reorganise variants_branches and branches
            new_branches = branches

            for (i in 4:length(branches)) {
                #Check if moved
                if (i != new_branch_order[i-3]) {
                    #New branches
                    new_branches[[i]] = branches[[new_branch_order[i-3]]]
                }
                #New variant branch assignments with '#' character
                variants_branches[
                    which(variants_branches[,"branch"] == new_branch_order[i-3] - 3),
                    "branch"
                ] = paste0("#", i-3)
            }

            #Assign new values
            branches = new_branches
            variants_branches[, "branch"] = as.numeric(sub("#", "", variants_branches[, "branch"]))
    } else {
        sample_order = 1
        names(sample_order) = samples
    }

# --------------------------------------------------
# Limiting to root and cutting it
# --------------------------------------------------

#Function for creating median VAFs for each variant
    medianVAFs = function(vars_tbl) {
        apply(
            vars_tbl[, endsWith(colnames(vars_tbl), ".SAC"), drop=F],
            1,
            function(sacs) {
                #Get read counts
                rcs = sacsToADs(sacs)
                #Compute median VAF
                median(rcs[2,] / colSums(rcs), na.rm=T)
            }
        )
    }

#Limit to root and split if at least one cutoff specified
    if (root_cutoffs[1] >= 0) {
        #Keep root only
        root_ind = which(variants_branches$branch == 1)
        other_branches_ind = which(variants_branches$branch > 1)

        #Create cutoffs for root as new branch names
        vafs = medianVAFs(variants[root_ind,, drop=F])

        if (root_cutoff_specific) {
            cutoff_count = length(root_cutoffs)
            cutoff_quantiles = quantile(vafs, unique( sort(c(0, 1 - root_cutoffs, 1)) ))
        } else {
            cutoff_count = root_cutoffs
            cutoff_quantiles = quantile(vafs, seq(0, 1, length.out=2 + root_cutoffs))
        }

        variant_groups = cut(
            vafs, cutoff_quantiles,
            labels = F, include.lowest = T
        )
        variant_groups = max(variant_groups) + 1 - variant_groups

        #Apply new branches
        variants_branches$branch[other_branches_ind] = variants_branches$branch[other_branches_ind] + cutoff_count
        variants_branches$branch[root_ind] = variant_groups
        branches = branches[c(
            1, 2, 3,
            rep(4, 1 + cutoff_count),
            if (length(branches) > 4) {5:length(branches)}
        )]
    }

# --------------------------------------------------
# Compute variant order by branches
# --------------------------------------------------

    #Order variants by branches
    vars_order = order(variants_branches[,"branch"])

    #Order low coverage and ambiguous variants (branches -2 and 0)
    if (length(samples) >= 2) {
        for (branch in c(-2, 0)) {
            variants_to_order = which(variants_branches[vars_order, "branch"] == branch)

            #Operation required if there are at least 2 ambiguous variants
            if (length(variants_to_order) >= 2) {
                #Compute order within branch
                new_order = do.call(
                    order,
                    c(
                        decreasing = T,
                        data.frame(
                            vars_tbl_rc_binary[
                                variants_branches[,"branch"] == branch,, drop=F
                            ][,sample_order, drop=F]
                        )
                    )
                )
                #Update order
                vars_order[variants_to_order] = vars_order[variants_to_order][new_order]
            }
        }
    }

# --------------------------------------------------
# Write results
# --------------------------------------------------

#Save relevant objects as RData, i.e. data required to join other variant features
#and draw sample trees
save(
    variants, variants_branches,
    samples, branches,
    sample_order, vars_order,
    min_rc_patient, min_rc_sample,
    sacsToADs, treeGetPath,
    file = out_file
)

# --------------------------------------------------
# End of file
# --------------------------------------------------



