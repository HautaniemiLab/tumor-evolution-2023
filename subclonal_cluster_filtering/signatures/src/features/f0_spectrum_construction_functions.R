#Collection of functions related to construction of mutational profiles (spectra
#or catalogues)
#List of contents:
# - Variant loading functions:
#   - readVarFile(file_names)
# - Spectra building functions:
#   - varVafTbl(vars_tbl)
#   - varBinaryVec(vars_tbl, ...)
#   - getMutTypes(variants, ...)
#   - buildCatalogue(variants, ref_sigs, ...)
#   - buildCataloguePatientSamples(vars_tbl, ref_sigs, ...)
#   - buildCatalogueCohortSamples(vars_tbl_list, ref_sigs, ...)

# --------------------------------------------------
# Packages
# --------------------------------------------------

#Parallelisation
suppressMessages(require(R.utils))
suppressMessages(require(parallel))

# --------------------------------------------------
# Variant loading functions
# --------------------------------------------------

#Function for reading variant files, excludes chrM
    readVarFile = function(file_names) {
        #Read files into list, omitting chrM variants
        out = mclapply(
            file_names,
            function(file_name) {
                tbl = read.table(file_name, header=T, as.is=T, check.names=F)
                tbl[tbl[,1] != "chrM",]
            }
        )

        #Fix sample/patient names and return
        names(out) = gsub("(.+/)|(_.+)|(\\..+)", "", file_names, perl=T)
        out
    }

# --------------------------------------------------
# Spectra building functions
# --------------------------------------------------

#Function for computing VAF table of each variant of each sample
    varVafTbl = function(vars_tbl) {
        #Get SAC column indices
        sac_indices = which(endsWith(colnames(vars_tbl), ".SAC"))

        #Compute VAF from SAC for each variant/sample pair
        #Make sure output is a table rather than a vector
        rbind(
            apply(
                vars_tbl[, sac_indices, drop=F],
                2,
                function(col) sapply(
                    col,
                    function(sac) {
                        if (is.na(sac)) return(NA)
                        read_counts = as.numeric( unlist(strsplit(sac, ",")) )
                        sum(read_counts[3:4]) / sum(read_counts)
                    }
                )
            )
        )
    }

#Function for computing binary vector given keys to grep like p,i,r or pOvaR1 (case
#insensitive).
#Options include 'invert' (inverts sample choice), 'exclusive' (variant is observed
#only in the selected samples), 'all' (variant must be observed in all of the selected
#samples). This is for sample-based selection of variants.
#Also weeds out any NA lines
    varBinaryVec = function(vars_tbl, keys=".", invert=F, exclusive=F, all=F) {
        #Compute VAF table from variants table
        tbl_vaf = varVafTbl(vars_tbl)

        #Get samples matching any of keys
        sel_idx = sort(unlist(sapply(
            keys,
            function(key) grep(
                key,
                colnames(tbl_vaf),
                ignore.case = T,
                invert = invert
            )
        )))

        #Find exclusive variants in a Boolean (vector), NAs result in F
        excl_binary = 1
        if (exclusive) {
            excl_binary = rowSums(
                tbl_vaf[,-sel_idx, drop=F]) == 0 &
                    !is.na( rowSums(tbl_vaf[,-sel_idx, drop=F] )
            )
        }

        #Final vector, make sure no NAs if all samples have to contain ALT variant
        if (all) {
            #Remove NA lines
            rowSums(tbl_vaf[,sel_idx, drop=F] > 0, na.rm=T) == length(sel_idx) &
                !is.na( rowSums(tbl_vaf[,sel_idx, drop=F]) ) &
                excl_binary
        } else {
            rowSums(tbl_vaf[,sel_idx, drop=F], na.rm=T) > 0 & excl_binary
        }
    }

#Function for returning vector of mutation types from variants table
    getMutTypes = function(variants, mut_class="SBS") {
        #Deal with empty input
        if (nrow(variants) == 0) return(NULL)
        stopifnot(mut_class %in% c("SBS", "DBS", "ID"))

        #SBS constructed from type and subtype columns, others within type
        if (mut_class == "SBS") {
            apply(
                variants[,c("type", "subtype"), drop=F],
                1,
                function(x) {
                    #Combine mutation type with trinucleotide context
                    bases = unlist(strsplit(x[2], ""))
                    paste0(bases[1], "[", x[1], "]", bases[3])
                }
            )
        } else {
            variants$type
        }
    }

#Function for building catalogues with built-in selection of varBinaryVec()
#Takes either a variant table or a list of them as input, and outputs a table for each
#list
    buildCatalogue = function(
        variants, ref_sigs, mut_class="SBS",
        keys=NULL, invert=F, exclusive=F, all=F
    ) {
        #Ensure a valid type is given
        stopifnot(mut_class %in% c("SBS", "DBS", "ID"))

        #Convert single input variant table to a list of table
        if (!is.null(dim(variants))) {
            variants = list(variants)
        }

        #Build catalogue for each table (e.g. patient) and combine them
        out = do.call(
            cbind,
            lapply(
                mclapply(
                    variants,
                    function(vars_tbl) {
                        #Handle NULL table (no variants) as 0 vector
                        if (nrow(vars_tbl) == 0) {
                            return(table(rownames(ref_sigs)) - 1)
                        }

                        #Skip selection if keys input is NULL
                        if (is.null(keys)) {
                            return(table(getMutTypes(vars_tbl, mut_class)))
                        }

                        #Selection if specified
                        sel_vec = varBinaryVec(vars_tbl, keys, invert, exclusive, all)
                        table(getMutTypes(vars_tbl[sel_vec,, drop=F], mut_class))
                    }
                ),
                "[",
                rownames(ref_sigs)
            )
        )
        rownames(out) = rownames(ref_sigs)

        #Set NAs to zero
        out[is.na(out)] = 0
        out
    }

    #Function for building catalogues from all single samples in input variants table.
    #Sample with variant VAF > 0 (i.e. at least one ALT read) will count the variant
    #in the catalogue for that sample.
    buildCataloguePatientSamples = function(vars_tbl, ref_sigs, mut_class="SBS") {
        #Ensure a valid type is given
        stopifnot(mut_class %in% c("SBS", "DBS", "ID"))

        #Compute VAF table from variants table
        tbl_vaf = varVafTbl(vars_tbl)

        #Binary table for ALT allele observed with VAF > 0, NAs considered false
        tbl_binary = tbl_vaf > 0
        tbl_binary[is.na(tbl_binary)] = F

        #Sample names
        if (!is.null(dim(tbl_binary))) {
            samples = colnames(tbl_binary)
        } else {
            samples = names(tbl_binary)
        }

        #Convert variants table to mutation type vector
        mut_types = getMutTypes(vars_tbl, mut_class)

        #Build catalogue for each sample and combine them
        out = do.call(
            cbind,
            lapply(
                mclapply(
                    samples,
                    function(sample) {
                        #Handle NULL table (no variants) as 0 vector
                        if (is.null(mut_types)) {
                            return(table(rownames(ref_sigs)) - 1)
                        }

                        table(mut_types[tbl_binary[,sample]])
                    }
                ),
                "[",
                rownames(ref_sigs)
            )
        )
        rownames(out) = rownames(ref_sigs)

        #Set sample names of output with '.SAC' removed
        colnames(out) = sub(".SAC", "", samples)

        #Set NAs to zero
        out[is.na(out)] = 0
        out
    }

    #Function for building catalogues from all single samples in a list of variant tables.
    #Concatenates each patient's catalogues into a single table and outputs an additional
    #metadata table with a patient column for each sample
    buildCatalogueCohortSamples = function(vars_tbl_list, ref_sigs, mut_class="SBS") {
        #Compute catalogues
        catalogues_list = mclapply(
            vars_tbl_list,
            function(vars_tbl) {
                #For each sample, create a catalogue with variants belonging to it
                buildCataloguePatientSamples(vars_tbl, ref_sigs, mut_class)
            }
        )

        #Patients metadata
        patients = rep(names(catalogues_list), sapply(catalogues_list, ncol))
        names(patients) = do.call(c, sapply(catalogues_list, colnames))

        #Output as list of metadata and catalogues
        list(
            "metadata" = data.frame(patient = patients),
            "spectra" = do.call(cbind, catalogues_list)
        )
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


