#Collection of functions related to signature attribution.
#List of contents:
# - Signature attribution functions:
#   - fitSigExposures(M, sigs, method=sigSSEOptim, obj_fun=selStatSSE, init_exps=NULL)
#   - fitSigFwd(M, sigs, start_sigs=c(), ...)
#   - fitSigBwd(M, sigs, start_sigs=colnames(sigs), ...)
# - Applied signature attribution functions:
#   - fitSigBF(M, sigs, start_sigs=colnames(sigs), ...)
#   - fitSigHierarchicalFwd(M, sigs, start_sigs=c(), ...)
# - Result manipulation functions:
#   - joinSigExposures(exposures_list, signatures)
#List of helper contents:
# - General comparison functions:
#   - Cosine similarity:
#     - cosSim(A,B)
#   - Generalised Kullback-Leibler divergence:
#     - gKLDiv(a,b)
#   - Sum of Squared Errors:
#     - sse(a,b)
#   - Harmonic mean:
#     - harmMean(x)
#   - AIC function:
#     - getAIC(log_likel, k, ...)
#   - BIC function:
#     - getBIC(log_likel, k, n, ...)
# - Selection statistic functions:
#   - Generalised Kullback-Leibler divergence
#     - selStatKL(M, signatures, exposures)
#   - Sum of Squared Errors:
#     - selStatSSE(M, signatures, exposures)
#   - Harmonic mean of accuracy (1 - cos.sim.) and F-norm
#     - selStatHarmMean(M, signatures, exposures)
#   - Accuracy (1 - cosine similarity)
#     - selStatAcc(M, signatures, exposures)
#   - Negative Gaussian log-likelihood
#     - normNLL(M, signatures, exposures)
#   - Negative Poisson log-likelihood
#     - poisNLL(M, signatures, exposures)
#   - BIC using Gaussian log-likelihood
#     - normLLBIC(M, W, H)
#   - BIC using Poisson log-likelihood
#     - poisLLBIC(M, W, H)
# - Optimisation functions:
#   - Method='SA'; Simulated annealing using GenSA:
#     - sigSAOptim(catalogue, signatures, init_exps, obj_fun)
#   - Method='PNNLS'; Least Squares optimisation using PNNLS from LSEI:
#     - sigSSEOptim(catalogue, signatures, init_exps, obj_fun)
#   - Method='NMF_KL'; KL-divergence using NMF H-update (sklearn or NMF):
#     - sigKLOptim(catalogue, signatures, init_exps, obj_fun)
# - Signature attribution rule functions:
#   - fitSigsRemoveByTSB(M, sigs, curr_sigs, M_tsb, tsb_rules, ...)
#   - fitSigsRemoveByMutCount(M, sigs, curr_sigs, mut_count_rules, ...)
#   - fitSigsAddConnected(M, sigs, curr_sigs, connected_sigs, ...)
#   - fitSigsAddForced = function(M, sigs, curr_sigs, forced_sigs, ...)
#   - fitSigApplyRules(M, sigs, curr_sigs, forced_sigs=NULL, connected_sigs=NULL, M_tsb=NULL, tsb_rules=NULL, mut_count_rules=NULL)

# --------------------------------------------------
# Packages
# --------------------------------------------------

#Optimisation packages
suppressMessages(require(lsei))
suppressMessages(require(GenSA))
#Not needed if sklearn is used
suppressMessages(require(NMF))
#Python sklearn via reticulate
suppressMessages(require(reticulate))
#use_python("~/virtual_py_env/bin/python", required=T)
sklearn_nmf = import("sklearn.decomposition")
#Parallelisation
suppressMessages(require(R.utils))
suppressMessages(require(parallel))

# --------------------------------------------------
#
# --------------------------------------------------

# --------------------------------------------------
# General comparison functions
# --------------------------------------------------

#Cosine similarity function
    cosSim = function(A,B) {
        A_num = as.numeric(A)
        B_num = as.numeric(B)
        A_num %*% B_num / sqrt( sum(A_num^2) * sum(B_num^2) )
    }

#Generalised Kullback-Leibler divergence
    gKLDiv = function(a,b) {
        A = as.numeric(a)
        B = as.numeric(b)
        sum(A * log(A/B) - A + B, na.rm=T)
    }

#SSE function
    sse = function(a,b) {
        sum((as.numeric(a) - as.numeric(b))^2)
    }

#Harmonic mean
    harmMean = function(x) {
        1 / mean(1/x)
    }

#Function for computing AIC
    getAIC = function(log_likel, k, ...) {
        2 * k - 2 * log_likel
    }

#Function for computing BIC
    getBIC = function(log_likel, k, n, ...) {
        log(n) * k - 2 * log_likel
    }

# --------------------------------------------------
# Model selection stat / objective functions
# --------------------------------------------------

#KL selection stat
    selStatKL = function(M, signatures, exposures) {
        #V ~ WH, proportions
        V = signatures %*% exposures
        gKLDiv(M, V)
    }

#SSE selection stat
    selStatSSE = function(M, signatures, exposures) {
        #V ~ WH, proportions
        V = signatures %*% exposures
        sum((M - V)^2)
    }

#Harmonic mean of accuracy and F-norm selection stat
    selStatHarmMean = function(M, signatures, exposures) {
        #V ~ WH
        V = signatures %*% exposures
        #Harmonic mean of accuracy and F-norm
        harmMean( c(1 - cosSim(M, V), norm(M - V, "F")) )
    }

#Accuracy selection stat
    selStatAcc = function(M, signatures, exposures) {
        #V ~ WH
        V = signatures %*% exposures
        #Accuracy
        1 - cosSim(M, V)
    }

#Function for computing Normal negative log-likelihood; input exposures assumed MLEs
    normNLL = function(M, signatures, exposures) {
        #V ~ WH
        V = signatures %*% exposures

        n = length(M)
        rss = sum((M-V)^2)
        #Gaussian NLL
        - (-n/2 * log(2*pi) - n/2 * log( rss / n ) - n/2)
    }

#Function for computing Poisson negative log-likelihood
    poisNLL = function(M, signatures, exposures) {
        #V ~ WH
        V = signatures %*% exposures
        #Poisson NLL; M given V
        - sum(dpois(M, V, T))
    }

#Gaussian log-likelihood BIC selection stat for attribution (note sum-to-total constraint)
    normLLBIC = function(M, W, H) {
        getBIC(-normNLL(M, W, H), length(H) - ncol(H), length(M))
    }

#Poisson log-likelihood BIC selection stat for attribution (note sum-to-total constraint)
    poisLLBIC = function(M, W, H) {
        getBIC(-poisNLL(M, W, H), length(H) - ncol(H), length(M))
    }

# --------------------------------------------------
# Signature attribution optimisation functions
# --------------------------------------------------

#Function for translating string to optimisation function
    optimMethod = function(method) {
        if (is.character(method)) {
            switch(method,
                SA = {method = sigSAOptim},
                PNNLS = {method = sigSSEOptim},
                NMF_KL = {method = sigKLOptim}
            )
        }
        method
    }

#GenSA objective optimisation function
    sigSAOptim = function(catalogue, signatures, init_exps, obj_fun) {
        #Variables
        n = colSums(catalogue) * 1.0
        k = ncol(signatures)

        #Force usage of correct function
        if (identical(sse, obj_fun)) obj_fun = selStatSSE
        if (identical(gKLDiv, obj_fun)) obj_fun = selStatKL

        #Use smoothed signatures; add minimum pseudocounts to ensure that optimisation algorithms work
        sigs_smooth = signatures + 1e-99
        for (i in 1:ncol(sigs_smooth)) {
            sigs_smooth[,i] = prop.table(sigs_smooth[,i])
        }

        #Handle multiple catalogues
        exposures = matrix(
            sapply(
                1:ncol(catalogue),
                function(i) {
                    #Temporary variables
                    catal = catalogue[,i,drop=F]
                    catal_size = n[i]
                    #Run simulated annealing and get result
                    GenSA(
                        init_exps[,i],
                        function(par) obj_fun(catal, sigs_smooth, prop.table(par)*catal_size),
                        lower = rep(0, k),
                        upper = rep(catal_size, k),
                        control = list("maxit" = 1000, "temperature" = 100, "trace.mat"=F)
                    )$par
                }
            ),
            ncol = ncol(catalogue)
        )
        #Set near-zeros to zero (likely caused by smoothing)
        exposures[exposures < 1e-1] = 0

        #Return result
        sweep(exposures, 2, colSums(catalogue) / margin.table(exposures, 2), "*")
    }

#LSEI PNNLS function for optimising SSE
    sigSSEOptim = function(catalogue, signatures, init_exps, obj_fun) {
        #Handle multiple catalogues
        matrix(
            sapply(
                1:ncol(catalogue),
                function(i) pnnls(signatures, catalogue[,i,drop=F], 0, sum(catalogue[,i]))$x
            ),
            ncol = ncol(catalogue)
        )
    }

if (T) {
#sklearn H-update function for optimising KL divergence
    sigKLOptim = function(catalogue, signatures, init_exps, obj_fun) {
        #Use smoothed signatures; add minimum pseudocounts to ensure that optimisation algorithms works
        sigs_smooth = signatures + 1e-99
        for (i in 1:ncol(sigs_smooth)) {
            sigs_smooth[,i] = prop.table(sigs_smooth[,i])
        }

        #Fit iteratively
        exposures = t(sklearn_nmf$non_negative_factorization(
            X = r_to_py(t(catalogue)),
            H = r_to_py(t(signatures)),
            W = r_to_py(t(init_exps)),
            init = 'custom',
            update_H = F,
            beta_loss = "kullback-leibler",
            solver = "mu",
            n_components = ncol(signatures),
            tol = 1e-12,
            max_iter = as.integer(25000)
        )[[1]])

        #Set near-zeros to zero (likely caused by smoothing)
        exposures[exposures < 1e-1] = 0

        #Return result
        sweep(exposures, 2, colSums(catalogue) / margin.table(exposures, 2), "*")
    }
} else {
#NMF H-update function for optimising KL divergence
#NOT IN USE AS SKLEARN IS FASTER
    sigKLOptim = function(catalogue, signatures, init_exps, obj_fun) {
        #Use smoothed signatures; add minimum pseudocounts to ensure that optimisation algorithms works
        sigs_smooth = signatures + 1e-99
        for (i in 1:ncol(sigs_smooth)) {
            sigs_smooth[,i] = prop.table(sigs_smooth[,i])
        }

        #Initialisation
        exposures = as.matrix(init_exps, ncol=ncol(catalogue), nrow=ncol(signatures))
        #Iterative matrix update
        for (i in 1:10000) {
            exposures = nmf_update.KL.h(catalogue, w=sigs_smooth, h=exposures)
        }

        #Set near-zeros to zero (likely caused by smoothing)
        exposures[exposures < 1e-1] = 0

        #Return result
        sweep(exposures, 2, colSums(catalogue) / margin.table(exposures, 2), "*")
    }
}

# --------------------------------------------------
# Signature attribution functions
# --------------------------------------------------

#Signature exposure fitting function (objective function optimisation) with GenSA
    fitSigExposures = function(M, sigs, method=sigSSEOptim, obj_fun=selStatSSE, init_exps=NULL) {
        #Make method refer to function if given as string name
        method = optimMethod(method)

        #Initialise init_exps
        if (is.null(init_exps)) {
            init_exps = matrix(runif(ncol(sigs) * ncol(M)), ncol=ncol(M))
            init_exps = sweep(init_exps, 2, colSums(M) / margin.table(init_exps, 2), "*")
        }

        #Run optimisation
        exposures = method(M, sigs, init_exps, obj_fun)

        #Return exposure vector with names
        rownames(exposures) = colnames(sigs)
        colnames(exposures) = colnames(M)
        exposures
    }

#Signature forward selection function
    fitSigFwd = function(M, sigs, start_sigs=c(), f_sel_stat=selStatAcc, f_stop_stat=selStatAcc, stop_thres=0.05, stop_n=ncol(sigs), method=sigSSEOptim, obj_fun=selStatSSE, keep_steps=F) {
        #Make method refer to function if given as string name
        method = optimMethod(method)

        #If more than 0 start_sigs, fit starting configuration and compute stop stat
        if (length(start_sigs) > 0) {
            #Starting with 1+ signatures
            prev_sigs = colnames(sigs[,start_sigs,drop=F])
            prev_exposures = fitSigExposures(M, sigs[,prev_sigs,drop=F], method, obj_fun)
            prev_stop_stat = f_stop_stat(M, sigs[,prev_sigs,drop=F], prev_exposures)
        } else {
            #Starting from zero signature; fit all signatures and compute selection stat
            sel_stats = sapply(
                1:ncol(sigs),
                function(i) {
                    f_sel_stat(M, sigs[,i,drop=F], fitSigExposures(
                        M,
                        sigs[,i,drop=F],
                        method,
                        obj_fun
                    ))
                }
            )

            #Get best signature and compute stop stat
            prev_sigs = colnames(sigs)[which(sel_stats == min(sel_stats))]
            prev_exposures = matrix(colSums(M), ncol=ncol(M))
            rownames(prev_exposures) = prev_sigs
            colnames(prev_exposures) = colnames(M)
            prev_stop_stat = f_stop_stat(M, sigs[,prev_sigs,drop=F], prev_exposures)
        }

        #Keep list of exposures if desired
        if (keep_steps) {exp_steps = list(prev_exposures)}

        #Forward selection; break if stop stat difference is smaller than threshold
        while( length(prev_sigs) < min(ncol(sigs),stop_n) ) {
            #Fit exposures with each extra signature
            next_exposures = mclapply(
                colnames(sigs)[!colnames(sigs) %in% prev_sigs],
                function(sig_label) {
                    fitSigExposures(
                        M,
                        sigs[,c(prev_sigs,sig_label),drop=F],
                        method,
                        obj_fun
                    )
                }
            )

            #Remove new signature if it returns zero and stop if no candidates remain
            next_exposures = next_exposures[sapply(next_exposures, function(x) sum(x[nrow(x),]) > 0)]
            if (length(next_exposures) == 0) break

            #Compute selection stat
            sel_stats = sapply(
                next_exposures,
                function(exposures) {
                    f_sel_stat(M, sigs[,rownames(exposures),drop=F], exposures)
                }
            )

            #Pick best signature
            curr_exposures = next_exposures[[ which(sel_stats == min(sel_stats))[1] ]]
            curr_sigs = rownames(curr_exposures)
            curr_stop_stat = f_stop_stat(M, sigs[,curr_sigs,drop=F], curr_exposures)

            #Update variables and continue if stopping threshold not reached
            if (prev_stop_stat - curr_stop_stat > stop_thres) {
                prev_exposures = curr_exposures
                prev_sigs = curr_sigs
                prev_stop_stat = curr_stop_stat

                #Keep list of exposures if desired
                if (keep_steps) {exp_steps[[length(exp_steps) + 1]] = prev_exposures}
            } else {
                #Add the rejected last step if desired
                if (keep_steps > 1) {exp_steps[[length(exp_steps) + 1]] = curr_exposures}
                break
            }
        }

        #Return value depends on if end result or all steps are required
        if (!keep_steps) {
            #Sort final exposures
            out = prev_exposures[colnames(sigs)[colnames(sigs) %in% prev_sigs],,drop=F]
        } else {
            #Sort all exposures
            out = lapply(
                exp_steps,
                function(exposures) {
                    exposures[colnames(sigs)[colnames(sigs) %in% rownames(exposures)],,drop=F]
                }
            )
        }

        out
    }

#Signature backward selection function
    fitSigBwd = function(M, sigs, start_sigs=colnames(sigs), f_sel_stat=selStatAcc, f_stop_stat=selStatAcc, stop_thres=0.01, stop_n=1, method=sigSSEOptim, obj_fun=selStatSSE, keep_steps=F) {
        #Do nothing if <1 signatures
        if (length(start_sigs) < 1) {
            stop("No start signatures")
        }

        #Make method refer to function if given as string name
        method = optimMethod(method)

        #Fit starting configuration and compute stop stat
        prev_sigs = colnames(sigs[,start_sigs,drop=F])
        prev_exposures = fitSigExposures(M, sigs[,prev_sigs,drop=F], method, obj_fun)
        prev_stop_stat = f_stop_stat(M, sigs[,prev_sigs,drop=F], prev_exposures)

        #Keep list of exposures if desired
        if (keep_steps) {exp_steps = list(prev_exposures)}

        #Get rid of 0 value signatures
        if (any(rowSums(prev_exposures) == 0)) {
            prev_sigs = names(which(rowSums(prev_exposures) > 0))
            prev_exposures = fitSigExposures(M, sigs[,prev_sigs,drop=F], method, obj_fun)
            prev_stop_stat = f_stop_stat(M, sigs[,prev_sigs,drop=F], prev_exposures)

            #Keep list of exposures if desired
            if (keep_steps) {exp_steps[[length(exp_steps) + 1]] = prev_exposures}
        }

        #Backward selection; break if stop stat difference is greater than threshold
        while( length(prev_sigs) > max(1,stop_n) ) {
            #Fit exposures with each extra signature
            next_exposures = mclapply(
                1:length(prev_sigs),
                function(i) {
                    fitSigExposures(
                        M,
                        sigs[,prev_sigs[-i],drop=F],
                        method,
                        obj_fun
                    )
                }
            )

            #Compute selection stat
            sel_stats = sapply(
                next_exposures,
                function(exposures) {
                    f_sel_stat(M, sigs[,rownames(exposures),drop=F], exposures)
                }
            )

            #Pick best signature
            curr_exposures = next_exposures[[ which(sel_stats == min(sel_stats))[1] ]]
            curr_sigs = rownames(curr_exposures)
            curr_stop_stat = f_stop_stat(M, sigs[,curr_sigs,drop=F], curr_exposures)

            #Update variables and continue if stopping threshold not reached
            if (curr_stop_stat - prev_stop_stat < stop_thres) {
                prev_exposures = curr_exposures
                prev_sigs = curr_sigs
                prev_stop_stat = curr_stop_stat

                #Keep list of exposures if desired
                if (keep_steps) {exp_steps[[length(exp_steps) + 1]] = prev_exposures}
            } else {
                #Add the rejected last step if desired
                if (keep_steps > 1) {exp_steps[[length(exp_steps) + 1]] = curr_exposures}
                break
            }
        }

        #Return value depends on if end result or all steps are required
        if (!keep_steps) {
            #Sort final exposures
            out = prev_exposures[colnames(sigs)[colnames(sigs) %in% prev_sigs],,drop=F]
        } else {
            #Sort all exposures
            out = lapply(
                exp_steps,
                function(exposures) {
                    exposures[colnames(sigs)[colnames(sigs) %in% rownames(exposures)],,drop=F]
                }
            )
        }

        out
    }

# --------------------------------------------------
# Signature attribution rule functions
# --------------------------------------------------

#Transcriptional strand bias rules
    fitSigsRemoveByTSB = function(M, sigs, curr_sigs, M_tsb, tsb_rules, method=sigSSEOptim, obj_fun=selStatSSE, no_fit=F) {
        #Select signatures present to test
        subset_rules = tsb_rules[rownames(tsb_rules) %in% curr_sigs,]
        tests = sapply(
            1:0,
            function(exon) {
                out = sapply(subset_rules[subset_rules$exon == exon, colnames(subset_rules) != "exon"], unique)
                out = out[sapply(out, length) > 1, drop=F]
                sapply(out, sum)
            }
        )

        #Tests to perform; combine over signatures
        tests = unique(do.call(
            rbind,
            apply(
                subset_rules,
                1,
                #Combine over mutation types
                function(row) do.call(
                    rbind,
                    lapply(
                        2:7,
                        function(i) if (row[i] != 0) {
                            #Output exome/genome, mutation type, test direction
                            c(
                                c("genome", "exome")[row[1]+1],
                                colnames(subset_rules)[i],
                                if (row[i] < 0) "less" else "greater"
                            )
                        }
                    )
                )
            )
        ))
        colnames(tests) = c("exon", "type", "alt")

        #Do the tests; loop over sample transcriptional strand information
        p_values = sapply(
            M_tsb,
            function(tsb_tbl) apply(
                tests,
                1,
                function(test) if (sum(tsb_tbl[test["type"],,test["exon"]]) > 0) {
                    binom.test(
                        tsb_tbl[test["type"],1,test["exon"]],
                        sum(tsb_tbl[test["type"],,test["exon"]]),
                        alternative = test["alt"]
                    )$p.value
                } else {
                    NA
                }
            )
        )

        #Apply FDR adjustment and determine if any sample passes the tests
        p_vals_adj = matrix(p.adjust(p_values, "fdr"), nrow = nrow(p_values))
        colnames(p_vals_adj) = colnames(p_values)
        test_results = cbind(tests, "pass"=(apply(p_vals_adj, 1, min, na.rm=T) < 0.01))

        #Determine and remove failing signatures
        sigs_to_remove = rownames(subset_rules)[
            rowSums(apply(
                tests[test_results[,"pass"] == F,],
                1,
                function(test) {
                    (subset_rules$exon == (test["exon"] == "exome")) & (
                        c("less","","greater")[subset_rules[,colnames(subset_rules) == test["type"]] + 2] == test["alt"]
                    )
                }
            )) > 0
        ]

        #Remove the signatures
        new_sigs = curr_sigs[! curr_sigs %in% sigs_to_remove]
        new_sigs = colnames(sigs)[colnames(sigs) %in% new_sigs]

        if (no_fit) return(new_sigs)

        #Fit with new signatures
        fitSigExposures(M, sigs[,new_sigs,drop=F], method=sigSSEOptim, obj_fun=selStatSSE)
    }

#Mutation count rules
    fitSigsRemoveByMutCount = function(M, sigs, curr_sigs, mut_count_rules, method=sigSSEOptim, obj_fun=selStatSSE, no_fit=F) {
        #Find signatures to remove
        sigs_to_remove = names(mut_count_rules)[
            rowSums(sapply(colSums(M), function(mc) mc > mut_count_rules)) == 0
        ]

        #Remove the signatures
        new_sigs = curr_sigs[! curr_sigs %in% sigs_to_remove]
        new_sigs = colnames(sigs)[colnames(sigs) %in% new_sigs]

        if (no_fit) return(new_sigs)

        #Fit with new signatures
        fitSigExposures(M, sigs[,new_sigs,drop=F], method=sigSSEOptim, obj_fun=selStatSSE)
    }

#Adding connected signatures
    fitSigsAddConnected = function(M, sigs, curr_sigs, connected_sigs, method=sigSSEOptim, obj_fun=selStatSSE, no_fit=F) {
        #Find connected signatures
        sigs_to_add = unlist(sapply(
            connected_sigs,
            function(sig_list) if (sum(curr_sigs %in% sig_list) > 0) sig_list else NULL
        ))

        #Restrict addition to signatures to consider
        sigs_to_add = sigs_to_add[sigs_to_add %in% colnames(sigs)]
        sigs_to_add = sigs_to_add[! sigs_to_add %in% curr_sigs]

        #Combine and sort signatures
        new_sigs = c(curr_sigs, sigs_to_add)
        new_sigs = colnames(sigs)[colnames(sigs) %in% new_sigs]

        if (no_fit) return(new_sigs)

        #Fit with new signatures
        exposures = fitSigExposures(M, sigs[,new_sigs,drop=F], method=sigSSEOptim, obj_fun=selStatSSE)

        #Remove new signatures with zero exposure
        if (any(rowSums(exposures)[sigs_to_add] == 0)) {
            remaining_sigs = names(which(rowSums(exposures) > 0))
            exposures = fitSigExposures(M, sigs[,remaining_sigs,drop=F], method, obj_fun)
        }

        exposures
    }

#Adding forced signatures
    fitSigsAddForced = function(M, sigs, curr_sigs, forced_sigs, method=sigSSEOptim, obj_fun=selStatSSE, no_fit=F) {
        #Restrict addition to signatures to consider
        sigs_to_add = forced_sigs[forced_sigs %in% colnames(sigs)]
        sigs_to_add = sigs_to_add[! sigs_to_add %in% curr_sigs]

        #Combine and sort signatures
        new_sigs = c(curr_sigs, sigs_to_add)
        new_sigs = colnames(sigs)[colnames(sigs) %in% new_sigs]
        if (no_fit) return(new_sigs)

        #Fit with new signatures
        exposures = fitSigExposures(M, sigs[,new_sigs,drop=F], method=sigSSEOptim, obj_fun=selStatSSE)

        #Remove new signatures with zero exposure
        if (any(rowSums(exposures)[sigs_to_add] == 0)) {
            remaining_sigs = names(which(rowSums(exposures) > 0))
            exposures = fitSigExposures(M, sigs[,remaining_sigs,drop=F], method, obj_fun)
        }

        exposures
    }

#Signature rule application function
    fitSigApplyRules = function(M, sigs, curr_sigs, forced_sigs=NULL, connected_sigs=NULL, M_tsb=NULL, tsb_rules=NULL, mut_count_rules=NULL) {
        new_sigs = curr_sigs

        #Remove signatures based on transcription bias
        if (!is.null(new_sigs) & !is.null(M_tsb) & !is.null(tsb_rules)) {
            new_sigs = fitSigsRemoveByTSB(M, sigs, new_sigs, M_tsb, tsb_rules, no_fit=T)
        }

        #Remove signatures based on mutation count
        if (!is.null(new_sigs) & !is.null(mut_count_rules)) {
            new_sigs = fitSigsRemoveByMutCount(M, sigs, new_sigs, mut_count_rules, no_fit=T)
        }

        #Add connected signatures
        if (!is.null(new_sigs) & !is.null(connected_sigs)) {
            new_sigs = fitSigsAddConnected(M, sigs, new_sigs, connected_sigs, no_fit=T)
        }

        #Add forced signatures
        if (!is.null(forced_sigs)) {
            new_sigs = fitSigsAddForced(M, sigs, new_sigs, forced_sigs, no_fit=T)
        }

        #Return resulting signatures
        if (is.null(new_sigs)) return(NULL)

        colnames(sigs)[colnames(sigs) %in% new_sigs]
    }

# --------------------------------------------------
# Applied signature attribution functions
# --------------------------------------------------

#Signature backward-forward selection function
#Includes options for signature rules (supplied in arguments), before backward search:
# - removal based on transcriptional strand bias of 6-type classification;
#   table with rownames as signature names, Boolean "exon" column signifying whether
#   to include intronic variants, and 6-type mutation columns (C>A, C>G, ...) with
#   (-1/)1 to test for (un)transcribed strand bias
# - removal based on mutation count;
#   vector of minimum counts with signature names (e.g. MMR signatures)
#as well as before and after each step:
# - addition of forced signatures;
#   vector of signature names (e.g. age signatures)
# - adding connected signatures;
#   list of vectors of signature names (e.g. APOBEC signatures)
    fitSigBF = function(
        M, sigs, start_sigs=colnames(sigs),
        f_sel_stat=selStatAcc, f_stop_stat=selStatAcc,
        bwd_thres=0.01, fwd_thres=0.05, bwd_stop_n=0, fwd_stop_n=ncol(sigs),
        forced_sigs=NULL, connected_sigs=NULL, M_tsb=NULL, tsb_rules=NULL, mut_count_rules=NULL,
        method=sigSSEOptim, obj_fun=selStatSSE,
        keep_steps=F
    ) {
        sigs_in = start_sigs

        #Apply signature rules
        sigs_in = fitSigApplyRules(M, sigs, sigs_in, forced_sigs, connected_sigs, M_tsb, tsb_rules, mut_count_rules)

        #Backward search only if there is at least one signature
        if (length(start_sigs) >= 1) {
            exposures = fitSigBwd(M, sigs, sigs_in, f_sel_stat, f_stop_stat, bwd_thres, bwd_stop_n, method, obj_fun, keep_steps)
            #Clear exposures if there is only one left to ensure that a poorly fitting signature does not remain
            if (!keep_steps) {
                sigs_in = rownames(exposures)
            } else {
                sigs_in = rownames(exposures[[length(exposures) - keep_steps + 1]])
            }
            if (length(sigs_in) == 1 & bwd_stop_n < 1) sigs_in = NULL
        } else {
            exposures = NULL
            sigs_in = NULL
        }

        #Apply signature rules (forced and connected only)
        sigs_in = fitSigApplyRules(M, sigs, sigs_in, forced_sigs, connected_sigs)

        #Forward search, slightly different input and processing if tracking steps
        if (!keep_steps) {
            exposures = fitSigFwd(M, sigs, sigs_in, f_sel_stat, f_stop_stat, fwd_thres, fwd_stop_n, method, obj_fun, keep_steps)
        } else {
            #Deal with duplicate fits if steps are kept:
            #At least one remaining signature
            if (length(sigs_in) > 1) {
                #Remove duplicate fit for tracked steps (keep_steps 1) when signature rules have no effect
                if (keep_steps == 1 & sum(rownames(exposures[[length(exposures)]]) %in% sigs_in) == length(sigs_in)) {
                    exposures = exposures[1:(length(exposures)-1)]
                #Add the final accepted step as duplicate fit if rejected steps are also kept (keep_steps 2)
                #and signature rules have an effect
                } else if (keep_steps == 2 & sum(rownames(exposures[[length(exposures)]]) %in% sigs_in) < length(sigs_in)) {
                    exposures = c(exposures, exposures[length(exposures)-1])
                }
            }
            #Apply forward selection
            exposures = c(
                exposures,
                fitSigFwd(M, sigs, sigs_in, f_sel_stat, f_stop_stat, fwd_thres, fwd_stop_n, method, obj_fun, keep_steps)
            )
        }

        #Apply signature rules (forced and connected only)
        if (!keep_steps) {
            sigs_in = rownames(exposures)
        } else {
            sigs_in = rownames(exposures[[length(exposures) - keep_steps + 1]])
        }
        sigs_in = fitSigApplyRules(M, sigs, sigs_in, forced_sigs, connected_sigs)

        #Fit the final signatures
        if (!keep_steps) {
            fitSigExposures(M, sigs[,sigs_in,drop=F], method, obj_fun)
        } else {
            #Remove duplicate fit for tracked steps (keep_steps 1) when signature rules have no effect
            if (keep_steps == 1 & sum(rownames(exposures[[length(exposures)]]) %in% sigs_in) == length(sigs_in)) {
                exposures = exposures[1:(length(exposures)-1)]
                #Add the final accepted step as duplicate fit if rejected steps are also kept (keep_steps 2)
                #and signature rules have an effect
            } else if (keep_steps == 2 & sum(rownames(exposures[[length(exposures)]]) %in% sigs_in) < length(sigs_in)) {
                exposures = c(exposures, exposures[length(exposures)-1])
            }
            c(exposures, list(fitSigExposures(M, sigs[,sigs_in,drop=F], method, obj_fun)))
        }
    }

#Dataset signature hierarchical forward selection
    fitSigHierarchicalFwd = function(
        M, sigs, start_sigs=c(),
        f_sel_stat=selStatAcc, f_stop_stat=selStatAcc,
        stop_thres=0.05, sigs_per_step=1000, min_remaining=2, min_remaining_prop=0.05,
        method=sigSSEOptim, obj_fun=selStatSSE,
        adj_stop_thres=F, log=F
    ) {
        #Prepare variables
        adj_thres = stop_thres
        M_poor = M
        n_sigs = 0
        exps_out = NULL
        exposures = NULL
        if (length(start_sigs) > 0) {
            exposures = matrix(1, nrow=length(start_sigs))
            rownames(exposures) = start_sigs
        }

        #Hierarchical fitting
        while(T) {
            #Stop threshold adjustment
            if (adj_stop_thres) {
                if (log) {
                    adj_thres = stop_thres + log(ncol(M)) - log(ncol(M_poor))
                } else {
                    adj_thres = stop_thres * ncol(M) / ncol(M_poor)
                }
            }

            #Fit signatures forward
            exposures = fitSigFwd(
                M_poor, sigs, rownames(exposures),
                f_sel_stat, f_stop_stat, adj_thres,
                min(n_sigs + sigs_per_step, ncol(sigs)),
                method, obj_fun
            )
            M_fit = sigs[, rownames(exposures), drop=F] %*% exposures

            #Compute cosine similarities of each sample and check threshold
            idx_M_poor = sapply(
                1:ncol(M_poor),
                function(i) cosSim(M_poor[,i], M_fit[,i])
            ) < .95

            #Count the number of poor fits and if any additional signature was fitted for stopping
            if (
                sum(idx_M_poor) == 0 |
                sum(idx_M_poor) < min_remaining |
                sum(idx_M_poor) / ncol(M) < min_remaining_prop |
                nrow(exposures) == n_sigs
            ) {
                #Bind exposures
                if (is.null(exps_out) | length(exps_out) == 0) {
                    exps_out = exposures
                } else {
                    missing_sigs = rownames(exposures)[! rownames(exposures) %in% rownames(exps_out)]
                    if (ncol(exps_out) > 0 & length(missing_sigs) > 0) {
                        zeroes = data.frame(
                            matrix(0, ncol=ncol(exps_out), nrow=length(missing_sigs)),
                            row.names = missing_sigs
                        )
                        colnames(zeroes) = colnames(exps_out)
                        exps_out = rbind(exps_out, zeroes)[rownames(exposures),]
                    }
                    exps_out = cbind(exps_out, exposures)
                }
                break
            }

            #Split well and poorly fitted catalogues; bind output exposures
            if (is.null(exps_out) | length(exps_out) == 0) {
                exps_out = exposures[, !idx_M_poor, drop=F]
            } else {
                missing_sigs = rownames(exposures)[! rownames(exposures) %in% rownames(exps_out)]
                if (ncol(exps_out) > 0 & length(missing_sigs) > 0) {
                    zeroes = data.frame(
                        matrix(0, ncol=ncol(exps_out), nrow=length(missing_sigs)),
                        row.names = missing_sigs
                    )
                    colnames(zeroes) = colnames(exps_out)
                    exps_out = rbind(exps_out, zeroes)[rownames(exposures),]
                }
                exps_out = cbind(exps_out, exposures[,!idx_M_poor, drop=F])
            }
            M_poor = M_poor[,idx_M_poor, drop=F]
            n_sigs = nrow(exposures)
        }

        #Set NAs to 0 and return sorted result
        exps_out[is.na(exps_out)] = 0
        as.matrix(exps_out[
            colnames(sigs)[colnames(sigs) %in% rownames(exps_out)],
            order(colnames(exps_out)),
            drop = F
        ])
    }

# --------------------------------------------------
# Result manipulation helper functions
# --------------------------------------------------

#Signature exposure list joining function
#Joins exposure tables with different signatures - missing contributions set to 0.
#Input signatures determines signature order in the output
    suppressMessages(require(plyr))
    joinSigExposures = function(exposures_list, signatures=NULL) {
        #Join the exposures by signatures as rows
        exps_joined = do.call(
            rbind.fill,
            lapply(
                exposures_list,
                function(x) {
                    data.frame(col_name = colnames(x), t(x))
                }
            )
        )

        #Restore column name and replace NAs with 0
        rownames(exps_joined) = exps_joined$col_name
        exps_joined = exps_joined[, !colnames(exps_joined) %in% "col_name", drop=F]
        exps_joined[is.na(exps_joined)] = 0

        #Sort signatures in reference signature order if given
        if (!is.null(signatures)) {
            exps_joined = exps_joined[
                ,
                colnames(ref_sigs)[
                    colnames(ref_sigs) %in% colnames(exps_joined)
                ],
                drop = F
            ]
        }

        #Output as matrix in original orientation
        as.matrix(t(exps_joined))
    }

# --------------------------------------------------
# End of file
# --------------------------------------------------


