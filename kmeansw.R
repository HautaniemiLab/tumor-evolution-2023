
#' K-means for weighted and missing data
#'
#' @param X Data matrix, rows are samples and columns are variables. NAs are
#'   converted to zero-weight elements.
#' @param W Weight matrix. Can be a constant (e.g. unity), a sample weight
#'   vector, a data element uncertainty weight (matrix).
#' @param k Number of clusters.
#' @param nstart Number of restarts (default: 1).
#' @param iter.max Maximum number of iterations (default: 1000).
#'
#' @return A list with the following members:
#'   $L Class labels for each of the samples.
#'   $Y Cluster centroids for each of the clusters.
#'   $sum.res Total model residual.
#'
kmeansw <- function(X, W = 1., k, nstart = 1L, iter.max = 1000L) {
	# get dimensions
	m <- nrow(X)
	n <- ncol(X)

	# prepare weights
	W <- W * matrix(1., m, n)
	stopifnot(identical(dim(W), c(m, n)))

	# convert NAs to zero weights
	W[is.na(X)] <- 0.
	X[is.na(X)] <- 0.

	# create labels
	L <- rep(NA_integer_, m)
	names(L) <- rownames(X)

	# create centroids
	Y <- matrix( NaN, k, n )
	rownames(Y) <- seq_len(k)
	colnames(Y) <- colnames(X)

	# set up 
	best = list( L = L, Y = Y, sum.res = Inf )
#write(sprintf('kmeansw(X = <%d-by-%d>, W, k = %d, nstart = %d, iter.max = %d)', m, n, k, nstart, iter.max), file = stderr())

	# run restart
	for (start in seq_len(nstart)) {
#write(sprintf('start %d..', start), file = stderr())
		# set up
		L[] <- sample(k, m, replace = T)
		sum.res <- Inf

		# optimize
		for (iter in seq_len(iter.max)) {
			#
			# TODO: here, I have independent variables (no cov)
			# should I estimate cov from data?
			#

			# infer centroids
			Y[] <- .kmeansw.bucket( L, W*X, L.max = k ) / .kmeansw.bucket( L, W, L.max = k )
			
			# relabel
			d <- rep( Inf, m )
			for (l in seq_len(k)) {
				R <- t( t(X) - Y[l, ] )
				d1 <- rowSums( W*(R*R) )
				mask <- d1 < d
				mask[is.na(mask)] <- F
				L[mask] <- l
				d[mask] <- d1[mask]
			}

			# get cost
			last.sum.res <- sum.res
			sum.res <- sum(d)
#write(sprintf('iter = %d, sum.res = %g', iter, sum(d)), file = stderr())

			# early exit
			if (!(sum.res < last.sum.res)) {
#write(sprintf('early exit at iter = %d', iter), file = stderr())
				break
			}

		}

		# get cost
		sum.res <- sum(d)
#write(sprintf('final sum.res = %g', sum.res), file = stderr())

		# keep best
		if (sum.res < best$sum.res)
			best <- list( L = L, Y = Y, sum.res = sum.res )
	}

	# return best
#write(sprintf('best sum.res = %g', sum.res), file = stderr())
	return (best)
}

#  bucket data
.kmeansw.bucket <- function(L, X, L.max = max(L), fill = NaN) {
	Y <- matrix( fill, L.max, ncol(X) )
	for (k in seq_len(L.max))
		Y[k, ] <- colSums( X[L == k, , drop = F] )
	return (Y)
}
