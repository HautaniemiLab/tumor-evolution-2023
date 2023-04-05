# author: Sanaz Jamalzadeh, Antti Häkkinen

# find models correspond to each patient from the phylogenetic tress as input
models <- read.table('path to selected mutation tress from pyclone as csv file', 
                     sep = '\t', header = T)

# 
# this assembles the results from the frequency tables
#   check if these are the implied models
#

# find files
files <- list.files('path to the generated trees',
	pattern = '*_cellular_freqs.csv', recursive = T, full.names = T)

results <- NULL
for (file in files) {
write(sprintf('processing \'%s\'..', file), file = stderr())
	patient <- sub( '^([^_]+\\d+)(_v2)?_vaf_(.*)_cellular_freqs\\.csv$', '\\1', basename(file) )

	freqs <- read.table(file, sep = '\t', header = T)
	freqs <- freqs[ freqs$model.num == models$model[grep(patient, models$patient)], ]

	unique.samples <- unique( freqs$sample.id )
	unique.clones <- seq_len(max( freqs$cloneID ))
	values <- matrix( 0., length(unique.clones), length(unique.samples) )
	rownames(values) <- unique.clones
	colnames(values) <- unique.samples
	values[ cbind( match( freqs$cloneID, unique.clones ), match( freqs$sample.id, unique.samples ) ) ] <- freqs$cell.freq / 100.

	stopifnot(all( 0 <= values & values <= 1. ))
	stopifnot(all( 1-.05 <= colSums(values) & colSums(values) <= 1+.05 ))

	p <- t( t(values) / colSums(values) )
	
	z <- p*log(p)
	z[!(p > 0.)] <- 0.
	hc <- colSums(-z)
	c <- exp(hc)
	
	q <- rowMeans(p)
	z <- p*log(q)
	z[!(p > 0.)] <- 0.
	hu <- colSums(-z)
	u <- exp(hu)
	
	aug <- as.data.frame(t(p))
	colnames(aug) <- sprintf('w%d', seq_len(ncol(aug)))
	
	results[[ length(results)+1L ]] <- data.frame( hc = hc, c = c, hu = hu, u = u, n = colSums( p > 0. ), aug )
}

source('toporbind.R')
results <- do.call(toporbind, results)

# dump
write.table( results, 'cellular_freqs.tsv', sep = '\t', row.names = T, col.names = NA)
