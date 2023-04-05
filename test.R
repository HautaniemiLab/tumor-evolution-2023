#' author: Antti Häkkinen
#!/usr/bin/env Rscript

set.seed(123L)  

dat <- read.table('path to cellular frequencies from discovery set',
	sep = '\t', header = T, row.names = 1L, stringsAsFactors = F)

meta <- read.table('path to metadata including clinical data if is needed',
	sep = ';', header = T, stringsAsFactors = F)

# excluding cell lines and other irrelevant samples
keep <- !grepl( '_CL\\d*', rownames(dat) ) &
	      !grepl( '_sph', rownames(dat) )

dat <- dat[keep, ]

dat$patient <- sub( '^(\\w+\\d+)_([pir])(.*)', '\\1', rownames(dat) )
dat$phase <- sub( '^(\\w+\\d+)_([pir])(.*)', '\\2', rownames(dat) )

# considering primary (diagnostic) samples only
	use.phases <- c('p')
#	use.phases <- 'i'
#	use.phases <- c('i', 'r')
#	use.phases <- 'r'

dat$hup <- rep( NA, nrow(dat) )
for (pat in unique(dat$patient)) {
	msk <- dat$patient == pat & dat$phase %in% use.phases
	V <- as.matrix( dat[ msk, grepl('^(w)(\\d+)$', colnames(dat)), drop = F ] )
	V <- V[, !is.na(colSums(V)), drop = F ]
	q <- colMeans(V)
	z <- t( t(V) * log(q) )
	z[!(V > 0.)] <- 0.
	dat[ msk, ]$hup <- -rowSums(z)
	if (nrow(V) == 1)
		dat[ msk, ]$hup <- NaN   # don't know, single sample
}

dat$dup <- dat$hup - dat$hc

avg.by.group <- function(x, group)
	return (rowsum( x, group) / rowsum( +!is.na(x), group))

	msk <- dat$phase %in% use.phases
avg <- data.frame( hcp = avg.by.group( dat$hc[msk], dat$patient[msk] ),
	hup = avg.by.group( dat$hu[msk], dat$patient[msk] ),
	dup = avg.by.group( dat$dup[msk], dat$patient[msk] ) )
avg$patient <- rownames(avg)

## new data
# X = avg$hcp
# W = avg$dup

colors <- c( 'no coexistence' = 'orange', 'competing' = 'red',
	'maintaining' = 'green', 'adaptive' = 'blue' )

	ig.pats <- c('H086')
	ex.pats <- c('H098', 'H152', 'H147', 'H233', 'H064') 

source('kmeansw.R')
X <- cbind( avg$hcp, avg$dup )
rownames(X) <- rownames( avg )
W <- 0*X + 1.
W[ avg$patient %in% ig.pats ] <- 1e-10
clu.colors <- c('blue', 'darkgreen', 'orange')
clu <- kmeansw( X, W, length(clu.colors), nstart = 1e3 )
X1 <- clu$Y[ clu$L, ]
X1[!is.na(X)] <- X[!is.na(X)]

clu.names <- rownames(cc)[ assign(cc) ]
clu.colors <- colors[clu.names]


