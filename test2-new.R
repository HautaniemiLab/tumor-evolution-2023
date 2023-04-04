#!/usr/bin/env Rscript

	set.seed(123L)  # deterministic clustering?

dat <- read.table('cellular_freqs-old.tsv',
	sep = '\t', header = T, row.names = 1L, stringsAsFactors = F)

meta <- read.table('Coexistence_analysis_from_trees.csv',
	sep = ';', header = T, stringsAsFactors = F)

keep <- !grepl( '_CL\\d*', rownames(dat) ) &
	      !grepl( '_sph', rownames(dat) )

dat <- dat[keep, ]

dat$patient <- sub( '^(\\w+\\d+)_([pir])(.*)', '\\1', rownames(dat) )
dat$phase <- sub( '^(\\w+\\d+)_([pir])(.*)', '\\2', rownames(dat) )

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

cc <- table( meta$COEXISTENCE_UNTREATED[ match(names(clu$L), meta$ID_old) ], clu$L )
		assign <- function(CC) {
			stopifnot( nrow(CC) >= ncol(CC) )
			perm <- rep( 0L, ncol(CC) )
			used <- rep( F, nrow(CC) )
			for (k in order(CC, decreasing = T)) {
				i <- ( (k-1) %%  nrow(CC) ) + 1L
				j <- ( (k-1) %/% nrow(CC) ) + 1L
				if (!used[i] && perm[j] == 0L) {
					perm[j] <- i
					used[i] <- T
				}
			}
			return (perm)
		}
 clu.names <- rownames(cc)[ assign(cc) ]
 clu.colors <- colors[clu.names]

 	dump <- meta
	if (file.exists('dump-new/dump.tsv'))
		dump <- read.table('dump-new/dump.tsv', sep = ';', header = T, stringsAsFactors = F)
	{
		idxs <- match( dump$ID_old, avg$pat )
		 nam <- paste( use.phases, collapse = '' )
print(nam)
		dump[sprintf('log.C.%s', nam)] <- avg$hcp[idxs]
		dump[sprintf('log.D.%s', nam)] <- avg$dup[idxs]
		dump[sprintf('clu.%s', nam)] <- clu.names[clu$L][idxs]
	write.table(dump, 'dump-new/dump.tsv', sep = ';', row.names = F, col.names = T)
	}

col <- clu.colors[clu$L]
 col[ avg$patient %in% ig.pats ] <- 'black'
pch <- rep('o', nrow(X1))
 pch[ rowSums( is.na(X) ) > 0 ] <- '*'

lim <- c( 1, 5.5 )
plot( exp(X1[,1]), exp(X1[,2]), log = 'xy', col = col, pch = pch,
	xlim = lim, ylim = lim, xlab = 'eff clonality', ylab = 'eff clonal divergence' )
	
	msk <- avg$patient %in% ex.pats
	points( exp(X1[,1])[msk], exp(X1[,2])[msk], col = col[msk], pch = 'x', cex = .75 )

	#points( exp(clu$Y[,1]), exp(clu$Y[,2]), col = clu.colors, pch = '+' )

	jx <- 0.025
	jy <- 0.025
text( exp( X1[,1]+jx ), exp( X1[,2]+jy ), avg$patient,
	col = col, cex = 0.5, srt = +45 )


