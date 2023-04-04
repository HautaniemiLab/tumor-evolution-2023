##############33 clustering for new data -- validation set

## ------------ run the beginning of clustering code for new data
set.seed(123L)  # deterministic clustering?

dat <- read.table('cellular_freqs-validation.tsv',
                  sep = '\t', header = T, row.names = 1L, stringsAsFactors = F)

meta <- read.table('Coexistence_analysis_from_trees.csv',
                   sep = ';', header = T, stringsAsFactors = F)

keep <- !grepl( '_CL\\d*', rownames(dat) ) &
  !grepl( '_sph', rownames(dat) )

dat <- dat[keep, ]

dat$patient <- sub( '^(\\w+\\d+)_([piro])(.*)', '\\1', rownames(dat) )
dat$phase <- sub( '^(\\w+\\d+)_([piro])(.*)', '\\2', rownames(dat) )

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
#avg$cp <- exp(avg$hcp)
#avg$dp <- exp(avg$dup)

#### ----- load the parameters from old clustering 

load('clu.p.RData')

#### ----  get X and W for new validation data 

new.X <- cbind( avg$hcp, avg$dup )
rownames(new.X) <- rownames( avg )

ig.pats <- c('H086')
new.W <- 0*new.X + 1.
new.W[ avg$patient %in% ig.pats ] <- 1e-10

old.clu <- clu

source('kmeansw_new.R')
new.clu <- kmeansw.relabel( old.clu, new.X, new.W )

names <- clu$names[new.clu$L]


result <- data.frame( patient = names(new.clu$L), label =new.clu$L,  cluster <- names )
colnames(result) <- c('patient', 'label', 'cluster')


write.table(result, file = 'clustering-relapse.csv', sep = '\t', row.names =T, col.names =NA)


##### write dump file 
dump <- avg
nam <- paste( use.phases, collapse = '' )

dump[sprintf('log.C.%s', nam)] <- avg$hcp
dump[sprintf('log.D.%s', nam)] <- avg$dup
dump[sprintf('clu.%s', nam)] <- clu$names[new.clu$L]
dump <- dump[,-c(1:4)]
write.table(dump, 'dump-new/dump-p.tsv', sep = ';', row.names = T, col.names = NA)


