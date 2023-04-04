
#' Merge data frame with topological variable sorting
#'
#' @param arg1 First data frame.
#' @param ... Rest of the data frames.
#' @param fill Value to use as a filler for missing columns. 
#'
#' @return Merged data frame.
#'
toporbind <- function(df1, ..., fill = NA) {
	# grab frames
	dfs <- list( df1, ... )

	# build rhs graph
	graph <- NULL
	for (i in seq_along( dfs )) {
		prev.name <- NULL
		for (name in colnames( dfs[[i]] )) {
			if (is.null( graph[[name]] ))
				graph[name] <- list(list( parent.count = 0L, children = c() ))

			if (!is.null( prev.name )) {
				graph[[name]]$parent.count <- graph[[name]]$parent.count + 1L
				graph[[prev.name]]$children <- c( graph[[prev.name]]$children, name )
			}

			prev.name <- name
		}
	}

	# sort that topologically
	queue <- NULL
	queue.p <- 0L
	for (name in names(graph))
		if (graph[[name]]$parent.count == 0L)
			queue <- c( queue, name )
	while (queue.p < length(queue)) {
		queue.p <- queue.p + 1L
		name <- queue[queue.p]

		for (next.name in graph[[name]]$children) {
			graph[[next.name]]$parent.count <- graph[[next.name]]$parent.count - 1L
			if (!(graph[[next.name]]$parent.count > 0L))
				queue <- c( queue, next.name )
		}

		graph[[name]] <- NULL
	}
	if (length(graph) > 0L)
		stop('topological order does not exist')

	# fill 
	for (i in seq_along(dfs))
		for (name in setdiff( queue, colnames(dfs[[i]]) ))
			dfs[[i]][, name] <- rep( fill, nrow(dfs[[i]]) )

	# permute
	for (i in seq_along(dfs)) {
		perm <- order(match( colnames(dfs[[i]]), queue ))
		dfs[[i]] <- dfs[[i]][, perm]
	}

	# merge
	return (do.call(rbind, dfs))
}
