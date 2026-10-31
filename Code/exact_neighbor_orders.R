# Binary neighbor shells defined by exact shortest-path distance.
# In particular, second neighbors exclude the focal sample and first neighbors.
exact_neighbor_orders <- function(adjacency, max_order) {
  stopifnot(is.matrix(adjacency), nrow(adjacency)==ncol(adjacency),
            all(is.finite(adjacency)), all(adjacency %in% c(0,1)),
            isSymmetric(adjacency), all(diag(adjacency)==0),
            max_order>=1, max_order==as.integer(max_order))
  n <- nrow(adjacency)
  seen <- diag(TRUE,n)
  frontier <- diag(1,n)
  shells <- vector('list',max_order)
  for(k in seq_len(max_order)) {
    frontier <- (frontier %*% adjacency > 0) & !seen
    diag(frontier) <- FALSE
    dimnames(frontier) <- dimnames(adjacency)
    shells[[k]] <- 1L*frontier
    seen <- seen | frontier
    stopifnot(isSymmetric(shells[[k]]), all(diag(shells[[k]])==0))
  }
  stopifnot(all(Reduce('+',shells)<=1))
  shells
}
