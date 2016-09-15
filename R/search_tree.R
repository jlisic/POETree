#' Search for k nearsest neighbors in a persistent expandable k-d tree
#' 
#' \code{search_tree} search a persistent expandable k-d tree.
#'
#' @param y A vector to search for in the k-d tree, or a matrix for multiple
#'  searches.
#' @param tree A K-d tree stored as a POETree object.
#' @param x The matrix associated with the POETree object.
#' @param k A scalar used to specify the number of neighbors to return from the 
#'   nearest neighbor search.
#' @return A vector of integers identifying nearest neighbors.
#'   
#' @examples 
#' x <- matrix(runif(20),10,2)
#' # find the four closest nearest neighbors
#' myTree <- create_tree(x)
#' y <- c(0,0)
#' neighbors <- search_tree(y,myTree,k=4)
#' @useDynLib POETree 
#' @export
search_tree <-
function(
  y,
  tree,
  x,
  k=5
) {

  if(is.matrix(y)) {
    p <- NCOL(y)
    m <- NROW(y)
  } else {
    p <- length(y)
    m <- 1
  }

  if( p > 1) {
    if( p != NCOL(x) ) stop("Dimensions do not match between x and y.") 
  } else {
    if( length(y) != NCOL(x) ) stop("Dimensions do not match between x and y.") 
  }

  y <- as.double(t(y))
  k <- as.integer(min(k,tree$tree[0,'size']))
  if( k < 1) stop( "k must be a positive integer.")

  neighbors <- as.integer(rep(-1,k*m))
  distances <- rep(Inf,k*m)

  .Call("R_POETree_search_tree", y, tree,t(x),k,p,neighbors,distances) 

  return( list(
               neighbors=matrix((neighbors+1),byrow=T,ncol=k), 
               distances=matrix(distances,byrow=T,ncol=k)
               ) 
  )
}



