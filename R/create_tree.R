#' Create a persistent expandable k-d tree
#' 
#' \code{create_tree} build a persistent expandable k-d tree from a vector or 
#' matrix. 
#'
#' @param x n vectors stored in an n by p matrix.  These vectors will be used
#'   to create the K-d tree.  It is assumed that all values are finite and 
#'   observed.
#' @param leafSize A scalar used to specify the number of points to store in the 
#'   leaf nodes. 
#' @param leafBuffer A scalar used to specify the number of empty index locations
#'   to retrain in the leaf nodes.  This is useful to allow for quicker insertions
#'   of new data points.  Note that empty points are bounded by leaf nodes, 
#'   therefore \code{leafBuffer} must be smallter than \code{leafSize}.
#' @return An S3 object cotaining the persistent tree.
#'   
#' @examples 
#' x <- matrix(runif(20),10,2)
#' myTree <- create_tree(x)
#' @useDynLib POETree 
#' @export
create_tree <-
function(
  x,
  leafSize=11,
  leafBuffer=0  
) {

  if(leafBuffer >= leafSize) stop("Invalid leafBuffer value.")

  # get data size
  n <- NROW(x)
  p <- NCOL(x) 

  # calculate the number of nodes
  Q <- 0 
  i <- 0
  while( 2^i*(leafSize -leafBuffer) < n ) {
    Q <- Q + 2^i 
    i = i + 1
  }
  Q <- as.integer(Q + 2^i)
  S <- as.integer(2^i)
  leafSize <- as.integer(leafSize)
  leafBuffer <- as.integer(leafBuffer)

  x_index  <- as.integer(rep(-1,Q)) 
  x_dim    <- as.integer(rep(0,Q))
  x_left   <- as.integer(rep(-1,Q)) 
  x_right  <- as.integer(rep(-1,Q))
  x_median <- as.double(rep(0,Q))
  x_size   <- as.integer(rep(0,Q))
  x_size[1] <- n

  x_leaf_index <- list()
  for( i in 1:S) x_leaf_index[[i]] <- as.integer(rep(-1,leafSize)) 
 
  .Call("R_POETree_build_tree",
    x,            # the data in column major form from R
    x_dim,        # dim in the tree
    x_index,      # final index for x 
    x_leaf_index, # the index tree  
    x_size,       # sum of subsequent nodes 
    x_left,       # left node index
    x_right,      # right node index
    x_median,     # splitting point
    n,            # number of obs
    p,            # dim of obs
    leafSize-leafBuffer   # leaf size
  ) 


  treeStructure <- data.frame(
    dim=x_dim,
    size=x_size,
    left=x_left,
    right=x_right,
    median=x_median,
    index=x_index)

  y <-  list( tree=treeStructure, index=x_leaf_index, type="k-d")
  attr(y,"class") <- "POETree"
  return(y)
}
