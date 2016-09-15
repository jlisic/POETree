#' test vector
#' 
#' \code{test_vector} this is a test.
#'
#' @param y input.
#' @return whatever.
#'   
#' @export 
test_vec <-
function(
  y
) {

  .Call("R_testsexp", y)
  
}
