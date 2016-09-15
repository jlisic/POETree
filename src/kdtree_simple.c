#include "kdtree_simple.h"

// function to find the minimal Euclidian distance 
void find_knn(
  double * x,           // vector of values
  double * y,           // what needs neighbors
  int * neighbors,      // closest neighbors
  double * distances,   // distances between y and neighbors
  int k,                // number of neighbors
  int p,                // length of y
  int size,             // size 
  int * index           // indexes in x of leaf node
  ) {

  int i,j,l;
  double *  z;
  double currentDist,tmp;
  

  for( i=0; i < size; i++) {

    if( index[i] < 0 ) continue;

    // get first index 
    z = x + index[i]*p; 

    // calculate L2 distance
    for( j=0, currentDist=0; j < p; j++) { 
      tmp = z[j] - y[j];
      tmp *= tmp;
      currentDist += tmp;
    }
    
    // if smaller than prior index update 
    // dist is ordered from smallest to biggest
    if( currentDist < distances[k-1] ) {

      for(l=k-2;l>=0;l--) { 
        // if the distance is bigger than current dist we will shift things down.
        if( currentDist < distances[l] ) {
          distances[l+1] = distances[l];
          neighbors[l+1] = neighbors[l];
        }
        else break;
      }
   
      distances[l+1] = currentDist;
      neighbors[l+1] = index[i];
    } 
  }

 return;
}




/* the kdtree consists of
 *
 * int    dim
 * int    index
 * int    size
 * int    left
 * int    right
 * double median 
 */ 
void search_tree(
      double * x,         // contents of tree
      double * y,         // what we are looking for
      int * x_dim,
      int * x_index,   // SEXP list to store final indexes
      int * x_size,    // length Q set of sizes from row  (Ptr)
      int * x_left,    // length Q set of lefts   (Ptr)
      int * x_right,   // length Q set of rights  (Ptr)
      double * x_median,  // length Q set of medians (Ptr)
      SEXP R_index,    // tree data
      int p,
      int * neighbors,   // neighbors
      double * distances, // distances for neighbors
      int k,           // number of neighbors
      int row           // current row
      ) {

  double distMin;

  // is there anything here?
  if (x_index[row] > -1) {
    find_knn( x, y, neighbors, distances, k, p, x_size[row], INTEGER(VECTOR_ELT(R_index,x_index[row])));
    return;
  }

  // if not let's look for leaves
  if( y[x_dim[row]] <= x_median[row] ) {

    // go to left
    search_tree(x,y,x_dim,x_index,x_size,x_left,x_right,x_median,R_index,p,neighbors,distances,k,x_left[row]);
  
    distMin = y[x_dim[row]] - x_median[row];
    distMin *= distMin;

    // check distance and go to right
    if(distMin < distances[k-1]) 
      search_tree(x,y,x_dim,x_index,x_size,x_left,x_right,x_median,R_index,p,neighbors,distances,k,x_right[row]);
  } else {

    // go to right 
    search_tree(x,y,x_dim,x_index,x_size,x_left,x_right,x_median,R_index,p,neighbors,distances,k,x_right[row]);
    
    distMin = y[x_dim[row]] - x_median[row];
    distMin *= distMin;
    // check distance and go to left 
    if(distMin < distances[k-1]) 
      search_tree(x,y,x_dim,x_index,x_size,x_left,x_right,x_median,R_index,p,neighbors,distances,k,x_left[row]);
  }  

  return;
}



/*************************************************************************************
Function to buildTree

*************************************************************************************/
int buildTree(
    SEXP x,     // the data in column major form from R
    SEXP x_dim,      // current dim
    int * index,  // length n set of indexes
    SEXP R_index,
    SEXP leaf_index, 
    SEXP x_size,   // length n set of sizes from row 
    SEXP x_left,   // length Q set of lefts 
    SEXP x_right,  // length Q set of rights
    SEXP x_median, // length Q set of medians
    int d,          // current dimension
    int n,          // total number of observations
    int p,          // length of each observation
    int m,          // max leaf node size
    int * q,       // current index for left/right/median etc... q in {1,...,Q}
    int * s
  ) {

  int i;
  int current_q = *q;
  int current_size = INTEGER(x_size)[current_q];
  int left_size;
  int right_size;
  int * left_index;
  int * right_index;
  SEXP current_leaf_index;

  /********** Check if we are done *********/
  if( current_size <= m ) {
    // record where we are storing the result
    
    INTEGER(R_index)[current_q] = *s;

    // record the leaf index 
    current_leaf_index = VECTOR_ELT(leaf_index,*s);
    for(i=0; i <current_size; i++) INTEGER(current_leaf_index)[i] = index[i]; 

    // free index
    free(index);
    
    //increment the storage pointer
    *s = *s + 1;

    return current_q; 
  }

  /********** Calculate the Median *********/
  double * y = calloc( current_size, sizeof(double) ); 
  double ** yPtr = calloc( current_size, sizeof(double * )); 

  // copy over the dim of interest to temporary arrays
  for(i=0;i<current_size;i++) {
    y[i] = REAL(x)[ n*d + index[i] ];  // dim of interest
    yPtr[i] = &(y[i]);                     // pointers to dim of interest
  }
  
  // get the size of left and right 
  left_size = current_size/2; 
  right_size = current_size - left_size;
  
  left_index = calloc( left_size, sizeof(int));
  right_index = calloc( right_size, sizeof(int));

  // save median
  REAL(x_median)[current_q] = quantile_quickSelectIndex(yPtr, left_size, current_size);


  // copy over the indexes on the left and right side of the split 
  for(i = 0; i < left_size; i++) left_index[i] = index[(int)(yPtr[i] - y)];

  for(i = left_size; i < current_size; i++) right_index[i-left_size] = index[(int)(yPtr[i] - y)];
  free(index);
 
  // free our temp arrays 
  free(yPtr);
  free(y);
 
  /****** save and update dim ******/

  // save dim
  INTEGER(x_dim)[current_q] = d;
  // change d
  d = (d +1) % p;

  // increment q 
  *q = *q + 1; 
 
  // to the left; 
  INTEGER(x_size)[*q] = left_size;
  INTEGER(x_left)[current_q] = buildTree(x, x_dim, left_index, R_index, leaf_index, x_size, x_left, x_right, 
      x_median, d, n, p, m, q, s);
  
  // increment q 
  *q = *q + 1; 
  
  // to the right
  INTEGER(x_size)[*q] = right_size;
  INTEGER(x_right)[current_q] = buildTree(x, x_dim, right_index, R_index, leaf_index, x_size, x_left, x_right, 
      x_median, d, n, p, m, q,s);

  return current_q;
} 





SEXP R_POETree_build_tree(
  SEXP R_x,
  SEXP R_x_dim,
  SEXP R_x_index,
  SEXP R_x_leaf_index,
  SEXP R_x_size,
  SEXP R_x_left,
  SEXP R_x_right,
  SEXP R_x_median, 
  SEXP R_n,
  SEXP R_p,
  SEXP R_m
  ) {
  

  int q = 0;
  int s = 0;

  int i;

  int n = INTEGER(R_n)[0];
  int p = INTEGER(R_p)[0];
  int m = INTEGER(R_m)[0];
  
  // this is freed in build tree
  int * index = (int *) calloc(n, sizeof(int));
  for(i=0; i < n;i++) index[i]= i;

  buildTree(
    R_x,            // the data in column major form from R (Ptr)
    R_x_dim,        // current dim (Ptr)
    index,          // length n set of indexes
    R_x_index,      // SEXP list to store final indexes
    R_x_leaf_index, // SEXP list to store final indexes
    R_x_size,       // length Q set of sizes from row  (Ptr)
    R_x_left,       // length Q set of lefts   (Ptr)
    R_x_right,      // length Q set of rights  (Ptr)
    R_x_median,     // length Q set of medians (Ptr)
    0,              // current dimension
    n,              // total number of observations
    p,              // length of each observation
    m,              // max leaf node size
    &q,             // current index for left/right/median etc... q in {1,...,Q} (Ptr)
    &s              // current index for leaf nodes, s in {1,...,S} (Ptr)  
    );

  return(R_NilValue);
} 





/* this is a wrapper function for the kd tree search */
/* this function returns a list of k nearest neighbors */
SEXP R_POETree_search_tree(
  SEXP R_y,
  SEXP R_tree,
  SEXP R_x,
  SEXP R_k,
  SEXP R_p,
  SEXP R_neighbors,
  SEXP R_distances
  ) {
    
  int i; 
  SEXP tree_names =  getAttrib(R_tree, R_NamesSymbol); // get the tree list names
  SEXP R_tree_structure = NULL;  // pointer to the tree structure
  SEXP R_index = NULL;           // pointer to tree contents 
  int k = INTEGER(R_k)[0];       // number of points to look for
  int p = INTEGER(R_p)[0];       // dim of each vector to look for
  int n = length(R_tree);        // number of elements in the tree
  int * x_dim;                   // vector of dims used to build tree
  int * x_index;                 // vector of indexes to the tree content list
  int * x_size;                  // vector of tree sizes 
  int * x_left;                  // vector of left child indexes
  int * x_right;                 // vector of right child indexes
  double * x_median;             // vector of medians (spliting points) 
  int * neighbors = INTEGER(R_neighbors);  // indexes of neighbors, initialized to -1
  double * distances = REAL(R_distances);  // distances to the neighbors, initialized to Inf
  double * y = REAL(R_y);                  // what we are searching for nn for
  int  m = length(R_y)/p;                  // number of y's to search for
  double * x = REAL(R_x);                    // the transpose of R_x

  // get index and tree structure
  for( i = 0; i <n; i++) {
    // get tree structure 
    if(strcmp(CHAR(STRING_ELT(tree_names, i)), "tree") == 0) {
      R_tree_structure = VECTOR_ELT(R_tree,i);
      //printf("found tree\n");
    // get the R_index
    } else if(strcmp(CHAR(STRING_ELT(tree_names, i)), "index") == 0) {
      R_index = VECTOR_ELT(R_tree,i);
      //printf("found index\n");
    }
  }

  // double check that this is a real tree!
  if( R_tree_structure == NULL || R_index == NULL ) {
    return(R_NilValue);
  }
  
  /* convert R_tree_structure back to C types */  
  x_dim    = INTEGER(VECTOR_ELT(R_tree_structure,0));
  x_size   = INTEGER(VECTOR_ELT(R_tree_structure,1));
  x_left   = INTEGER(VECTOR_ELT(R_tree_structure,2));
  x_right  = INTEGER(VECTOR_ELT(R_tree_structure,3));
  x_median = REAL(VECTOR_ELT(R_tree_structure,4)); 
  x_index  = INTEGER(VECTOR_ELT(R_tree_structure,5));


  for( i=0; i<m; i++) {

    
    search_tree(
      x,
      y + i*p,   // what we are looking for
      x_dim,
      x_index,   // SEXP list to store final indexes
      x_size,    // length Q set of sizes from row  (Ptr)
      x_left,    // length Q set of lefts   (Ptr)
      x_right,   // length Q set of rights  (Ptr)
      x_median,  // length Q set of medians (Ptr)
      R_index,   // tree data
      p,
      neighbors, // neighbors
      distances, // distances
      k,         // number of neighbors
      0
    );
      
    neighbors = neighbors + k;
    distances = distances + k;
  }

  return(R_NilValue);
} 




// test function
SEXP R_testsexp( SEXP x ) {

  int i,j,m;
  int n = length(x);
  SEXP y;

  printf("type = %d length = %d\n", TYPEOF(x), n );


  if( TYPEOF(x) == 19 ) {

    for( i=0; i < n; i ++) { 
      y=VECTOR_ELT(x,i);
      printf("i = %d, type = %d, length = %d\n", i, TYPEOF(y), length(y)  );
      m = length(y);
      for(j = 0; j < m; j++) INTEGER(y)[j] = 12; 
    }

  }

  return(R_NilValue);
}



