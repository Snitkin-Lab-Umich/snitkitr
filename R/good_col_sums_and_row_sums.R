#' Good Column Sums - wont throw an error if you input a vector!
#' @description like base R colSums() but allows you to input a vector.
#' Thus if you enter a vector of a single row, col_sums will return that vector 
#' (instead of throwing an error).
#' If you enter a vector of a single column, col_sums will sum() that vector 
#' (instead of throwing an error)
#' 
#' @param vector_or_mat - a vector or a matrix 
#' @param row_or_col - a flag indicating if you are entering a row (1) or a column (2)
#'
#' @return column sums
#' @export
#'
#' @examples
col_sums <- function(vector_or_mat, row_or_col){
  
  if (row_or_col == 1){ # inputting a row
    
    if(is.null(dim(vector_or_mat))){
      vector_or_mat
    }else{
      colSums(vector_or_mat)
    }
    
  }else if (row_or_col == 2){ #inputting a column
    
    if(is.null(dim(vector_or_mat))){
      sum(vector_or_mat)
    }else{
      colSums(vector_or_mat)
    }
    
  }else {
    'Error - row_or_col must have value of 1 or 2'
  }
  
  
  
}


#' Good Row Sums - wont throw an error if you input a vector!
#'
#' @param vector_or_mat - a vector or a matrix 
#' @param row_or_col - a flag indicating if you are entering a row (1) or a column (2)
#' @description like base R rowSums() but allows you to input a vector.
#' Thus if you enter a vector of a single column, col_sums will return that vector 
#' (instead of throwing an error).
#' If you enter a vector of a single row, col_sums will sum() that vector 
#' (instead of throwing an error)
#' 
#' @return row sums 
#' @export
#'
#' @examples
row_sums <- function(vector_or_mat, row_or_col){
  if (row_or_col == 1){ #inputting a row
    
    if(is.null(dim(vector_or_mat))){
      sum(vector_or_mat)
    }else if(nrow(vector_or_mat) == 0){ # should add to everything - if its 2 dimensions but empty
      NA
    }else{
      rowSums(vector_or_mat)
    }
    
  }else if(row_or_col == 2){ #inputting a column
    
    if(is.null(dim(vector_or_mat))){
      vector_or_mat
    }else{
      rowSums(vector_or_mat)
    }
    
  }else{
    
    'Error - row_or_col must have value of 1 or 2'
  
  }
  
}
