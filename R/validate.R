#' Tests if an input object has the specified class.
#'
#' @param obj Any R object.
#' @param current_class String. Name of the expected class of the R object.
#'
#' @return is_this_class: Logical.
#' @noRd
#' @export
is_this_class <- function(obj, current_class){
  if (length(current_class) != 1) {
    stop("Current_class must have a length of 1")
  }
  if (!methods::is(current_class, "character")) {
    stop("Current_class is expected to be a string describing a class")
  }
  r_classes <- c("character",
                 "numeric",
                 "integer",
                 "logical",
                 "complex",
                 "phylo",
                 "DNAbin",
                 "phyDat",
                 "matrix",
                 "data.frame",
                 "factor",
                 "vcfR")
  if (!(current_class %in% r_classes)) {
    stop("current_class is expected to be a R class")
  }

  is_this_class <- methods::is(obj, current_class)
  return(is_this_class)
}

#' Checks if an object is of the expected R class.
#'
#' Doesn't return anything. Gives error if the object is not of the expected R
#' class.
#'
#' @param obj Any R object.
#' @param current_class Character string. Name of R class
#' @noRd
#' @export
check_is_this_class <- function(obj, current_class){
  class_log <- is_this_class(obj, current_class)
  if (class_log != TRUE) {
    stop(paste("Input must be a", current_class))
  }
}

#' Check that the input tree is actually a 'phylo' object.
#'
#' Doesn't return anything. Gives an error message if the object is not a
#' 'phylo' object.
#'
#' @param tree Phylogenetic tree.
#' @noRd
#' @export
check_is_tree <- function(tree){
  if (!is_this_class(tree, "phylo")) {
    stop("Input requires either a path to a tree file or an ape phylo object")
  }
}

#' Check that the tree has a root.
#'
#' Doesn't return anything. Gives an error message if the object is not a
#' rooted tree.
#'
#' @param tree Phylogenetic tree.
#' @noRd
#' @export
check_tree_is_rooted <- function(tree){
  check_is_tree(tree)
  if (!ape::is.rooted(tree)) {
    stop("Tree must be rooted.")
  }
}

#' Confirm that the tree and variant matrix contain exactly the same samples
#'
#' Doesn't return anything. Gives error if the two inputs do not match.
#'
#' @param tip_labels Character. Vector of tree$tip.labels.
#' @param colnames_mat Character. Vector of column names from variant matrix.
#' @noRd
#' @export
check_setequal_tree_mat <- function(tip_labels, colnames_mat){
  if (!setequal(tip_labels, colnames_mat)) {
    stop("Tree and variant matrix sample names do not match.")
  }
}
