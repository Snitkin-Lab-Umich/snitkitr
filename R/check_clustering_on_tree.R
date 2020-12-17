# check clustering of outcome and source

#' Reverse list structure
#'
#' @param ls list you want to reverse
#'
#' @return reversed list
#' @export
#'
#' @examples
#' #Reference with example: https://stackoverflow.com/questions/15263146/revert-list-structure
reverse_list_str <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- future.apply::future_lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  future.apply::future_apply(do.call(rbind, x), 2, as.list) 
}


#' Get largest pure subtrees
#'
#' @param subtrs Subtrees created using ape::subtrees to look for clustering on. Should include all isolates of interest.
#' @param isolate_labels Named vector of labels by which pure clusters are defined. Names must be equivalent to tree tip label names.
#' @param control_labels Named vector of labels known to cluster. Names must be equivalent to tree tip label names. This controls for clustering by requiring that the pure clusters must contain multiple of the control labels. 
#' @param bootstrap Bootstrap support to use to filter unconfident tree edges (keeps > bootstrap; NULL = keep all; default: 90).
#'
#' @return list containing the largest pure subtree that each isolate belongs to, the index of that subtree, and the edges in that subtree.
#' @export
#'
#' @examples
get_largest_subtree <- function(subtrs, isolate_labels, control_labels=NULL, bootstrap = 90){
  largest_st_info = future.apply::future_lapply(names(isolate_labels), function(i){
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. 
    #CLUSTERS ARE DEFINED AS:
    # 1) HAVE ONLY A SINGLE EPI LABEL, 2) HAVE BOOTSTRAP SUPPORT GREATER THAN 90, 3) INCLUDE MORE THAN ONE CONTROL LABEL
    sts = future.apply::future_sapply(subtrs, FUN = function(st){ 
      i_in_subtree = sum(grepl(i, st$tip.label, perl = TRUE)) > 0 # isolate is in subtree
      one_label = length(unique(isolate_labels[intersect(st$tip.label, names(isolate_labels))])) == 1 # only one label in subtree
      good_bootstrap = rep(TRUE, length(st$node.label[[1]]))
      if(!is.null(bootstrap)){
        good_bootstrap = !is.na(as.numeric(st$node.label[[1]])) && as.numeric(st$node.label[[1]]) > bootstrap # bootstrap support > 90
      }
      multiple_control = ifelse(is.null(control_labels), TRUE, # always true if not controlling for another variable
                                length(unique(control_labels[intersect(st$tip.label, names(control_labels))])) > 1) # more than one control label in subtree
      if(i_in_subtree && one_label && good_bootstrap && multiple_control){ 
        length(intersect(names(isolate_labels), st$tip.label))
      }else{
        0
      }
    })
    # GET THE LARGEST SUBTREE
    largest_st = max(sts)
    # GET THE INDEX OF THE LARGEST SUBTREE
    largest_st_i = which.max(sts)
    # GET EDGES BELONGING TO SUBTREES
    largest_st_edges = ape::which.edge(subtrs[[1]], subtrs[[largest_st_i]]$tip.label)
    return(list(largest_st=largest_st,
                largest_st_i=largest_st_i,
                largest_st_edges=largest_st_edges))
  })# end future_apply
  # reverse lists
  largest_st_info = reverse_list_str(largest_st_info)
  return(largest_st_info)
}#end get_largest_subtree

#' Check for clustering on tree
#' @description This function is used to test for clustering of a certain epi factor on a phylogenetic tree. 
#'
#' @param tree Tree to look for clustering on.
#' @param isolate_labels Named vector of labels by which pure clusters are defined. Names must be equivalent to tree tip label names.
#' @param nperm Number of permutations to perform.
#' @param control_labels Named vector of labels known to cluster. Names must be equivalent to tree tip label names. This controls for clustering by requiring that the pure clusters must contain multiple of the control labels. 
#' @param plot_path Path to output histogram. No plot created if NA. 
#'
#' @return Vector of real (1st element) and permuted (all other elements) pure cluster sizes.
#' @export
#'
#' @examples
#' tree = ape::rtree(50)
#' labs = sample(c(0,1),50,replace=TRUE)
#' names(labs) = tree$tip.label
#' check_tree_clustering(tree,labs,nperm=10)
check_tree_clustering = function(tree,isolate_labels, nperm = 1000, control_labels=NULL, plot_path='cluster.pdf'){
  subtrs = ape::subtrees(tree)
  # find largest subtrees
  if(is.null(control_labels)){
    largest_subtr <- get_largest_subtree(subtrs, isolate_labels)
  }else{
    largest_subtr <- get_largest_subtree(subtrs, isolate_labels, control_labels=control_labels)
  }
  #GET CLUSTERS FOR RANDOMIZED LABELS
  in_cluster = numeric(nperm + 1)
  in_cluster[1] = sum(largest_subtr[[2]] > 1 & largest_subtr[[1]] > 1)
  rand_labels = isolate_labels
  in_cluster_perm = future.apply::future_sapply(1:nperm, function(r){
    names(rand_labels) = sample(names(isolate_labels))
    if(is.null(control_labels)){
      largest_st_rand_list <- get_largest_subtree(subtrs, rand_labels)
    }else{
      largest_st_rand_list <- get_largest_subtree_loc(subtrs, rand_labels, control_labels=control_labels)
    }
    sum(largest_st_rand_list[[2]] > 1 & largest_st_rand_list[[1]] > 1)    
  })
  
  in_cluster[2:(nperm+1)] = in_cluster_perm
  
  if(!is.na(plot_path)){
    
    # get p value
    p = (1 + sum(in_cluster[2: length(in_cluster)] >= in_cluster[1]))/ (1 + length(in_cluster)-1)
    
    # plot histogram
    pdf(plot_path, height = 8 , width = 8)
    h = hist(in_cluster[2:length(in_cluster)], 20, xlim = c(0, max(in_cluster + 2)), 
             main = paste('Null distribution of number of isolates in pure clusters ', '\n', "(p = ", round(p,3), ")", sep = ""), 
             xlab = "Number of isolates in a pure cluster", col = "lightgray")
    par(new = TRUE)
    points(in_cluster[1],0, col = "black", bg = "red", cex = 2, pch = 23, 
           xlim = c(0, max(in_cluster + 2)), ylim = c(0,max(h$counts) + 5))
    dev.off()
  }
  
  return(in_cluster)
}






