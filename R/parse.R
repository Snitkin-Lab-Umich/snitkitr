#' Remove variant matrix rows with bugs in the annotations
#' @description Removes rows from matrix in which row.names have bugs. Bugs include
#' row annotations with 1) warnings, 2) incorrect number of pipes (should be a
#' multiple of 9), 3) the string CHR_END (we want this to be the last gene in the
#' annotated genome instead of CHR_END), 4) no strand or locus tag information or 5) rows
#' annotated with "None". Eventually can get rid of this when all the bugs are fixed,
#' or can keep as a sanity check to make sure we aren't seeing these bugs in the
#' row annotation. You should run varmat_code and varmat_allele through this function separately
#' and should expect the same rows to be removed.
#'
#' @param varmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return returns a varmat (class = data.frame) with rows with bugs in the
#' annotation removed. Also writes a file called YEAR_MONTH_DATE_rows_removed_from_varmatNAME_due_to_bugs.txt
#' logging the row names of the removed rows.
#'
#' @export
#' @noRd
remove_rows_with_bugs <- function(varmat){
  # Intialize a filename to log the removed rows
  filename = paste0(Sys.Date(), '_rows_removed_from_', deparse(substitute(varmat)), 'due_to_bugs.txt')

  # 1. Remove rows with warnings in the row annotation
  rows_with_warnings = grep('WARNING', row.names(varmat))
  write(row.names(varmat)[rows_with_warnings], file = filename, append = TRUE)
  if (length(rows_with_warnings) > 0) {
    varmat = varmat[-rows_with_warnings,]
  }

  # 2. Remove rows with the incorrect number of pipes in the row annotation
  # number of |
  num_pipes = stringr::str_count(row.names(varmat), '\\|')
  table(num_pipes) # remove not intervals of 9

  num_semicolon = stringr::str_count(row.names(varmat), ';')
  table(num_semicolon)

  write(row.names(varmat)[(num_pipes/(num_semicolon - 1)) %% 9 != 0], file = filename, append = TRUE)

  # only keep rows with correct number of pipes
  varmat = varmat[(num_pipes/(num_semicolon - 1)) %% 9 == 0,] # must be a multiple of 9

  # 3. Remove rows that still have 'CHR_END' in the row annotation
  rows_with_chr_end = grep('CHR_END', row.names(varmat))
  write(row.names(varmat)[rows_with_chr_end], file = filename, append = TRUE)

  if (length(rows_with_chr_end) > 0) {
    varmat = varmat[-rows_with_chr_end, ]
  }

  # 4. Remove rows with not enough locus tag information - have to wait for Ali
  # to convert all reported gene symbols to locus_tags
  # split_annotations <- strsplit(row.names(varmat), ";")
  # sapply(split_annotations, function(split){
  #   annot = split[1]
  #   locus_tag = gsub('^.+locus_tag=','',annot) %>% gsub(' Strand .*$','',.)
  # })

  #5. Remove rows with no strand information or locus tag information in the
  # row annotations
  locus_tag = unname(sapply(row.names(varmat), function(row){
    gsub('^.+locus_tag=', '', row) %>% gsub(' Strand .*$', '', .)
  }))
  no_locus_tag_listed = grep('NULL', locus_tag)

  no_strand_info_listed = grep('No Strand Information found', row.names(varmat))

  remove_bc_lack_of_info = union(no_locus_tag_listed, no_strand_info_listed)
  write(row.names(varmat)[remove_bc_lack_of_info], file = filename, append = TRUE)
  if (length(remove_bc_lack_of_info) > 0) {
    varmat = varmat[-remove_bc_lack_of_info, ]
  }

  #6. Remove rows with "None" annotation
  remove_bc_none_annotation = grep('None', row.names(varmat))
  write(row.names(varmat)[remove_bc_none_annotation], file = filename, append = TRUE)
  if (length(remove_bc_none_annotation) > 0) {
    varmat = varmat[-remove_bc_none_annotation, ]
  }
  return(varmat)
}

#' Remove rows with no variants or that are completely masked
#' @description Removes rows that have no variants (the same allele in every
#' sample +/- N and dash) or rows that are completely masked (all Ns).
#'
#' @param varmat_code - data.frame where the rows are variants (numeric description
#' variants: numbers ranging from -4 to 3), the columns are genomes, and the
#' row.names are annotations
#' @param varmat_allele - data.frame where the rows are variants (character description
#' variants: A,C,T,G,N,-), the columns are genomes, and the row.names are annotations
#'
#' @return Returns a list with the following elements (in order):
#' 1. varmat_code (class = data.frame) with non-variant rows removed. All rows have at least
#' one sample with a variant.
#' 2. varmat_allele (class = data.frame) with non-variant rows removed. All rows have at least
#' one sample with a variant.
#' 1 and 2 should be the same dimensions and have the same row names.
#' Also writes a file called YEAR_MONTH_DATE_rows_removed_because_no_variants
#' logging the row names of the removed rows.
#' @export
remove_rows_with_no_variants_or_completely_masked <- function(varmat_code, varmat_allele){
  rows_with_one_allele_or_all_Ns_or_dashes = apply(varmat_allele, 1, function(row){
    length(unique(row)) == 1
  })

  file = paste0(Sys.Date(), '_rows_removed_because_no_variants')
  write(x = row.names(varmat_code)[rows_with_one_allele_or_all_Ns_or_dashes], file = file)

  varmat_code_rows_removed = varmat_code[!rows_with_one_allele_or_all_Ns_or_dashes,]
  varmat_allele_rows_removed = varmat_allele[!rows_with_one_allele_or_all_Ns_or_dashes,]

  return(list(varmat_code_rows_removed, varmat_allele_rows_removed))
}

#' Split rows that have multiple annotations
#' @description Rows that have X number of annotations are replicated X number
#'   of times. Multiple annotations can be due to 1) multiallelic sites which
#'   will result in an annotation for each individual allele and 2) SNPs in
#'   overlapping genes which will result in an annotation for each gene or 3) a
#'   combination of 1 and 2. The row names will be changed to have one
#'   annotation per row (that is, multiallelic SNPs are represented as biallelic
#'   sites and each snp in a gene that overlaps with another gene will be
#'   represented on a single line). The contents of the data.frame is replicated
#'   -- that is, the contents of the replicated rows are NOT changed. You should
#'   run varmat_code and varmat_allele through this function separately and
#'   should expect the data.frames to have the same dimensions and same
#'   duplicated rows.
#'
#' @param varmat - data.frame where the rows are variants, the columns are
#'   genomes, and the row.names are annotations
#' @param snp_parser_log - logical. TRUE when handling snp information. FALSE
#'   when handling indel information.
#'
#' @return Returns a list with the following elements (in order): 1.
#'   rows_with_multiple_annots_log - a logical vector with length of
#'   nrow(varmat_added) indicating which rows once had multiple annotations
#'   (that is, were split from one row into multiple rows) 2.
#'   rows_with_mult_var_allele_log -  a logical vector with length of
#'   nrow(varmat_added) indicating which rows once had multiple annotations in
#'   the form of  multiallelic sites (that is, were split from multiallelic
#'   sites to biallelic sites) 3. rows_with_overlapping_genes_log -> -  a
#'   logical vector with length of nrow(varmat_added) indicating which rows once
#'   had multiple annotations in the form of overlapping genes (that is, were
#'   split from a SNP in multiple genes to each gene being represented on a
#'   single line) 4. split_rows_flag - an integer vector indicating which rows
#'   were split from from a row with multiple annotations (For example if varmat
#'   had 4 rows: 1, 2, 3, 4 with row 2 having 3 annotations and row 4 having 2
#'   annotations, the vector would be 1 2 2 2 3 4 4). 5. varmat_added - a
#'   data.frame where the rows are variants, the columns are genomes, and the
#'   row.names are SPLIT annotations (each overlapping gene and multiallelic
#'   site represented as a single line).
#' @export
#'
split_rows_with_multiple_annots <- function(varmat, snp_parser_log){

  # Returns a vector of numbers
  num_dividers <- sapply(1:nrow(varmat), function(x) lengths(regmatches(row.names(varmat)[x], gregexpr(";[A,C,G,T]", row.names(varmat)[x]))))

  # Returns a vector of numbers
  rows_with_multiple_annotations <- c(1:nrow(varmat))[num_dividers >= 1 & stringr::str_count(row.names(varmat), '\\|') > 9]

  # Get rows with multallelic sites
  if (snp_parser_log) {
    rows_with_multi_allelic_sites = grep('^.+> [A,C,T,G],[A,C,T,G]', row.names(varmat))
  } else {
    rows_with_multi_allelic_sites = grep('^.+> [A,C,T,G]+,[A,C,T,G]+', row.names(varmat))
  }

  # Get SNVs present in overlapping genes
  split_annotations <- strsplit(row.names(varmat)[rows_with_multiple_annotations], ";")

  num_genes_per_site = sapply(split_annotations, function(annots){
    unique(sapply(2:length(annots), function(i){
      unlist(stringr::str_split(annots[i], '[|]'))[4]
    }))
  })
  rows_with_overlapping_genes = rows_with_multiple_annotations[sapply(num_genes_per_site, length) > 1]

  # Duplicate rows with multiallelic sites
  row_indices = 1:nrow(varmat)

  varmat_added = varmat[rep(row_indices, num_dividers),]

  print("I")
  print(row.names(varmat_added)[which(grepl("Coding Indel at 440983", row.names(varmat_added)))])

  # When rows are duplicated .1, .2, .3, etc are added to the end
  # (depending on how many times they were duplicated)
  # Remove to make the duplicated rows have the exact same name
  #names_of_rows = row.names(varmat_added) %>% gsub(';\\.[0-9].*$', ';', .)

  split_rows_flag = rep(row_indices, num_dividers)

  dup = unique(split_rows_flag[duplicated(split_rows_flag)]) # rows that were duplicated

  split_annotations <- strsplit(row.names(varmat_added)[split_rows_flag %in% dup], ";")

  # FIX ANNOTS OF SNP MAT ADDED - RELIES ON THE .1, .2, .3, ... etc flag
  row.names(varmat_added)[split_rows_flag %in% dup] =  sapply(split_annotations, function(r){
    # r is the vector of each list element
    last_vector_entry <- r[length(r)]
    num_vector_entry <- length(r)

    if (num_vector_entry == 3) {
      # This is the case for the first list element (No .#; .1 means the second entry)
      paste(r[1], r[2], sep = ';')
    } else if (num_vector_entry > 3 & !grepl(pattern = "[.][0-9]+", last_vector_entry)) {
      # Neither the first list element nor a .1, .2, .3, etc...
      paste(r[1], r[2], sep = ';')
    } else {
      # .1, .2, .3, etc....
      index = as.numeric(gsub('\\.','', last_vector_entry))
      paste(r[1], r[index + 2], sep = ';')
    }
  })
  print("II")
  print(row.names(varmat_added)[which(grepl("Coding Indel at 440983", row.names(varmat_added)))])


  rows_with_multiple_annots_log = split_rows_flag %in% rows_with_multiple_annotations
  rows_with_mult_var_allele_log = split_rows_flag %in% rows_with_multi_allelic_sites
  rows_with_overlapping_genes_log = split_rows_flag %in% rows_with_overlapping_genes

  # FIX ANNOTS OF SNP MAT ADDED - ROWS WITH MULT VAR ALLELE
  if (snp_parser_log) {
    row.names(varmat_added)[rows_with_mult_var_allele_log] = sapply(row.names(varmat_added)[rows_with_mult_var_allele_log], function(r){
      if (grepl('> [A,C,T,G]+,[A,C,T,G]+.*functional=', r)) {
        var =  gsub('^.*Strand Information:', '', r) %>% gsub('\\|.*$', '', .) %>% substr(.,nchar(.),nchar(.))
        gsub('> [A,C,T,G]+,[A,C,T,G]+.*functional=', paste('>', var, 'functional='), r)
      }
    })
  } else {
    row.names(varmat_added)[rows_with_mult_var_allele_log] = sapply(row.names(varmat_added)[rows_with_mult_var_allele_log], function(r){
      if (grepl('> [A,C,T,G]+,[A,C,T,G]+.*functional=', r)) {
        var =  gsub('^.*Strand Information:', '', r) %>% gsub('\\|.*$', '', .) %>% gsub(".*;", "", .)
        gsub('> [A,C,T,G]+,[A,C,T,G]+.*functional=', paste('>', var, 'functional='), r)
      }
    })
  }

  print("III")
  print(row.names(varmat_added)[which(grepl("Coding Indel at 440983", row.names(varmat_added)))])


  return(list(rows_with_multiple_annots_log,
              rows_with_mult_var_allele_log,
              rows_with_overlapping_genes_log,
              split_rows_flag,
              varmat_added))
}

#' Remove any rows with multiple annotations
#' @description - Bypass dealing with rows with multiple annotations (due to
#' overlapping genes or multiallelic sites) by removing them from the data.frame.
#' Useful especially as we are testing this function and the functionality to deal
#' with sites with multiple annotations is not ready.
#'
#' @param varmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return Returns a varmat (class = data.frame) with all rows with multiple annotations
#' removed. Also writes a file called YEAR_MONTH_DATE_rows_with_multiple_annots_removed
#' indicating which rows were removed.
#'
#' @export
remove_rows_with_multiple_annots <- function(varmat){
  # IDENTIFY ROWS WITH MULTIPLE ANNOTATIONS
  num_dividers <- sapply(1:nrow(varmat), function(x) lengths(regmatches(row.names(varmat)[x], gregexpr(";[A,C,G,T]", row.names(varmat)[x]))))
  rows_with_multiple_annotations <- c(1:nrow(varmat))[num_dividers >= 2 & stringr::str_count(row.names(varmat), '\\|') > 9]

  # SAVE TO LOG FILE
  log_file = paste0(Sys.Date(), '_rows_with_multiple_annots_removed')
  write('The following rows with multiple annotations were removed:', log_file)
  write(row.names(varmat)[rows_with_multiple_annotations], log_file, append = TRUE)

  # REMOVE ROWS WITH MULTIPLE ANNOTATIONS
  if (length(rows_with_multiple_annotations) > 0) {
    varmat = varmat[-rows_with_multiple_annotations,]
  }

  return(varmat)
}

#' Update code matrix such that alternative alleles are 0s
#' @description Input the split matrix where rows that once had multiple
#'   annotations on single line are now represented on multiple lines. For the
#'   sites that once were multiallelic sites and are now represented as
#'   biallelic, thus function will change the contents of varmat_code such that
#'   the alternative allele(s) are 0. For example, T -> G, C is split into two
#'   lines: T -> G and T -> C. In the code matrix, turn the codes correspoding
#'   to the allele C in the row T -> G to 0 and the codes corresponding to the
#'   allele G in the row T -> C to 0.
#'
#' @param varmat_code_split - data.frame where the rows are variants (numeric
#'   description variants: numbers ranging from -4 to 3), the columns are
#'   genomes, and the row.names are annotations, and each line has a single
#'   annotation
#' @param varmat_allele_split - data.frame where the rows are variants
#'   (character description variants: A,C,T,G,N,-), the columns are genomes, and
#'   the row.names are annotations, and each line has a single annotation
#' @param ref - character vector length nrow(varmat_code_split) =
#'   nrow(varmat_allele_split) indicating the reference allele in terms of the
#'   positive strand
#' @param var - character vector length nrow(varmat_code_split) =
#'   nrow(varmat_allele_split) indicating the variant allele in terms of the
#'   positive strand
#' @param rows_with_mult_var_allele_log - logical vector length
#'   nrow(varmat_code_split) = nrow(varmat_allele_split) indicating which rows
#'   are multiallelic sites
#'
#' @return - varmat_code - data.frame where the rows are variants (numeric
#'   description variants: numbers ranging from -4 to 3), the columns are
#'   genomes, and the row.names are annotations, and each line has a single
#'   annotation where the alternative/minor allele in a
#'   biallelic-represrentation of a multiallelic site is now 0
#' @export
remove_alt_allele_code_from_split_rows <- function(varmat_code_split, varmat_allele_split, ref, var, rows_with_mult_var_allele_log){

  # UPDATE CODE MATRIX:
  index_mult_var = (1:length(rows_with_mult_var_allele_log))[rows_with_mult_var_allele_log]

  for (i in index_mult_var) {
    # CHANGE THE ALT ALLELE TO 0 IN BIALLELIC REPRESENTATION OF MULTIALLELIC POSITION
    # CHANGE !REF, !VAR, !N, or !-  TO 0
    # T > C,G
    # T > C:   T C G N -; 0 1 1 0 0
    # T > G:   T C G N -; 0 1 1 0 0
    # INTO
    # T > C:   T C G N -; 0 1 0 0 0
    # T > G:   T C G N -; 0 0 1 0 0

    varmat_code_split[i,!(as.character(varmat_allele_split[i,]) %in% c(as.character(var[i]), as.character(ref[i]), 'N', '-'))] = 0
  }
  return(varmat_code_split)
}

#' Root tree on outgroup
#' @description Root tree based on outgroup. If outgroup is null and the tree
#'   isn't rooted, midpoint root tree. If tree is rooted, return tree as is.
#'
#' @param tree phylogenetic tree or file path of tree
#' @param outgroup tip name of outgroup in phylogeny. If NULL, midpoint root if
#'   not rooted
#'
#' @return rooted tree without outgroup
#' @export
root_tree_og = function(tree, outgroup = NULL){
  if (is.character(tree)) {
    # LOAD IN TREE
    tree = ape::read.tree(tree)
  }
  # IF NO OUTGROUP AND TREE IS UNROOTED
  if (is.null(outgroup) & !ape::is.rooted(tree)) {
    # MIDPOINT ROOT TREE
    tree = phytools::midpoint.root(tree)
  } else if (!is.null(outgroup)) {
    # ROOT TREE ON OUTGROUP
    tree = ape::root(tree, outgroup)
    tree = ape::drop.tip(tree, outgroup)
  }
  return(tree)
}

#' Get ancestral state of alleles
#' @description Rereference alleles based on rooted tree
#'
#' @param tree rooted tree
#' @param mat allele matrix (rows are variants, columns are samples)
#'
#' @return matrix of most likely ancestral allele for each row in allele matrix and probability that that is the ancestral state
#' @export
#'
#' @examples
get_anc_alleles = function(tree,mat){
  future::plan(future::multiprocess)

  if (sum(!(tree$tip.label %in% colnames(mat))) > 0) {
    stop('Some samples in tree are not in allele matrix.')
  }

  if (sum(!(colnames(mat) %in% tree$tip.label)) > 0) {
    stop('Some samples in allele matrix are not in tree.')
  }

  if (!ape::is.rooted(tree)) {
    stop('Tree must be rooted.')
  }

  # ORDER MATRIX TO MATCH TREE TIP LABELS
  mat = mat[ ,tree$tip.label]

  if (sum(tree$edge.length == 0) > 0) {
    warning('All zero branch lengths changed to small non-zero number to be able to perform ancestral reconstruction.')
    # Change any edge lengths that are zero to a very small number (so ancestral reconstruction doesn't break)
    tree$edge.length[tree$edge.length == 0] = min(tree$edge.length[tree$edge.length > 0]) / 1000
  }

  # Get ancestral state of root; 1st column = var absent (0), 2nd column = var present (1)
  ar_all = t(future.apply::future_apply(mat, 1, function(tip_states){
    tip_state = unique(tip_states)
    if (length(tip_state) > 1) {
      ar = ape::ace(x = tip_states,phy = tree, type = 'discrete')
      states = ar$lik.anc[1, ]
      tip_state = names(states)[which.max(states)]
      prob = states[which.max(states)]
      c(tip_state, prob)
    } else {
      c(tip_states, 1)
    }
  }))
  return(ar_all)
}


#' Load matrix from path if needed
#'
#' @param mat - loaded in data.frame of varmat or character string of a path to a varmat
#' @description Loads variant matrix from path if not already loaded
#'
#' @return variant matrix
#' @export
load_if_path = function(mat){
  if (is.character(mat)) {
    mat = read.table(mat,
                     header = TRUE,
                     stringsAsFactors = FALSE,
                     sep = "\t",
                     quote = "",
                     row.names = 1)
  }
  return(mat)
}

#' Remove unknown ancestral states
#' @description Remove rows from variant matrix where the ancestral state is unknown (- or N)
#'
#' @param varmat_code
#' @param varmat_allele
#' @param annots
#'
#' @return
#' @export
remove_unknown_anc = function(varmat_code, varmat_allele, annots){
  unknown = annots$anc %in% c('-', 'N')
  removed = rownames(varmat_code)[unknown]
  filename = paste0(Sys.Date(), '_rows_removed_because_unknown_ancestral_state.txt')
  write.table(removed,
              file = filename,
              sep = '\n',
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  return(list(varmat_code = varmat_code[!unknown, ],
              varmat_allele = varmat_allele[!unknown, ],
              annots = annots[!unknown, ]))
}

