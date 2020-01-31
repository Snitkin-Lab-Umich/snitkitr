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
#'
#' @examples

remove_rows_with_bugs <- function(varmat){
  library(magrittr)
  library(Biostrings)
  library(stringr)

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
  num_pipes = str_count(row.names(varmat), '\\|')
  table(num_pipes) # remove not intervals of 9

  num_semicolon = str_count(row.names(varmat), ';')
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
}#end remove_rows_with_bugs
