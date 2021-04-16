check_anc_rerefencing <- function(parsed_indel, row) {
  # with these inputs:
  # remove_multi_annots = TRUE,
  # return_binary_matrix = TRUE,
  # ref_to_anc = TRUE,
  # keep_conf_only = TRUE,
  # mat_suffix = "_R1_001.fastq.gz",
  # ref_to_maj = FALSE,

  # grab the inferred anc allele and then check that in the bin mat the positions that correspond to that allele are 0s
  anc_allele <- unname(unlist(parsed_indel$bin$annots$anc))[row]
  ref_allele <- as.character((parsed_indel$allele$annots$ref[row]))
  anc_index <- which(unname(unlist(parsed_indel$allele$mat[row, ])) == anc_allele)
  not_anc_index <- c(1:nrow(parsed_indel$bin$mat))[!c(1:nrow(parsed_indel$bin$mat)) %in% anc_index]

  original_encoding_for_ancestral_allele <- unlist(unname(
    parsed_indel$code$mat[row,][which(parsed_indel$allele$mat[row,] == anc_allele)[1]]))
  # if ancestral allele equals reference allele then:
  if (original_encoding_for_ancestral_allele != 0) {# ancestral allele initial endocded as a not 0
    expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[anc_index], rep(0, length(anc_index)))
    expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[not_anc_index], rep(1, length(not_anc_index)))
  } else {# ancestral allele initially encoded as a 0
    expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[anc_index], rep(1, length(anc_index)))
    expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[not_anc_index], rep(0, length(not_anc_index)))
  }
}

check_maj_rerefencing <- function(parsed_indel, row, col_index) {
  # grab the major allele and then check that in the bin mat the positions that correspond to that allele are 0s
  most_common_occurance <- unname(unlist(parsed_indel$allele$mat[row, ])) %>% table(.) %>% unname() %>% max()
  most_common_occurance_tbl_log <-  unname(unlist(parsed_indel$allele$mat[row, ])) %>% table() == most_common_occurance
  maj <- names(most_common_occurance_tbl_log)[most_common_occurance_tbl_log == TRUE]
  maj_index <- which(unname(unlist(parsed_indel$allele$mat[row, ])) == maj)
  min_index <- c(col_index)[!c(col_index) %in% maj_index]
  expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[maj_index], rep(0, length(maj_index)))
  expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[min_index], rep(1, length(min_index)))
}

check_match_ref_genome <- function(parsed_indel, row) {
  # grab the ref allele and then check that in the bin mat the positions that correspond to that allele are 0s
  ref_allele <- unname(unlist(parsed_indel$bin$annots$ref))[row]
  ref_index <- which(unname(unlist(parsed_indel$allele$mat[row, ])) == ref_allele)
  not_ref_index <- c(1:nrow(parsed_indel$bin$mat))[!c(1:nrow(parsed_indel$bin$mat)) %in% ref_index]
  expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[ref_index], rep(0, length(ref_index)))
  expect_identical(unname(unlist(parsed_indel$bin$mat[row, ]))[not_ref_index], rep(1, length(not_ref_index)))
}
