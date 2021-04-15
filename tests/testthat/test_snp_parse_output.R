test_that("parsed_snp bin mat expected given ref_to_maj", {
  col_index <- c(1:3)
  parsed_snp <- parse_snps(varmat_code = snitkitr::cdif_snp_code_mat[1:200, col_index],
                               varmat_allele = snitkitr::cdif_snp_allele_mat[1:200, col_index],
                               tree = NULL,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = FALSE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = TRUE)

  check_maj_rerefencing(parsed_snp, 1, col_index)
  check_maj_rerefencing(parsed_snp, 2, col_index)
  check_maj_rerefencing(parsed_snp, 3, col_index)
  check_maj_rerefencing(parsed_snp, 4, col_index)
  check_maj_rerefencing(parsed_snp, 5, col_index)
  check_maj_rerefencing(parsed_snp, 6, col_index)
  check_maj_rerefencing(parsed_snp, 7, col_index)
  check_maj_rerefencing(parsed_snp, 8, col_index)
  check_maj_rerefencing(parsed_snp, 9, col_index)
  check_maj_rerefencing(parsed_snp, 10, col_index)
})

test_that("parsed_snp bin mat expected given ref_to_anc", {
  parsed_snp <- parse_snps(varmat_code = snitkitr::cdif_snp_code_mat[1:100, ],
                               varmat_allele = snitkitr::cdif_snp_allele_mat[1:100, ],
                               tree = snitkitr::cdif_tree,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = TRUE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = FALSE)

  check_anc_rerefencing(parsed_snp, 1)
  check_anc_rerefencing(parsed_snp, 2)
  check_anc_rerefencing(parsed_snp, 3)
  check_anc_rerefencing(parsed_snp, 4)
  check_anc_rerefencing(parsed_snp, 5)
  check_anc_rerefencing(parsed_snp, 6)
  check_anc_rerefencing(parsed_snp, 7)
  check_anc_rerefencing(parsed_snp, 8)
  check_anc_rerefencing(parsed_snp, 9)
  check_anc_rerefencing(parsed_snp, 10)
})


test_that("parsed_snp bin mat expected given reference to reference genome", {
  parsed_snp <- parse_snps(varmat_code = snitkitr::cdif_snp_code_mat[1:100, ],
                               varmat_allele = snitkitr::cdif_snp_allele_mat[1:100, ],
                               tree = NULL,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = FALSE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = FALSE)

  check_match_ref_genome(parsed_snp, 1)
  check_match_ref_genome(parsed_snp, 2)
  check_match_ref_genome(parsed_snp, 3)
  check_match_ref_genome(parsed_snp, 4)
  check_match_ref_genome(parsed_snp, 5)
  check_match_ref_genome(parsed_snp, 6)
  check_match_ref_genome(parsed_snp, 7)
  check_match_ref_genome(parsed_snp, 8)
  check_match_ref_genome(parsed_snp, 9)
  check_match_ref_genome(parsed_snp, 10)
})
