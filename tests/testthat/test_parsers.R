test_that("Parse_snps returns expected output 1", {
  parsed_snp_cdif_test1 <- parse_snps(varmat_code = snitkitr::cdif_snp_code_mat,
                                  varmat_allele = snitkitr::cdif_snp_allele_mat,
                                  tree = NULL,
                                  og = NULL,
                                  remove_multi_annots = FALSE,
                                  return_binary_matrix = TRUE,
                                  ref_to_anc = FALSE,
                                  keep_conf_only = TRUE,
                                  mat_suffix = "_R1_001.fastq.gz")
  expect_true(identical(parsed_snp_cdif_test1, snitkitr::parsed_snp_cdif1))
  expect_identical(parsed_snp_cdif_test1, snitkitr::parsed_snp_cdif1)
  expect_equal(parsed_snp_cdif_test1, snitkitr::parsed_snp_cdif1)

  print(identical(parsed_snp_cdif_test1, snitkitr::parsed_snp_cdif1))


})


# test_that("Parse_snps returns expected output 2", {
#   parsed_snp_cdif_test2 <- parse_snps(varmat_code = cdif_snp_code_mat,
#                              varmat_allele = cdif_snp_allele_mat,
#                              tree = NULL,
#                              og = NULL,
#                              remove_multi_annots = TRUE,
#                              return_binary_matrix = TRUE,
#                              ref_to_anc = FALSE,
#                              keep_conf_only = TRUE,
#                              mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test2, parsed_snp_cdif2)
# })
#
# test_that("Parse_snps returns expected output 3", {
#   parsed_snp_cdif_test3 <- parse_snps(varmat_code = cdif_snp_code_mat,
#                              varmat_allele = cdif_snp_allele_mat,
#                              tree = NULL,
#                              og = NULL,
#                              remove_multi_annots = FALSE,
#                              return_binary_matrix = FALSE,
#                              ref_to_anc = FALSE,
#                              keep_conf_only = TRUE,
#                              mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test3, parsed_snp_cdif3)
# })
#
# test_that("Parse_snps returns expected output 4", {
#   parsed_snp_cdif_test4 <- parse_snps(varmat_code = cdif_snp_code_mat,
#                              varmat_allele = cdif_snp_allele_mat,
#                              tree = NULL,
#                              og = NULL,
#                              remove_multi_annots = FALSE,
#                              return_binary_matrix = FALSE,
#                              ref_to_anc = FALSE,
#                              keep_conf_only = FALSE,
#                              mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test4, parsed_snp_cdif4)
# })
#
# test_that("Parse_snps returns expected output 5", {
#     parsed_snp_cdif_test5 <- parse_snps(varmat_code = cdif_snp_code_mat,
#                              varmat_allele = cdif_snp_allele_mat,
#                              tree = cdif_tree,
#                              og = NULL,
#                              remove_multi_annots = FALSE,
#                              return_binary_matrix = TRUE,
#                              ref_to_anc = TRUE,
#                              keep_conf_only = TRUE,
#                              mat_suffix = "_R1_001.fastq.gz")
#     expect_identical(parsed_snp_cdif_test5, parsed_snp_cdif5)
# })
# test_that("Parse_snps returns expected output 6", {
#   parsed_snp_cdif_test6 <- parse_snps(varmat_code = cdif_snp_code_mat,
#                              varmat_allele = cdif_snp_allele_mat,
#                              tree = cdif_tree,
#                              og = NULL,
#                              remove_multi_annots = TRUE,
#                              return_binary_matrix = TRUE,
#                              ref_to_anc = TRUE,
#                              keep_conf_only = TRUE,
#                              mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test6, parsed_snp_cdif6)
# })
#
# test_that("Parse_snps returns expected output 7", {
#   parsed_snp_cdif_test7 <- parse_snps(varmat_code = cdif_snp_code_mat,
#                              varmat_allele = cdif_snp_allele_mat,
#                              tree = cdif_tree,
#                              og = NULL,
#                              remove_multi_annots = TRUE,
#                              return_binary_matrix = TRUE,
#                              ref_to_anc = FALSE,
#                              keep_conf_only = TRUE,
#                              mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test7, parsed_snp_cdif7)
# })
#
# test_that("Parse_snps returns expected output 8", {
#   parsed_snp_cdif_test8 <- parse_snps(varmat_code = cdif_snp_code_mat,
#                                   varmat_allele = cdif_snp_allele_mat,
#                                   tree = cdif_tree,
#                                   og = NULL,
#                                   remove_multi_annots = FALSE,
#                                   return_binary_matrix = FALSE,
#                                   ref_to_anc = FALSE,
#                                   keep_conf_only = FALSE,
#                                   mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test8, parsed_snp_cdif8)
# })
#
# test_that("Parse_indels returns expected output 1", {
#     parsed_snp_cdif_test1 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                              varmat_allele = cdif_indel_allele_mat,
#                              tree = NULL,
#                              og = NULL,
#                              remove_multi_annots = FALSE,
#                              return_binary_matrix = FALSE,
#                              ref_to_anc = FALSE,
#                              keep_conf_only = TRUE,
#                              mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test1, parsed_snp_cdif1)
# })

# test_that("Parse_indels returns expected output 3", {
#     parsed_snp_cdif_test3 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                              varmat_allele = cdif_indel_allele_mat,
#                              tree = NULL,
#                              og = NULL,
#                              remove_multi_annots = FALSE,
#                              return_binary_matrix = FALSE,
#                              ref_to_anc = FALSE,
#                              keep_conf_only = TRUE,
#                              mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test3, parsed_snp_cdif3)
# })
#
# test_that("Parse_indels returns expected output 4", {
#   parsed_snp_cdif_test4 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                                     varmat_allele = cdif_indel_allele_mat,
#                                     tree = NULL,
#                                     og = NULL,
#                                     remove_multi_annots = FALSE,
#                                     return_binary_matrix = FALSE,
#                                     ref_to_anc = FALSE,
#                                     keep_conf_only = FALSE,
#                                     mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test4, parsed_snp_cdif4)
# })
#
# test_that("Parse_indels returns expected output 5", {
#   # With tree
#   parsed_snp_cdif_test5 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                                     varmat_allele = cdif_indel_allele_mat,
#                                     tree = cdif_tree,
#                                     og = NULL,
#                                     remove_multi_annots = FALSE,
#                                     return_binary_matrix = TRUE,
#                                     ref_to_anc = TRUE,
#                                     keep_conf_only = TRUE,
#                                     mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test5, parsed_snp_cdif5)
# })
#
# test_that("Parse_indels returns expected output 8", {
#   parsed_snp_cdif_test8 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                                     varmat_allele = cdif_indel_allele_mat,
#                                     tree = cdif_tree,
#                                     og = NULL,
#                                     remove_multi_annots = FALSE,
#                                     return_binary_matrix = FALSE,
#                                     ref_to_anc = FALSE,
#                                     keep_conf_only = FALSE,
#                                     mat_suffix = "_R1_001.fastq.gz")
#   expect_identical(parsed_snp_cdif_test8, parsed_snp_cdif8)
# })
#
