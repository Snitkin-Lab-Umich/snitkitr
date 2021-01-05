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
  expect_false(identical(parsed_snp_cdif_test1, snitkitr::parsed_snp_cdif1))
  # These should no longer be identical beacuse the saved version references to major

  parsed_snp_cdif_test1b <- parse_snps(varmat_code = snitkitr::cdif_snp_code_mat,
             varmat_allele = snitkitr::cdif_snp_allele_mat,
             tree = NULL,
             og = NULL,
             remove_multi_annots = FALSE,
             return_binary_matrix = TRUE,
             ref_to_anc = FALSE,
             keep_conf_only = TRUE,
             mat_suffix = "_R1_001.fastq.gz",
             ref_to_maj = TRUE)

  expect_identical(parsed_snp_cdif_test1b, snitkitr::parsed_snp_cdif1)
  # But this should pass now because I'm setting ref_to_maj = TRUE

  # the ref_to_maj=TRUE version should have no and yes in reref, while
  # ref_to_maj=FALSE (default) shoudl have just no in reref based on how it's
  # currently coded
  expect_length(table(parsed_snp_cdif_test1$bin$annots$reref), 1)
  expect_length(table(parsed_snp_cdif1$bin$annots$reref), 2)

})


test_that("Parse_snps returns expected output 2", {
  parsed_snp_cdif_test2 <- parse_snps(varmat_code = cdif_snp_code_mat,
                             varmat_allele = cdif_snp_allele_mat,
                             tree = NULL,
                             og = NULL,
                             remove_multi_annots = TRUE,
                             return_binary_matrix = TRUE,
                             ref_to_anc = FALSE,
                             keep_conf_only = TRUE,
                             mat_suffix = "_R1_001.fastq.gz")
  expect_false(identical(parsed_snp_cdif_test2, parsed_snp_cdif2))
  # These should no longer be identical beacuse the saved version references to major while the just run version references to reference genome allele

  parsed_snp_cdif_test2b <- parse_snps(varmat_code = snitkitr::cdif_snp_code_mat,
                                       varmat_allele = snitkitr::cdif_snp_allele_mat,
                                       tree = NULL,
                                       og = NULL,
                                       remove_multi_annots = TRUE,
                                       return_binary_matrix = TRUE,
                                       ref_to_anc = FALSE,
                                       keep_conf_only = TRUE,
                                       mat_suffix = "_R1_001.fastq.gz",
                                       ref_to_maj = TRUE)

  expect_identical(parsed_snp_cdif_test2b, snitkitr::parsed_snp_cdif2)
  # But this should pass now because I'm setting ref_to_maj = TRUE

  # the ref_to_maj=TRUE version should have no and yes in reref, while
  # ref_to_maj=FALSE (default) shoudl have just no in reref based on how it's
  # currently coded
  expect_length(table(parsed_snp_cdif_test2$bin$annots$reref), 1)
  expect_length(table(parsed_snp_cdif2$bin$annots$reref), 2)

})

test_that("Parse_snps returns expected output 3", {
  parsed_snp_cdif_test3 <- parse_snps(varmat_code = cdif_snp_code_mat,
                             varmat_allele = cdif_snp_allele_mat,
                             tree = NULL,
                             og = NULL,
                             remove_multi_annots = FALSE,
                             return_binary_matrix = FALSE,
                             ref_to_anc = FALSE,
                             keep_conf_only = TRUE,
                             mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_snp_cdif_test3, parsed_snp_cdif3)
})

test_that("Parse_snps returns expected output 4", {
  parsed_snp_cdif_test4 <- parse_snps(varmat_code = cdif_snp_code_mat,
                             varmat_allele = cdif_snp_allele_mat,
                             tree = NULL,
                             og = NULL,
                             remove_multi_annots = FALSE,
                             return_binary_matrix = FALSE,
                             ref_to_anc = FALSE,
                             keep_conf_only = FALSE,
                             mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_snp_cdif_test4, parsed_snp_cdif4)
})

test_that("Parse_snps returns expected output 5", {
    parsed_snp_cdif_test5 <- parse_snps(varmat_code = cdif_snp_code_mat,
                             varmat_allele = cdif_snp_allele_mat,
                             tree = cdif_tree,
                             og = NULL,
                             remove_multi_annots = FALSE,
                             return_binary_matrix = TRUE,
                             ref_to_anc = TRUE,
                             keep_conf_only = TRUE,
                             mat_suffix = "_R1_001.fastq.gz")
    expect_identical(parsed_snp_cdif_test5, parsed_snp_cdif5)
})
test_that("Parse_snps returns expected output 6", {
  parsed_snp_cdif_test6 <- parse_snps(varmat_code = cdif_snp_code_mat,
                             varmat_allele = cdif_snp_allele_mat,
                             tree = cdif_tree,
                             og = NULL,
                             remove_multi_annots = TRUE,
                             return_binary_matrix = TRUE,
                             ref_to_anc = TRUE,
                             keep_conf_only = TRUE,
                             mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_snp_cdif_test6, parsed_snp_cdif6)
})

test_that("Parse_snps returns expected output 7", {
  parsed_snp_cdif_test7 <- parse_snps(varmat_code = cdif_snp_code_mat,
                             varmat_allele = cdif_snp_allele_mat,
                             tree = cdif_tree,
                             og = NULL,
                             remove_multi_annots = TRUE,
                             return_binary_matrix = TRUE,
                             ref_to_anc = FALSE,
                             keep_conf_only = TRUE,
                             mat_suffix = "_R1_001.fastq.gz")
  expect_false(identical(parsed_snp_cdif_test7, parsed_snp_cdif7))

  parsed_snp_cdif_test7b <- parse_snps(varmat_code = snitkitr::cdif_snp_code_mat,
                                       varmat_allele = snitkitr::cdif_snp_allele_mat,
                                       tree = cdif_tree,
                                       og = NULL,
                                       remove_multi_annots = TRUE,
                                       return_binary_matrix = TRUE,
                                       ref_to_anc = FALSE,
                                       keep_conf_only = TRUE,
                                       mat_suffix = "_R1_001.fastq.gz",
                                       ref_to_maj = TRUE)

  expect_identical(parsed_snp_cdif_test7b, snitkitr::parsed_snp_cdif7)
  # But this should pass now because I'm setting ref_to_maj = TRUE

  # the ref_to_maj=TRUE version should have no and yes in reref, while
  # ref_to_maj=FALSE (default) shoudl have just no in reref based on how it's
  # currently coded
  expect_length(table(parsed_snp_cdif_test7$bin$annots$reref), 1)
  expect_length(table(parsed_snp_cdif7$bin$annots$reref), 2)


})

test_that("Parse_snps returns expected output 8", {
  parsed_snp_cdif_test8 <- parse_snps(varmat_code = cdif_snp_code_mat,
                                  varmat_allele = cdif_snp_allele_mat,
                                  tree = cdif_tree,
                                  og = NULL,
                                  remove_multi_annots = FALSE,
                                  return_binary_matrix = FALSE,
                                  ref_to_anc = FALSE,
                                  keep_conf_only = FALSE,
                                  mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_snp_cdif_test8, parsed_snp_cdif8)
})

test_that("Parse_indels returns expected output 1", {
    parsed_indel_cdif_test1 <- parse_indels(varmat_code = cdif_indel_code_mat,
                             varmat_allele = cdif_indel_allele_mat,
                             tree = NULL,
                             og = NULL,
                             remove_multi_annots = FALSE,
                             return_binary_matrix = TRUE,
                             ref_to_anc = FALSE,
                             keep_conf_only = TRUE,
                             mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_indel_cdif_test1, parsed_indel_cdif1)
})

test_that("Parse_indels returns expected output 3", {
    parsed_indel_cdif_test3 <- parse_indels(varmat_code = cdif_indel_code_mat,
                                          varmat_allele = cdif_indel_allele_mat,
                                          tree = NULL,
                                          og = NULL,
                                          remove_multi_annots = FALSE,
                                          return_binary_matrix = FALSE,
                                          ref_to_anc = FALSE,
                                          keep_conf_only = TRUE,
                                          mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_indel_cdif_test3, parsed_indel_cdif3)
})

test_that("Parse_indels returns expected output 4", {
  parsed_indel_cdif_test4 <- parse_indels(varmat_code = cdif_indel_code_mat,
                                    varmat_allele = cdif_indel_allele_mat,
                                    tree = NULL,
                                    og = NULL,
                                    remove_multi_annots = FALSE,
                                    return_binary_matrix = FALSE,
                                    ref_to_anc = FALSE,
                                    keep_conf_only = FALSE,
                                    mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_indel_cdif_test4, parsed_indel_cdif4)
})

test_that("Parse_indels returns expected output 5", {
  # With tree
  parsed_indel_cdif_test5 <- parse_indels(varmat_code = cdif_indel_code_mat,
                                    varmat_allele = cdif_indel_allele_mat,
                                    tree = cdif_tree,
                                    og = NULL,
                                    remove_multi_annots = FALSE,
                                    return_binary_matrix = TRUE,
                                    ref_to_anc = TRUE,
                                    keep_conf_only = TRUE,
                                    mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_indel_cdif_test5, parsed_indel_cdif5)
})

test_that("Parse_indels returns expected output 8", {
  parsed_indel_cdif_test8 <- parse_indels(varmat_code = cdif_indel_code_mat,
                                    varmat_allele = cdif_indel_allele_mat,
                                    tree = cdif_tree,
                                    og = NULL,
                                    remove_multi_annots = FALSE,
                                    return_binary_matrix = FALSE,
                                    ref_to_anc = FALSE,
                                    keep_conf_only = FALSE,
                                    mat_suffix = "_R1_001.fastq.gz")
  expect_identical(parsed_indel_cdif_test8, parsed_indel_cdif8)
})

