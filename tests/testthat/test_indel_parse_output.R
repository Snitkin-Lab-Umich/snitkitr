test_that("parsed_indel matrices are all the same size & have the same row names & number of annotation entries", {
  parsed_indel <- parse_indels(varmat_code = snitkitr::cdif_indel_code_mat[1:100, 1:3],
                               varmat_allele = snitkitr::cdif_indel_allele_mat[1:100, 1:3],
                               tree = NULL,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = FALSE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = TRUE,
                               parallelization = "multisession")
  # dimensions
  expect_identical(dim(parsed_indel$bin$mat), dim(parsed_indel$code$mat))
  expect_identical(dim(parsed_indel$bin$mat), dim(parsed_indel$allele$mat))

  # rownames
  expect_identical(row.names(parsed_indel$bin$mat), row.names(parsed_indel$code$mat))
  expect_identical(row.names(parsed_indel$bin$mat), row.names(parsed_indel$allele$mat))

  # nrow
  expect_identical(nrow(parsed_indel$bin$annots), nrow(parsed_indel$code$annots))
  expect_identical(nrow(parsed_indel$bin$annots), nrow(parsed_indel$allele$annots))

  # nrow annots == nrow mat
  expect_identical(nrow(parsed_indel$bin$annots), nrow(parsed_indel$code$mat))
  expect_identical(nrow(parsed_indel$bin$annots), nrow(parsed_indel$allele$mat))
})

test_that("parsed_indel annot are expected parsings of the row names", {
  parsed_indel <- parse_indels(varmat_code = snitkitr::cdif_indel_code_mat[1:100, 1:3],
                               varmat_allele = snitkitr::cdif_indel_allele_mat[1:100, 1:3],
                               tree = NULL,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = FALSE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = TRUE,
                               parallelization = "multisession")
  expect_identical(parsed_indel$allele$annots, parsed_indel$code$annots)
  expect_false(identical(parsed_indel$bin$annots, parsed_indel$code$annots))

  # row.names(parsed_indel$allele$mat)[1]
  # "Coding Indel at 3009 > T functional=NULL_NULL_NULL locus_tag=CD630_00030
  #  Strand Information: CD630_00030=+;T|frameshift_variant|HIGH|CD630_00030|
  #  c.196delA|p.Val67fs|196/207|66/68|null or hypothetical protein|putative
  #  RNA-binding mediating protein;"
  # t(parsed_indel$allele$annots[1, ] )
  # label                           "Coding Indel"
  expect_true(parsed_indel$allele$annots$label[1] == "Coding Indel")
  # pos                             "3009"
  expect_true(parsed_indel$allele$annots$pos[1] == 3009)
  # phage                           "NULL"
  # repeated_region                 "NULL"
  # masked                          "NULL"
  expect_true(parsed_indel$allele$annots$phage[1] == "NULL")
  expect_true(parsed_indel$allele$annots$repeated_region[1] == "NULL")
  expect_true(parsed_indel$allele$annots$masked[1] == "NULL")
  # locus_tag                       "CD630_00030"
  expect_true(parsed_indel$allele$annots$locus_tag[1] == "CD630_00030")
  # strand_info                     "CD630_00030=+"
  expect_true(parsed_indel$allele$annots$strand_info[1] == "CD630_00030=+")
  # strand                          "+"
  expect_true(parsed_indel$allele$annots$strand[1] == "+")
  # ref                             "TA"
  # var                             "T"
  # variant_type                    "frameshift_variant"
  expect_true(parsed_indel$allele$annots$variant_type[1] == "frameshift_variant")
  # snpeff_impact                   "HIGH"
  expect_true(parsed_indel$allele$annots$snpeff_impact[1] == "HIGH")
  # nuc_pos_in_gene                 "196"
  expect_true(parsed_indel$allele$annots$nuc_pos_in_gene[1] == "196")
  # aa_pos_in_gene                  "66"
  expect_true(parsed_indel$allele$annots$aa_pos_in_gene[1] == "66")
  # gene_length_in_bp               "207"
  expect_true(parsed_indel$allele$annots$gene_length_in_bp[1] == "207")
  # annotation_1                    "null or hypothetical protein"
  expect_true(parsed_indel$allele$annots$annotation_1[1] == "null or hypothetical protein")
  # annotation_2                    "putative RNA-binding mediating protein;"
  expect_true(parsed_indel$allele$annots$annotation_2[1] == "putative RNA-binding mediating protein;")
  # ig_gene1                        ""
  # ig_gene2                        ""
  # intergenic                      "FALSE"
  # raw_rownames                    "Coding Indel at 3009 > T functional=NULL_NULL_NULL locus_tag=CD630_00030 Strand Information: CD630_00030=+;T|frameshift_variant|HIGH|CD630_00030|c.196delA|p.Val67fs|196/207|66/68|null or hypothetical protein|putative RNA-binding mediating protein;"
  # indel_type                      "del"
  expect_true(parsed_indel$allele$annots$indel_type[1] == "del")
  # indel_nuc                       "A"
  expect_true(parsed_indel$allele$annots$indel_nuc[1] == "A")
  # rows_with_multiple_annots_log   "FALSE"
  # ^ Yes, should be false because original matrix has a rowname identical to this rowname
  expect_true(parsed_indel$allele$annots$raw_rownames[1] %in%
                row.names(snitkitr::cdif_indel_code_mat[1:100,]))
  expect_true(parsed_indel$allele$annots$rows_with_multiple_annots_log[1] == "FALSE")
  # rows_with_mult_var_allele_log   "FALSE"
  # rows_with_overlapping_genes_log "FALSE"
  # split_rows_flag                 "1"
  # maj                             "TA"
  most_common_occurance <- unname(unlist(parsed_indel$allele$mat[1, ])) %>% table(.) %>% unname() %>% max()
  most_common_occurance_tbl_log <-  unname(unlist(parsed_indel$allele$mat[1, ])) %>% table() == most_common_occurance
  maj <- names(most_common_occurance_tbl_log)[most_common_occurance_tbl_log == TRUE]
  expect_identical(maj, parsed_indel$allele$annots$maj[1] %>% as.character())

  # Double check some of the trickier stuff for an intergenic indel now
  # row.names(parsed_indel$allele$mat)[2]
  # "Non-Coding Indel at 278620 > A functional=NULL_NULL_NULL
  # locus_tag=CD630_02140-CD630_02150 Strand Information:
  # CD630_02140=-/CD630_02150=+/;A|intergenic_region|MODIFIER|
  # CD630_02140-CD630_02150|n.278621_278624delATAT||||null or hypothetical
  # protein,null or hypothetical protein|uncharacterised protein,Transposase
  # IS200/IS605-like OrfA,;"
  expect_true(parsed_indel$allele$annots$locus_tag[2] == "CD630_02140-CD630_02150")
  expect_true(parsed_indel$allele$annots$strand_info[2] == "CD630_02140=-/CD630_02150=+/")
  expect_true(parsed_indel$allele$annots$strand[2] == "-/+")
  expect_true(parsed_indel$allele$annots$intergenic[2] == "TRUE")
  expect_true(parsed_indel$allele$annots$ig_gene1[2] == "CD630_02140")
  expect_true(parsed_indel$allele$annots$ig_gene2[2] == "CD630_02150")
  expect_true(parsed_indel$allele$annots$variant_type[2] == "intergenic_region")
})

test_that("parsed_indel bin mat expected given ref_to_maj", {
  col_index <- c(1:3)
  parsed_indel <- parse_indels(varmat_code = snitkitr::cdif_indel_code_mat[1:200, col_index],
                               varmat_allele = snitkitr::cdif_indel_allele_mat[1:200, col_index],
                               tree = NULL,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = FALSE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = TRUE,
                               parallelization = "multisession")

  check_maj_rerefencing(parsed_indel, 1, col_index)
  check_maj_rerefencing(parsed_indel, 2, col_index)
  check_maj_rerefencing(parsed_indel, 3, col_index)
  check_maj_rerefencing(parsed_indel, 4, col_index)
  check_maj_rerefencing(parsed_indel, 5, col_index)
  check_maj_rerefencing(parsed_indel, 6, col_index)
  check_maj_rerefencing(parsed_indel, 7, col_index)
  check_maj_rerefencing(parsed_indel, 8, col_index)
  check_maj_rerefencing(parsed_indel, 9, col_index)
  check_maj_rerefencing(parsed_indel, 10, col_index)
})

test_that("parsed_indel bin mat expected given ref_to_anc", {
  parsed_indel <- parse_indels(varmat_code = snitkitr::cdif_indel_code_mat[1:100, ],
                               varmat_allele = snitkitr::cdif_indel_allele_mat[1:100, ],
                               tree = snitkitr::cdif_tree,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = TRUE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = FALSE,
                               parallelization = "multisession")

  check_anc_rerefencing(parsed_indel, 1)
  check_anc_rerefencing(parsed_indel, 2)
  check_anc_rerefencing(parsed_indel, 3)
  check_anc_rerefencing(parsed_indel, 4)
  check_anc_rerefencing(parsed_indel, 5)
  check_anc_rerefencing(parsed_indel, 6)
  check_anc_rerefencing(parsed_indel, 7)
  check_anc_rerefencing(parsed_indel, 8)
  check_anc_rerefencing(parsed_indel, 9)
  check_anc_rerefencing(parsed_indel, 10)
})


test_that("parsed_indel bin mat expected given reference to reference genome", {
  parsed_indel <- parse_indels(varmat_code = snitkitr::cdif_indel_code_mat[1:100, ],
                               varmat_allele = snitkitr::cdif_indel_allele_mat[1:100, ],
                               tree = NULL,
                               og = NULL,
                               remove_multi_annots = TRUE,
                               return_binary_matrix = TRUE,
                               ref_to_anc = FALSE,
                               keep_conf_only = TRUE,
                               mat_suffix = "_R1_001.fastq.gz",
                               ref_to_maj = FALSE,
                               parallelization = "multisession")

  check_match_ref_genome(parsed_indel, 1)
  check_match_ref_genome(parsed_indel, 2)
  check_match_ref_genome(parsed_indel, 3)
  check_match_ref_genome(parsed_indel, 4)
  check_match_ref_genome(parsed_indel, 5)
  check_match_ref_genome(parsed_indel, 6)
  check_match_ref_genome(parsed_indel, 7)
  check_match_ref_genome(parsed_indel, 8)
  check_match_ref_genome(parsed_indel, 9)
  check_match_ref_genome(parsed_indel, 10)
})


