# # # # These matrices were copied from:
# # # # /nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/Project_propensity_score_match/consensus/2020_12_07_16_25_20_core_results/data_matrix/matrices
# # # # SNP_matrix_code.tsv and SNP_matrix_allele.tsv
# # # # Indel_matrix_code.tsv and Indel_matrix_allele.tsv
# # # #
# # code_mat <- read.table("data/SNP_matrix_code.tsv",
# #                  header = TRUE,
# #                  stringsAsFactors = FALSE,
# #                  sep = "\t",
# #                  quote = "",
# #                  row.names = 1)
# # # #
# # allele_mat <- read.table("data/SNP_matrix_allele.tsv",
# #                          header = TRUE,
# #                          stringsAsFactors = FALSE,
# #                          sep = "\t",
# #                          quote = "",
# #                          row.names = 1)
# # nrow(code_mat) # 247564
# # nrow(allele_mat) # 247564
# #
# # # Make dataset a lot smaller by making it 10000 rows and 25 columnes
# # rows_to_keep <- floor(seq(1, nrow(code_mat), nrow(code_mat) / 1000))
# # cols_to_keep <- floor(seq(1, ncol(code_mat), ncol(code_mat) / 25))
# #
# # cdif_snp_code_mat <- code_mat[rows_to_keep, cols_to_keep]
# # cdif_snp_allele_mat <- allele_mat[rows_to_keep, cols_to_keep]
# # #
# # dim(cdif_snp_code_mat) # 1000 25
# # dim(cdif_snp_allele_mat) # 1000 25
# # #
# # identical(ncol(cdif_snp_code_mat), ncol(cdif_snp_allele_mat)) # TRUE
# # identical(colnames(cdif_snp_code_mat), colnames(cdif_snp_allele_mat)) # TRUE
# # identical(nrow(cdif_snp_code_mat), nrow(cdif_snp_allele_mat)) # TRUE
# # # #
# # save(cdif_snp_code_mat, file = "data/cdif_snp_code_matrix.RData")
# # save(cdif_snp_allele_mat, file = "data/cdif_snp_allele_matrix.RData")
# # #
# # # # Save a tree
# # tree <- ape::read.tree("data/severity_tree_rooted.tree")
# # ape::is.rooted(tree)
# # #
# # tips_to_keep <- gsub("_R1_001.fastq.gz", "", colnames(cdif_snp_code_mat))
# # tips_to_drop <- tree$tip.label[!tree$tip.label %in% tips_to_keep]
# # #
# # cdif_tree <- ape::drop.tip(tree, tips_to_drop)
# # save(cdif_tree, file = "data/cdif_tree.RData")
#
#
# # All of this parsing was done prior to making any changes to the parser on 1/4/2021.
# # No tree
# # parsed_snp_cdif1 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = NULL,
# #                            og = NULL,
# #                            remove_multi_annots = FALSE,
# #                            return_binary_matrix = TRUE,
# #                            ref_to_anc = FALSE,
# #                            keep_conf_only = TRUE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif1, file = "data/cdif_snp_parsed1.RData")
# #
# # parsed_snp_cdif2 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = NULL,
# #                            og = NULL,
# #                            remove_multi_annots = TRUE,
# #                            return_binary_matrix = TRUE,
# #                            ref_to_anc = FALSE,
# #                            keep_conf_only = TRUE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif2, file = "data/cdif_snp_parsed2.RData")
# # #
# # parsed_snp_cdif3 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = NULL,
# #                            og = NULL,
# #                            remove_multi_annots = FALSE,
# #                            return_binary_matrix = FALSE,
# #                            ref_to_anc = FALSE,
# #                            keep_conf_only = TRUE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif3, file = "data/cdif_snp_parsed3.RData")
# # #
# # parsed_snp_cdif4 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = NULL,
# #                            og = NULL,
# #                            remove_multi_annots = FALSE,
# #                            return_binary_matrix = FALSE,
# #                            ref_to_anc = FALSE,
# #                            keep_conf_only = FALSE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif4, file = "data/cdif_snp_parsed4.RData")
# # #
# # # # With tree
# # parsed_snp_cdif5 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = cdif_tree,
# #                            og = NULL,
# #                            remove_multi_annots = FALSE,
# #                            return_binary_matrix = TRUE,
# #                            ref_to_anc = TRUE,
# #                            keep_conf_only = TRUE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif5, file = "data/cdif_snp_parsed5.RData")
# # #
# # #
# # parsed_snp_cdif6 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = cdif_tree,
# #                            og = NULL,
# #                            remove_multi_annots = TRUE,
# #                            return_binary_matrix = TRUE,
# #                            ref_to_anc = TRUE,
# #                            keep_conf_only = TRUE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif6, file = "data/cdif_snp_parsed6.RData")
# # #
# # parsed_snp_cdif7 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = cdif_tree,
# #                            og = NULL,
# #                            remove_multi_annots = TRUE,
# #                            return_binary_matrix = TRUE,
# #                            ref_to_anc = FALSE,
# #                            keep_conf_only = TRUE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif7, file = "data/cdif_snp_parsed7.RData")
# # #
# # parsed_snp_cdif8 <- parse_snps(varmat_code = cdif_snp_code_mat,
# #                            varmat_allele = cdif_snp_allele_mat,
# #                            tree = cdif_tree,
# #                            og = NULL,
# #                            remove_multi_annots = FALSE,
# #                            return_binary_matrix = FALSE,
# #                            ref_to_anc = FALSE,
# #                            keep_conf_only = FALSE,
# #                            mat_suffix = "_R1_001.fastq.gz")
# # save(parsed_snp_cdif8, file = "data/cdif_snp_parsed8.RData")
# # ----
# #
# # # Indels -----
# # indel_code_mat <- read.table("data/Indel_matrix_code.tsv",
# #                   header = TRUE,
# #                   stringsAsFactors = FALSE,
# #                   sep = "\t",
# #                   quote = "",
# #                   row.names = 1)
# # #
# # indel_allele_mat <- read.table("data/Indel_matrix_allele.tsv",
# #                           header = TRUE,
# #                           stringsAsFactors = FALSE,
# #                           sep = "\t",
# #                           quote = "",
# #                           row.names = 1)
# # #
# # nrow(indel_code_mat) # 12401
# # nrow(indel_allele_mat) # 12401
# # #
# # # # Make dataset a lot smaller by making it 1000 rows and 25 columnes
# # rows_to_keep <- floor(seq(1, nrow(indel_code_mat), nrow(indel_code_mat) / 1000))
# # cols_to_keep <- floor(seq(1, ncol(indel_code_mat), ncol(indel_code_mat) / 25))
# # #
# # cdif_indel_code_mat <- indel_code_mat[rows_to_keep, cols_to_keep]
# # cdif_indel_allele_mat <- indel_allele_mat[rows_to_keep, cols_to_keep]
# # #
# # dim(cdif_indel_code_mat) # 1000 25
# # dim(cdif_indel_allele_mat) # 1000 25
# # #
# # identical(ncol(cdif_indel_code_mat), ncol(cdif_indel_allele_mat)) # TRUE
# # identical(colnames(cdif_indel_code_mat), colnames(cdif_indel_allele_mat)) # TRUE
# # identical(nrow(cdif_indel_code_mat), nrow(cdif_indel_allele_mat)) # TRUE
# #
# # save(cdif_indel_code_mat, file = "data/cdif_indel_code_matrix.RData")
# # save(cdif_indel_allele_mat, file = "data/cdif_indel_allele_matrix.RData")
# #
#
# load("data/cdif_tree.RData")
# load("data/cdif_indel_allele_matrix.RData")
# load("data/cdif_indel_code_matrix.RData")
# # # No tree
# parsed_indel_cdif1 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = NULL,
#                            og = NULL,
#                            remove_multi_annots = FALSE,
#                            return_binary_matrix = TRUE,
#                            ref_to_anc = FALSE,
#                            keep_conf_only = TRUE,
#                            mat_suffix = "_R1_001.fastq.gz")
# save(parsed_indel_cdif1, file = "data/cdif_indel_parsed1.RData")
#
# parsed_indel_cdif2 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = NULL,
#                            og = NULL,
#                            remove_multi_annots = TRUE,
#                            return_binary_matrix = TRUE,
#                            ref_to_anc = FALSE,
#                            keep_conf_only = TRUE,
#                            mat_suffix = "_R1_001.fastq.gz")
# save(parsed_indel_cdif2, file = "data/cdif_indel_parsed2.RData")
# # Can't save because errors out
#
# parsed_indel_cdif3 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = NULL,
#                            og = NULL,
#                            remove_multi_annots = FALSE,
#                            return_binary_matrix = FALSE,
#                            ref_to_anc = FALSE,
#                            keep_conf_only = TRUE,
#                            mat_suffix = "_R1_001.fastq.gz")
# save(parsed_indel_cdif3, file = "data/cdif_indel_parsed3.RData")
#
# parsed_indel_cdif4 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = NULL,
#                            og = NULL,
#                            remove_multi_annots = FALSE,
#                            return_binary_matrix = FALSE,
#                            ref_to_anc = FALSE,
#                            keep_conf_only = FALSE,
#                            mat_suffix = "_R1_001.fastq.gz")
# save(parsed_indel_cdif4, file = "data/cdif_indel_parsed4.RData")
#
# # With tree
# parsed_indel_cdif5 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = cdif_tree,
#                            og = NULL,
#                            remove_multi_annots = FALSE,
#                            return_binary_matrix = TRUE,
#                            ref_to_anc = TRUE,
#                            keep_conf_only = TRUE,
#                            mat_suffix = "_R1_001.fastq.gz")
# save(parsed_indel_cdif5, file = "data/cdif_indel_parsed5.RData")
#
# parsed_indel_cdif6 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = cdif_tree,
#                            og = NULL,
#                            remove_multi_annots = TRUE,
#                            return_binary_matrix = TRUE,
#                            ref_to_anc = TRUE,
#                            keep_conf_only = TRUE,
#                            mat_suffix = "_R1_001.fastq.gz")
# # Couldnt' save because parser didn't work
# save(parsed_indel_cdif6, file = "data/cdif_indel_parsed6.RData")
#
# parsed_indel_cdif7 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = cdif_tree,
#                            og = NULL,
#                            remove_multi_annots = TRUE,
#                            return_binary_matrix = TRUE,
#                            ref_to_anc = FALSE,
#                            keep_conf_only = TRUE,
#                            mat_suffix = "_R1_001.fastq.gz")
# # Couldnt' save because parser didn't work
# save(parsed_indel_cdif7, file = "data/cdif_indel_parsed7.RData")
# #
# parsed_indel_cdif8 <- parse_indels(varmat_code = cdif_indel_code_mat,
#                            varmat_allele = cdif_indel_allele_mat,
#                            tree = cdif_tree,
#                            og = NULL,
#                            remove_multi_annots = FALSE,
#                            return_binary_matrix = FALSE,
#                            ref_to_anc = FALSE,
#                            keep_conf_only = FALSE,
#                            mat_suffix = "_R1_001.fastq.gz")
# save(parsed_indel_cdif8, file = "data/cdif_indel_parsed8.RData")
#
#
