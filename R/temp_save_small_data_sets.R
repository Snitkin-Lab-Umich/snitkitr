# These matrices were copied from:
# /nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/Project_propensity_score_match/consensus/2020_12_07_16_25_20_core_results/data_matrix/matrices
# SNP_matrix_code.tsv and SNP_matrix_allele.tsv

code_mat <- read.table("data/SNP_matrix_code.tsv",
                 header = TRUE,
                 stringsAsFactors = FALSE,
                 sep = "\t",
                 quote = "",
                 row.names = 1)

allele_mat <- read.table("data/SNP_matrix_allele.tsv",
                         header = TRUE,
                         stringsAsFactors = FALSE,
                         sep = "\t",
                         quote = "",
                         row.names = 1)

nrow(code_mat)
nrow(allele_mat)

# Make dataset a lot smaller by making it 1/10th the size in each dimension.
rows_to_keep <- seq(1, nrow(code_mat), nrow(code_mat) / 10)
cols_to_keep <- seq(1, ncol(code_mat), ncol(code_mat / 10))

cdif_snp_code_mat <- cdif_code_mat[rows_to_keep, cols_to_keep]
cdif_snp_allele_mat <- allele_mat[rows_to_keep, cols_to_keep]

identical(ncol(cdif_snp_code_mat), ncol(cdif_snp_allele_mat))
identical(colnames(cdif_snp_code_mat), colnames(cdif_snp_allele_mat))
identical(nrow(cdif_snp_code_mat), nrow(cdif_snp_allele_mat))

save(cdif_snp_code_mat, file = "data/cdif_snp_code_matrix.RData")
save(cdif_snp_allele_mat, file = "data/cdif_snp_allele_matrix.RData")
