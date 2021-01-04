#' SNP Code matrix for the Cdif trehalose/severity data set
#'
#' A subset SNP code matrix to use as an example
#' @format A matrix with 10,000 rows and 25 columns. Each row name corresponds to a
#'   variant. Each column corresponds to a sample.
#' \describe{
#'   \item{snp_code_mat}{The number in the code matrix indicates confidence and
#'   quality metrics. This particular data matrix is a subset of a matrix from:
#'   /nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/Project_propensity_score_match/consensus/2020_12_07_16_25_20_core_results/data_matrix/matrices}
#'  }
"cdif_snp_code_mat"

#' SNP Allele matrix for the Cdif trehalose/severity data set
#'
#' A subset SNP allele matrix to use as an example
#' @format A matrix with 10,000 rows and 25 columns. Each row name corresponds to a
#'   variant. Each column corresponds to a sample.
#' \describe{
#'   \item{snp_allele_mat}{The letter in the allele matrix indicates the
#'   variant in that sample. This particular data matrix is a subset of a matrix
#'   from: /nfs/turbo/umms-esnitkin/Project_Cdiff/Sequence_data/Project_propensity_score_match/consensus/2020_12_07_16_25_20_core_results/data_matrix/matrices}
#'  }
"cdif_snp_allele_mat"


#' Tree for the Cdif trehalose/severity data set
#'
#' A subset tree to use as an example
#' @format Ape tree format
#' \describe{
#'   \item{cdif_tree}{This particular tree  is a subset of a tree from:
#'   /nfs/turbo/umms-esnitkin/Project_Cdiff/Analysis/Propensity_paper/clinical_cdifficile_trehalose_variants/data/inputs/input_tree.tree }
#'  }
"cdif_tree"

# Parsed SNP matrices ----
#' Parsed SNP matrices 1
#' @format list
"parsed_snp_cdif1"

#' Parsed SNP matrices 2
#' @format list
"parsed_snp_cdif2"

#' Parsed SNP matrices 3
#' @format list
"parsed_snp_cdif3"

#' Parsed SNP matrices 4
#' @format list
"parsed_snp_cdif4"

#' Parsed SNP matrices 5
#' @format list
"parsed_snp_cdif5"

#' Parsed SNP matrices 6
#' @format list
"parsed_snp_cdif6"

#' Parsed SNP matrices 7
#' @format list
"parsed_snp_cdif7"

#' Parsed SNP matrices 8
#' @format list
"parsed_snp_cdif8"

# Parsed indel matrices -----

#' Parsed indel matrices 1
#' @format list
"parsed_indel_cdif1"

#' Parsed indel matrices 2
#' @format list
# "parsed_indel_cdif2"

#' Parsed indel matrices 3
#' @format list
"parsed_indel_cdif3"

#' Parsed indel matrices 4
#' @format list
"parsed_indel_cdif4"

#' Parsed indel matrices 5
#' @format list
"parsed_indel_cdif5"

#' Parsed indel matrices 6
#' @format list
# "parsed_indel_cdif6"

#' Parsed indel matrices 7
#' @format list
# "parsed_indel_cdif7"

#' Parsed indel matrices 8
#' @format list
#' \describe{
#'   \item{parsed_indel_cdif8}{
#'   parse_indels(varmat_code = cdif_indel_code_mat,
#'                              varmat_allele = cdif_indel_allele_mat,
#'                              tree = cdif_tree,
#'                              og = NULL,
#'                              remove_multi_annots = FALSE,
#'                              return_binary_matrix = FALSE,
#'                              ref_to_anc = FALSE,
#'                              keep_conf_only = FALSE,
#'                              mat_suffix = "_R1_001.fastq.gz")}
#'  }
"parsed_indel_cdif8"

