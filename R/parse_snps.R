#' Grab information from each annotation
#' @description Parse annotations from row names.
#' @param varmat - data.frame where the rows are variants, the columns are
#'   genomes, and the row.names are annotations
#'
#' @return data.frame of length nrow(varmat) containing the following columns:
#'   label - string indicating "Coding SNP or Non-Coding SNP" pos - position of
#'   variant in the reference genome phage -  the word NULL or the word PHAGE
#'   repeated_region -  the word NULL or the word REPEAT masked -  the word NULL
#'   or the word MASKED locus_tag - locus tag from gff file strand_info -
#'   returns strand information (currently incorrect until Ali fixes a bug)
#'   strand - returns + or - depending if gene is on positive or negative strand
#'   (determining this myself) ref - allele in the reference genome var -
#'   variant allele in terms of the positive strand aa_change - amino acid
#'   change in the form p.Tyr95Tyr or p.Glu75Ala variant_type - snpeff
#'   determined variant type (e.g. missense_variant, intergenic_region see
#'   snpeff manual for more info:
#'   http://snpeff.sourceforge.net/SnpEff_manual.html) snpeff_impact - snpeff
#'   determined impacted (e.g. LOW, MODERATE, MODIFIER, HIGH) nuc_pos_in_gene -
#'   variant position relative to the gene length aa_pos_in_gene - position of
#'   the codon-containing-variant relative to the gene length gene_length_in_bp
#'   - gene length in nucleotide base pairs annotation_1 - information about
#'   variant annotation_2 - information about variant ig_gene1 - if intergenic
#'   variant, first intergenic locus tag surrounding the variant ig_gene2 - if
#'   intergenic variant, second intergenic locus tag surrounding the variant
#'   intergenic - logical indicating if the variant is in an intergenic region
#'
#' @export
get_snp_info_from_annotations <- function(varmat){
  # GET REF AND VAR ALLELE, STRAND
  label <- pos <- phage <- repeated_region <- masked <- locus_tag <-
    strand_info <- ref <- var <- aa_change <- variant_type <- snpeff_impact <-
    nuc_pos_in_gene <- aa_pos_in_gene <-  gene_length_in_bp <- annotation_1 <-
    annotation_2 <- strand <- ig_gene1 <- ig_gene2 <- intergenic <-
    rep(NA, nrow(varmat))

  for (i in 1:nrow(varmat)) {
    row <- row.names(varmat)[i]
    split_row <- unlist(stringr::str_split(row, pattern = "[|]"))

    # Coding SNP, Non coding SNP - for checking
    label[i] <- gsub(" at [1-9].*$", "", split_row[1])

    # Position in genome
    pos[i] <- gsub("^.*at ", "", split_row[1]) %>% gsub(" > [A,C,T,G].*$", "", .)

    # PHAGE, REPEAT, MASK
    functional_temp <- gsub("^.*functional=", "", split_row[1]) %>%
      gsub(" locus_tag.*$", "", .) %>%
      stringr::str_split(., "_") %>%
      unlist
    phage[i] <- functional_temp[1]
    repeated_region[i] <- functional_temp[2]
    masked[i] <- functional_temp[3]

    # STRAND INFO
    strand_info[i] <- gsub("^.*Strand Information: ", "", split_row[1]) %>%
      gsub(";.*$", "", .)

    # VARIANT TYPE
    variant_type[i] <- split_row[2]

    # SNPEFF IMPACT
    snpeff_impact[i] <- split_row[3]

    # LOCUS TAG (OR GENE SYMBOL UNTIL ERROR IS FIXED)
    locus_tag[i] <- split_row[4]

    # REF AND VAR - in terms of the positive strand
    var_1 <- substr(split_row[1], nchar(split_row[1]), nchar(split_row[1]))
    var_2 <- substr(split_row[5], nchar(split_row[5]), nchar(split_row[5]))

    var[i] <- var_1

    ref_temp <- substr(split_row[5], nchar(split_row[5]) - 2, nchar(split_row[5]) - 2)

    if (var_1 != var_2) {
      ref[i] <- as.character(Biostrings::complement(Biostrings::DNAString(ref_temp)))
      strand[i] <- "-"
    } else {
      ref[i] <- ref_temp
      strand[i] <- "+"
    }

    # AMINO ACID CHANGE
    aa_change[i] <- split_row[6]

    # GENE LENGTH AND POSITION OF MUTATION IN RELATION TO THE GENE
    nuc_pos_in_gene[i] <-
      (stringr::str_split(split_row[7], "/") %>% unlist())[1]
    gene_length_in_bp[i] <-
      (stringr::str_split(split_row[7], "/") %>% unlist())[2]
    aa_pos_in_gene[i] <-
      (stringr::str_split(split_row[8], "/") %>% unlist())[1]

    # ANNOTATIONS
    annotation_1[i] <- split_row[9]
    annotation_2[i] <- split_row[10]

    # INTERGENIC REGIONS
    ig_gene1[i] <- ""
    ig_gene2[i] <- ""
    intergenic[i] <- FALSE
    if (variant_type[i] == "intergenic_region") {
      ig_gene1[i] <- (stringr::str_split(split_row[4], "-") %>% unlist())[1]
      ig_gene2[i] <- (stringr::str_split(split_row[4], "-") %>% unlist())[2]
      intergenic[i] <- TRUE
    }
  }

  annotations <- data.frame(label,
                           pos,
                           phage,
                           repeated_region,
                           masked,
                           locus_tag,
                           strand_info,
                           strand,
                           ref,
                           var,
                           aa_change,
                           variant_type,
                           snpeff_impact,
                           nuc_pos_in_gene,
                           aa_pos_in_gene,
                           gene_length_in_bp,
                           annotation_1,
                           annotation_2,
                           ig_gene1,
                           ig_gene2,
                           intergenic)
  return(annotations)
}

#' Parse SNP matrix from Ali's pipeline
#' @description Input matrices generated from internal (Ali's) variant calling
#'   pipeline. Always returns parsed annotation info. In addition, you have the
#'   option to: 1. split rows with multiple annotations (snps in overlapping
#'   genes, multiallelic snps) 2. Re-reference to the ancestral allele at that
#'   position (instead of to the reference genome) 3. Simplify the code matrix -
#'   which contains numbers from -4 to 3 indicating different information about
#'   the variants - to a binary matrix indicating simple presence/absence of a
#'   SNP at that site.
#' @param varmat_code Loaded data.frame or path to the varmat_code file
#'   generated from internal variant calling pipeline
#' @param varmat_allele Loaded data.frame or path to the varmat_allele file
#'   generated from internal variant calling pipeline
#' @param tree Optional: path to tree file or loaded in tree (class = phylo)
#' @param og Optional: character string of the name of the outgroup (has to
#'   match what it is called in the tree)
#' @param remove_multi_annots Logical flag indicating if you want to remove
#'   rows with multiple annotations - alternative is to split rows with mutliple
#'   annotations (default = FALSE)
#' @param return_binary_matrix Logical flag indicating if you want to return a
#'   binary matrix (default = TRUE)
#' @param ref_to_anc Whether to reference to the ancestral allele to create the binary marix (default = TRUE)
#' @param keep_conf_only Logical flag indicating if only confident variants should be kept (1's in Ali's pipeline, otherwise 3's are also kept) (default = TRUE)
#' @param mat_suffix Suffix to remove from code and allele matrices so the names match with the tree tip labels.
#'
#' @return list of allele mat, code mat, binary mat and corresponding parsed
#'   annotations. output will depend on arguments to the function.
#' @export
parse_snps <- function(varmat_code,
                       varmat_allele,
                       tree = NULL,
                       og = NULL,
                       remove_multi_annots = FALSE,
                       return_binary_matrix = TRUE,
                       ref_to_anc = TRUE,
                       keep_conf_only = TRUE,
                       mat_suffix = '_R1_001.fastq.gz|_R1.fastq.gz|_1.fastq.gz',
                       ref_to_maj = FALSE){
  parse_snp_or_indel(varmat_code = varmat_code,
                     varmat_allele = varmat_allele,
                     mat_type = "SNP",
                     tree = tree,
                     og = og,
                     remove_multi_annots = remove_multi_annots,
                     return_binary_matrix = return_binary_matrix,
                     ref_to_anc = ref_to_anc,
                     keep_conf_only = keep_conf_only,
                     mat_suffix = mat_suffix,
                     ref_to_maj = ref_to_maj)
}
