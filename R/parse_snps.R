#-----------------------------------------------------------------------------
# SUB FUNCTIONS
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# remove_rows_with_bugs
#-----------------------------------------------------------------------------

#' remove_rows_with_bugs
#' @description Removes rows from matrix in which row.names have bugs. Bugs include
#' row annotations with 1) warnings, 2) incorrect number of pipes (should be a
#' multiple of 9), 3) the string CHR_END (we want this to be the last gene in the
#' annotated genome instead of CHR_END), 4) no strand or locus tag information or 5) rows
#' annotated with "None". Eventually can get rid of this when all the bugs are fixed,
#' or can keep as a sanity check to make sure we aren't seeing these bugs in the
#' row annotation. You should run snpmat_code and snpmat_allele through this function separately
#' and should expect the same rows to be removed.
#'
#' @param snpmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return returns a snpmat (class = data.frame) with rows with bugs in the
#' annotation removed. Also writes a file called YEAR_MONTH_DATE_rows_removed_from_SNPMATNAME_due_to_bugs.txt
#' logging the row names of the removed rows.
#'
#' @export
#'
#' @examples

remove_rows_with_bugs <- function(snpmat){
  library(magrittr)
  library(Biostrings)
  library(stringr)

  # Intialize a filename to log the removed rows
  filename = paste0(Sys.Date(), '_rows_removed_from_', deparse(substitute(snpmat)), 'due_to_bugs.txt')

  # 1. Remove rows with warnings in the row annotation
  rows_with_warnings = grep('WARNING', row.names(snpmat))
  write(row.names(snpmat)[rows_with_warnings], file = filename, append = TRUE)
  if (length(rows_with_warnings) > 0){
    snpmat = snpmat[-rows_with_warnings,]
  }

  # 2. Remove rows with the incorrect number of pipes in the row annotation
  # number of |
  num_pipes = str_count(row.names(snpmat), '\\|')
  table(num_pipes) # remove not intervals of 9

  num_semicolon = str_count(row.names(snpmat), ';')
  table(num_semicolon)

  write(row.names(snpmat)[(num_pipes/(num_semicolon-1))%%9 != 0], file = filename, append = TRUE)

  # only keep rows with correct number of pipes
  snpmat = snpmat[(num_pipes/(num_semicolon-1))%%9 == 0,] # must be a multiple of 9

  # 3. Remove rows that still have 'CHR_END' in the row annotation
  rows_with_chr_end = grep('CHR_END', row.names(snpmat))
  write(row.names(snpmat)[rows_with_chr_end], file = filename, append = TRUE)

  if (length(rows_with_chr_end) > 0){
    snpmat = snpmat[-rows_with_chr_end,]
  }

  # 4. Remove rows with not enough locus tag information - have to wait for Ali
  # to convert all reported gene symbols to locus_tags
  # split_annotations <- strsplit(row.names(snpmat), ";")
  # sapply(split_annotations, function(split){
  #   annot = split[1]
  #   locus_tag = gsub('^.+locus_tag=','',annot) %>% gsub(' Strand .*$','',.)
  # })

  #5. Remove rows with no strand information or locus tag information in the
  # row annotations
  locus_tag = unname(sapply(row.names(snpmat), function(row){
    gsub('^.+locus_tag=','',row) %>% gsub(' Strand .*$','',.)
  }))
  no_locus_tag_listed = grep('NULL', locus_tag)

  no_strand_info_listed = grep('No Strand Information found', row.names(snpmat))

  remove_bc_lack_of_info = union(no_locus_tag_listed, no_strand_info_listed)
  write(row.names(snpmat)[remove_bc_lack_of_info], file = filename, append = TRUE)
  if (length(remove_bc_lack_of_info) > 0){
    snpmat = snpmat[-remove_bc_lack_of_info,]
  }

  #6. Remove rows with "None" annotation
  remove_bc_none_annotation = grep('None', row.names(snpmat))
  write(row.names(snpmat)[remove_bc_none_annotation], file = filename, append = TRUE)
  if (length(remove_bc_none_annotation) > 0){
    snpmat = snpmat[-remove_bc_none_annotation,]
  }
  return(snpmat)
}#end remove_rows_with_bugs

#-----------------------------------------------------------------------------
# remove_rows_with_no_variants_or_completely_masked
#-----------------------------------------------------------------------------

#' remove_rows_with_no_variants_or_completely_masked
#' @description Removes rows that have no variants (the same allele in every
#' sample +/- N and dash) or rows that are completely masked (all Ns).
#'
#' @param snpmat_code - data.frame where the rows are variants (numeric description
#' variants: numbers ranging from -4 to 3), the columns are genomes, and the
#' row.names are annotations
#' @param snpmat_allele - data.frame where the rows are variants (character description
#' variants: A,C,T,G,N,-), the columns are genomes, and the row.names are annotations
#'
#' @return Returns a list with the following elements (in order):
#' 1. snpmat_code (class = data.frame) with non-variant rows removed. All rows have at least
#' one sample with a variant.
#' 2. snpmat_allele (class = data.frame) with non-variant rows removed. All rows have at least
#' one sample with a variant.
#' 1 and 2 should be the same dimensions and have the same row names.
#' Also writes a file called YEAR_MONTH_DATE_rows_removed_because_no_variants
#' logging the row names of the removed rows.
#' @export
#'
#' @examples

remove_rows_with_no_variants_or_completely_masked <- function(snpmat_code, snpmat_allele){
  rows_with_one_allele_or_all_Ns_or_dashes = apply(snpmat_allele, 1, function(row){
    length(unique(row)) == 1
  })

  file = paste0(Sys.Date(), '_rows_removed_because_no_variants')
  write(x = row.names(snpmat_code)[rows_with_one_allele_or_all_Ns_or_dashes], file = file)

  snpmat_code_rows_removed = snpmat_code[!rows_with_one_allele_or_all_Ns_or_dashes,]
  snpmat_allele_rows_removed = snpmat_allele[!rows_with_one_allele_or_all_Ns_or_dashes,]

  return(list(snpmat_code_rows_removed, snpmat_allele_rows_removed))
}# end remove_rows_with_no_variants_or_completely_masked

#-----------------------------------------------------------------------------
# split_rows_with_multiple_annots
#-----------------------------------------------------------------------------

#' split_rows_with_multiple_annots
#' @description Rows that have X number of annotations are replicated X number of times.
#' Multiple annotations can be due to 1) multiallelic sites which will result in an
#' annotation for each individual allele and 2) SNPs in overlapping genes which will result in
#' an annotation for each gene or 3) a combination of 1 and 2. The row names will be
#' changed to have one annotation per row (that is, multiallelic SNPs are represented
#' as biallelic sites and each snp in a gene that overlaps with another gene will
#' be represented on a single line). The contents of the data.frame is replicated --
#' that is, the contents of the replicated rows are NOT changed. You should run
#' snpmat_code and snpmat_allele through this function separately and should expect
#' the data.frames to have the same dimensions and same duplicated rows.
#'
#' @param snpmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return Returns a list with the following elements (in order):
#' 1. rows_with_multiple_annots_log - a logical vector with length of
#' nrow(snpmat_added) indicating which rows once had multiple annotations (that is,
#' were split from one row into multiple rows)
#' 2. rows_with_mult_var_allele_log -  a logical vector with length of
#' nrow(snpmat_added) indicating which rows once had multiple annotations in the
#' form of  multiallelic sites (that is, were split from multiallelic sites to
#' biallelic sites)
#' 3. rows_with_overlapping_genes_log -> -  a logical vector with length of
#' nrow(snpmat_added) indicating which rows once had multiple annotations in the
#' form of overlapping genes (that is, were split from a SNP in multiple genes
#' to each gene being represented on a single line)
#' 4. split_rows_flag - an integer vector indicating which rows were split from
#' from a row with multiple annotations (For example if snpmat had 4 rows: 1, 2, 3, 4
#' with row 2 having 3 annotations and row 4 having 2 annotations, the vector would
#' be 1 2 2 2 3 4 4).
#' 5. snpmat_added - a data.frame where the rows are variants, the columns are genomes,
#' and the row.names are SPLIT annotations (each overlapping gene and multiallelic site
#' represented as a single line).
#' @export
#'
#' @examples

split_rows_with_multiple_annots <- function(snpmat){

  num_dividers <- sapply(1:nrow(snpmat), function(x) lengths(regmatches(row.names(snpmat)[x], gregexpr(";[A,C,G,T]", row.names(snpmat)[x]))))

  rows_with_multiple_annotations <- c(1:nrow(snpmat))[num_dividers >= 1 & str_count(row.names(snpmat), '\\|') > 9]

  # Get rows with multallelic sites
  rows_with_multi_allelic_sites = grep('^.+> [A,C,T,G],[A,C,T,G]', row.names(snpmat))

  # Get SNVs present in overlapping genes
  split_annotations <- strsplit(row.names(snpmat)[rows_with_multiple_annotations], ";")

  num_genes_per_site = sapply(split_annotations, function(annots){
    unique(sapply(2:length(annots), function(i){
      unlist(str_split(annots[i], '[|]'))[4]
    }))
  })
  rows_with_overlapping_genes = rows_with_multiple_annotations[sapply(num_genes_per_site, length) > 1]

  # Duplicate rows with multiallelic sites
  row_indices = 1:nrow(snpmat)

  snpmat_added = snpmat[rep(row_indices, num_dividers),]

  # When rows are duplicated .1, .2, .3, etc are added to the end
  # (depending on how many times they were duplicated)
  # Remove to make the duplicated rows have the exact same name
  #names_of_rows = row.names(snpmat_added) %>% gsub(';\\.[0-9].*$', ';', .)

  split_rows_flag = rep(row_indices, num_dividers)

    #as.numeric(factor(names_of_rows, levels = unique(names_of_rows))) # could be because multiple annotations, could be because different variants

  # print(rep(row_indices, num_dividers)[5990:6010])
  # print(split_rows_flag[5990:6010])
  # print(rownames(snpmat)[5952])

  dup = unique(split_rows_flag[duplicated(split_rows_flag)]) # rows that were duplicated

  split_annotations <- strsplit(row.names(snpmat_added)[split_rows_flag %in% dup], ";")


  # FIX ANNOTS OF SNP MAT ADDED - RELIES ON THE .1, .2, .3, ... etc flag
  row.names(snpmat_added)[split_rows_flag %in% dup] =  sapply(split_annotations, function(r){
    if(length(r) == 3){
      paste(r[1], r[2], sep = ';')
    }else if(length(r) > 3 & length(str_split(r[length(r)], '')[[1]]) > 2){
      paste(r[1], r[2], sep = ';')
    }else{
      index = as.numeric(gsub('\\.','',r[length(r)]))
      paste(r[1], r[index+2], sep = ';')
    }
  })

  rows_with_multiple_annots_log = split_rows_flag %in% rows_with_multiple_annotations
  rows_with_mult_var_allele_log = split_rows_flag %in% rows_with_multi_allelic_sites
  rows_with_overlapping_genes_log = split_rows_flag %in% rows_with_overlapping_genes

  # FIX ANNOTS OF SNP MAT ADDED - ROWS WITH MULT VAR ALLELE
  row.names(snpmat_added)[rows_with_mult_var_allele_log] = sapply(row.names(snpmat_added)[rows_with_mult_var_allele_log], function(r){
    if(grepl('> [A,C,T,G],[A,C,T,G].*functional=', r)){
      var =  gsub('^.*Strand Information:','',r) %>% gsub('\\|.*$', '', .) %>% substr(.,nchar(.),nchar(.))
      gsub('> [A,C,T,G],[A,C,T,G].*functional=', paste('>', var, 'functional='), r)
    }
  })



  return(list(rows_with_multiple_annots_log,
              rows_with_mult_var_allele_log,
              rows_with_overlapping_genes_log,
              split_rows_flag,
              snpmat_added))


} # end split_rows_with_multiple_annots

#-----------------------------------------------------------------------------
# remove_rows_with_multiple_annots
#-----------------------------------------------------------------------------

#' remove_rows_with_multiple_annots
#' @description - Bypass dealing with rows with multiple annotations (due to
#' overlapping genes or multiallelic sites) by removing them from the data.frame.
#' Useful especially as we are testing this function and the functionality to deal
#' with sites with multiple annotations is not ready.
#'
#' @param snpmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return Returns a snpmat (class = data.frame) with all rows with multiple annotations
#' removed. Also writes a file called YEAR_MONTH_DATE_rows_with_multiple_annots_removed
#' indicating which rows were removed.
#'
#' @export
#'
#' @examples

remove_rows_with_multiple_annots <- function(snpmat){
  # IDENTIFY ROWS WITH MULTIPLE ANNOTATIONS
  num_dividers <- sapply(1:nrow(snpmat), function(x) lengths(regmatches(row.names(snpmat)[x], gregexpr(";[A,C,G,T]", row.names(snpmat)[x]))))
  rows_with_multiple_annotations <- c(1:nrow(snpmat))[num_dividers >= 2 & str_count(row.names(snpmat), '\\|') > 9]

  # SAVE TO LOG FILE
  log_file = paste0(Sys.Date(), '_rows_with_multiple_annots_removed')
  write('The following rows with multiple annotations were removed:', log_file)
  write(row.names(snpmat)[rows_with_multiple_annotations], log_file, append = TRUE)

  # REMOVE ROWS WITH MULTIPLE ANNOTATIONS
  if(length(rows_with_multiple_annotations)>0){
    snpmat = snpmat[-rows_with_multiple_annotations,]
  }

  return(snpmat)
} # end remove_rows_with_multiple_annots

#-----------------------------------------------------------------------------
# get_info_from_annotations
#-----------------------------------------------------------------------------

#' get_info_from_annotations
#' @description Parse annotations from row names.
#' @param snpmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return data.frame of length nrow(snpmat) containing the following columns:
#' label - string indicating "Coding SNP or Non-Coding SNP"
#' pos - position of variant in the reference genome
#' phage -  the word NULL or the word PHAGE
#' repeated_region -  the word NULL or the word REPEAT
#' masked -  the word NULL or the word MASKED
#' locus_tag - locus tag from gff file
#' strand_info - returns strand information (currently incorrect until Ali fixes a bug)
#' strand - returns + or - depending if gene is on positive or negative strand
#' (determining this myself)
#' ref - allele in the reference genome
#' var - variant allele in terms of the positive strand
#' aa_change - amino acid change in the form p.Tyr95Tyr or p.Glu75Ala
#' variant_type - snpeff determined variant type (e.g. missense_variant, intergenic_region
#' see snpeff manual for more info: http://snpeff.sourceforge.net/SnpEff_manual.html)
#' snpeff_impact - snpeff determined impacted (e.g. LOW, MODERATE, MODIFIER, HIGH)
#' nuc_pos_in_gene - variant position relative to the gene length
#' aa_pos_in_gene - position of the codon-containing-variant relative to the gene length
#' gene_length_in_bp - gene length in nucleotide base pairs
#' annotation_1 - information about variant
#' annotation_2 - information about variant
#' ig_gene1 - if intergenic variant, first intergenic locus tag surrounding the variant
#' ig_gene2 - if intergenic variant, second intergenic locus tag surrounding the variant
#' intergenic - logical indicating if the variant is in an intergenic region
#'
#' @export
#'
#' @examples

get_info_from_annotations <- function(snpmat){
  library(magrittr)
  library(Biostrings)
  library(stringr)

  # GET REF AND VAR ALLELE, STRAND
  label = rep(NA, nrow(snpmat))

  pos = rep(NA, nrow(snpmat))

  phage = rep(NA, nrow(snpmat))
  repeated_region = rep(NA, nrow(snpmat))
  masked= rep(NA, nrow(snpmat))

  locus_tag = rep(NA, nrow(snpmat))

  strand_info = rep(NA, nrow(snpmat))

  ref = rep(NA, nrow(snpmat))
  var = rep(NA, nrow(snpmat))

  aa_change = rep(NA, nrow(snpmat))

  variant_type = rep(NA, nrow(snpmat))

  snpeff_impact = rep(NA, nrow(snpmat))

  nuc_pos_in_gene =  rep(NA, nrow(snpmat))
  aa_pos_in_gene = rep(NA, nrow(snpmat))

  gene_length_in_bp = rep(NA, nrow(snpmat))

  annotation_1 = rep(NA, nrow(snpmat))
  annotation_2 = rep(NA, nrow(snpmat))

  strand = rep(NA, nrow(snpmat))

  ig_gene1 = rep(NA, nrow(snpmat))
  ig_gene2 = rep(NA, nrow(snpmat))
  intergenic = rep(NA, nrow(snpmat))


  for (i in 1:nrow(snpmat)){

    row = row.names(snpmat)[i]

    split_row = unlist(str_split(row, pattern = '[|]'))

    # Coding SNP, Non coding SNP - for checking
    label[i] = gsub(' at [1-9].*$','', split_row[1])

    # Position in genome
    pos[i] = gsub('^.*at ','', split_row[1]) %>% gsub(' > [A,C,T,G].*$', '', .)

    # PHAGE, REPEAT, MASK
    functional_temp = gsub('^.*functional=','',  split_row[1]) %>% gsub(' locus_tag.*$','',.) %>% str_split(., '_') %>% unlist
    phage[i] = functional_temp[1]
    repeated_region[i] = functional_temp[2]
    masked[i] = functional_temp[3]

    # STRAND INFO
    strand_info[i] = gsub('^.*Strand Information: ','',split_row[1]) %>% gsub(';.*$','',.)

    # VARIANT TYPE
    variant_type[i] = split_row[2]

    # SNPEFF IMPACT
    snpeff_impact[i] = split_row[3]

    # LOCUS TAG (OR GENE SYMBOL UNTIL ERROR IS FIXED)
    locus_tag[i] = split_row[4]

    # REF AND VAR - in terms of the positive strand
    var_1 = substr(split_row[1], nchar(split_row[1]), nchar(split_row[1]))
    var_2 = substr(split_row[5], nchar(split_row[5]), nchar(split_row[5]))

    var[i] = var_1

    ref_temp = substr(split_row[5], nchar(split_row[5])-2, nchar(split_row[5])-2)

    if(var_1 != var_2){
      ref[i] = as.character(complement(DNAString(ref_temp)))
      strand[i] = '-'
    }else{
      ref[i] = ref_temp
      strand[i] = '+'
    }

    # AMINO ACID CHANGE
    aa_change[i] = split_row[6]

    # GENE LENGTH AND POSITION OF MUTATION IN RELATION TO THE GENE
    nuc_pos_in_gene[i] = (str_split(split_row[7], '/') %>% unlist())[1]
    gene_length_in_bp[i] = (str_split(split_row[7], '/') %>% unlist())[2]
    aa_pos_in_gene[i] = (str_split(split_row[8], '/') %>% unlist())[1]

    # ANNOTATIONS
    annotation_1[i] = split_row[9]
    annotation_2[i] = split_row[10]

    # INTERGENIC REGIONS
    if(variant_type[i] == 'intergenic_region'){
      ig_gene1[i] = (str_split(split_row[4], '-') %>% unlist())[1]
      ig_gene2[i] = (str_split(split_row[4], '-') %>% unlist())[2]
      intergenic[i] = TRUE
    }else{
      ig_gene1[i] = ''
      ig_gene2[i] = ''
      intergenic[i] = FALSE
    }

  }

  annotations = data.frame(label, pos, phage, repeated_region, masked, locus_tag, strand_info, strand,
                           ref, var, aa_change, variant_type, snpeff_impact, nuc_pos_in_gene,
                           aa_pos_in_gene, gene_length_in_bp, annotation_1, annotation_2,
                           ig_gene1, ig_gene2, intergenic)
  return(annotations)
}# end get_info_from_annotations

#-----------------------------------------------------------------------------
# remove_alt_allele_code_from_split_rows
#-----------------------------------------------------------------------------

#' remove_alt_allele_code_from_split_rows
#' @description Input the split matrix where rows that once had multiple annotations on
#' single line are now represented on multiple lines. For the sites that once were
#' multiallelic sites and are now represented as biallelic, thus function will
#' change the contents of snpmat_code such that the alternative allele(s) are 0.
#' For example, T -> G, C is split into two lines: T -> G and T -> C. In the code matrix,
#' turn the codes correspoding to the allele C in the row T -> G to 0 and the codes
#' corresponding to the allele G in the row T -> C to 0.
#'
#' @param snpmat_code_split - data.frame where the rows are variants (numeric description
#' variants: numbers ranging from -4 to 3), the columns are genomes,
#' and the row.names are annotations, and each line has a single annotation
#' @param snpmat_allele_split - data.frame where the rows are variants (character description
#' variants: A,C,T,G,N,-), the columns are genomes, and the row.names are annotations,
#' and each line has a single annotation
#' @param ref - character vector length nrow(snpmat_code_split) = nrow(snpmat_allele_split) indicating
#' the reference allele in terms of the positive strand
#' @param var - character vector length nrow(snpmat_code_split) = nrow(snpmat_allele_split) indicating
#' the variant allele in terms of the positive strand
#' @param rows_with_mult_var_allele_log - logical vector length nrow(snpmat_code_split) = nrow(snpmat_allele_split)
#' indicating which rows are multiallelic sites
#'
#' @return - snpmat_code - data.frame where the rows are variants (numeric description
#' variants: numbers ranging from -4 to 3), the columns are genomes,
#' and the row.names are annotations, and each line has a single annotation where
#' the alternative/minor allele in a biallelic-represrentation of a multiallelic site is now 0
#' @export
#'
#' @examples

remove_alt_allele_code_from_split_rows <- function(snpmat_code_split, snpmat_allele_split, ref, var, rows_with_mult_var_allele_log){

  # UPDATE CODE MATRIX:
  index_mult_var = (1:length(rows_with_mult_var_allele_log))[rows_with_mult_var_allele_log]

  for (i in index_mult_var){

    # CHANGE THE ALT ALLELE TO 0 IN BIALLELIC REPRESENTATION OF MULTIALLELIC POSITION
    # CHANGE !REF, !VAR, !N, or !-  TO 0
    # T > C,G
    # T > C:   T C G N -; 0 1 1 0 0
    # T > G:   T C G N -; 0 1 1 0 0
    # INTO
    # T > C:   T C G N -; 0 1 0 0 0
    # T > G:   T C G N -; 0 0 1 0 0

    snpmat_code_split[i,!(as.character(snpmat_allele_split[i,]) %in% c(as.character(var[i]), as.character(ref[i]), 'N', '-'))] = 0

  }

  return(snpmat_code_split)

} # end remove_alt_allele_code_from_split_rows


#' Root tree on outgroup
#' @description Root tree based on outgroup. If outgroup is null and the tree isn't rooted, midpoint root tree. If tree is rooted, return tree as is.
#'
#' @param tree phylogenetic tree or file path of tree
#' @param outgroup tip name of outgroup in phylogeny. If NULL, midpoint root if not rooted
#'
#' @return rooted tree without outgroup
#' @export
#'
#' @examples
#' tree = rcoal(100)
#' is.rooted(tree)
#' root_tree_og(tree)
#' is.rooted(tree)
root_tree_og = function(tree,outgroup=NULL){
  library(phytools)
  library(ape)

  if(is.character(tree)){
    # LOAD IN TREE
    tree = read.tree(tree)
  }
  # IF NO OUTGROUP AND TREE IS UNROOTED
  if(is.null(outgroup) & !is.rooted(tree)){
    # MIDPOINT ROOT TREE
    tree = midpoint.root(tree)
  }else if(!is.null(outgroup)){
    # ROOT TREE ON OUTGROUP
    tree = root(tree,outgroup)
    tree = drop.tip(tree,outgroup)
  }
  return(tree)
}

#' Get ancestral state of alleles
#' @description Rereference alleles based on rooted tree
#'
#' @param tree rooted tree
#' @param mat allele matrix (rows are variants, columns are samples)
#'
#' @return matrix of most likely ancestral allele for each row in allele matrix and probability that that is the ancestral state
#' @export
#'
#' @examples
get_anc_alleles = function(tree,mat){
  library(future.apply)
  plan(multiprocess)

  if(sum(!(tree$tip.label %in% colnames(mat))) > 0){
    stop('Some samples in tree are not in allele matrix.')
  }

  if(sum(!(colnames(mat) %in% tree$tip.label)) > 0){
    stop('Some samples in allele matrix are not in tree.')
  }

  if(!is.rooted(tree)){
    stop('Tree must be rooted.')
  }

  # ORDER MATRIX TO MATCH TREE TIP LABELS
  mat = mat[,tree$tip.label]

  if(sum(tree$edge.length == 0) > 0){
    warning('All zero branch lengths changed to small non-zero number to be able to perform ancestral reconstruction.')
    # Change any edge lengths that are zero to a very small number (so ancestral reconstruction doesn't break)
    tree$edge.length[tree$edge.length == 0] = min(tree$edge.length[tree$edge.length > 0])/1000
  }

  # Get ancestral state of root; 1st column = var absent (0), 2nd column = var present (1)
  ar_all = t(future_apply(mat,1,function(tip_states){
    tip_state = unique(tip_states)
    if(length(tip_state) > 1){
      ar = ace(x = tip_states,phy = tree,type = 'discrete')
      states = ar$lik.anc[1,]
      tip_state = names(states)[which.max(states)]
      prob = states[which.max(states)]
      c(tip_state,prob)
    }else{
      c(tip_states,1)
    }
  }))

  return(ar_all)

}

#' Load matrix from path if needed
#'
#' @param mat - loaded in data.frame of snpmat or character string of a path to a snpmat
#' @description Loads variant matrix from path if not already loaded
#'
#' @return variant matrix
#' @export
#'
#' @examples
load_if_path = function(mat){
  if(is.character(mat)){
    mat = read.table(mat,
               header = TRUE,
               stringsAsFactors = FALSE,
               sep = "\t",
               quote = "",
               row.names = 1)
  }
  return(mat)
}

#' Remove unknown ancestral states
#' @description Remove rows from variant matrix where the ancestral state is unknown (- or N)
#'
#' @param snpmat_code
#' @param snpmat_allele
#' @param annots
#'
#' @return
#' @export
#'
#' @examples
remove_unknown_anc = function(snpmat_code, snpmat_allele, annots){
  unknown = annots$anc %in% c('-','N')
  removed = rownames(snpmat_code)[unknown]
  filename = paste0(Sys.Date(), '_rows_removed_because_unknown_ancestral_state.txt')
  write.table(removed,file=filename,sep='\n',quote=F,row.names=F,col.names=F)
  return(list(snpmat_code=snpmat_code[!unknown,],
              snpmat_allele=snpmat_allele[!unknown,],
              annots=annots[!unknown,]))
}

#-----------------------------------------------------------------------------
# MAIN FUNCTION
#-----------------------------------------------------------------------------

#' parse_snps
#' @description Input matrices generated from internal (Ali's) variant calling pipeline.
#' Always returns parsed annotation info. In addition, you have the option to:
#' 1. split rows with multiple annotations (snps in overlapping genes, multiallelic snps)
#' 2. Re-reference to the ancestral allele at that position (instead of to the reference genome)
#' 3. Simplify the code matrix - which contains numbers from -4 to 3 indicating
#' different information about the variants - to a binary matrix indicating
#' simple presence/absence of a SNP at that site.
#' @param snpmat_code - loaded data.frame or path to the snpmat_code file generated
#' from internal variant calling pipeline
#' @param snpmat_allele - loaded data.frame or path to the snpmat_allele file generated
#' from internal variant calling pipeline
#' @param tree - optional: path to tree file or loaded in tree (class = phylo)
#' @param og - optional: character string of the name of the outgroup (has to match what
#' it is called in the tree)
#' @param remove_multi_annots - logical flag indicating if you want to remove rows
#' with multiple annotations - alternative is to split rows with mutliple annotations (default = FALSE)
#' @param return_binary_matrix - logical flag indicating if you want to return a binary matrix (default = FALSE)
#'
#' @return list of allele mat, code mat, binary mat and corresponding parsed annotations.
#' output will depend on arguments to the function.
#' @export
#'
#' @examples

parse_snps <- function(snpmat_code, snpmat_allele, tree=NULL, og = NULL, remove_multi_annots = FALSE, return_binary_matrix = FALSE, ref_to_anc = T){

  if(is.null(tree) & return_binary_matrix & ref_to_anc){
    stop('Tree file required when returning a binary matrix.')
  }

  #-----------------------------------------------------------------------------
  # READ IN SNPMAT CODE AND SNPMAT ALLELE
  #-----------------------------------------------------------------------------

  snpmat_code <- load_if_path(snpmat_code)

  snpmat_allele <- load_if_path(snpmat_allele)

  # add semicolons to the end of the row names that don't have semicolons
  row.names(snpmat_code)[!grepl(';$', row.names(snpmat_code))] = paste0(row.names(snpmat_code)[!grepl(';$', row.names(snpmat_code))], ';')
  row.names(snpmat_allele)[!grepl(';$', row.names(snpmat_allele))] = paste0(row.names(snpmat_allele)[!grepl(';$', row.names(snpmat_allele))], ';')



  #-----------------------------------------------------------------------------
  # REMOVE BUGS
  #-----------------------------------------------------------------------------
  snpmat_code = remove_rows_with_bugs(snpmat_code)
  snpmat_allele = remove_rows_with_bugs(snpmat_allele)

  #-----------------------------------------------------------------------------
  # REMOVE LINES WITH NO VARIANTS - NO VARIANT OR ALL MASKED
  #-----------------------------------------------------------------------------
  snpmats = remove_rows_with_no_variants_or_completely_masked(snpmat_code, snpmat_allele)
  snpmat_code = snpmats[[1]]
  snpmat_allele = snpmats[[2]]

  #-----------------------------------------------------------------------------
  # EITHER (1) REMOVE ROWS WITH MULTIPLE ANNOTATIONS OR (2) SPLIT ROWS WITH
  # MULTIPLE ANNOTATIONS - DEPENDING ON VALUE OF REMOVE_MULTI_ANNOTS FLAG
  # (TRUE/FALSE)
  #-----------------------------------------------------------------------------
  if(remove_multi_annots){
    # REMOVE ROWS WITH MULTIPLE ANNOTATIONS
    snpmat_code = remove_rows_with_multiple_annots(snpmat_code)
    snpmat_allele = remove_rows_with_multiple_annots(snpmat_allele)

    #-----------------------------------------------------------------------------
    # FIND ANCESTRAL STATE OF EACH ALLELE
    #-----------------------------------------------------------------------------

    major_alleles = get_major_alleles(data.matrix(snpmat_allele))

    if(return_binary_matrix){
      # REROOT TREE
      tree = root_tree_og(tree)
      # GET ANCESTRAL ALLELE FOR EACH VARIANT
      if(ref_to_anc){
        alleles = get_anc_alleles(tree, snpmat_allele)
      }else{
        # REFERENCE TO MAJOR ALLELE
        alleles = major_alleles
      }

    }

    split_rows_flag = 1:nrow(snpmat_allele)

    rows_with_multiple_annots_log = rep(FALSE, nrow(snpmat_allele))
    rows_with_mult_var_allele_log = rep(FALSE, nrow(snpmat_allele))
    rows_with_overlapping_genes_log = rep(FALSE, nrow(snpmat_allele))

    # GET ANNOTATIONS
    annots = cbind(get_info_from_annotations(snpmat_code), rows_with_multiple_annots_log,
                   rows_with_mult_var_allele_log, rows_with_overlapping_genes_log,
                   split_rows_flag)

    annots$maj = major_alleles

  }else{
      #-----------------------------------------------------------------------------
    # FIND ANCESTRAL STATE OF EACH ALLELE
    #-----------------------------------------------------------------------------

    major_alleles = get_major_alleles(snpmat_allele)

    if(return_binary_matrix){
      if(ref_to_anc){
        # REROOT TREE
        tree = root_tree_og(tree)

        # GET ANCESTRAL ALLELE FOR EACH VARIANT
        alleles = get_anc_alleles(tree, snpmat_allele)

      }else{
        # REFERENCE TO MAJOR ALLELE
        alleles = major_alleles
      }

    }

  # RAW ROWNAMES
  raw_rownames = row.names(snpmat_code)

  # SPLIT MATRICES
  snpmat_code_split_list = split_rows_with_multiple_annots(snpmat_code)
  snpmat_allele_split_list = split_rows_with_multiple_annots(snpmat_allele)

  snpmat_code = snpmat_code_split_list[[5]]
  snpmat_allele = snpmat_allele_split_list[[5]]

  rows_with_multiple_annots_log = snpmat_code_split_list[[1]]
  rows_with_mult_var_allele_log = snpmat_code_split_list[[2]]
  rows_with_overlapping_genes_log = snpmat_code_split_list[[3]]
  split_rows_flag = snpmat_code_split_list[[4]]

  if(return_binary_matrix){
    alleles = alleles[split_rows_flag,]
  }

  major_alleles = major_alleles[split_rows_flag]

  # EXPAND RAW ROW NAMES
  raw_rownames = raw_rownames[split_rows_flag]

  # GET ANNOTATIONS
  annots = cbind(get_info_from_annotations(snpmat_code), rows_with_multiple_annots_log,
                 rows_with_mult_var_allele_log, rows_with_overlapping_genes_log,
                 split_rows_flag,maj=major_alleles, raw_rownames = raw_rownames)

  # CHANGE SNPMAT CODE TO REFLECT BIALLELIC REPRESENTATION OF A MULTIALLELIC SITE
  snpmat_code = remove_alt_allele_code_from_split_rows(snpmat_code,
                                         snpmat_allele,
                                         annots$ref,
                                         annots$var,
                                         rows_with_mult_var_allele_log)

  }

  #annots$maj = major_alleles

  if(return_binary_matrix){
    if(ref_to_anc){
      # ADD ANCESTRAL ALLELE INFO TO ANNOTATIONS
      annots$anc = alleles[,1]
      annots$anc_prob = alleles[,2]

      # remove sites with unknown ancestor
      snpmats = remove_unknown_anc(snpmat_code, snpmat_allele, annots)
      snpmat_code = snpmats$snpmat_code
      snpmat_allele = snpmats$snpmat_allele
      annots = snpmats$annots

      # MAKE BINARY MATRIX
      snpmat_bin = snpmat_code
      to_keep = !(rowSums(snpmat_bin ==  2) > 0 |
                    rowSums(snpmat_bin == -2) > 0 |
                    rowSums(snpmat_bin == -3) > 0 |
                    rowSums(snpmat_bin == -4) > 0)
      snpmat_bin = snpmat_bin[to_keep,]
      snpmat_bin[snpmat_bin == 3] = 1
      snpmat_bin[snpmat_bin == -1] = 0

      annots_bin = annots[to_keep,]

      snpmat_bin_reref = data.frame(t(sapply(1:nrow(snpmat_bin), function(x){
        if(annots_bin$ref[x] == annots_bin$anc[x]){
          unlist(snpmat_bin[x,])
        }else if(!annots_bin$rows_with_mult_var_allele_log[x]){
          unlist(as.numeric(!snpmat_bin[x,]))
        }else if(annots_bin$var[x] == annots_bin$anc[x]){
          unlist(snpmat_bin[x,])
        }else{
          unlist(rep(NA,ncol(snpmat_bin)))
        }
      })))

      reref = sapply(1:nrow(snpmat_bin), function(x){
        if(annots_bin$ref[x] == annots_bin$anc[x]){
          'no'
        }else if(!annots_bin$rows_with_mult_var_allele_log[x]){
          'yes'
        }else if(annots_bin$var[x] == annots_bin$anc[x]){
          'no'
        }else{
          'complicated'
        }
      })
    }else{




    # MAKE BINARY MATRIX
    snpmat_bin = snpmat_code
    to_keep = !(rowSums(snpmat_bin ==  2) > 0 |
                  rowSums(snpmat_bin == -2) > 0 |
                  rowSums(snpmat_bin == -3) > 0 |
                  rowSums(snpmat_bin == -4) > 0)
    snpmat_bin = snpmat_bin[to_keep,]
    snpmat_bin[snpmat_bin == 3] = 1
    snpmat_bin[snpmat_bin == -1] = 0

    annots_bin = annots[to_keep,]

    snpmat_bin_reref = data.frame(t(sapply(1:nrow(snpmat_bin), function(x){
      if(annots_bin$ref[x] == annots_bin$maj[x]){
        unlist(snpmat_bin[x,])
      }else if(!annots_bin$rows_with_mult_var_allele_log[x]){
        unlist(as.numeric(!snpmat_bin[x,]))
      }else if(annots_bin$var[x] == annots_bin$maj[x]){
        unlist(snpmat_bin[x,])
      }else{
        unlist(rep(NA,ncol(snpmat_bin)))
      }
    })))

    reref = sapply(1:nrow(snpmat_bin), function(x){
      if(annots_bin$ref[x] == annots_bin$maj[x]){
        'no'
      }else if(!annots_bin$rows_with_mult_var_allele_log[x]){
        'yes'
      }else if(annots_bin$var[x] == annots_bin$maj[x]){
        'no'
      }else{
        'complicated'
      }
    })
    }

    parsed = list(code=list(mat=snpmat_code,annots=annots),
         allele=list(mat=snpmat_allele,annots=annots),
         bin=list(mat=snpmat_bin_reref,annots=cbind(annots_bin,reref=reref)))
    save(parsed, file = 'SNP_parsed.RData')
    return(parsed)
  }

  parsed = list(code=list(mat=snpmat_code, annots=annots),
       allele=list(mat=snpmat_allele, annots=annots))
  save(parsed, file = 'SNP_parsed.RData')
  return(parsed)
}# end parse_snps
