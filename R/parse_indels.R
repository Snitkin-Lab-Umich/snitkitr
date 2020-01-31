#-----------------------------------------------------------------------------
# SUB FUNCTIONS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# remove_rows_with_multiple_annots
#-----------------------------------------------------------------------------

#' remove_rows_with_multiple_annots
#' @description - Bypass dealing with rows with multiple annotations (due to
#' overlapping genes or multiallelic sites) by removing them from the data.frame.
#' Useful especially as we are testing this function and the functionality to deal
#' with sites with multiple annotations is not ready.
#'
#' @param indelmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return Returns a indelmat (class = data.frame) with all rows with multiple annotations
#' removed. Also writes a file called YEAR_MONTH_DATE_rows_with_multiple_annots_removed
#' indicating which rows were removed.
#'
#' @export
#'
#' @examples

remove_rows_with_multiple_annots <- function(indelmat){
  # IDENTIFY ROWS WITH MULTIPLE ANNOTATIONS
  num_dividers <- sapply(1:nrow(indelmat), function(x) lengths(regmatches(row.names(indelmat)[x], gregexpr(";[A,C,G,T]", row.names(indelmat)[x]))))
  rows_with_multiple_annotations <- c(1:nrow(indelmat))[num_dividers >= 2 & str_count(row.names(indelmat), '\\|') > 9]

  # SAVE TO LOG FILE
  log_file = paste0(Sys.Date(), '_rows_with_multiple_annots_removed')
  write('The following rows with multiple annotations were removed:', log_file)
  write(row.names(indelmat)[rows_with_multiple_annotations], log_file, append = TRUE)

  # REMOVE ROWS WITH MULTIPLE ANNOTATIONS
  if(length(rows_with_multiple_annotations)>0){
    indelmat = indelmat[-rows_with_multiple_annotations,]
  }

  return(indelmat)
} # end remove_rows_with_multiple_annots

#-----------------------------------------------------------------------------
# get_info_from_annotations
#-----------------------------------------------------------------------------

#' get_info_from_annotations
#' @description Parse annotations from row names.
#' @param indelmat - data.frame where the rows are variants, the columns are genomes,
#' and the row.names are annotations
#'
#' @return data.frame of length nrow(indelmat) containing the following columns:
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

get_info_from_annotations <- function(indelmat){
  library(magrittr)
  library(Biostrings)
  library(stringr)

  # GET REF AND VAR ALLELE, STRAND
  label = rep(NA, nrow(indelmat))

  pos = rep(NA, nrow(indelmat))

  phage = rep(NA, nrow(indelmat))
  repeated_region = rep(NA, nrow(indelmat))
  masked= rep(NA, nrow(indelmat))

  locus_tag = rep(NA, nrow(indelmat))

  strand_info = rep(NA, nrow(indelmat))

  ref = rep(NA, nrow(indelmat))
  var = rep(NA, nrow(indelmat))

  aa_change = rep(NA, nrow(indelmat))

  variant_type = rep(NA, nrow(indelmat))

  snpeff_impact = rep(NA, nrow(indelmat))

  nuc_pos_in_gene =  rep(NA, nrow(indelmat))
  aa_pos_in_gene = rep(NA, nrow(indelmat))

  gene_length_in_bp = rep(NA, nrow(indelmat))

  annotation_1 = rep(NA, nrow(indelmat))
  annotation_2 = rep(NA, nrow(indelmat))

  strand = rep(NA, nrow(indelmat))

  ig_gene1 = rep(NA, nrow(indelmat))
  ig_gene2 = rep(NA, nrow(indelmat))
  intergenic = rep(NA, nrow(indelmat))


  for (i in 1:nrow(indelmat)){

    row = row.names(indelmat)[i]

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

    # STRAND INFORMATION
    strand[i] = gsub('.*=|;.*','',split_row[1])

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
                           ig_gene1, ig_gene2, intergenic,full_annots=rownames(indelmat))
  return(annotations)
}# end get_info_from_annotations

#-----------------------------------------------------------------------------
# remove_alt_allele_code_from_split_rows
#-----------------------------------------------------------------------------

#' remove_alt_allele_code_from_split_rows
#' @description Input the split matrix where rows that once had multiple annotations on
#' single line are now represented on multiple lines. For the sites that once were
#' multiallelic sites and are now represented as biallelic, thus function will
#' change the contents of indelmat_code such that the alternative allele(s) are 0.
#' For example, T -> G, C is split into two lines: T -> G and T -> C. In the code matrix,
#' turn the codes correspoding to the allele C in the row T -> G to 0 and the codes
#' corresponding to the allele G in the row T -> C to 0.
#'
#' @param indelmat_code_split - data.frame where the rows are variants (numeric description
#' variants: numbers ranging from -4 to 3), the columns are genomes,
#' and the row.names are annotations, and each line has a single annotation
#' @param indelmat_allele_split - data.frame where the rows are variants (character description
#' variants: A,C,T,G,N,-), the columns are genomes, and the row.names are annotations,
#' and each line has a single annotation
#' @param ref - character vector length nrow(indelmat_code_split) = nrow(indelmat_allele_split) indicating
#' the reference allele in terms of the positive strand
#' @param var - character vector length nrow(indelmat_code_split) = nrow(indelmat_allele_split) indicating
#' the variant allele in terms of the positive strand
#' @param rows_with_mult_var_allele_log - logical vector length nrow(indelmat_code_split) = nrow(indelmat_allele_split)
#' indicating which rows are multiallelic sites
#'
#' @return - indelmat_code - data.frame where the rows are variants (numeric description
#' variants: numbers ranging from -4 to 3), the columns are genomes,
#' and the row.names are annotations, and each line has a single annotation where
#' the alternative/minor allele in a biallelic-represrentation of a multiallelic site is now 0
#' @export
#'
#' @examples

remove_alt_allele_code_from_split_rows <- function(indelmat_code_split, indelmat_allele_split, ref, var, rows_with_mult_var_allele_log){

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

    indelmat_code_split[i,!(as.character(indelmat_allele_split[i,]) %in% c(as.character(var[i]), as.character(ref[i]), 'N', '-'))] = 0

  }

  return(indelmat_code_split)

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
#' @param mat - loaded in data.frame of indelmat or character string of a path to a indelmat
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
#' @param indelmat_code
#' @param indelmat_allele
#' @param annots
#'
#' @return
#' @export
#'
#' @examples
remove_unknown_anc = function(indelmat_code, indelmat_allele, annots){
  unknown = annots$anc %in% c('-','N')
  removed = rownames(indelmat_code)[unknown]
  filename = paste0(Sys.Date(), '_rows_removed_because_unknown_ancestral_state.txt')
  write.table(removed,file=filename,sep='\n',quote=F,row.names=F,col.names=F)
  return(list(indelmat_code=indelmat_code[!unknown,],
              indelmat_allele=indelmat_allele[!unknown,],
              annots=annots[!unknown,]))
}

#-----------------------------------------------------------------------------
# MAIN FUNCTION
#-----------------------------------------------------------------------------

#' parse_indels
#' @description Input matrices generated from internal (Ali's) variant calling pipeline.
#' Always returns parsed annotation info. In addition, you have the option to:
#' 1. split rows with multiple annotations (snps in overlapping genes, multiallelic snps)
#' 2. Re-reference to the ancestral allele at that position (instead of to the reference genome)
#' 3. Simplify the code matrix - which contains numbers from -4 to 3 indicating
#' different information about the variants - to a binary matrix indicating
#' simple presence/absence of a SNP at that site.
#' @param indelmat_code - loaded data.frame or path to the indelmat_code file generated
#' from internal variant calling pipeline
#' @param indelmat_allele - loaded data.frame or path to the indelmat_allele file generated
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

parse_indels <- function(indelmat_code, indelmat_allele, tree=NULL, og = NULL, remove_multi_annots = FALSE, return_binary_matrix = FALSE, ref_to_anc = T){

  if(is.null(tree) & return_binary_matrix){
    stop('Tree file required when returning a binary matrix.')
  }

  #-----------------------------------------------------------------------------
  # READ IN indelmat CODE AND indelmat ALLELE
  #-----------------------------------------------------------------------------

  indelmat_code <- load_if_path(indelmat_code)

  indelmat_allele <- load_if_path(indelmat_allele)

  # add semicolons to the end of the row names that don't have semicolons
  row.names(indelmat_code)[!grepl(';$', row.names(indelmat_code))] = paste0(row.names(indelmat_code)[!grepl(';$', row.names(indelmat_code))], ';')
  row.names(indelmat_allele)[!grepl(';$', row.names(indelmat_allele))] = paste0(row.names(indelmat_allele)[!grepl(';$', row.names(indelmat_allele))], ';')



  #-----------------------------------------------------------------------------
  # REMOVE BUGS
  #-----------------------------------------------------------------------------
  indelmat_code = remove_rows_with_bugs(indelmat_code)
  indelmat_allele = remove_rows_with_bugs(indelmat_allele)

  #-----------------------------------------------------------------------------
  # REMOVE LINES WITH NO VARIANTS - NO VARIANT OR ALL MASKED
  #-----------------------------------------------------------------------------
  indelmats = remove_rows_with_no_variants_or_completely_masked(indelmat_code, indelmat_allele)
  indelmat_code = indelmats[[1]]
  indelmat_allele = indelmats[[2]]

  #-----------------------------------------------------------------------------
  # EITHER (1) REMOVE ROWS WITH MULTIPLE ANNOTATIONS OR (2) SPLIT ROWS WITH
  # MULTIPLE ANNOTATIONS - DEPENDING ON VALUE OF REMOVE_MULTI_ANNOTS FLAG
  # (TRUE/FALSE)
  #-----------------------------------------------------------------------------
  if(remove_multi_annots){
    # REMOVE ROWS WITH MULTIPLE ANNOTATIONS
    indelmat_code = remove_rows_with_multiple_annots(indelmat_code)
    indelmat_allele = remove_rows_with_multiple_annots(indelmat_allele)

    #-----------------------------------------------------------------------------
    # FIND ANCESTRAL STATE OF EACH ALLELE
    #-----------------------------------------------------------------------------

    if(return_binary_matrix){
      # REROOT TREE
      tree = root_tree_og(tree)
      # GET ANCESTRAL ALLELE FOR EACH VARIANT
      if(ref_to_anc){
        alleles = get_anc_alleles(tree, indelmat_allele)
      }else{
        # REFERENCE TO MAJOR ALLELE
        alleles = get_major_alleles(data.matrix(indelmat_allele))
      }

    }


    split_rows_flag = 1:nrow(indelmat_allele)

    rows_with_multiple_annots_log = rep(FALSE, nrow(indelmat_allele))
    rows_with_mult_var_allele_log = rep(FALSE, nrow(indelmat_allele))
    rows_with_overlapping_genes_log = rep(FALSE, nrow(indelmat_allele))

    # GET ANNOTATIONS
    annots = cbind(get_info_from_annotations(indelmat_code), rows_with_multiple_annots_log,
                   rows_with_mult_var_allele_log, rows_with_overlapping_genes_log,
                   split_rows_flag)

  }else{

    #-----------------------------------------------------------------------------
    # FIND ANCESTRAL STATE OF EACH ALLELE
    #-----------------------------------------------------------------------------

    if(return_binary_matrix){
      # REROOT TREE
      tree = root_tree_og(tree)
      if(ref_to_anc){
        # GET ANCESTRAL ALLELE FOR EACH VARIANT
        alleles = get_anc_alleles(tree, indelmat_allele)
      }else{
        # REFERENCE TO MAJOR ALLELE
        alleles = get_major_alleles(indelmat_allele)
      }

    }

    # SPLIT MATRICES
    indelmat_code_split_list = split_rows_with_multiple_annots(indelmat_code)
    indelmat_allele_split_list = split_rows_with_multiple_annots(indelmat_allele)

    indelmat_code = indelmat_code_split_list[[5]]
    indelmat_allele = indelmat_allele_split_list[[5]]

    rows_with_multiple_annots_log = indelmat_code_split_list[[1]]
    rows_with_mult_var_allele_log = indelmat_code_split_list[[2]]
    rows_with_overlapping_genes_log = indelmat_code_split_list[[3]]
    split_rows_flag = indelmat_code_split_list[[4]]

    if(return_binary_matrix){
      alleles = alleles[split_rows_flag,]
    }

    # GET ANNOTATIONS
    annots = cbind(get_info_from_annotations(indelmat_code), rows_with_multiple_annots_log,
                   rows_with_mult_var_allele_log, rows_with_overlapping_genes_log,
                   split_rows_flag)

    # CHANGE indelmat CODE TO REFLECT BIALLELIC REPRESENTATION OF A MULTIALLELIC SITE
    indelmat_code = remove_alt_allele_code_from_split_rows(indelmat_code,
                                                         indelmat_allele,
                                                         annots$ref,
                                                         annots$var,
                                                         rows_with_mult_var_allele_log)

  }

  if(return_binary_matrix){
    if(ref_to_anc){
      # ADD ANCESTRAL ALLELE INFO TO ANNOTATIONS
      annots$anc = alleles[,1]
      annots$anc_prob = alleles[,2]
    }else{
      annots$maj = alleles
    }

    # remove sites with unknown ancestor
    if(ref_to_anc){
      indelmats = remove_unknown_anc(indelmat_code, indelmat_allele, annots)
      indelmat_code = indelmats$indelmat_code
      indelmat_allele = indelmats$indelmat_allele
      annots = indelmats$annots
    }

    # MAKE BINARY MATRIX
    indelmat_bin = indelmat_code
    to_keep = !(rowSums(indelmat_bin ==  2) > 0 |
                  rowSums(indelmat_bin == -2) > 0 |
                  rowSums(indelmat_bin == -3) > 0 |
                  rowSums(indelmat_bin == -4) > 0)
    indelmat_bin = indelmat_bin[to_keep,]
    indelmat_bin[indelmat_bin == 3] = 1
    indelmat_bin[indelmat_bin == -1] = NA

    annots_bin = annots[to_keep,]

    if(ref_to_anc){
      indelmat_bin_reref = data.frame(t(sapply(1:nrow(indelmat_bin), function(x){
        if(annots_bin$ref[x] == annots_bin$anc[x]){
          unlist(indelmat_bin[x,])
        }else if(!annots_bin$rows_with_mult_var_allele_log[x]){
          unlist(as.numeric(!indelmat_bin[x,]))
        }else if(annots_bin$var[x] == annots_bin$anc[x]){
          unlist(indelmat_bin[x,])
        }else{
          unlist(rep(NA,ncol(indelmat_bin)))
        }
      })))

      reref = sapply(1:nrow(indelmat_bin), function(x){
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
      #indelmat_bin = indelmat_allele[to_keep,]
      names_indelmat_bin = names(indelmat_bin)
      indelmat_bin_reref = data.frame(t(sapply(1:nrow(indelmat_bin), function(x){
        if(sum(indelmat_bin[x,]==1,na.rm=T) > sum(indelmat_bin[x,]==0,na.rm=T)){
          return(as.numeric(indelmat_bin[x,]==0))
        }
        return(as.numeric(indelmat_bin[x,]))
      })))
      names(indelmat_bin_reref) = names_indelmat_bin
      reref = rep(NA,nrow(indelmat_bin))
    }

    return(list(code=list(mat=indelmat_code,annots=annots),
                allele=list(mat=indelmat_allele,annots=annots),
                bin=list(mat=indelmat_bin_reref,annots=cbind(annots_bin,reref=reref))))
  }


  return(list(code=list(mat=indelmat_code, annots=annots),
              allele=list(mat=indelmat_allele, annots=annots)))
}# end parse_snps
