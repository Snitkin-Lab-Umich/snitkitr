#-----------------------------------------------------------------------------
# SUB FUNCTIONS
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# remove_rows_with_bugs
#-----------------------------------------------------------------------------

#' remove_rows_with_bugs
#'
#' @param snpmat 
#'
#' @return snpmat with rows removed, writes a file logging what was removed
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
  write(row.names(snpmat)[num_pipes%%9 != 0], file = filename, append = TRUE)
  
  # only keep rows with correct number of pipes 
  snpmat = snpmat[num_pipes%%9 == 0,] # must be a multiple of 9 
  
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
#'
#' @param snpmat_code 
#' @param snpmat_allele 
#'
#' @return
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
#'
#' @param snpmat 
#'
#' @return
#' list including: 
#' 1. rows_with_multiple_annots_log -> length of snpmat_added
#' 2. rows_with_mult_var_allele_log-> length of snpmat_added
#' 3. rows_with_overlapping_genes_log -> length of snpmat_added
#' 4. split_rows_flag -> length of snpmat_added
#' 5. snpmat_added
#' Also writes a file to make sure all of the multiple annottaions are accounted for 
#' by either multialleic sites or overlapping genes 
#' @export
#'
#' @examples

split_rows_with_multiple_annots <- function(snpmat){
  
  num_dividers <- sapply(1:nrow(snpmat), function(x) lengths(regmatches(row.names(snpmat)[x], gregexpr(";[A,C,G,T]", row.names(snpmat)[x]))))
  
  rows_with_multiple_annotations <- c(1:nrow(snpmat))[num_dividers >= 2 & str_count(row.names(snpmat), '\\|') > 9]
  
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
  names_of_rows = row.names(snpmat_added) %>% gsub(';\\.[0-9].*$', ';', .)
  
  split_rows_flag = as.numeric(factor(names_of_rows, levels = unique(names_of_rows))) # could be because multiple annotations, could be because different variants 
  
  dup = unique(split_rows_flag[duplicated(split_rows_flag)]) # rows that were duplicated
  
  split_annotations <- strsplit(row.names(snpmat_added)[split_rows_flag %in% dup], ";")
  
  
  
  # FIX ANNOTS OF SNP MAT ADDED - RELIES ON THE .1 flag 
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
    if(grepl('> [A,C,T,G],[A,C,T,G].*functional', r)){
      var =  gsub('^.*Strand Information:','',r) %>% gsub('\\|.*$', '', .) %>% substr(.,nchar(.),nchar(.))
      gsub('> [A,C,T,G],[A,C,T,G].*functional', paste('>', var, 'functional'), r)
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
#'
#' @param snpmat 
#'
#' @return
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
#'
#' @param snpmat 
#'
#' @return
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
#'
#' @param snpmat_code_split 
#' @param snpmat_allele_split 
#' @param ref 
#' @param var 
#' @param rows_with_mult_var_allele_log 
#'
#' @return
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

#-----------------------------------------------------------------------------
# MAIN FUNCTION
#-----------------------------------------------------------------------------

#' parse_snps
#'
#' @param path_to_snpmat_code 
#' @param path_to_snpmat_allele 
#' @param remove_multi_annots - logical flag indicating if you want to remove rows
#' with multiple annotations - alternative is to split rows with mutliple annotations
#'
#' @return
#' @export
#'
#' @examples

parse_snps <- function(path_to_snpmat_code, path_to_snpmat_allele, remove_multi_annots = FALSE){
  
  #-----------------------------------------------------------------------------
  # READ IN SNPMAT CODE AND SNPMAT ALLELE 
  #-----------------------------------------------------------------------------
  snpmat_code <-   read.table(path_to_snpmat_code,
                              header = TRUE,
                              stringsAsFactors = FALSE,
                              sep = "\t",
                              quote = "", 
                              row.names = 1)
  

  snpmat_allele <-   read.table(path_to_snpmat_allele,
                                header = TRUE,
                                stringsAsFactors = FALSE,
                                sep = "\t",
                                quote = "", 
                                row.names = 1)
  
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
    
    # GET ANNOTATIONS
    annots = get_info_from_annotations(snpmat_code)
    
    

  }else{
  # SPLIT MATRICES
  snpmat_code_split_list = split_rows_with_multiple_annots(snpmat_code)
  snpmat_allele_split_list = split_rows_with_multiple_annots(snpmat_allele)
  
  snpmat_code = snpmat_code_split_list[[5]]
  snpmat_allele = snpmat_allele_split_list[[5]]
  
  rows_with_multiple_annots_log = snpmat_code_split_list[[1]]
  rows_with_mult_var_allele_log = snpmat_code_split_list[[2]]
  rows_with_overlapping_genes_log = snpmat_code_split_list[[3]]
  split_rows_flag = snpmat_code_split_list[[4]]
  
  # GET ANNOTATIONS
  annots = cbind(get_info_from_annotations(snpmat_code), rows_with_multiple_annots_log, 
                 rows_with_mult_var_allele_log, rows_with_overlapping_genes_log, 
                 split_rows_flag)
  
  # CHANGE SNPMAT CODE TO REFLECT BIALLELIC REPRESENTATION OF A MULTIALLELIC SITE
  snpmat_code = remove_alt_allele_code_from_split_rows(snpmat_code, 
                                         snpmat_allele, 
                                         annots$ref, 
                                         annots$var, 
                                         rows_with_mult_var_allele_log)
  
  }
  
  return(list(list(snpmat_code, annots), list(snpmat_allele, annots)))
}# end parse_snps
