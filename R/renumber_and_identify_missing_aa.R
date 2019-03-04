#' Renumber and identify missing amino acids in pairwise alignment results
#'
#' This function allows you to renumber and identify missing amino acids in pairwise alignment results
#' @param pwa_results Result from function pairwiseAlignment()
#' @param start_difference Difference between starting amino acid in subject and pattern alignments
#' @keywords pairwise alignment
#' @author Emma Gail
#' @export
#' @examples
#' renumbered_missing_amino_acids <- renumber_and_identify_missing_aa(pwa_strings,start_difference)


renumber_and_identify_missing_aa <- function(pwa_strings,start_difference){
  subject_string_vector <- pwa_strings$subject_string
  pattern_string_vector <- pwa_strings$pattern_string

  if(length(subject_string_vector) == 1){
    subject_string_vector <- strsplit(subject_string_vector,split='')[[1]]
    pattern_string_vector <- strsplit(pattern_string_vector,split='')[[1]]
  }

  missing_animo_acids_pattern <- c()
  missing_animo_acids_subject <- c()
  other_mutations <- c()
  new_pdb_numbering <- c()
  maa_pattern_numbering <- c()
  for(vector_num in 1:length(subject_string_vector)){
    if(pattern_string_vector[vector_num] != subject_string_vector[vector_num]){
      if(pattern_string_vector[vector_num] == '-'){
        missing_animo_acids_pattern <- c(missing_animo_acids_pattern,vector_num)
        new_resno <- vector_num + start_difference
        maa_pattern_numbering <- c(maa_pattern_numbering,new_resno)
      }
      else if(subject_string_vector[vector_num] == '-'){
        missing_animo_acids_subject <- c(missing_animo_acids_subject,vector_num)
      } else {
        other_mutations <- c(other_mutations,vector_num)
        new_resno <- vector_num + start_difference
        new_pdb_numbering <- c(new_pdb_numbering, new_resno)
      }
      #else == pattern and subject match
    } else {
      new_resno <- vector_num + start_difference
      new_pdb_numbering <- c(new_pdb_numbering, new_resno)
    }
  }

  missing_animo_acids_list <- list(pattern=missing_animo_acids_pattern,
                                   subject=missing_animo_acids_subject,
                                   mutations=other_mutations,
                                   new_pdb_numbering=new_pdb_numbering,
                                   maa_pattern_numbering=maa_pattern_numbering)

  return(missing_animo_acids_list)

}
