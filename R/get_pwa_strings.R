#' Get string from pairwise alignment results
#'
#' This function allows you to return strings from pairwise alignment results
#' @param pwa_results Result from function pairwiseAlignment()
#' @param result_as_vector If TRUE, returns string as a vector
#' @keywords pairwise alignment
#' @export
#' @examples
#' pwa_strings <- get_pwa_strings(pwa_results)


get_pwa_strings <- function(pwa_results,result_as_vector = FALSE){
  subject_string <- toString(subject(pwa_results))
  pattern_string <- toString(pattern(pwa_results))

  if(result_as_vector == FALSE){
    subject_string <- strsplit(subject_string,split='')[[1]]
    pattern_string <- strsplit(pattern_string,split='')[[1]]
  }

  pwa_strings <- list(subject_string=subject_string,
                      pattern_string=pattern_string)

  return(pwa_strings)
}
