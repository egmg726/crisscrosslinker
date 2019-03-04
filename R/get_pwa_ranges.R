#' Get ranges from pairwise alignment results
#'
#' This function allows you to return ranges (start and end) from pairwise alignment results
#' @param pwa_results Result from function pairwiseAlignment()
#' @keywords pairwise alignment
#' @export
#' @examples
#' pwa_ranges <- get_pwa_ranges(pwa_results)


get_pwa_ranges <- function(pwa_results){

  start_pattern <- start(pattern(pwa_results))
  end_pattern <- end(pattern(pwa_results))

  start_subject <- start(subject(pwa_results))
  end_subject <- end(subject(pwa_results))

  start_difference <- start_subject - start_pattern

  pwa_ranges <- list(start_pattern = start_pattern,
                     end_pattern = end_pattern,
                     start_subject = start_subject,
                     end_subject = end_subject,
                     start_difference = start_difference)

  return(pwa_ranges)

}
