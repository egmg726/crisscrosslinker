#'RBD: Get Binding Sequence - Casetello et. al 2016 Method
#'
#'This function gets the binding sequence using a similar method to Castello et. al 2016
#'
#'@param trypticPeptide trypticPeptide, peptide detected by RBDmap/MS method
#'@param enzyme enzyme used: LysC or ArgC
#'@param sourceSequence sourceSequence, used as the reference to locate the binding sequence
#'@param cleave_offset Cleave offset for the sequence. Defaults to 4
#'
#'@export

rbd.getBindingSeq3 <- function(trypticPeptide,enzyme,sourceSequence,
                               cleave_offset = 4){
  #include_ambiguous = FALSE, proteolytic_fragments = FALSE){

  if(enzyme == 'LysC'){
    cleave_aa <- 'K'
  } else if(enzyme == 'ArgC'){
    cleave_aa <- 'R'
  } else {
    warning('Unknown protease detected --> returning NULL')
    return(NULL)
  }

  # if(is.na(protID)){
  #   return(NULL)
  # }

  uniprotSeq <- sourceSequence

  if(is.na(uniprotSeq)){
    return(NULL)
  }

  list_of_cleaves <- str_locate_all(uniprotSeq,cleave_aa)[[1]][,1]

  #CleavageSites(ProtSeq = uniprotSeq, cleavagePattern = ".{1}.{1}.{1}[K]{1}.{1}.{1}")

  cleavagePattern <- paste0(".{1}.{1}.{1}[",cleave_aa,"]{1}.{1}.{1}")

  cleavage_sites <- gregexpr(pattern=cleavagePattern,text=as.character(uniprotSeq))[[1]]
  #cleavage_sites
  start_and_end_patterns <- str_locate_all(uniprotSeq,trypticPeptide)[[1]][1,]
  start_pattern <- as.numeric(start_and_end_patterns[1])
  end_pattern <- as.numeric(start_and_end_patterns[2])

  uniprotSeq_split <- strsplit(uniprotSeq,split='')[[1]]
  new_start_pattern <- tail(cleavage_sites[cleavage_sites < start_pattern],n=1)+cleave_offset
  if(length(new_start_pattern) == 0){
    new_start_pattern <- 1
  }
  if(new_start_pattern > start_pattern){
    new_start_pattern <- tail(cleavage_sites[cleavage_sites < start_pattern],n=2)[1]+cleave_offset
  }
  #new_end_pattern <- head(cleavage_sites[cleavage_sites > end_pattern],n=1)+(cleave_offset-1)
  new_end_pattern <- head(list_of_cleaves[list_of_cleaves >= end_pattern],n=1)
  # if(uniprotSeq_split[end_pattern] == cleave_aa){
  #   new_end_pattern <- end_pattern
  # } else {
  #   new_end_pattern <- head(cleavage_sites[cleavage_sites >= end_pattern],n=1)+cleave_offset-1
  # }
  if(length(new_end_pattern) == 0){
    new_end_pattern <- length(uniprotSeq_split)
  }

  #(cleavage_sites[cleavage_sites < start_pattern],n=1)

  #

  bindSeq <- paste0(uniprotSeq_split[new_start_pattern:new_end_pattern],collapse='')

  return(bindSeq)

} #end function rbd.getBindingSeq2
