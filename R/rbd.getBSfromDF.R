#'RBD: Get binding sequences from Dataframe
#'
#'Retrieves the binding sequences based on the information given in the dataframe
#'and includes them as part of the output.
#'
#'@param rbd.df Data.frame containing the columns protID, trypticPeptide, and enzyme.
#'@param fasta Name of fasta file or loaded fasta file by seqinr::read.fasta(). Defaults to NULL.
#'@param cleave_offset Number of amino acids (AAs) that have to be between the AAs cut by the enzyme. Defaults to 4 as defined by the original RBDmap experiment.
#'@return The original rbd.df with 3 columns added: the proteolytic fragment, the fragment start, and fragment end.
#'@author Emma Gail
#'@export
rbd.getBSfromDF <- function(rbd.df, fasta = NULL, cleave_offset = 4){

  #first check to make sure it has the 3 columns
  if(FALSE %in% (c('ProtID','enzyme','trypticPeptide') %in% colnames(rbd.df))){
    warning('rbd.df must contain the following columns: ProtID, enzyme, and trypticPeptide')
    return(NULL)
  } #end if(FALSE %in% (c('ProtID','enzyme','trypticPeptide') %in% colnames(rbd.df))){

  # binding_site_start_list <- c()
  # binding_site_end_list <- c()
  # binding_sequence_list <- c()
  # db_list <- c()
  # db_id_list <- c()
  # source_sequence_list <- c()

  proteolyticFragment_list <- c()
  fragmentStart_list <- c()
  fragmentStop_list <- c()

  for(row_num in 1:nrow(rbd.df)){
    rbd_row <- rbd.df[row_num,]
    trypticPeptide <- as.character(rbd_row$trypticPeptide)
    enzyme <- as.character(rbd_row$enzyme)
    protID <- as.character(rbd_row$ProtID)



    if(!is.null(fasta)){

      #check if file or list
      if(typeof(fasta) == 'character'){
        fasta <- seqinr::read.fasta(fasta)
      } #end if(typeof(fasta) == 'character'){

      if(protID %in% names(fasta)){
        sourceSequence <- toupper(paste0(fasta[[protID]],collapse = ''))
      } else {
        sourceSequence <- toupper(paste0(fasta[[1]],collapse = ''))
      }


    } else {#end if(!is.null(fasta)){

      #tryCatch here?
      sourceSequence <- uniprot.fasta(uniprot.id = protID)

    }#end if(!is.null(fasta)){

    #uniprotSeq <- uniprot.fasta(uniprot_id = protID)

    bindSeq <- rbd.getBindingSeq3(trypticPeptide = trypticPeptide,
                                  enzyme = enzyme, sourceSequence = sourceSequence,
                                  cleave_offset = cleave_offset)

    start_end <- str_locate_all(sourceSequence,bindSeq)[[1]]
    startp <- start_end[1]
    endp <- start_end[2]

    proteolyticFragment_list <- c(proteolyticFragment_list,bindSeq)
    fragmentStart_list <- c(fragmentStart_list,startp)
    fragmentStop_list <- c(fragmentStop_list,endp)

  } #end for(row_num in 1:nrow(rbd.df)){

  #have basically the same table but just add the following columns:
  #proteolyticFragment
  #fragmentStart
  #fragmentEnd

  #the output has to be the same/similar to bs_output
  #return(bs_output)

  rbd.df <- cbind(rbd.df,data.frame(list(proteolyticFragment=proteolyticFragment_list,
                                         fragmentStart=fragmentStart_list,
                                         fragmentStop=fragmentStop_list)))

  return(rbd.df)
} #end function rbd.getBSfromDF
