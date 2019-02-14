#'Fetch FASTA sequence from Uniprot
#'
#'This function fetches the FASTA sequence from UniProt
#'
#'@param uniprot_id uniprot_id
#'@param download_fasta download_fasta
#'@param console_message console_message
#'@param fasta_directory fasta_directory
#'
#'@export
uniprot.fasta <- function(uniprot_id, download_fasta = FALSE, console_message = TRUE,
                          fasta_directory=getwd()){

  uniprot_fasta <- strsplit(getURL(paste0('https://www.uniprot.org/uniprot/',uniprot_id,'.fasta')),'\n')[[1]]

  if(TRUE %in% grepl('</html>',uniprot_fasta)){
    #if TRUE --> there is some kind of error in here because it's just a generic error page
    warning('Error page detected --> NA substituted')
    uniprot_fasta_seq <- NA
  } else {
    if(console_message == TRUE){
      cat(paste('Getting sequence from Uniprot:',uniprot_fasta[1],'\n'))
    }
    uniprot_fasta_seq <- paste0(uniprot_fasta[2:length(uniprot_fasta)],collapse = '')
  }

  #cat(paste('Getting sequence from Uniprot:',uniprot_fasta[1],'\n'))
  if(length(uniprot_fasta) == 0){
    #if TRUE --> blank
    #go through loops again?
    #can put within a while loop

    warning(paste('Protein name not found on Uniprot: ',uniprot_id))
    uniprot_fasta_seq <- NA
  } else {#end if(length(uniprot_fasta) == 0)
    uniprot_fasta_seq <- paste0(uniprot_fasta[2:length(uniprot_fasta)],collapse = '')
  } #end else to #end if(length(uniprot_fasta) == 0)

  #can have a while loop here so that if  uniprot_fasta_seq is NA
  #can try to go through the menu loops
  # while(is.na(uniprot_fasta_seq)){
  #
  #   #can have option for re-running the sequence
  #
  #
  # }


  if(download_fasta == TRUE){
    write(uniprot_fasta,paste0(fasta_directory,'/',uniprot_id,'.fasta'))
  } #end if(download_fasta == TRUE){


  return(uniprot_fasta_seq)

} #end function uniprot.fasta
