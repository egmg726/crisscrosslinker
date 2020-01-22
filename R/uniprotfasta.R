#'Fetch FASTA sequence from Uniprot
#'
#'This function fetches the FASTA sequence directly from the UniProt .fasta page
#'
#'@param uniprot.id Uniprot ID of sequence wanted in the form "ID-isoform". Leave out isoform if using the canonical sequence.
#'@param download.fasta If TRUE, will download the FASTA file to your chosen directory. Defaults to FALSE.
#'@param console.message If TRUE, will display the name of the FASTA file when being loaded into R in the console. Defaults to TRUE.
#'@param fasta.directory Directory chosen for your downloaded FASTA file. Defaults to current directory.
#'@author Emma Gail
#'
#'@export
uniprot.fasta <- function(uniprot.id, download.fasta = FALSE, console.message = TRUE,
                          fasta.directory=getwd()){

  # uniprot_fasta <- strsplit(getURL(paste0('https://www.uniprot.org/uniprot/',uniprot.id,'.fasta')),'\n')[[1]]
  # 
  # if(TRUE %in% grepl('</html>',uniprot_fasta)){
  #   #if TRUE --> there is some kind of error in here because it's just a generic error page
  #   warning('Error page detected --> NA substituted')
  #   uniprot_fasta_seq <- NA
  # } else {
  #   if(console.message == TRUE){
  #     cat(paste('Getting sequence from Uniprot:',uniprot_fasta[1],'\n'))
  #   }
  #   uniprot_fasta_seq <- paste0(uniprot_fasta[2:length(uniprot_fasta)],collapse = '')
  # }
  
  uniprot_fasta <- tryCatch(readLines(curl(paste0('https://www.uniprot.org/uniprot/',uniprot.id,'.fasta'))),
                            error = function(err){
                              cat('Error page detected --> NA substituted\n')
                              return(NA)
                            })

  if(!is.na(uniprot_fasta)){
    uniprot_fasta_seq <- paste0(uniprot_fasta[2:length(uniprot_fasta)],collapse = '')
  }
  
  #cat(paste('Getting sequence from Uniprot:',uniprot_fasta[1],'\n'))
  # if(length(uniprot_fasta) == 0){
  #   #if TRUE --> blank
  #   #go through loops again?
  #   #can put within a while loop
  # 
  #   warning(paste('Protein name not found on Uniprot: ',uniprot.id))
  #   uniprot_fasta_seq <- NA
  # } else {#end if(length(uniprot_fasta) == 0)
  #   uniprot_fasta_seq <- paste0(uniprot_fasta[2:length(uniprot_fasta)],collapse = '')
  # } #end else to #end if(length(uniprot_fasta) == 0)

  #can have a while loop here so that if  uniprot_fasta_seq is NA
  #can try to go through the menu loops
  # while(is.na(uniprot_fasta_seq)){
  #
  #   #can have option for re-running the sequence
  #
  #
  # }


  if(download.fasta == TRUE){
    write(uniprot_fasta,paste0(fasta.directory,'/',uniprot.id,'.fasta'))
  } #end if(download.fasta == TRUE){


  return(uniprot_fasta_seq)

} #end function uniprot.fasta
