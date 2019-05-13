#'UniProt PDBmap
#'
#'This function is able to connect with the RCSB API to be able to match a residue number of either a PDB chain or UniProt canonical sequence to its match.
#'
#'@param pdb_id Valid ID of a PDB accession number within RCSB
#'@param resno Residue number(s) to match against the PDB/UniProt number.
#'@param chain Chain ID within the PDB to match. Needed if matching PDB residue number to UniProt ID. Defaults to NULL
#'@param uniprot_id UniProt ID of protein. Needed if matching UniProt sequence to residue within PDB
#'@param output Choose between 'pdb' or 'uniprot' for which ID the output residue number should correspond to
#'
#'@author Emma Gail
#'
#'@export

uniprot.PDBmap <- function(pdb_id,resno,chain=NULL,uniprot_id=NULL,output=c('pdb','uniprot')){

  uniprot_match <- c()

  use_pdb <- FALSE
  use_uniprot <- FALSE
  pdb_ids <- pdb_id
  resnos <- resno

  if(!is.null(chain)){
    use_pdb <- TRUE
    pdbmap_id <- paste0(pdb_id,'.',chain)
  } else if(!is.null(uniprot_id)){
    use_uniprot <- TRUE
    pdbmap_id <- uniprot_id
  } else {
    warning('Chain or UniProt ID must be selected')
    return(NULL)
  }

  if((output != 'pdb') && (output != 'uniprot')){
    warning('PDB or UniProt must be selected for output')
    return(NULL)
  } #end if((output != 'pdb') || (output != 'uniprot')){

  # if((use_uniprot == TRUE) && (output == 'pdb')){
  #   uniprot_match <- data.frame()
  # }

  if(length(pdb_ids) == 1){
    pdb_uniprot_mapping <- GET(paste0('https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=',pdb_ids))
    xml_data <- xmlParse(pdb_uniprot_mapping)
    xml_data <- xmlToList(xml_data)
  }

  for(resno_index in 1:length(resnos)){
    resno <- resnos[resno_index]
    match_made <- FALSE

    if(length(pdb_ids) > 1){
      #if more than 1 use the index of the current
      pdb_id <- pdb_ids[resno_index]
      pdb_uniprot_mapping <- GET(paste0('https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=',pdb_id))
      xml_data <- xmlParse(pdb_uniprot_mapping)
      xml_data <- xmlToList(xml_data)
    } #end if(length(pdb_ids) > 1){

    alignment_indices <- (1:length(names(xml_data)))[names(xml_data) == 'alignment']

    for(alignment_index in alignment_indices){
      block_data <- xml_data[[alignment_index]][names(xml_data[[alignment_index]]) == 'block']
      for(block_index in 1:length(block_data)){
        segment_data <- block_data[[block_index]][names(block_data[[block_index]]) == 'segment']
        segment_df <- t(data.frame(segment_data))
        if(nrow(segment_df) == 2){
          if(pdbmap_id %in% segment_df[,'intObjectId']){

            #cat(pdbmap_id)
            if(use_pdb == TRUE){
              uniprot_line <- segment_df[segment_df[,"intObjectId"] != pdbmap_id,] #the uniprot line
              pdb_line <- segment_df[segment_df[,"intObjectId"] == pdbmap_id,] #the PDB line
            } else if(use_uniprot == TRUE){
              uniprot_line <- segment_df[segment_df[,"intObjectId"] == pdbmap_id,] #the uniprot line
              pdb_line <- segment_df[segment_df[,"intObjectId"] != pdbmap_id,] #the PDB line
            }
            # uniprot_line <- segment_df[segment_df[,"intObjectId"] != pdbmap_id,] #the uniprot line
            # pdb_line <- segment_df[segment_df[,"intObjectId"] == pdbmap_id,] #the PDB line --> should check to make sure the PDB ID is in there

            #grepl(pdb_info$pdb_id,as.character(pdb_line[names(pdb_line) == "intObjectId"]))

            uniprot_id <- as.character(uniprot_line[names(uniprot_line) == "intObjectId"])
            pdb_line_id <- as.character(pdb_line[names(pdb_line) == "intObjectId"])
            pdb_line_chain <- strsplit(pdb_line_id,'\\.')[[1]][2]
            #can double check here to make sure the chain matches the chain originally given

            pdb_line_start <- as.numeric(pdb_line[names(pdb_line) == "start"])
            pdb_line_end <- as.numeric(pdb_line[names(pdb_line) == "end"])

            uniprot_line_start <- as.numeric(uniprot_line[names(pdb_line) == "start"])
            uniprot_line_end <- as.numeric(uniprot_line[names(pdb_line) == "end"])


            #get the number from the range


            #should see if the number actually falls between the two (since some of these are split up)
            #should it be resno? Or doesn't matter if within a function?

            if(output == 'pdb'){

              #cat(uniprot_line_start:uniprot_line_end)
              if(!(resno %in% uniprot_line_start:uniprot_line_end)){
                #uniprot_match <- c(uniprot_match,NA)
                next #move on since not within the two
              } #end

              #if it is within the two
              resno_index <- match(resno,uniprot_line_start:uniprot_line_end)
              uniprot_match <- c(uniprot_match,(pdb_line_start:pdb_line_end)[resno_index])

              if(use_uniprot == TRUE){
                #rbind(uniprot_match,data.frame(resno=)
                uniprot_match <- c(uniprot_match,pdb_line_chain)
                #uniprot_match <- c(uniprot_match,paste0((pdb_line_start:pdb_line_end)[resno_index],'_',pdb_line_chain))
              }

              match_made <- TRUE

            } else if(output == 'uniprot'){

              if(!(resno %in% pdb_line_start:pdb_line_end)){
                #uniprot_match <- c(uniprot_match,NA)
                next #move on since not within the two
              } #end

              #if it is within the two
              resno_index <- match(resno,pdb_line_start:pdb_line_end)
              uniprot_match <- c(uniprot_match,(uniprot_line_start:uniprot_line_end)[resno_index])
              match_made <- TRUE
            }



          } #end if(uniprot_id %in% segment_df[,'intObjectId'])

        } else {
          warning('More than 2 rows found in block.segment in pdb_uniprot_mapping --> check data\n')
        } #end else to if(nrow(segment_df) == 2)
      } #end for(block_index in 1:length(block_data))
    } #end for(alignment_index in alignment_indices)

    if(match_made == FALSE){
      uniprot_match <- c(uniprot_match,NA)
    }

  } #end for(resno_index in 1:length(resnos)){

  #within the for loop?

  return(uniprot_match)

} #end function uniprot.PDBmap
