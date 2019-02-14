#should be able to use multiple resnos
uniprot.PDBmap <- function(pdb_id,chain,resno,output=c('pdb','uniprot')){

  uniprot_match <- NA
  pdbmap_id <- paste0(pdb_id,'.',chain)

  pdb_uniprot_mapping <- GET(paste0('https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=',pdb_id))
  xml_data <- xmlParse(pdb_uniprot_mapping)
  xml_data <- xmlToList(xml_data)

  alignment_indices <- (1:length(names(xml_data)))[names(xml_data) == 'alignment']

  for(alignment_index in alignment_indices){
    block_data <- xml_data[[alignment_index]][names(xml_data[[alignment_index]]) == 'block']
    for(block_index in 1:length(block_data)){
      segment_data <- block_data[[block_index]][names(block_data[[block_index]]) == 'segment']
      segment_df <- t(data.frame(segment_data))
      if(nrow(segment_df) == 2){
        if(pdbmap_id %in% segment_df[,'intObjectId']){
          uniprot_line <- segment_df[segment_df[,"intObjectId"] != pdbmap_id,] #the uniprot line
          pdb_line <- segment_df[segment_df[,"intObjectId"] == pdbmap_id,] #the PDB line --> should check to make sure the PDB ID is in there

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

            if(!(resno %in% uniprot_line_start:uniprot_line_end)){
              next #move on since not within the two
            } #end

            #if it is within the two
            resno_index <- match(resno,uniprot_line_start:uniprot_line_end)
            uniprot_match <- (pdb_line_start:pdb_line_end)[resno_index]


          } else if(output == 'uniprot'){

            if(!(resno %in% pdb_line_start:pdb_line_end)){
              next #move on since not within the two
            } #end

            #if it is within the two
            resno_index <- match(resno,pdb_line_start:pdb_line_end)
            uniprot_match <- (uniprot_line_start:uniprot_line_end)[resno_index]

          }



        } #end if(uniprot_id %in% segment_df[,'intObjectId'])

      } else {
        warning('More than 2 rows found in block.segment in pdb_uniprot_mapping --> check data\n')
      } #end else to if(nrow(segment_df) == 2)
    } #end for(block_index in 1:length(block_data))
  } #end for(alignment_index in alignment_indices)

  return(uniprot_match)

} #end function uniprot.PDBmap
