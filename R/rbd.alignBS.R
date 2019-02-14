#'RBD: Align Binding Sequences to Uniprot and PDB
#'
#'This function aligns the
#'
#'@param rbd.df rbd.df
#'@param alignIDs alignIDs
#'@param alignTo alignTo
#'@param uniprot2pdb Boolean, defaults to TRUE. This will align to Uniprot sequence before PDB and using the uniprot.PDBmap function to improve alignment accuracy. Highly recommended for PDB alignments.
#'@param allowpartialBS Boolean, defaults to FALSE.
#'@export


rbd.alignBS <- function(rbd.df,alignIDs,
                        alignTo=c('pdb','uniprot','fasta'),
                        uniprot2pdb=TRUE,allowPartialBS=FALSE){
  #make it so that if uniprot2pdb is true, will match to the uniprot then to the pdb
  #for improved accuracy

  #protID <- as.character(unique(rbd.df2$ProtID))
  binding_site_start_list <- c()
  binding_site_end_list <- c()
  binding_sequence_list <- c()
  db_list <- c()
  db_id_list <- c()
  source_sequence_list <- c()

  #do the for loop up here

  for(row_num in 1:nrow(rbd.df)){

    rbddf_row <- rbd.df[row_num,]
    fragmentStart <- as.character(rbddf_row$fragmentStart)
    fragmentStop <- as.character(rbddf_row$fragmentStop)
    proteolyticFragment <- as.character(rbddf_row$proteolyticFragment)
    protID <- as.character(rbddf_row$ProtID)

    binding_sequence_list <- c(binding_sequence_list,proteolyticFragment)

    if(alignTo == 'pdb'){
      alignID <- as.character(alignIDs[alignIDs$protID==protID,'pdbID'])

      if((uniprot2pdb == TRUE)){

        #if TRUE --> will align to uniprotseq first
        alignID_split <- strsplit(alignID,'_')[[1]]
        pdb_id <- alignID_split[1]
        chain <- alignID_split[2]
        uniprot_id <- alignID <- as.character(alignIDs[alignIDs$protID==protID,'uniprotID'])
        if(!(uniprot_id == protID)){
          #if TRUE --> they do not match so can't just use fragmentStart and fragmentStop
          #need to do an alignment with the uniprot sequence
          #check
          #do str_locate_all
          if(!file.exists(paste0(uniprot_id,'.fasta'))){
            uniprotSeq <- uniprot.fasta(uniprot_id = uniprot_id)
          } else {
            uniprotSeq <- paste0(toupper(read.fasta(paste0(uniprot_id,'.fasta'))[[1]]),collapse='')
          }
          start_end <- t(str_locate_all(uniprotSeq,proteolyticFragment)[[1]])[,1]
          fragmentStart <- start_end[1]
          fragmentStop <- start_end[2]
        }
        #fragmentStart and fragmentStop need to be matched to uniprot first
        #can check if protID == uniprotID
        #or if protID is a valid Uniprot ID by checking with the uniprot.fasta() function

        pdb_positions <- c(uniprot.PDBmap(pdb_id,chain=chain,resno=fragmentStart,output='pdb'),
                           uniprot.PDBmap(pdb_id,chain=chain,resno=fragmentStop,output='pdb'))

        #will need if(NA %in% pdb_positions) too
        #will exclude those matches if allowPartialBS == FALSE
        pdb_read <- check_download_read_pdb(pdb_id) #add ID to list of IDs that need to be loaded?
        resno_and_resid <- quick_resno_and_resid(pdb_read = pdb_read, chain = chain)
        pdb_seq <- paste0(resno_and_resid$resid,collapse='')

        #add the

      } #end if((uniprot2pdb == TRUE)){

      if(NA %in% pdb_positions){

        db <- "UniProt"
        db_list <- c(db_list,db)
        db_id_list <- c(db_id_list,protID)
        binding_site_start_list <- c(binding_site_start_list,fragmentStart)
        binding_site_end_list <- c(binding_site_end_list,fragmentStop)
        source_sequence_list <- c(source_sequence_list,uniprot.fasta(protID))

      } else { #end if(NA %in% pdb_positions){
        #there is a PDB position for both of the segments
        #can then save those positions into a table

        db <- "PDB"
        db_list <- c(db_list,db)
        db_id_list <- c(db_id_list,paste0(pdb_id,'_',chain))
        binding_site_start_list <- c(binding_site_start_list,pdb_positions[1])
        binding_site_end_list <- c(binding_site_end_list,pdb_positions[2])
        #pdb_seq <- pdb.annotate(pdb_id)[paste0(pdb_id,'_',chain),]$sequence
        pdb_read <- read.pdb2(pdb_id) #add ID to list of IDs that need to be loaded?
        resno_and_resid <- quick_resno_and_resid(pdb_read = pdb_read, chain = chain)
        pdb_seq <- paste0(resno_and_resid$resid,collapse='')
        source_sequence_list <- c(source_sequence_list,pdb_seq)

      } #end else to if(NA %in% pdb_positions){
    } else if(alignTo == 'uniprot'){ #end (alignTo == 'pdb')
      alignID <- as.character(alignIDs[alignIDs$protID==protID,'uniprotID'])
      #will need to align to the uniprot file
      #check if there is a fasta file matching the uniprot ID
      #if not --> get the uniprot

    }



  } #end for(row_num in 1:nrow(rbd.df)){





  bs_output <- data.frame(list(binding_site_start = binding_site_start_list,
                               binding_site_end = binding_site_end_list,
                               binding_sequence = binding_sequence_list,
                               db = db_list,
                               db_id = db_id_list,
                               source_sequence = source_sequence_list))

  return(bs_output)
} #end function rbd.alignBS()
