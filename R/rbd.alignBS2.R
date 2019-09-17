#'RBD: Align Binding Sequences to Uniprot and PDB
#'
#'This function aligns the binding sequences 
#'
#'@param bs.output Output from rbd.getBSfromIET()
#'@param alignIDs alignIDs containing columns of protein names, UniProt and/or PDB identifiers
#'@param alignTo What the sequences should be aligned to. Choose between 'pdb' and 'uniprot'
#'@param uniprot2pdb Boolean, defaults to TRUE. This will align to Uniprot sequence before PDB and using the uniprot.PDBmap function to improve alignment accuracy. Highly recommended for PDB alignments.
#'@param allowpartialBS Boolean, defaults to FALSE. will allow partial matches of binding sequences to PDB/UniProt sequences
#'
#'@author Emma Gail
#'
#'@export
rbd.alignBS2 <- function(bs.output,alignIDs,
                        alignTo=c('pdb','uniprot','fasta'),
                        uniprot2pdb=TRUE,allowPartialBS=FALSE){

  
  #filtered_input_eluate_table
  
  #need to store the Uniprot sequence for future use so it doesn't get it multiple times
  
  
  stored_uniprots <- list()
  
  #protID <- as.character(unique(rbd.df2$ProtID))
  binding_site_start_list <- c()
  binding_site_end_list <- c()
  binding_sequence_list <- c()
  db_list <- c()
  db_id_list <- c()
  source_sequence_list <- c()
  
  #do the for loop up here
  
  pdb_reads <- list()
  
  for(row_num in 1:nrow(bs.output)){
    
    rbddf_row <- bs.output[row_num,]
    fragmentStart <- as.character(rbddf_row$binding_site_start)
    fragmentStop <- as.character(rbddf_row$binding_site_end)
    proteolyticFragment <- as.character(rbddf_row$binding_sequence)
    protID <- as.character(rbddf_row$pro_name)
    
    binding_sequence_list <- c(binding_sequence_list,proteolyticFragment)
    
    if(alignTo == 'pdb'){
      alignID <- as.character(alignIDs[alignIDs$protID==protID,'pdbID'])
      
      #there is no false option for this --> need to fix
      if((uniprot2pdb == TRUE)){
        
        
        
        #if TRUE --> will align to uniprotseq first
        alignID_split <- strsplit(alignID,'_')[[1]]
        pdb_id <- alignID_split[1]
        chain <- alignID_split[2]
        # if(is.na(chain)){
        #   chain <- NULL
        # }
        uniprot_id <- as.character(alignIDs[alignIDs$protID==protID,'uniprotID'])
        uniprot_id <- strsplit(uniprot_id,'-')[[1]][1]
        if(!(uniprot_id == protID)){
          #if TRUE --> they do not match so can't just use fragmentStart and fragmentStop
          #need to do an alignment with the uniprot sequence
          #check
          #do str_locate_all
          if(!file.exists(paste0(uniprot_id,'.fasta'))){
            
            if(!(uniprot_id %in% names(stored_uniprots))){
              uniprotSeq <- uniprot.fasta(uniprot_id)
              stored_uniprots[[uniprot_id]] <- uniprotSeq
            } else {
              uniprotSeq <- stored_uniprots[[uniprot_id]]
            }
            
          } else { #end if(!file.exists(paste0(uniprot_id,'.fasta'))){
            uniprotSeq <- paste0(toupper(read.fasta(paste0(uniprot_id,'.fasta'))[[1]]),collapse='')
          } #end else to if(!file.exists(paste0(uniprot_id,'.fasta'))){
          
          if(grepl(proteolyticFragment,uniprotSeq)){
            start_end <- t(str_locate_all(uniprotSeq,proteolyticFragment)[[1]])[,1]
            fragmentStart <- start_end[1]
            fragmentStop <- start_end[2]
          } else {
            fragmentStart <- NA
            fragmentStop <- NA
          }
          

        } #end if(!(uniprot_id == protID)){
        
        #fragmentStart and fragmentStop need to be matched to uniprot first
        #can check if protID == uniprotID
        #or if protID is a valid Uniprot ID by checking with the uniprot.fasta() function
        
        
        
        # pdb_positions <- c(uniprot.PDBmap(pdb_id,chain=chain,resno=fragmentStart,output='pdb'),
        #                    uniprot.PDBmap(pdb_id,chain=chain,resno=fragmentStop,output='pdb'))
        
        
        if(!is.na((fragmentStart))){
          if(!is.na(chain)){
            pdb_positions <- uniprot.PDBmap(pdb_id,chain=chain,resno=fragmentStart:fragmentStop,output='pdb')
          } else {
            pdb_positions <- uniprot.PDBmap(pdb_id,uniprot_id = uniprot_id,resno=fragmentStart:fragmentStop,output='pdb')
            #pdb_positions <- na.omit(pdb_positions)
            #pdb_positions[rep(c(TRUE,FALSE),length(pdb_positions)/2)]
            unique_chains <- unique(pdb_positions[rep(c(FALSE,TRUE),length(pdb_positions)/2)])
            if(length(unique_chains) == 1){
              chain <- unique_chains
            } else {
              warning('More than 1 chain detected')
            }
            # is.sequential(as.numeric(pdb_positions[rep(c(TRUE,FALSE),length(pdb_positions)/2)]))
            pdb_positions <- (pdb_positions[rep(c(TRUE,FALSE),length(pdb_positions)/2)])
          }
        } else {
          pdb_positions <- NA
        }
       
        
        
        
        
        if(NA %in% pdb_positions){
          #if there is NA --> not a perfect match
          warning('NA within uniprot.PDBmap range')
          pdb_positions <- c(pdb_positions[1],pdb_positions[length(pdb_positions)])
        } else {
          pdb_positions <- c(pdb_positions[1],pdb_positions[length(pdb_positions)])
        }
        
        #will need if(NA %in% pdb_positions) too
        #will exclude those matches if allowPartialBS == FALSE
        if(pdb_id %in% names(pdb_reads)){
          pdb_read <- pdb_reads[[pdb_id]]
        } else {
          pdb_read <- read.pdb2(pdb_id, verbose = FALSE)
          pdb_reads[[pdb_id]] <- pdb_read
        }
        #add ID to list of IDs that need to be loaded?
        
        if(!is.na(chain)){
          resno_and_resid <- quick_resno_and_resid(pdb_read = pdb_read, chain = chain)
          pdb_seq <- paste0(resno_and_resid$resid,collapse='')
        }

        
        #add the
        
      } #end if((uniprot2pdb == TRUE)){
      
      if(NA %in% pdb_positions){
        
        db <- "UniProt"
        db_list <- c(db_list,db)
        db_id_list <- c(db_id_list,uniprot_id)
        binding_site_start_list <- c(binding_site_start_list,fragmentStart)
        binding_site_end_list <- c(binding_site_end_list,fragmentStop)
        source_sequence_list <- c(source_sequence_list,uniprot.fasta(uniprot_id))
        
      } else { #end if(NA %in% pdb_positions){
        #there is a PDB position for both of the segments
        #can then save those positions into a table
        
        
        db <- "PDB"
        db_list <- c(db_list,db)
        db_id_list <- c(db_id_list,paste0(pdb_id,'_',chain))
        binding_site_start_list <- c(binding_site_start_list,pdb_positions[1])
        binding_site_end_list <- c(binding_site_end_list,pdb_positions[2])
        #pdb_seq <- pdb.annotate(pdb_id)[paste0(pdb_id,'_',chain),]$sequence
        #pdb_read <- read.pdb2(pdb_id) #add ID to list of IDs that need to be loaded?
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
} #end function rbd.alignBS2()
