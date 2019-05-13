#'Get Binding Sequences from Input/Eluate Table
#'
#'This function gets the binding sequences from the input eluate table
#'
#'@param filtered_input_eluate_table Output from rbd.makeIETable()
#'@param fasta_file FASTA file from seqinr::read.fasta()
#'@param align_to What source to align the BS to: "PDB","FASTA","UniProt"
#'@param protein_to_uniprot_id Table containing protein_id,uniprot_id,pdb_id
#'@param protein_dict protein_dict
#'@param include_ambiguous include_ambiguous
#'@export

rbd.getBSfromIET <- function(filtered_input_eluate_table,
                             fasta_file,
                             align_to = 'PDB',
                             protein_to_uniprot_id = NULL,
                             protein_dict = NULL,
                             include_ambiguous = FALSE){

  #filter the table beforehand as well for QC??
  filtered_input_eluate_table <- filtered_input_eluate_table[filtered_input_eluate_table$eluate > 0,]

  #assign protein dict to protein to uniprot id before
  #do a search and replace in the future to get rid of it later??



  #need
  bs_output_mega <- data.frame()
  protein_names_already_run <- list()



  #num <- 2
  for(num in 1:nrow(filtered_input_eluate_table)){

    protein_name <- filtered_input_eluate_table$protein[num]
    # if(protein_name == 'SUZ12_Q15022'){ #need to make changes later to make sure SUZ12 can be accounted for later
    #   next
    # }

    filtered_input_eluate_table_row <- filtered_input_eluate_table[num,]
    input_sequence <- filtered_input_eluate_table_row$sequence
    protein_name <- filtered_input_eluate_table_row$protein
    fasta_seq <- fasta_file[[protein_name]]

    #check if this has a '_' in it before splitting it?
    #need to get last
    protease_split <- strsplit(filtered_input_eluate_table_row$category,'_')[[1]]
    protease <- protease_split[length(protease_split)]

    bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                   protease = protease,
                                   sequence_for_alignment=toupper(paste0(fasta_seq,collapse='')),
                                   protein_name = protein_name,
                                   database_name = 'FASTA',
                                   database_id = protein_name,
                                   include_ambiguous = include_ambiguous)


    binding_seq <- as.character(bs_output$binding_sequence)

    #make a list of sequences that have already been run through the
    #loop so that it doesn't occur multiple times




    #can have an align_to variable


    #database_name <- align_to
    #if database_name is not changed before rbd.

    if(align_to == 'PDB'){

      #see if the protein has already been run (same as before)
      #add in protein_dict functionality

      if(!is.null(protein_dict)){
        #protein_name
        pdb_id <- as.character(protein_dict[protein_dict$protein_id == protein_name,'pdb_id'])
        #need to split pdb_id into pdb_id and chain before going into pdb_info
        pdb_split <- strsplit(pdb_id,'_')[[1]]
        pdb_id <- pdb_split[1]
        chain <- pdb_split[2]
        pdb_read <- check_download_read_pdb(pdb_id = pdb_id)
        resno_and_resid <- quick_resno_and_resid(pdb_read = pdb_read, chain = chain)

        pdb_info <- list(pdb_id=pdb_id,
                         chain=chain,
                         resno=resno_and_resid$resno,
                         resid=resno_and_resid$resid)

        #if not NULL --> skip the protein_names already run?
        #will need to get pdb_info from the data
        #may just need the id, resid, and resno
        #do pairwise alignment anyway with the PDB file?
        #check_download with strsplit(pdb_id,'_)[[1]][1]?

      } else { #end if(!is.null(protein_dict)){

        #no protein dictionary entered --> display the menus to get pdb_info that way

        if(!(protein_name %in% names(protein_names_already_run))){
          pdb_info <- display_preferred_pdb_structure_menu(protein_name = protein_name, fasta_file)
          protein_names_already_run[[protein_name]] <- pdb_info
        } else { #end if(!(protein_name %in% names(protein_names_already_run)))
          pdb_info <- protein_names_already_run[[protein_name]]

        } #end else to if(!(protein_name %in% names(protein_names_already_run)))


      } #end else to if(!is.null(protein_dict)){

      #if alignment can't be done --> need to change database_name


      #if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
      if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),binding_seq)[[1]]) == 0){

        #try a pairwise alignment
        # if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) > 0){
        #   #the input sequence is found in the
        #
        # }


        #pwa_results <- pairwiseAlignment(paste0(pdb_info$resid,collapse=''),binding_seq)
        #score(pwa_results) --> do by the number of dashes?



        cat('Binding sequence found in FASTA not found in PDB file\n')
        #bs_output_db <- rbd.menuDBSearch(input_sequence,fasta_file,protein_name,pdb_info = pdb_info)

        #bs_output <- rbd.menuDBSearch('NEFISEYCGEIISQDEADRR',fasta_file,protein_name)
        #next
      } else { #end if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
        #the input sequence really has shown up in the pdb_info

        bs_s_e <- str_locate_all(paste0(pdb_info$resid,collapse=''),binding_seq)[[1]]
        bs_output$binding_site_start <- pdb_info$resno[bs_s_e[1]]
        bs_output$binding_site_end <- pdb_info$resno[bs_s_e[2]]
        ms_s_e <- str_locate_all(paste0(pdb_info$resid,collapse=''),as.character(bs_output$ms2_peptide_seq))[[1]]
        if(length(ms_s_e) == 0){
          bs_output$ms2_start <- NA
          bs_output$ms2_end <- NA
        } else {
          bs_output$ms2_start <- pdb_info$resno[ms_s_e[1]]
          bs_output$ms2_end <- pdb_info$resno[ms_s_e[2]]
        }

        bs_output$source_sequence <- paste0(pdb_info$resid,collapse='')

        bs_output$db <- 'PDB'
        bs_output$db_id <- paste0(pdb_info$pdb_id,'_',pdb_info$chain)

        # bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
        #                                protease = protease,
        #                                sequence_for_alignment=paste0(pdb_info$resid,collapse=''),
        #                                protein_name = protein_name,
        #                                database_name = 'PDB',
        #                                database_id = paste0(pdb_info$pdb_id,'_',pdb_info$chain))
      } #end else to if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0)

    } else if(align_to == 'UniProt'){ #end if(align_to == 'PDB')

      #how to get the UniProt ID --> blast.pdb with the swissprot database?
      #will need to update the menu function
      #check if there's a

      if(is.null(protein_to_uniprot_id)){
        #if it is NULL --> need to search the database
        #should compile a list for future use? or it will be part of the bs_output so no need?

        if(!(protein_name %in% names(protein_names_already_run))){
          pdb_info <- go_through_menu_loops(fasta_sequence_vector = fasta_seq, database='swissprot')
          protein_names_already_run[[protein_name]] <- pdb_info
        } else { #end if(!(protein_name %in% names(protein_names_already_run)))
          pdb_info <- protein_names_already_run[[protein_name]]

        } #end else to if(!(protein_name %in% names(protein_names_already_run)))


        isoform_num <- substrRight(pdb_info$pdb_id,2)
        if((startsWith(isoform_num,'.')) & (!is.na(as.numeric(substrRight(pdb_info$pdb_id,1))))){
          if(as.numeric(substrRight(pdb_info$pdb_id,1)) > 1){
            cat('Different isoform than 1 detected but using isoform 1 from Uniprot. Recommend downloading FASTA file for analysis.')
          }
        }

        uniprot_id <- strsplit(pdb_info$pdb_id,'\\.')[[1]][1]

        uniprot_info <- uniprot(pdb_info$pdb_id)
        uniprot_seq <- uniprot_info$sequence

      } else { #end if(is.null(protein_to_uniprot_id))
        #something is here and should match the FASTA files to those in the file

        #will search for the files ending in .fa or .fasta

        uniprot_id <- as.character(protein_to_uniprot_id[protein_to_uniprot_id$protein_id == protein_name,'uniprot_id'])

        if(paste0(uniprot_id,'.fa') %in% list.files()){
          uniprot_seq <- paste0(toupper(seqinr::read.fasta(paste0(uniprot_id,'.fa'))),collapse='')
        } else if(paste0(uniprot_id,'.fasta') %in% list.files()){
          uniprot_seq <- paste0(toupper(seqinr::read.fasta(paste0(uniprot_id,'.fasta'))),collapse='')
        } else {
          #does not exist in current directory --> fetching from UniProt
          cat("\nCan't find FASTA file in current directory --> fetching from UniProt\n")
          uniprot_seq <- uniprot(uniprot_id)$sequence
        } #end else to if(paste0(uniprot_id,'.fa') %in% list.files()){




      } #end else to if(is.null(protein_to_uniprot_id))

      #return(uniprot_seq)

      #grepl(uniprot_info$accession,pdb_info$pdb_id)

      #return(length(str_locate_all(uniprot_seq,binding_seq)[[1]]))

      #if(length(str_locate_all(uniprot_seq,input_sequence)[[1]]) == 0){
      if(length(str_locate_all(uniprot_seq,binding_seq)[[1]]) == 0){
        #bs_output <- rbd.menuDBSearch(input_sequence,fasta_file,protein_name,pdb_info = pdb_info)
        cat('Binding sequence not found in chosen UniProt sequence')

      } else { #end if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
        # bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
        #                                protease = protease,
        #                                sequence_for_alignment=uniprot_seq,
        #                                protein_name = protein_name,
        #                                database_name = 'UniProt',
        #                                database_id = uniprot_id)

        bs_s_e <- str_locate_all(uniprot_seq,binding_seq)[[1]]
        #return(bs_s_e)

        #return(bs_output)
        bs_output$binding_site_start <- bs_s_e[1]
        bs_output$binding_site_end <- bs_s_e[2]
        ms_s_e <- str_locate_all(uniprot_seq,as.character(bs_output$ms2_peptide_seq))[[1]]
        #return(ms_s_e)
        bs_output$ms2_start <- ms_s_e[1]
        bs_output$ms2_end <- ms_s_e[2]
        bs_output$source_sequence <- uniprot_seq

        bs_output$db <- 'UniProt'
        bs_output$db_id <- uniprot_id


        #return(bs_output)

      } #end else to if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){


    } else if(align_to == 'FASTA') { #end else if(align_to == 'UniProt')

      #if chosen --> get the FASTA sequence and set as the sequence for rbd.getBindingSeq
      bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                     protease = protease,
                                     sequence_for_alignment=toupper(paste0(fasta_file[[protein_name]],collapse='')),
                                     protein_name = protein_name,
                                     database_name = 'FASTA',
                                     database_id = protein_name,
                                     include_ambiguous = include_ambiguous)


      # if(include_ambiguous == TRUE){
      #   if(nchar(as.character(bs_output$binding_sequence)) == 0){
      #
      #     #levels(bs_output$binding_sequence) <- c(levels(bs_output$binding_sequence),as.character(bs_output$ms2_peptide_seq))
      #     bs_output$binding_sequence <- as.character(bs_output$ms2_peptide_seq)
      #     bs_output$binding_site_start <- as.character(bs_output$ms2_start)
      #     bs_output$binding_site_end <- as.character(bs_output$ms2_end)
      #
      #   } #end if(nchar(as.character(bs_output$binding_sequence))){
      # } #end if(include_ambigious == TRUE){


    } #end else if(align_to == 'FASTA')

    #can then run rbd.getBindingSeq out here

    bs_output_mega <- rbind(bs_output_mega,bs_output)

    #return(bs_output_mega)
  } #end for(num in 1:nrow(filtered_input_eluate_table)){

  return(bs_output_mega)

} #end function rbd.getBSfromIET

