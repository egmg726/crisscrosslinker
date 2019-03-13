#'PPI: Match PDB
#'
#'This function matches combined data to PDB files.
#'
#'@param xlink_df data.frame() as created by ppi.combineData()
#'@param fasta_file FASTA file as loaded by seqinr::read.fasta()
#'@param pdb_numbering pdb_numbering
#'@param pdb_directory pdb_directory
#'@param pdb_match_vector pdb_match_vector
#'@param csv_pdb_input_file csv_pdb_input_file
#'@param fasta_names_to_generate_all_2d_structures fasta_names_to_generate_all_2d_structures
#'@param file_type_2d file_type_2d
#'@param protein_alternative_names_dict protein_alternative_names_dict
#'@export

ppi.matchPDB <- function(xlink_df, fasta_file, pdb_numbering = FALSE,
                         pdb_directory = NULL, pdb_match_vector = NULL,
                         csv_pdb_input_file = NULL,
                         fasta_names_to_generate_all_2d_structures = NULL,
                         protein_alternative_names_dict = NULL,
                         download_pdb = TRUE,
                         file_type_2d = 'single_dot'){

  #dist_xlink_list <- c()
  #pdb1_xlink_list <- c()
  #pdb2_xlink_list <- c()

  if(!is.null(pdb_directory)){
    experiment_directory <- pdb_directory
  } else {
    experiment_directory <- getwd()
  }

  xlink_list <- list()

  list_of_protein_names <- names(fasta_file)

  if(!is.null(protein_alternative_names_dict)){

    protein_alternative_names_dict <- read.csv(protein_alternative_names_dict)

    for(pand_name in names(protein_alternative_names_dict)){
      #add entire column to list

      pand_list <- as.character(protein_alternative_names_dict[,pand_name])
      list_of_protein_names <- c(list_of_protein_names,pand_list)

    }

    #add to list of protein names
  }

  list_of_protein_names <- unique(list_of_protein_names)

  current_directory <- getwd()


  #Sys.setenv(TZ='Australia/Melbourne')
  #create a system log file that has sys time on it instead?
  #what about the %03d?
  #new_directory_title <- paste('make_diff_analysis Run on',Sys.time())
  #new_directory_title <- gsub(':','-',new_directory_title)

  #new_directory <- paste(current_directory,'/',new_directory_title,sep = '')
  #dir.create(new_directory)

  pdb_suffix <- file_type_2d

  # if(file_type_2d == 'single dot'){
  #   pdb_suffix <- 'single_dot'
  # } else {
  #   pdb_suffix <- file_type_2d
  # }



  if(is.null(csv_pdb_input_file)){
    if(pdb_numbering == TRUE){
      if(is.null(pdb_match_vector)){
        pdb_match_vector <- ppi.alignPDB(fasta_file = fasta_file)
      }
    } else {
      #pdb_numbering == FALSE --> do what has been done before
      generated_pdb_lists <- generate_pdb_lists_and_files_from_fasta(fasta_file,pdb_suffix,
                                                                     fasta_names_to_generate_all_2d_structures)

      #output this as a csv file?
      #if option is selected yes --> can output as csv file
      make_start_end_pdb_df_output2(generated_pdb_lists = generated_pdb_lists)
      write.csv(make_start_end_pdb_df_output2(generated_pdb_lists = generated_pdb_lists),'pdb_start_stop.csv')
    }

  } else { #if someone has indicated a PDB file to be used

    generated_pdb_lists <- generate_pdb_lists_from_pdb_csv(csv_pdb_input_file)
    #need same input as list of all pdb files
  }

  #return(generated_pdb_lists)


  #if pdb_numbering == TRUE --> will need a substitute for these or otherwise ignore them
  if(pdb_numbering == FALSE){
    list_of_all_pdb_files <- generated_pdb_lists$all_pdb_files
    list_of_start_and_end_pdbs <- generated_pdb_lists$start_end_chain
  } else {
    #end if(pdb_numbering == FALSE)
    list_of_all_pdb_files <- pdb_match_vector
  }

  #xlink_df <- xlink.df
  xlink_df$dist <- rep(NA,nrow(xlink_df))
  xlink_df$pdb1 <- rep(NA,nrow(xlink_df))
  xlink_df$pdb2 <- rep(NA,nrow(xlink_df))

  stored_pdbs <- list()

  #-----xlink_df cycle starts---------

  for(row_num in 1:nrow(xlink_df)){
    xlink_row <- xlink_df[row_num,]
    otps_list <- c()
    xyz_coord_list <- list()
    for(pro_num in 1:2){
      protein_name <- as.character(xlink_row[[paste0('pro_name',pro_num)]])
      protein_pos <- as.character(xlink_row[[paste0('pro_pos',pro_num)]])
      if(pdb_numbering == TRUE){
        protein_pos_match <- match(protein_pos,pdb_match_vector[[protein_name]]$fasta)
        pdb_protein_pos <- as.numeric(pdb_match_vector[[protein_name]]$pdb[protein_pos_match])
        xlink_row[[paste0('pro_pos',pro_num)]] <- pdb_protein_pos
        #protein_split_list[[index_num]][2] <- pdb_protein_pos
        #xlink_row
        on_this_pdb_structure <- pdb_match_vector[[protein_name]]$chain[protein_pos_match]
        xlink_df[row_num,paste0('pdb',pro_num)] <- on_this_pdb_structure
        otps_list <- c(otps_list,on_this_pdb_structure)
      } else {
        on_this_pdb_structure <- get_pdb_structure_from_protein_pos(protein_name,
                                                                    protein_pos,
                                                                    list_of_all_pdb_files,
                                                                    list_of_start_and_end_pdbs)



        if(!is.null(on_this_pdb_structure)){
          otps_list <- c(otps_list,on_this_pdb_structure)
          xlink_df[row_num,paste0('pdb',pro_num)] <- on_this_pdb_structure
        } #end if(!is.null(on_this_pdb_structure)){
      } #end else to if(pdb_numbering == TRUE)
    } #end for(pro_num in 1:2){

    pdbs_are_same <- TRUE
    if(pdb_numbering == TRUE){
      pdbs_are_same <- (length(unique(unlist(strsplit(otps_list,'_'))[c(T,F)])) == 1)
    }


    if((length(otps_list) != 2) || ('-' %in% otps_list) || (pdbs_are_same == FALSE)){

      cat('Skipping this crosslinking site, adding to list of XL sites not listed\n')
      next
      #cat('Skipping this crosslinking site, adding to list of XL sites not listed')
      #add crosslinking site to list of not included crosslinking sites?
      #list or dataframe?

      #not_included_xl_sites <- c(not_included_xl_sites,seq_xlink2)

    } #end if(length(otps_list) != 2)

    #should change where these are?
    #pdb1_xlink_list <- c(pdb1_xlink_list,otps_list[1])
    #pdb2_xlink_list <- c(pdb2_xlink_list,otps_list[2])

    xlink_df[row_num,paste0('pdb',pro_num)] <- on_this_pdb_structure

    # if(row_num == 5){
    #   return(otps_list)
    # }


    if((!(FALSE %in% grepl('___',otps_list))) || pdb_numbering == TRUE){ #make sure that it is all TRUEs to go forward with dist calculation
      for(otps_index in 1:length(otps_list)){ #need to load the protein position as well
        #otps_index <- 1
        otps <- otps_list[otps_index]

        if(pdb_numbering == TRUE){

          #can check the temporary directory for the PDB file
          #see if pdb_id.pdb exists within the tempdir
          #if yes, can
          otps_split <- strsplit(otps,'_')[[1]]

          # if((row_num == 5) && (otps_index == 2)){
          #   return(otps_split)
          # }

          pdb_id <- otps_split[1]
          chain <- otps_split[2]




          if(is.na(pdb_id)){
            break
          }

          #print(pdb_id)
          if(paste0(pdb_id,'.pdb') %in% list.files(tempdir())){

            pdb_read <- read.pdb2(paste0(tempdir(),'/',pdb_id,'.pdb'))

          } else {

            if(!(pdb_id %in% names(stored_pdbs))){
              pdb_read <- read.pdb2(pdb_id)
              stored_pdbs[[pdb_id]] <- pdb_read
            } else {
              pdb_read <- stored_pdbs[[pdb_id]]
            }


          } #end if(paste0(pdb_id,'.pdb') %in% list.files(tempdir()))


          #once the pdb has been read --> get just the chain?

          pdb_read$atom <- pdb_read$atom[pdb_read$atom$chain == chain,]

          #use the read.pdb function to be able to

        } else { #end if(pdb_numbering == TRUE)

          if(!is.null(pdb_directory)){
            #paste the directory to otps
            if(!endsWith('/',pdb_directory)){
              pdb_directory <- paste(pdb_directory,'/',sep='')
            }

            pdb_path <- paste(pdb_directory,otps,sep='')
            pdb_read <- read.pdb2(pdb_path)

          } else {
            pdb_read <- read.pdb2(otps)
          }
        } #end  else to if(pdb_numbering == TRUE)

        # if((row_num == 5) && (otps_index == 2)){
        #   return(pdb_read)
        # }

        #xyz_coords <- get_xyz_coordinates_pdb(pdb_read)
        protein_name <- as.character(xlink_row[[paste0('pro_name',otps_index)]])
        protein_pos <- as.character(xlink_row[[paste0('pro_pos',otps_index)]])
        #if pdb_numbering == TRUE, this may need to be changed slightly
        #pos_match <- match(protein_pos,pdb_read$atom$resno)



        # if((row_num == 5) && (otps_index == 2)){
        #   return(pos_match)
        # }

        # if(!is.na(pos_match)){ #checking if there is a true match
        #   atom_match <- pdb_read$atom$elety[pos_match]
        #   while(atom_match != 'CA'){ #can choose any atom here, can change to N later on
        #     #make 'CA' as a variable when this is made into a function
        #     #with 'CA' as potential default in function
        #     pos_match <- pos_match + 1
        #     atom_match <- pdb_read$atom$elety[pos_match]
        #   }
        #   if(pdb_read$atom$resno[pos_match] == protein_pos){
        #
        #     #should anything be done in here?
        #     pos_match <- pos_match
        #
        #   } else {
        #     #what to do if it has been overshot
        #     cat('No match to chosen atom, reverting to first match')
        #     pos_match <- match(protein_pos,pdb_read$atom$resno)
        #
        #   }
        # } else {
        #   #if it is NA and there is no match
        #   warning("No match\n")
        # }
        #
        #

        # xyz_matches <- c()
        # #make list of the xyz coordinates
        # for(xyz_name in names(xyz_coords)){
        #   #xyz1 <- xyz_coords[[xyz_name]][pos_match]
        #   xyz_matches <- c(xyz_matches,xyz1)
        # }
        # xyz_coord_list[[otps_index]] <- xyz_matches

        pdb_read_chain <- pdb_read$atom[pdb_read$atom$resno == protein_pos,]
        xyz_matches <- as.numeric(pdb_read_chain[pdb_read_chain$elety == 'CA',][c('x','y','z')])
        xyz_coord_list[[otps_index]] <- xyz_matches

      } #end of otps_index in 1:length(opts_list) for loop
      #calculate distance here
      #double-check formula before using in the final version of distance

      # if(row_num == 5){
      #   return(xyz_coord_list)
      # }

      if(length(xyz_coord_list) == 2){
        dist_xyz <- dist.xyz(xyz_coord_list[[1]],xyz_coord_list[[2]])[1,1]
        #dist_xlink_list <- c(dist_xlink_list,dist_xyz)
        xlink_df[row_num,'dist'] <- dist_xyz
      } else {
        #dist_xlink_list <- c(dist_xlink_list,NA)
      }
    } else { #one or more PDB structures is not "real" so distance should not be calculated
      #dist_xlink_list <- c(dist_xlink_list,NA)
    }
  } #end for(row_num in 1:nrow(xlink_df))

  # if(length(dist_xlink_list) != nrow(xlink_df)){
  #   warning('Number of distances does not match original data.frame length\nReturning NULL')
  #   return(NULL)
  # }

  #xlink_df$dist <- dist_xlink_list
  #xlink_df$pdb1 <- pdb1_xlink_list
  #xlink_df$pdb2 <- pdb2_xlink_list

  return(xlink_df)

} #end ppi.matchPDB function

