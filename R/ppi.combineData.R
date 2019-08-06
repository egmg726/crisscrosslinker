#'PPI: Combine Data
#'
#'This function combines data and lets you set specifications for how you would like to filter
#'the data.
#'
#'@param xlink_mega_list List made by ppi.loadData()
#'@param fasta_file FASTA file as loaded by seqinr::read.fasta()
#'@param freq_cutoff Frequency cutoff, will discard values below this number. Defaults to NA.
#'@param score_cutoff Score cutoff, will discard values above this number. Defaults to NA.
#'@param cutoff_cond Cutoff conditions. If you have values for both freq_cutoff and score_cutoff, this will determine if it will use an "and" or "or" statement.
#'@param protein_alternative_names_dict data.frame() or string representing file name. If you have the same proteins represented by different names (NOT RECOMMENDED), you can use this to make sure they will all be combined successfully. Defaults to NULL.
#'@param category_color_input_file data.frame() or string representing file name. If you have specific colors you would like to use, use this field.
#'@param console_messages Will display console messages about if a match has been made. Primarily for debugging purposes.
#'@export

ppi.combineData <- function(xlink_mega_list, fasta_file, freq_cutoff = NA,
                            score_cutoff = NA, cutoff_cond = c('and','or'),
                            protein_alternative_names_dict = NULL,
                            category_color_input_file = NULL, console_messages = FALSE){

  if((!is.null(category_color_input_file)) && (typeof(category_color_input_file) == 'character')){
    category_color_input_file <- read.csv(category_color_input_file)
  }

  if(typeof(protein_alternative_names_dict) == 'character'){
    protein_alternative_names_dict <- read.csv(protein_alternative_names_dict)
  }

  pymol_reds <- c('red','tv_red','raspberry','darksalmon','salmon','deepsalmon',
                  'warmpink','firebrick','ruby','chocolate','brown')
  pymol_greens <- c('green','tv_green','chartreuse','splitpea','smudge',
                    'palegreen','limegreen','lime','limon','forest')
  pymol_blues <- c('blue','tv_blue','marine','slate','lightblue','skyblue',
                   'purpleblue','deepblue','density')
  pymol_yellows <- c('yellow','tv_yellow','paleyellow','yelloworange',
                     'limon','wheat','sand')
  pymol_magentas <- c('magenta','lightmagenta','hotpink','pink','lightpink',
                      'dirtyviolet','violet','violetpurple','purple',
                      'deeppurple')
  pymol_cyans <- c('cyan','palecyan','aquamarine','greencyan','teal','deepteal',
                   'lightteal')
  pymol_oranges <- c('orange','tv_orange','brightorange','lightorange','yelloworange',
                     'olive','deepolive')

  pymol_colors <- c()
  for(index_num in 1:7){
    pymol_colors <- c(pymol_colors,pymol_reds[index_num])
    pymol_colors <- c(pymol_colors,pymol_blues[index_num])
    pymol_colors <- c(pymol_colors,pymol_yellows[index_num])
    pymol_colors <- c(pymol_colors,pymol_greens[index_num])
    pymol_colors <- c(pymol_colors,pymol_oranges[index_num])
    pymol_colors <- c(pymol_colors,pymol_magentas[index_num])
    pymol_colors <- c(pymol_colors,pymol_cyans[index_num])

  }

  #add in new variable for frequency_color_list to make the colors for ppi.pymol
  frequency_color_list <- pymol_colors
  seq_xlink_list <- c()
  pro_xlink_list <- c()
  dist_xlink_list <- c()
  freq_xlink_list <- c()
  files_xlink_list <- c()
  site1_xlink_list <- c()
  site2_xlink_list <- c()
  pdb1_xlink_list <- c()
  pdb2_xlink_list <- c()
  freq_color_xlink_list <- c()
  pro_pos1_xlink_list <- c()
  pro_pos2_xlink_list <- c()
  pro_name1_xlink_list <- c()
  pro_name2_xlink_list <- c()
  pep_seq1_xlink_list <- c()
  pep_seq2_xlink_list <- c()
  pep_pos1_xlink_list <- c()
  pep_pos2_xlink_list <- c()
  score_xlink_list <- c()

  xlink_df_list <- list() #list or dataframe?
  list_of_protein_names <- names(fasta_file)
  not_included_xl_sites <- c()
  #should be able to accept the list of protein names or the fasta file depending on the type
  #is the actual fasta file necessary for this loop or just the protein names

  #should also make a new directory??
  #can put all of the new files in the new directory
  #as well as the names of the files used in a text file
  #should have a better name for the directory that is being made

  #can make the name for the directory
  #keep track of the time it took to run the whole program?

  #file_names <- names(xlink_list)
  for(file_name in names(xlink_mega_list)){

    xlink_list <- xlink_mega_list[[file_name]]
    #once the data is loaded
    #will need to add to a list and/or dataframe

    #need another loop to go through the
    for(xlink_index in 1:length(xlink_list)){

      protein_split_list <- NULL
      seq_xlink <- NULL
      if(xlink_index > length(xlink_list)){
        break
      }
      seq_xlink <- xlink_list[[xlink_index]]$sequence_xlink
      if((length(xlink_list[[xlink_index]]$proteins_xlink) != 0)){
        seq_xlink <- strsplit(seq_xlink,':0')[[1]]
      }
      if((length(xlink_list[[xlink_index]]$proteins_xlink) != 0) && grepl("/",xlink_list[[xlink_index]]$proteins_xlink)){
        xlink_list[[xlink_index]]$proteins_xlink <- strsplit(xlink_list[[xlink_index]]$proteins_xlink,'/')[[1]]
      }
      for(protein_xlink_in_list in xlink_list[[xlink_index]]$proteins_xlink){
        protein_split_list <- extract_protein_name_and_peptide_num(split_proteins = protein_xlink_in_list,
                                                                   list_of_protein_names = list_of_protein_names)
        if(!is.null(protein_split_list)){
          pro_output <- protein_split_list
          new_protein_list <- c()
          for(pro_index1 in 1:length(protein_split_list)){
            protein_name <- protein_split_list[[pro_index1]][1]
            if(!is.null(protein_alternative_names_dict)){
              protein_name <- check_and_change_protein_name(protein_name,protein_alternative_names_dict,list_of_protein_names,fasta_file)
            }
            protein_split_list[[pro_index1]][1] <- protein_name
            new_protein_list <- c(new_protein_list,protein_name)
          }
          pro_xlink <- paste(new_protein_list[1],'(',pro_output[[1]][2],')-',new_protein_list[2],'(',pro_output[[2]][2],')',sep='')
          break
        }
      } #end for(protein_xlink_in_list in xlink_list[[xlink_index]]$proteins_xlink)

      if(is.null(protein_split_list)){
        nixl_site <- paste(xlink_list[[xlink_index]]$proteins_xlink[1],' - ',file_name)
        not_included_xl_sites <- c(not_included_xl_sites,nixl_site)
        next
      }

      score_xlink <- min(as.numeric(xlink_list[[xlink_index]]$score_xlink))

      if(is.null(pro_xlink) || is.null(seq_xlink)){
        next
      }

      #print(xlink_index)

      seq_index <- match(seq_xlink,seq_xlink_list) #if == NA, means that it has not been added to the list yet
      pro_index <- match(pro_xlink,pro_xlink_list)

      if(length(unique(is.na(c(seq_index,pro_index)))) == 2){
        #both TRUE and FALSE exist in the list
        #if both are true --> probably should re-add to the list

        seq_index <- NA
        pro_index <- NA

        #should there be a note added?

      }


      if(seq_index != pro_index && !is.na(seq_index) && !is.na(pro_index)){
        if(pro_xlink_list[seq_index] == pro_xlink){
          pro_index <- seq_index
        } else if(seq_xlink_list[pro_index] == seq_xlink){
          seq_index <- pro_index
        } else {
          if(console_messages == TRUE){
            cat(paste('Index error for',seq_xlink,'and',pro_xlink))
          }
          warning((paste('Index error for',seq_xlink,'and',pro_xlink)))
        }
      }

      if((seq_index == pro_index) && !is.na(seq_index)){ #make sure that the indeces match and are real numbers
        #indiciates that they have already been added to the lists respectively

        #check if the file has not already been checked?

        if(!grepl(file_name,files_xlink_list[seq_index])){

          #need to change this so it will change it to the right frequency/color based on
          #file_name

          if(!is.null(category_color_input_file)){
            #if it is not null, will need to change the color based on the
            new_color <- as.character(category_color_input_file[category_color_input_file[['file_name']] == file_name,][['color']])
            #then check if the old color matches the new color

            if(new_color != freq_color_xlink_list[seq_index]){
              #colors don't match --> need to update the color

              #update the frequency +1 same as before
              freq_xlink_list[seq_index] <- freq_xlink_list[seq_index] + 1
              freq_color_xlink_list[seq_index] <- 'grey'


            }

            #if it does --> do not change the frequency
            #it doesn't --> change the color to "grey" and update the frequency

          } else {

            freq_xlink_list[seq_index] <- freq_xlink_list[seq_index] + 1
            freq_color_xlink_list[seq_index] <- frequency_color_list[freq_xlink_list[seq_index]]
            if(score_xlink_list[seq_index] > score_xlink){
              score_xlink_list[seq_index] <- score_xlink
            }

          }



          #will using a comma confuse the csv file? --> using '+' for now

          #will also have to have a contingency here for the file_name name if a list and
          #not a character



          files_xlink_list[seq_index] <- paste(files_xlink_list[seq_index],'+',file_name,sep='')


        }

      } else { #need to double check to make sure that both are not NA?
        #have not been added to the list
        #seq_xlink_list <- c(seq_xlink_list,seq_xlink)
        #pro_xlink_list <- c(pro_xlink_list,pro_xlink)
        #need to keep track of index or if everything works well, that will be taken care of already?

        seq_xlink2 <- seq_xlink

        sequence_xlink_list <- split_sequences2(xlink_list[[xlink_index]]$sequence_xlink, proteins = FALSE)
        #protein_split_list <- split_sequences2(xlink_list[[xlink_index]]$proteins_xlink[1], proteins = TRUE)

        for(index_num in 1:length(sequence_xlink_list)){ #go through peptides (1 and 2)
          seq_xlink <- sequence_xlink_list[[index_num]]
          pro_split <- protein_split_list[[index_num]]
          protein_name <- pro_split[1]
          protein_pos <- as.numeric(pro_split[2])
          peptide_seq <- seq_xlink[1]
          peptide_pos <- as.numeric(seq_xlink[2])


          if(!(protein_name %in% names(list_of_protein_names))){ #does it have a PDB file?
            #if it does not have a PDB sturcture, check if it has a fasta file
            #first check if there's a fasta file or check the dictionary?

            #check if it has a valid number?



            protein_name_has_been_changed <- FALSE

            #then check here if it has an entry in the protein dictionary
            if(!is.null(protein_alternative_names_dict)){

              #protein_alternative_names_dict <- read.csv(protein_alternative_names_dict)

              #copy the code here for the protein alternative names dict
              #change the name of variable?

              protein_original_names <- as.character(protein_alternative_names_dict$original_name)

              #cat(protein_original_names)

              if(!(protein_name %in% protein_original_names)){

                #check if the protein_name has an alternative

                #should maybe go through the rows

                for(row_num in 1:nrow(protein_alternative_names_dict)){


                  pnames_a_row <- protein_alternative_names_dict[row_num,]

                  pnames_as_list <- as.character(t(pnames_a_row))

                  if(protein_name %in% pnames_as_list){
                    #change the protein name to the first column
                    #check to make sure that it's named

                    if(console_messages == TRUE){
                    cat(paste(protein_name,'changed to',as.character(protein_alternative_names_dict[row_num,'original_name']),
                              'in output\n'))
                    }
                    protein_name <- as.character(protein_alternative_names_dict[row_num,'original_name'])
                    #activate boolean?
                    #what to do if protein name does not show up in this list and does not have a
                    #fasta file?

                    protein_name_has_been_changed <- TRUE

                    #add protein to a list that will be outputted at the end of the loop

                    break

                  } #end if(protein_name %in% pnames_as_list)

                  #turn row into list
                  #see if protein name is in the list
                  #if it is in that list, make the protein name the first name in the list


                  #if it is not in the list,
                  #skip the protein (the crosslink entirely) and add the protein to a list of
                  #not included proteins
                  #should probably include all of the crosslinking information
                  #possibly indicate "null" or something for the protein structure and everything
                } #for(row_num in 1:nrow(pnames_alternatives))

                #if it does, rename the protein to the previous

                #if it does not, put up a message to the
              }



              #if protein name has been changed -->
              #protein_name_has_been_changed <- TRUE



              #if it's not in the protein dictionary, move to see if it has an entry
              #in the fasta file (same as the else if that comes afterwards)
              #can also change to just an if and add a boolean so it will only
              #go that if if the name has not been changed

            }

            if((protein_name %in% names(fasta_file)) && (protein_name_has_been_changed == FALSE)){
              #if it's not in the protein dictionary, check if it is in the fasta file

              if(console_messages == TRUE){
                cat('Is in fasta file but protein name has not been changed yet\n')
              }


              #ask if the user wants to create a separate PDB file for it

              #menu loops function here

            } else {
              #the user has not included it in neither the fasta file nor the dictionary
              #add the whole crosslinking site to a dataframe to be
              #outputted if the xl site does not meet the requirements

              if(console_messages == TRUE){
              cat('Not in fasta file or in the dictionary\n')
              }

            }

            #should also check protein_pos?

            #can move the is.null if statement within this if statement
            #since it might be helpful to check beforehand if there is a PDB file for
            #the protein

          } #end if(!(protein_name %in% names(list_of_protein_names)))



        } #end for(index_num in 1:length(sequence_xlink_list))

        for(index_num in 1:length(sequence_xlink_list)){ #go through peptides (1 and 2)
          seq_xlink <- sequence_xlink_list[[index_num]]
          pro_split <- protein_split_list[[index_num]]
          protein_name <- pro_split[1]
          protein_pos <- as.numeric(pro_split[2])
          peptide_seq <- seq_xlink[1]
          peptide_pos <- as.numeric(seq_xlink[2])


          if(!(protein_name %in% names(list_of_protein_names))){ #does it have a PDB file?
            #if it does not have a PDB sturcture, check if it has a fasta file
            #first check if there's a fasta file or check the dictionary?

            #check if it has a valid number?



            protein_name_has_been_changed <- FALSE

            #then check here if it has an entry in the protein dictionary
            if(!is.null(protein_alternative_names_dict)){

              #protein_alternative_names_dict <- read.csv(protein_alternative_names_dict)

              #copy the code here for the protein alternative names dict
              #change the name of variable?

              protein_original_names <- as.character(protein_alternative_names_dict$original_name)

              #cat(protein_original_names)

              if(!(protein_name %in% protein_original_names)){

                #check if the protein_name has an alternative

                #should maybe go through the rows

                for(row_num in 1:nrow(protein_alternative_names_dict)){


                  pnames_a_row <- protein_alternative_names_dict[row_num,]

                  pnames_as_list <- as.character(t(pnames_a_row))

                  if(protein_name %in% pnames_as_list){
                    #change the protein name to the first column
                    #check to make sure that it's named

                    if(console_messages == TRUE){
                    cat(paste(protein_name,'changed to',as.character(protein_alternative_names_dict[row_num,'original_name']),
                              'in output\n'))

                    }
                    protein_name <- as.character(protein_alternative_names_dict[row_num,'original_name'])
                    #activate boolean?
                    #what to do if protein name does not show up in this list and does not have a
                    #fasta file?

                    protein_name_has_been_changed <- TRUE

                    #add protein to a list that will be outputted at the end of the loop

                    break

                  } #end if(protein_name %in% pnames_as_list)

                  #turn row into list
                  #see if protein name is in the list
                  #if it is in that list, make the protein name the first name in the list


                  #if it is not in the list,
                  #skip the protein (the crosslink entirely) and add the protein to a list of
                  #not included proteins
                  #should probably include all of the crosslinking information
                  #possibly indicate "null" or something for the protein structure and everything
                } #for(row_num in 1:nrow(pnames_alternatives))

                #if it does, rename the protein to the previous

                #if it does not, put up a message to the
              }



              #if protein name has been changed -->
              #protein_name_has_been_changed <- TRUE



              #if it's not in the protein dictionary, move to see if it has an entry
              #in the fasta file (same as the else if that comes afterwards)
              #can also change to just an if and add a boolean so it will only
              #go that if if the name has not been changed

            }

            if((protein_name %in% names(fasta_file)) && (protein_name_has_been_changed == FALSE)){
              #if it's not in the protein dictionary, check if it is in the fasta file


              if(console_messages == TRUE){
              cat('Is in fasta file but protein name has not been changed yet\n')
              }

              #ask if the user wants to create a separate PDB file for it

              #menu loops function here

            } else {
              #the user has not included it in neither the fasta file nor the dictionary
              #add the whole crosslinking site to a dataframe to be
              #outputted if the xl site does not meet the requirements

              if(console_messages == TRUE){
              cat('Not in fasta file or in the dictionary\n')
              }

            }

            #should also check protein_pos?

            #can move the is.null if statement within this if statement
            #since it might be helpful to check beforehand if there is a PDB file for
            #the protein

          } #end if(!(protein_name %in% names(list_of_protein_names)))


          #if it does not have a PDB file but it does have a fasta file
          #ask if the user does not if they want to run it in the fasta

          #first check if the protein really does need to be changed
          #does it match a PDB file and/or a name in a fasta file?


          #if there is no

          #make sure



          #will also have to have protein conversion list
          #put function here changing the protein name

          #add any protein name/seq link not used
          #maybe dataframe with the seq_xlink and the protein name
          #should include the whole one (each in their own column?)


          #does this eject a NULL when there is an error?
          #can skip over this if pdb_numbering is TRUE

          #return(protein_split_list)

          #return(list(all_pdb=list_of_protein_names,start_and_end=list_of_start_and_end_pdbs))

          #add protein positions to another list?

          #cat('On this PDB Structure\n')
          #print(on_this_pdb_structure)



        } #end for(index_num in 1:length(sequence_xlink_list))

        #end if(length(otps_list) != 2)

        #cat('OTPS List')
        #print(otps_list)





        if(!is.null(category_color_input_file)){
          #need to make the color match the category
          new_color <- as.character(category_color_input_file[category_color_input_file[['file_name']] == file_name,][['color']])
          freq_color_xlink_list <- c(freq_color_xlink_list,new_color)

        } else {

          freq_color_xlink_list <- c(freq_color_xlink_list,frequency_color_list[1])

        }

        freq_xlink_list <- c(freq_xlink_list,1)

        if(typeof(file_name) == 'list'){
          #if list --> need to have some kind of character or something to replace it instead
          #just have a dummy name for now?
        }
        files_xlink_list <- c(files_xlink_list,file_name)

        #add each one to a list?

        pro_pos1 <- protein_split_list[[c(1,2)]]
        pro_pos2 <- protein_split_list[[c(2,2)]]
        pro_name1 <- protein_split_list[[c(1,1)]]
        pro_name2 <- protein_split_list[[c(2,1)]]

        pep_seq1 <- sequence_xlink_list[[c(1,1)]]
        pep_seq2 <- sequence_xlink_list[[c(2,1)]]
        pep_pos1 <- sequence_xlink_list[[c(1,2)]]
        pep_pos2 <- sequence_xlink_list[[c(2,2)]]

        pro_pos1_xlink_list <- c(pro_pos1_xlink_list,pro_pos1)
        pro_pos2_xlink_list <- c(pro_pos2_xlink_list,pro_pos2)
        pro_name1_xlink_list <- c(pro_name1_xlink_list,pro_name1)
        pro_name2_xlink_list <- c(pro_name2_xlink_list,pro_name2)
        pep_seq1_xlink_list <- c(pep_seq1_xlink_list,pep_seq1)
        pep_seq2_xlink_list <- c(pep_seq2_xlink_list,pep_seq2)
        pep_pos1_xlink_list <- c(pep_pos1_xlink_list,pep_pos1)
        pep_pos2_xlink_list <- c(pep_pos2_xlink_list,pep_pos2)
        score_xlink_list <- c(score_xlink_list,score_xlink)
        seq_xlink_list <- c(seq_xlink_list,seq_xlink2)
        pro_xlink_list <- c(pro_xlink_list,pro_xlink)


      } #end MEGA else





    } #end for(xlink_index in 1:length(xlink_list)){

  } #end for(file_name in names(xlink_list)){

  #optional input: freq_cutoff and score_cutoff that will be used to determine
  #if any data will b excluded from the final

  #step 1: set up the loop that will used for the rest of the code


  xlink_mega_list <- list(seq=seq_xlink_list,
                          pro=pro_xlink_list,
                          #dist=dist_xlink_list,
                          freq=freq_xlink_list,
                          freq_color=freq_color_xlink_list,
                          files=files_xlink_list,
                          #pdb1=pdb1_xlink_list,
                          #pdb2=pdb2_xlink_list,
                          pro_pos1=pro_pos1_xlink_list,
                          pro_pos2=pro_pos2_xlink_list,
                          pro_name1=pro_name1_xlink_list,
                          pro_name2=pro_name2_xlink_list,
                          pep_seq1=pep_seq1_xlink_list,
                          pep_seq2=pep_seq2_xlink_list,
                          pep_pos1=pep_pos1_xlink_list,
                          pep_pos2=pep_pos2_xlink_list,
                          score=score_xlink_list)

  #return(xlink_mega_list)

  #remove the .pdb from the file names in pdb1 and pdb2 afterwards?
  xlink_mega_df <- data.frame(xlink_mega_list)

  #add in the cutoffs here
  #can also exclude color if wished
  #freq_cutoff = NA, score_cutoff = NA

  #cutoff_cond
  #making conditional statements?


  #should check if both are NA
  #can change this to an and statement and then go from there
  #if both --> else freq_cutoff --> else score_cutoff --> else output saying nothing filtered
  if(!is.na(freq_cutoff) && !is.na(score_cutoff)){

    freq_bool <- NULL
    score_bool <- NULL

    if(!is.na(freq_cutoff)){
      #also be able to customize the sign?
      #make a separate function if this complicated?
      freq_bool <- xlink_mega_df$freq >= freq_cutoff
    }

    if(!is.na(freq_cutoff)){
      #also be able to customize the sign?
      #make a separate function if this complicated?
      score_bool <- xlink_mega_df$score <= score_cutoff
    }

    if(cutoff_cond == 'or'){
      xlink_mega_df <- xlink_mega_df[(freq_bool | score_bool),]
    } else if(cutoff_cond == 'and'){
      xlink_mega_df <- xlink_mega_df[(freq_bool & score_bool),]
    } else {
      warning('Invalid condition statement: select and/or')
    }


  } else if(!is.na(freq_cutoff)){

    freq_bool <- xlink_mega_df$freq >= freq_cutoff
    xlink_mega_df <- xlink_mega_df[(freq_bool),]

  } else if(!is.na(score_cutoff)){

    score_bool <- xlink_mega_df$score <= score_cutoff
    xlink_mega_df <- xlink_mega_df[(score_bool),]

  } else {
    if(console_messages == TRUE){
    cat('Nothing filtered\n')
    }
  }



  #new export option: automatic

  return(xlink_mega_df)

  #return xlink_df
} #end function ppi.combineData
