#'Write PyMOL file for protein-protein interaction data
#'
#'Generates a file for PyMOL based on the XL sites aligned to PDB structures.
#'
#'@param xlink.df.pdb data.frame output from ppi.matchPDB() or ppi.matchPDB2()
#'@param list_of_start_and_end_pdbs list_of_start_and_end_pdbs
#'@param show_only_real_structures show_only_real_structures
#'@param write_file If TRUE, will write the produced .pml file.
#'@param colors A list of colors as R names or hexcodes or a palette defined by RColorBrewer/viridis. If left blank, will use the colors that are in the column freq_color.
#'@param color_by If colors is not NULL, choose either 'freq' or 'dist' for which variable should be colored.
#'@param file.name File name for the exported PyMOL file. Should end in '.pml'
#'@export

ppi.pymol <- function(xlink.df.pdb,list_of_start_and_end_pdbs = NULL,show_only_real_structures = NULL, write_file = TRUE,
                      colors = NULL, color_by = 'freq', experiment_directory = NULL,
                      pdb_numbering = FALSE,file.name='xlink.pml'){
  #rename the columns to the original column names for easier integration to the rest of
  #the function

  if(is.null(experiment_directory)){
    experiment_directory <- getwd()
  }

  if(!is.null(colors)){
    custom.color <- TRUE
  } else {
    custom.color <- FALSE
  }

  #if there is no start and end PDBs can just get the info from the PDB file itself?
  #sd_pdb <- read.pdb2('Spc98-yeast_missing_sequence_3_single_dot.pdb')
  #if(endsWith(pdb_name, 'single_dot.pdb'))
  #sd_pdb$atom$resno #this will be the position for the making of the pymol file

  #colnames(xlink.df.pdb)

  if('X' %in% colnames(xlink.df.pdb)){
    xlink.df.pdb$X <- NULL
  }

  og_xlink_colnames <- c("seq","pro","dist","freq","freq_color","files","pdb1","pdb2",
                         "pro_pos1","pro_pos2","pro_name1","pro_name2","pep_seq1","pep_seq2",
                         "pep_pos1","pep_pos2","score")

  if('Protein.Position.1' %in% colnames(xlink.df.pdb)){
    colnames(xlink.df.pdb) <- og_xlink_colnames
  }

  pymol_lines <- c()

  #if custom color is TRUE
  #would vars be frequency?
  #xlink.df.pdb <- xl_df2

  if(color_by == 'dist'){
    #do the cut here
    #have another variable for the cutoff for dist?
    #this will also need to be true even if custom.color is FALSE?

    #xlink.df.pdb.pdb$freq_color <-

    #will need to choose the colors by the custom.color
    #if custom.color is false can just choose a random palette?

    #can customize
    #dist.int <- c('start','end','int')


    xlink.df.pdb$freq_color <- cut(xlink.df.pdb$dist, c(-Inf,10,15,20,Inf), c("blue", "green", "orange","red"))

    #xlink.df.pdb.pdb

    #can just have a differe

    pymol_cc <- color.pymol(1:4, colors = c("blue", "green", "orange","red"))
    pymol_lines <- c(pymol_lines,pymol_cc$set_colors)

  } #end if(color_by == 'dist')


  if((custom.color == TRUE) && color_by != 'dist'){


      pymol_cc <- color.pymol(sort(unique(xlink.df.pdb[[color_by]])), colors = colors)


    #pymol_cc$vars
    #pymol_cc$color_names

    pymol_lines <- c(pymol_lines,pymol_cc$set_colors)

    levels(xlink.df.pdb$freq_color) <- c(levels(xlink.df.pdb$freq_color),pymol_cc$color_names)

    #xlink.df.pdb$freq

    for(index_num in 1:length(pymol_cc$vars)){

      var_name <- pymol_cc$vars[index_num]

      #will rename the freq_color by whatever the var name is


      xlink.df.pdb[xlink.df.pdb[[color_by]] == var_name, 'freq_color'] <- pymol_cc$color_names[index_num]

      #within(xlink.df.pdb, freq_color[color_by == var_name] <- pymol_cc$color_names[index_num])

      #need to select and change the variable name based on the value


    } #end for(index_num in 1:length(pymol_cc$vars))

  } else {#end if custom color == TRUE
    #custom_color == FALSE

    #make the PyMOL legend still even if they're using the original colors
    #just input those colors into the colors/vars

  }
  #should still produce a legend even if they have not opted for a custom color


  #would have to replace the values in freq_color


  is_pdb_match_vector <- FALSE
  pdbs_in_df <- unique(c(levels(xlink.df.pdb$pdb1),levels(xlink.df.pdb$pdb2)))
  if(is.null(pdbs_in_df)){
    pdbs_in_df <- unique(c((xlink.df.pdb$pdb1),(xlink.df.pdb$pdb2)))
  }

  #return(unique(c((xlink.df.pdb$pdb1),(xlink.df.pdb$pdb2))))

  if(!(grepl('\\.pdb',paste0(pdbs_in_df,'|',collapse = '')))){
    #If TRUE --> it's a pdb_match_vector
    is_pdb_match_vector <- TRUE
    pdbs_in_df <- pdbs_in_df[pdbs_in_df != '-']
    #return(pdbs_in_df)
    pdbs_in_df <- unique(unlist(strsplit(pdbs_in_df,'_'))[c(T,F)])
    for(pdb_name in pdbs_in_df){
      py_line <- paste0('fetch ',pdb_name,",async=0")
      pymol_lines <- c(pymol_lines,py_line)
    }
  } else { #end if(!(grepl('\\.pdb',paste0(pdbs_in_df,'|',collapse = '')))){
    for(pdb_name in pdbs_in_df){
      #load pdb files

      if(!is.null(show_only_real_structures)){
        pdb_in_sors <- FALSE
        for(sors in show_only_real_structures){
          #noob solution to current problem
          if(grepl(sors,pdb_name)){
            pdb_in_sors <- TRUE
          }


        } #end first for(sors in show_only_real_structures)


        #will also need to include option in case real PDB structures are being used
        #won't have .pdb at the end?
        #will split the ID and the chain to be able to get the

        for(sors in show_only_real_structures){
          #include the line of code if the protein is in list and is real
          #or or if it does not show up in the list at all
          #or statement needs to be able to



          if((grepl(sors,pdb_name) && grepl('___',pdb_name))){

            #do the code here
            py_line <- paste('load ',experiment_directory,'/',pdb_name,
                             sep='')
            pymol_lines <- c(pymol_lines,py_line)

          } #end if((grepl(sors,pdb_name) && grepl('___',pdb_name)) || !grepl(sors,pdb_name))
        } #end for(sors in show_only_real_structures)

        if(!pdb_in_sors){

          py_line <- paste('load ',experiment_directory,'/',pdb_name,
                           sep='')
          pymol_lines <- c(pymol_lines,py_line)


        }

      } else {#end if(!is.null(show_only_real_structures))

        py_line <- paste('load ',experiment_directory,'/',pdb_name,
                         sep='')
        pymol_lines <- c(pymol_lines,py_line)


      }

      # py_line <- paste('load ',experiment_directory,'/',pdb_name,
      #                  sep='')
      # pymol_lines <- c(pymol_lines,py_line)


      #if this boolean is turned on,
      #color color_name, protein_name

      #should each of the PDBs be colored based on the protein they come from
      #can make this a user option

    } #end for(pdb_name in pdbs_in_df)
  } #end else to if(!(grepl('\\.pdb',paste0(pdbs_in_df,'|',collapse = '')))){


  #loop through the pdbs
  #pdb_name <- pdbs_in_df[1]


  #can just have one option
  #if show_surface, color_gray == TRUE
  #py_line <- 'show surface'
  #pymol_lines <- c(pymol_lines,py_line)
  py_line <- 'color gray'
  pymol_lines <- c(pymol_lines,py_line)

  mega_distance_count <- 0
  for(row_num in 1:nrow(xlink.df.pdb)){

    if(!endsWith(as.character(xlink.df.pdb$pdb1[row_num]),'pdb')){
      #this means it is a pdb_match_vector
      pdb1 <- strsplit(as.character(xlink.df.pdb$pdb1[row_num]),'_')[[1]][1]
      pdb2 <- strsplit(as.character(xlink.df.pdb$pdb2[row_num]),'_')[[1]][1]
    } else {
      pdb1 <- strsplit(as.character(xlink.df.pdb$pdb1[row_num]),'.pdb')[[1]]
      pdb2 <- strsplit(as.character(xlink.df.pdb$pdb2[row_num]),'.pdb')[[1]]
    }


    #if grepl '___' get last character for the chain
    #otherwise use 'A'

    pro_pos1 <- as.character(xlink.df.pdb$pro_pos1[row_num])
    pro_pos2 <- as.character(xlink.df.pdb$pro_pos2[row_num])

    #if statement to select the right chain


    #can do the check here in the grepl to do the check if it a single dot PDB structure, linear or curved

    #should make a list of pdbs so that this is not repeated twice?
    #can make output as a list that is then used for the next part of the function

    pdb_names <- c(pdb1,pdb2)
    pro_pos_list <- c(pro_pos1,pro_pos2)
    #file_type_2d <- 'single atom'

    #store new variables in a list?
    #check if list is empty or store the original values inside if not single atom?

    selection_name_list <- list()

    pdb_num_count <- 0

    for(pdb_num in pdb_names){
      pdb_num_count <- pdb_num_count + 1
      chain_name <- paste('chain',as.character(pdb_num_count),sep='')
      pro_pos_name <- paste('pro_pos',pdb_num_count,sep='')
      pdb_num_name <- paste('pdb',pdb_num_count,sep='')
      selection_name_list[[pdb_num_name]] <- pdb_num

      if(grepl('___',pdb_num)){ #if 'real' PDB structure use last character since that is the chain
        chain <- substrRight(pdb_num,1)
        selection_name_list[[chain_name]] <- chain

      } else if(is_pdb_match_vector == TRUE){
        chain <- strsplit(as.character(xlink.df.pdb[[paste0('pdb',pdb_num_count)]][row_num]),'_')[[1]][2]
        selection_name_list[[chain_name]] <- chain
      }else {

        #should it be in this part
        #should have an or statement or excessive?
        #do we need file_type --> unless the file is being created within the function
        #if(file_type_2d == 'single atom' && grepl('single_atom',pdb_num)){
        if(grepl('single_dot',pdb_num)){


          if(is.null(list_of_start_and_end_pdbs)){
            sd_pdb <- read.pdb2(paste0(pdb_num,'.pdb'))
            pro_pos <- sd_pdb$atom$resno
          } else {
            pro_pos <- list_of_start_and_end_pdbs[[paste(pdb_num,'.pdb',sep='')]][1]
          }

          #transform the protein position for the PyMol script so it is the first in the list

          selection_name_list[[pro_pos_name]] <- pro_pos


        } #series of if else statements for each of the other options or just else?

        chain <- 'A'
        selection_name_list[[chain_name]] <- chain
      }

      if(is.null(selection_name_list[[pro_pos_name]])){
        #if has not been changed in the single atom loop, will establish the original
        #protein position here
        selection_name_list[[pro_pos_name]] <- pro_pos_list[pdb_num_count]

      }


    } #end loop for(pdb_num in pdb_names)


    #getting the right selection names here
    #assigning the derived values to the selection names

    #pull out by the number?
    #then pull out by name of the variable ie 'chain'


    #create selection names

    completed_selection_names <- c()
    amino_acid_selection_names <- c()

    sors_pdb_check <- list(is_in_list=c(),
                           is_it_real=c())

    both_pdbs_pass_sors <- TRUE
    for(num in 1:2){
      selection_names_filtered <- names(selection_name_list)[grepl(as.character(num),names(selection_name_list))]
      sn_chain <- selection_names_filtered[grepl('chain',selection_names_filtered)]
      sn_chain <- selection_name_list[[sn_chain]]
      sn_pro_pos <- selection_names_filtered[grepl('pro_pos',selection_names_filtered)]
      sn_pro_pos <- selection_name_list[[sn_pro_pos]]
      sn_pdb <- selection_names_filtered[grepl('pdb',selection_names_filtered)]
      sn_pdb <- selection_name_list[[sn_pdb]]

      selection_name <- paste('/',sn_pdb,'//',sn_chain,'/',sn_pro_pos,'/CA',sep='')
      aa_select_name <- paste('/',sn_pdb,'//',sn_chain,'/',sn_pro_pos, sep='')

      completed_selection_names <- c(completed_selection_names,selection_name)
      amino_acid_selection_names <- c(amino_acid_selection_names,aa_select_name)

      #can do the process in here for determining whether or not the xl events should be
      #in the PyMOL output file
      #only one that does not meet the criteria should be excluded
      #start with true before loop and then if one does not meet critertia, then
      #make variable == FALSE

      if(!is.null(show_only_real_structures)){

        #check sn_pdb here
        #there should be a way to grepl a whole list
        for(sors in show_only_real_structures){

          #grepl check sn_pdb
          #another boolean?

          # sors_pdb_check <- list(is_in_list=c(),
          #                        is_it_real=c())

          if(grepl(sors,sn_pdb)){

            sors_pdb_check$is_in_list <- c(sors_pdb_check$is_in_list,TRUE)
            is_in_list <- TRUE

          } else {

            sors_pdb_check$is_in_list <- c(sors_pdb_check$is_in_list,FALSE)
            is_in_list <- FALSE
          }

          if(grepl('___',sn_pdb)){

            sors_pdb_check$is_it_real <- c(sors_pdb_check$is_it_real,TRUE)
            is_it_real <- TRUE

          } else {

            sors_pdb_check$is_it_real <- c(sors_pdb_check$is_it_real,FALSE)
            is_it_real <- FALSE
          }

          if(!((is_in_list && is_it_real) || !is_in_list)){
            both_pdbs_pass_sors <- FALSE
          }



          #similanteously do other check
          #or just store booleans in two lists
          #1st list: whether or not shows up in sors list
          #2nd list: is it a real pdb structure

          #either 1+2 == TRUE or 1==FALSE

        }


      } #end if(!is.null(show_only_real_structures))


    } #end for(num in 1:2)

    #enclose all of the following with pymol additions into the if/for/if loop

    #need to change these loops to accomodate for the fact that there are two selection names
    #here

    if(!is.null(show_only_real_structures)){

      #check list here
      #sors_pdb_check$is_in_list
      #sors_pdb_check$is_in_real

      #need to account for the fact that there are 2 structures being compared
      # pdbs_in_sors <- FALSE
      # for(sors in show_only_real_structures){
      #   if(grepl(sors,cs_name) && grepl('___',cs_name)){
      #     pdbs_in_sors <- TRUE
      #
      #   }
      #
      #   #can make two T/F vectors and then compare them to see if they're both T/F
      #
      # for(cs_name in completed_selection_names){
      #   for(sors in show_only_real_structures){
      #     if(grepl(sors,cs_name) && grepl('___',cs_name)){
      #       pdbs_in_sors <- TRUE
      #
      #     }
      #
      #   }
      # }


      #include the line of code if the protein is in list and is real
      #or or if it does not show up in the list at all
      if(both_pdbs_pass_sors){

        #do the code here
        freq_color <- as.character(xlink.df.pdb$freq_color[row_num])
        py_line <- paste('color ',freq_color,', ',amino_acid_selection_names[1],sep='')
        pymol_lines <- c(pymol_lines,py_line)
        py_line <- paste('color ',freq_color,', ',amino_acid_selection_names[2],sep='')
        pymol_lines <- c(pymol_lines,py_line)
        #accounting for peptides that show up in more than one XL event

        mega_distance_count <- mega_distance_count + 1
        py_line <- paste('distance ',completed_selection_names[1],', ',completed_selection_names[2], sep='')
        pymol_lines <- c(pymol_lines,py_line)
        #keep numerical list of distances?


        if(nchar(mega_distance_count) == 1){
          dist_name <- paste('dist0',mega_distance_count,sep='')
        } else {
          dist_name <- paste('dist',mega_distance_count,sep='')
        }

        freq_color <- as.character(xlink.df.pdb$freq_color[row_num])
        py_line <- paste('color ',freq_color,', ',dist_name,sep='')
        pymol_lines <- c(pymol_lines,py_line)

      } #end if((grepl(sors,pdb_name) && grepl('___',pdb_name)) || !grepl(sors,pdb_name))

    } else {#end if(!is.null(show_only_real_structures))

      freq_color <- as.character(xlink.df.pdb$freq_color[row_num])
      py_line <- paste('color ',freq_color,', ',amino_acid_selection_names[1],sep='')
      pymol_lines <- c(pymol_lines,py_line)
      py_line <- paste('color ',freq_color,', ',amino_acid_selection_names[2],sep='')
      pymol_lines <- c(pymol_lines,py_line)
      #accounting for peptides that show up in more than one XL event

      mega_distance_count <- mega_distance_count + 1
      py_line <- paste('distance ',completed_selection_names[1],', ',completed_selection_names[2], sep='')
      pymol_lines <- c(pymol_lines,py_line)
      #keep numerical list of distances?


      if(nchar(mega_distance_count) == 1){
        dist_name <- paste('dist0',mega_distance_count,sep='')
      } else {
        dist_name <- paste('dist',mega_distance_count,sep='')
      }

      freq_color <- as.character(xlink.df.pdb$freq_color[row_num])
      py_line <- paste('color ',freq_color,', ',dist_name,sep='')
      pymol_lines <- c(pymol_lines,py_line)

    }
    #color the amino acids
    # freq_color <- as.character(xlink.df.pdb$freq_color[row_num])
    # py_line <- paste('color ',freq_color,', ',amino_acid_selection_names[1],sep='')
    # pymol_lines <- c(pymol_lines,py_line)
    # py_line <- paste('color ',freq_color,', ',amino_acid_selection_names[2],sep='')
    # pymol_lines <- c(pymol_lines,py_line)
    # #accounting for peptides that show up in more than one XL event
    #
    # py_line <- paste('distance ',completed_selection_names[1],', ',completed_selection_names[2], sep='')
    # pymol_lines <- c(pymol_lines,py_line)
    # #keep numerical list of distances?
    #
    #
    # if(nchar(row_num) == 1){
    #   dist_name <- paste('dist0',row_num,sep='')
    # } else {
    #   dist_name <- paste('dist',row_num,sep='')
    # }
    #
    # freq_color <- as.character(xlink.df.pdb$freq_color[row_num])
    # py_line <- paste('color ',freq_color,', ',dist_name,sep='')
    # pymol_lines <- c(pymol_lines,py_line)
    #

    #can rename each of the crosslinking sites to xlink01 or something instead of dist01

    #make distances for each of them
    #add to list

    #draw lines between

    #use freq_color to color the distance line


  }

  py_line <- paste('hide labels',sep='')
  pymol_lines <- c(pymol_lines,py_line)


  if(write_file == TRUE){
    write(paste(pymol_lines,collapse = '\n'),file.name)
  }

  return(pymol_lines)

  #need to have 2 options: write and output a list that can be used
} #end function ppi.pymol()


