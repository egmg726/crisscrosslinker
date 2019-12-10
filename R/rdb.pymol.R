#'Write PyMol file for RBDmap Data
#'
#'This function makes a PyMOL file (.pml) that can be used to visualize data obtained from rdb.getBindingSeq()
#'
#'@param bs_output Output from rbd.getBindingSeq()
#'@param color_by The variable that will be used to determine color. Must a column in bs_output or use "freq" to color by frequency of amino acid position in binding sequences. Defaults to "binding_sequence".
#'@param colors Colors to be used for color.pymol() function.
#'@param file.name If write.file = T, name of the output file. Defaults to 'rbd_pymol.pml'
#'@param write.file Boolean. If TRUE, will output file with name of file.name. If FALSE, will return file as list.
#'@param experiment.dir Directory of experiment, where PDB files are. If left as NULL, will use current directory.
#'@param gray0 Boolean. If TRUE, will shift color scheme for frequency analysis and label 0s as "gray" in PyMOL.
#'@param assembly Assembly number for PyMOL output that corresponds to RCSB database. Works for PyMOL 1.8 and above.
#'@export

rbd.pymol <- function(bs_output, color_by = 'binding_sequence',
                      colors = NULL, file.name = 'rbd_pymol.pml',
                      write.file = TRUE, experiment.dir = NULL,
                      gray0 = FALSE, heatmap = TRUE, assembly = 0,
                      fetch = TRUE){

  #what is color_by is frequency?

  if(is.null(experiment.dir)){
    experiment.dir <- getwd()
  }

  #freq_vector
  pymol_lines <- c()

  if(is.null(colnames(bs_output))){
    bs_output <- data.frame(bs_output)
  }

  bs_output <- bs_output[bs_output$db == 'PDB',]
  #load the PDB files up here?

  split_pdbs <- unlist(strsplit(as.character(bs_output$db_id),'_'))
  unique_pdbs <- unique(split_pdbs[nchar(split_pdbs) == 4])

  py_line <- paste0('set assembly, ', assembly)
  pymol_lines <- c(pymol_lines,py_line)

  for(u_pdb in unique_pdbs){
    #py_line <- paste('fetch',u_pdb)
    if(fetch == TRUE){
      py_line <- paste0('fetch ',u_pdb,', async=0')
    } else {
      py_line <- paste0('load ',experiment.dir,'/',u_pdb,'.pdb')
    }

    pymol_lines <- c(pymol_lines,py_line)
  }



  if(color_by == 'freq'){

    #rename the names give to the bs_freqVector --> name_by?

    #bs_output$source_sequence <- paste0(resno_and_resid$resid,collapse='')
    # bs_output <- rbd.getBindingSeq(ms_sequence = 'ERTEILNQEWK',
    #                   protease = 'ArgC',
    #                   sequence_for_alignment =paste0(resno_and_resid$resid,collapse=''),
    #                   protein_name = 'EZH2',
    #                   database_name = 'PDB',
    #                   database_id = '6C23_C')
    # bs_output <- data.frame(bs_output)


    #need to add in colors for freqVector
    bs_freqVector <- rbd.freqVector(bs_output = bs_output, name_by = 'db_id',heatmap = FALSE, db_selection = 'PDB')

    #return(bs_freqVector)
    bsfv_vars <- min(unlist(bs_freqVector)):max(unlist(bs_freqVector))
    #explicitly make the bar??

    #return(bsfv_vars)
    #return(list(vars=bsfv_vars,colors = colors,gray0 = gray0))


    if(length(colors) != length(bsfv_vars)){
      colors <- colorRampPalette(colors = colors)(length(bsfv_vars))
    }
    pymol_cc <- color.pymol(bsfv_vars, colors = colors, png.name = 'freqvector_legend_test.svg', gray0 = gray0, print.legend = write.file)


    #return(pymol_cc)

    #return(pymol_cc)
    bs_freqVector2 <- rbd.freqVector(bs_output = bs_output, name_by = 'db_id',heatmap = heatmap, db_selection = 'PDB', colors = pymol_cc$hexcodes)

    #return(bs_freqVector2)

    #pymol_lines <- c()
    pymol_lines <- c(pymol_lines,pymol_cc$set_colors)

    #can alter this so the first color set that is equal to 0 will be the 0
    if(gray0 == TRUE){
      pymol_lines <- c(pymol_lines,'color gray')
    } else {#if gray0 == FALSE
      if(length(paste0('color ',pymol_cc$color_names[pymol_cc$vars == 0])) > 0){
        pymol_lines <- c(pymol_lines,paste0('color ',pymol_cc$color_names[pymol_cc$vars == 0]))
      } else {
        cat('No var with value 0 detected --> coloring gray')
        pymol_lines <- c(pymol_lines,'color gray')
      }

    } #end else to if(gray0 == TRUE){

    #should have this changed
    pymol_lines <- c(pymol_lines,'show surface')
    #pymol_cc$color_names[pymol_cc$vars == 0]


    #return(pymol_lines)

    for(pdb_name in names(bs_freqVector)){
      freq_vector <- bs_freqVector[[pdb_name]]

      #freq_vector

      pdb_id <- strsplit(pdb_name,'_')[[1]][1]
      chain <- strsplit(pdb_name,'_')[[1]][2]

      pdb_read <- check_download_read_pdb(pdb_id) #add ID to list of IDs that need to be loaded?
      #or can just remove repeating lines from final PyMOL file

      resno_and_resid <- quick_resno_and_resid(pdb_read = pdb_read, chain = chain)

      if(length(resno_and_resid$resno) == length(freq_vector)){
        #the lengths match --> probably the right sequence

        #get the resno for each
        #make sure that the color matches to the pymol_cc
        #
        #what about making grays 0s?

        #for loop to go through vars and make them different colors?
        #then go through the entire list?
        #pymol_cc$vars

        for(index_num in 1:length(freq_vector)){
          aa_freq <- freq_vector[index_num]
          match_index <- match(aa_freq,pymol_cc$vars)
          #pymol_cc$vars[match_index]
          color_name <- pymol_cc$color_names[match_index]
          resno_index <- resno_and_resid$resno[index_num]

          py_line <- paste0('color ',color_name,', /',pdb_id,'//',chain,'/',as.character(resno_index))
          pymol_lines <- c(pymol_lines,py_line)

          #construct the string here and make into a line

          #use the match index to get the right colors


          #match the aa_freq

          #use the index_num to get the right numbering from resno

        } #end for(index_num in 1:length(freq_vector))

      } else { #end if(length(resno_and_resid$resno) == length(freq_vector))

        warning("The lengths don't match --> try a different sequence")
        #add specifics about what sequence it is?
        return(list(resno=resno_and_resid$resno,freq_vector=freq_vector))
      } #end else to if(length(resno_and_resid$resno) == length(freq_vector))




    }  #end  for(pdb_name in names(bs_freqVector))


    #keeping track of the top frequency number


  } else { #end if(color_by == 'freq')
    #do what's below
    #assume that it's the name of a column in bs_output

    pymol_lines <- c(pymol_lines,'show surface')

    # if((is.numeric(bs_output[[color_by]]) == FALSE) && (gray0 == FALSE)){
    #   cat('color_by not numeric --> gray0 is now true\n')
    #   gray0 <- TRUE
    # }

    pymol_cc <- color.pymol(sort(unique(bs_output[[color_by]])), colors = colors, gray0 = gray0)
    pymol_lines <- c(pymol_lines,pymol_cc$set_colors)

    if(gray0 == TRUE){
      pymol_lines <- c(pymol_lines,'color gray')
    } else {#if gray0 == FALSE
      if(length(pymol_cc$color_names[pymol_cc$vars == 0]) > 0){
        pymol_lines <- c(pymol_lines,paste0('color ',pymol_cc$color_names[pymol_cc$vars == 0]))
      } else {
        cat('No var with value 0 detected --> coloring gray')
        pymol_lines <- c(pymol_lines,'color gray')
      }

    } #end else to if(gray0 == TRUE){

    bs_output$color <- rep(NA,nrow(bs_output))
    levels(bs_output$color) <- c(levels(bs_output$color),pymol_cc$color_names)

    for(index_num in 1:length(pymol_cc$vars)){
      var_name <- pymol_cc$vars[index_num]
      bs_output[bs_output[[color_by]] == var_name, 'color'] <- pymol_cc$color_names[index_num]
    } #end for(index_num in 1:length(pymol_cc$vars))



    #need to then go through each of the rows of the bs_output

    for(row_num in 1:nrow(bs_output)){
      bso_row <- bs_output[row_num,]
      #write each of the lines and add to pymol_lines

      pdb_id <- strsplit(as.character(bso_row$db_id),'_')[[1]][1]
      chain <- strsplit(as.character(bso_row$db_id),'_')[[1]][2]

      py_line <- paste0('color ',bso_row$color,', /',pdb_id,'//',chain,'/',bso_row$binding_site_start,'-',bso_row$binding_site_end)
      pymol_lines <- c(pymol_lines,py_line)

    } #end for(row_num in 1:nrow(bs_output))


  } #end else { #end if(color_by == 'freq')


  if(write.file == TRUE){
    write(paste(unique(pymol_lines),collapse = '\n'),file.name)
  } else {
    return(pymol_lines)
  }


} #end function rbd.pymol

