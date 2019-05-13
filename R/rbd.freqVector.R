#'Get Frequency Vector for RBDmap Data
#'
#'This function gets the frequency vector of the binding sequences detected
#'
#'@param bs_output data.frame() output from rbd.getBSfromIET() or rbd.getBSfromDF()
#'@param name_by Name of column in bs_output that will correspond to where the binding sites will be measured
#'@param heatmap If TRUE, a heatmap will be produced. If FALSE, the data.frame used for the heatmap will be returned
#'@param db_selection Database to be used for the heatmap (e.g. FASTA, PDB, or UniProt)
#'@param colors Colors to be used for the heatmap, if TRUE. Can be a list of hexcodes/standard R colors or RColorBrewer palette
#'@param save_plot If TRUE, will save the plot to your working directory
#'@author Emma Gail
#'
#'@export

rbd.freqVector <- function(bs_output, name_by = 'pro_name', heatmap = TRUE, db_selection = NULL, colors = c('#d0d0d0','#2f5ac6','#e50000'),save_plot=TRUE){

  #will only work if the data is all Uniprot or FASTA
  #will be a list
  #if((length(levels(bs_output$db))) == 1 && ((levels(bs_output$db) == 'Uniprot') || (levels(bs_output$db) == 'FASTA'))


  #could theoretically do just PDB as well as long as length(levels()) == 1
  if(is.null(db_selection)){

    if(length(levels(bs_output$db)) > 1){
      #if(!((levels(bs_output$db) == 'Uniprot') || (levels(bs_output$db) == 'FASTA'))){

      #if the levels do not match either of these
      #function cannot proceed
      warning('All rows in database must belong to only 1 value')


      menu_selection <- menu(choices = levels(bs_output$db), title = 'Choose a database type for your numbering')

      db_selection <- levels(bs_output$db)[menu_selection] #new selection --> new level

      #can also put up menu option that will ask if they want to filter by just one level
      #with the levels as the options for the menu selection
      #return(NULL)

    } else {

      if(is.null(levels(bs_output$db))){
        db_selection <- unique(bs_output$db)
      } else {
        db_selection <- levels(bs_output$db)
      }


    }


  } #end if(!is.null(db_selection))


  #return(db_selection)

  bs_output <- bs_output[bs_output$db == db_selection,]

  #return(bs_output)

  list_of_aa_vectors <- list()

  aa_mega_df <- data.frame()


  list_of_protein_names <- levels(factor(bs_output[[name_by]]))
  if(is.null(list_of_protein_names)){
    list_of_protein_names <- unique(bs_output[[name_by]])
  }

  #return(list_of_protein_names)

  for(protein_level in list_of_protein_names){

    bs_o_df_rows <- bs_output[bs_output[[name_by]] == protein_level,]
    rownames(bs_o_df_rows) <- 1:nrow(bs_o_df_rows)



    #Get the sequence based on the
    # if(db_selection == 'PDB'){
    #   #check_download_read_pdb() #get ID from db_id
    #   #get the chain from db_id
    #   #get the sequence
    #   #
    #
    # } else if(db_selection == 'UniProt'){
    #   #get the sequence from Uniprot
    #
    #
    #
    # } else if(db_selection == 'FASTA'){
    #   #get the FASTA sequence that corresponds to the
    #
    # } else {
    #   #unknown database selection
    #   warning('Unknown database selection')
    # }

    selection_seq <- unique(as.character(bs_o_df_rows$source_sequence)) #the sequence to be used for the vector

    #return(selection_seq)

    #check if protein_to_uniprot_id is NULL
    #will have to find the Uniprot sequence with blast or the PDB ID if available
    #the Uniprot ID might already be in there?

    #Uniprot not needed here then
    #uniprot_id <- protein_to_uniprot_id[protein_to_uniprot_id$protein_id == protein_level,'uniprot_id']
    #uniprot_fasta <- seqinr::read.fasta(paste(fasta_file_directory,'/',uniprot_id,'.fasta',sep=''))

    #uniprot_fasta_seq <- toupper(uniprot_fasta[[names(uniprot_fasta)]])
    #make vector of 0s?
    #use the indices as a count
    aa_count_vector <- NULL
    aa_count_vector <- rep(0,nchar(selection_seq))

    #return(aa_count_vector)

    #aa_count_vector[1:6] <- aa_count_vector[1:6] + 1

    for(row_num in 1:nrow(bs_o_df_rows)){

      bs_o_df_row <- bs_o_df_rows[row_num,]

      #bs_start <- as.numeric(bs_o_df_row$binding_site_start)
      #bs_end <- as.numeric(bs_o_df_row$binding_site_end)
      bs_seq <- as.character(bs_o_df_row$binding_sequence)

      bse <- str_locate_all(selection_seq,bs_seq)[[1]]
      #need to check the length of this

      bs_start <- bse[1] #make it so it's relative to the original? #get the numeric vector?
      bs_end <- bse[2]

      #see about the original sequence if

      if(is.na(bs_start) || is.na(bs_end)){
        next
      }
      aa_count_vector[bs_start:bs_end] <- aa_count_vector[bs_start:bs_end] + 1

    }

    #once the vector is done --> need to add to list

    list_of_aa_vectors[[protein_level]] <- aa_count_vector

    aa_df <- data.frame(xx=1:length(list_of_aa_vectors[[protein_level]]),
                        Count=list_of_aa_vectors[[protein_level]],
                        zz=rep(protein_level,length(list_of_aa_vectors[[protein_level]])))

    aa_mega_df <- rbind(aa_mega_df,data.frame(aa_df))



  }

  #return(aa_mega_df)

  #make this heatmap for all of the input and eluate as well
  #4 heatmaps --> + and - UV  for both input and eluate

  if(heatmap == TRUE){

    #use colors variable


    #color_palette1 <- colorRampPalette(c('#d0d0d0','#2f5ac6','#e50000'),3)
    #color_palette1 <- colorRampPalette(c('#d0d0d0','#2f5ac6','#e50000'))

    #hm.palette <- colorRampPalette(rev(brewer.pal(4, 'Spectral')), space='Lab')
    #hm.palette <- colorRampPalette(c('#d0d0d0','#2f5ac6','#e50000'), space='Lab')

    viridis_colors <- c("magma","inferno","plasma","viridis","cividis")

    if((colors %in% rownames(brewer.pal.info)) || (colors %in% viridis_colors)){
      freq_vars <- sort(unique(unlist(bs_freqVector)))
      freq_colors <- color.pymol(freq_vars,colors=colors)
      hm.palette <- colorRampPalette(colors = freq_colors$hexcodes)
    } else {
      hm.palette <- colorRampPalette(colors = colors)
    }






    #startsWith(colors,'#')

    #color.pymol(,colors = colors)

    #color.pymol()

    #colors <- c('#d0d0d0','#2f5ac6','#e50000')


    #color.pymol(1:leng)


    #colorRampPalette('Blues',space='Lab')


    #hm.palette(10)

    #unique(aa_mega_df$Count)

    aa_plot <- ggplot(aa_mega_df, aes(xx,zz)) + geom_tile(aes(fill=Count)) + xlab('AA Position in Protein Sequence') + ylab('Protein Name') + coord_fixed(ratio = 50) +
      theme(
        panel.background = element_rect(fill = "#ffffff",
                                        colour = "#000000",
                                        size = 0.5, linetype = "solid"),
        legend.key = element_rect(fill = "#eaf2ff")
      ) + scale_fill_gradientn(colours = hm.palette(100))

    print(aa_plot)

    if(save_plot == TRUE){
      ggsave('uvxlms_heatmap_all_proteins.svg')
    }




  } #end if(heatmap == TRUE)




  return(list_of_aa_vectors)

} #end rbd.freq_vector function

