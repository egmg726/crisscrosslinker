#----Make PDB menu options 2-----

#'Make PDB Menu Options
#'
#'This function creates a menu from a pdb_hit_table that is retreived by bio3d's blast.pdb function
#'
#'@param pdb_hit_table Hit table retrieved from bio3d::blast.pdb()
#'@param pdb If it searches BLAST PDB. Defaults to true
#'@export

#need to make it so that the FASTA file name is included in the menu
blast.menuOptions <- function(pdb_hit_table, database = 'pdb'){

  pdb_menu_mega_list <- list()
  #booleans as part of list to indicate whether or not the list needs
  #a statement saying to go back and forth through the menu options

  pdb <- TRUE
  #if(datbase == 'pdb')

  num_rows <- nrow(pdb_hit_table)
  #if num_rows <= 10 should take away option to have line in list that
  #will give you the option to show next
  #option in case list is less than

  more_than_ten_rows <- TRUE
  if(num_rows <= 10){
    more_than_ten_rows <- FALSE
  }
  #make menu for if there is less than 10 rows


  row_num <- 0
  for(first_row in seq(1,num_rows,10)){
    row_num <- row_num + 1

    if(first_row+10 <= num_rows){

      rows_range <- first_row:(first_row+9)
      sub_hit_table <- pdb_hit_table[rows_range,]

      if(pdb == TRUE){
        hit_id <- sub_hit_table$subjectids

      } else {
        hit_id <- strsplit(sub_hit_table$subject_ids[row_num],'_NA')[[1]]
      }

      identity_ht <- sub_hit_table$identity
      alignment_length <- sub_hit_table$alignmentlength
      positives <- sub_hit_table$positives

      menu_options_list <- c()
      for(num in 1:length(hit_id)){

        menu_option <- paste(hit_id[num],
                             "|",
                             identity_ht[num],
                             "|",
                             alignment_length[num],
                             "|",
                             positives[num])

        menu_options_list <- c(menu_options_list,menu_option)
      }

      #can customize the next page option for PDB or NCBI hits
      if(more_than_ten_rows == TRUE){
        next_page_option <- "Go to next page of hits"
        menu_options_list <- c(menu_options_list,next_page_option)
      }

      if(first_row != 1){
        previous_page_option <- "Go back to previous page of hits"
        menu_options_list <- c(menu_options_list,previous_page_option)
      }

      #should there be a function to go back one page?
      #if first_row != 1 --> add option to go to previous page of hits

      if(pdb == TRUE){
        make_2d_option <- "More options"
        menu_options_list <- c(menu_options_list,make_2d_option)
      }

      #first index in mega list is list of booleans that correspond to
      #?

      #do something similar to num_row_range
      #add to pdb_menu_mega_list
      #iterating through mega list should be in the main script

      pdb_menu_mega_list[[row_num]] <- menu_options_list

    } else {
      #if the +10 will cause it go outside of the indeces of the hit table

      #need for loop that will be contained within
      #first_row and length(pdb_hit_table)
      sub_hit_table <- pdb_hit_table[first_row:num_rows,]

      if(pdb == TRUE){
        hit_id <- sub_hit_table$subjectids

      } else {
        hit_id <- strsplit(sub_hit_table$pdb.id[num],'_NA')[[1]]
      }

      alignment_length <- sub_hit_table$alignmentlength
      positives <- sub_hit_table$positives

      menu_options_list <- c()
      for(num in 1:length(hit_id)){

        menu_option <- paste(hit_id[num],
                             "|",
                             alignment_length[num],
                             "|",
                             positives[num])

        menu_options_list <- c(menu_options_list,menu_option)
      }

      #option to go to the first
      if(more_than_ten_rows == TRUE){
        next_page_option <- "Go back to first page of hits"
        menu_options_list <- c(menu_options_list,next_page_option)
        previous_page_option <- "Go back to previous page of hits"
        menu_options_list <- c(menu_options_list,previous_page_option)
      }



      if(pdb == TRUE){
        make_2d_option <- "More options"
        menu_options_list <- c(menu_options_list,make_2d_option)
      }

      pdb_menu_mega_list[[row_num]] <- menu_options_list

    }
  }

  return(pdb_menu_mega_list)

} #end function
