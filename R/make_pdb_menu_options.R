#----Make PDB menu options-----

#'Make PDB Menu Options
#'
#'This function makes the PDB menu options
#'
#'@param pdb_hit_table Hit table from blast.pdb()
#'@param pdb Boolean. If TRUE, uses pdb databases.
#'@export

make_pdb_menu_options <- function(pdb_hit_table, pdb = TRUE){

  pdb_menu_mega_list <- list()

  #add booleans to tell if menu should show the next 10 results
  #go back option that will take you back to previous menu?
  #should previous menu then be saved in another variable?
  #all menus be saved in a master list and users are simply scrolling
  #between mega list of menus?
  #unnamed indices of list so they can be retrieved by index + 1 and
  #index - 1




  num_rows <- nrow(pdb_hit_table)


  for(first_row in seq(1,num_rows,10)){

    first_row <- 1
    if(first_row+10 <= num_rows){

      rows_range <- first_row:(first_row+9)
      sub_hit_table <- pdb_hit_table[rows_range,]



      #do something similar to num_row_range
      #add to pdb_menu_mega_list
      #iterating through mega list should be in the main script

    } else {
      #if the +10 will cause it go outside of the indeces of the hit table

    }
  }
  #need to make condition for last row which may be less than 10


  more_rows_option <- FALSE
  num_row_range <- 1:num_rows

  #add number of hits as part of results
  #one output

  # if(num_rows > 10){
  #   num_row_range <- 1:10
  #   more_rows_option <- TRUE
  # } else {
  #   num_row_range <- 1:num_rows
  #   more_rows_option <- FALSE
  # }

  blasted_pdb_menu_options <- c()


  for(num in num_row_range){

    if(pdb == TRUE){
      hit_id <- pdb_hit_table$pdb.id[num]

    } else {
      hit_id <- strsplit(pdb_hit_table$pdb.id[num],'_NA')[[1]]
    }


    menu_option <- paste(hit_id,
                         "|",
                         pdb_hit_table$alignmentlength[num],
                         "|",
                         pdb_hit_table$positives[num])

    blasted_pdb_menu_options <- c(blasted_pdb_menu_options,menu_option)
  }

  #accounting for row number for more options?
  #make menu options part of a list where there are booleans that will
  #trigger certain options

  if(more_rows_option == TRUE){
    more_rows_choice <- "See more PDB matches"
    blasted_pdb_menu_options <- c(blasted_pdb_menu_options,more_rows_choice)

  }

  #c() with strings? pdb pdb pdb pdb more_rows make_2d
  if(pdb == TRUE){
    twod_option <- "Do not use any PDB structure -- generate 2D structure from sequence"
    blasted_pdb_menu_options <- c(blasted_pdb_menu_options,twod_option)
  }

  #manually input PDB ID and chain?



  return(blasted_pdb_menu_options)

}
