#'BLAST sequence and go through menu loops
#'
#'This function will BLAST the given sequence and give the user an interactive menu to select the best option.
#'@param sequence Sequence as a vector or string to be used to BLAST against
#'@param database Menu selected for blast.pdb(). Options: "pdb","nr", or "swissprot". Defaults to "pdb".
#'@param time_out Integer indicating time in seconds the blast.pdb() will search before time out occurs. Defaults to
#'@export

blast.menu <- function(sequence, database = 'pdb', time_out = NULL){

  fasta_sequence_vector <- sequence
  ff_sequence <- toupper(paste(fasta_sequence_vector,collapse=''))

  # if(pdb == TRUE){
  #   selected_database <- 'pdb'
  # } else {
  #   selected_database <- 'nr'
  # }

  ff_blast <- tryCatch(blast.pdb(ff_sequence, database = database, time.out = time_out),
                       error = function(err){
                         cat('No results for found for this sequence\n')
                         return(NULL)
                       })

  if(is.null(ff_blast)){
    return(NULL)
  }

  #would need to update the blast.menuOptions so that it was also not a true/false
  #if including swissprot as an option
  ff_menu_options <- blast.menuOptions(pdb_hit_table = ff_blast$hit.tbl, database = database)
  menu_index <- 1

  pdb_selected <- FALSE
  go_to_more_options <- FALSE
  use_searched_pdb <- FALSE

  while(pdb_selected == FALSE){

    menu_displayed <- ff_menu_options[[menu_index]]
    blasted_title <- paste("Select best matching alignment for your needs using your keyboard.",
                           "\n\nTop PDB hits (page ",as.character(menu_index),"/",
                           as.character(length(ff_menu_options)),")",":\n",
                           "pdb.id | identity | alignment length | % positives",sep='')
    menu_selection <- menu(menu_displayed, title = blasted_title)


    if(menu_index == 1){ #first menu
      #if there is only one menu list
      if(menu_index == length(ff_menu_options)){

        pdb_selected <- TRUE
        if(menu_selection == length(menu_displayed)){
          pdb_selected <- FALSE
          go_to_more_options <- TRUE
        }

      } else { #else if there is more than one

        if(menu_selection == (length(menu_displayed))){
          #pdb_selected <- TRUE
          go_to_more_options <- TRUE
        } else if(menu_selection == (length(menu_displayed)-1)){
          menu_index <- menu_index + 1
          #next
        } else {
          pdb_selected <- TRUE
        }

      }


    } else if(menu_index < length(ff_menu_options)){
      #if index is less than the length of the menu options
      #but not 1

      if(menu_selection == length(menu_displayed)){
        go_to_more_options <- TRUE
      } else if(menu_selection == (length(menu_displayed) - 1)){
        menu_index <- menu_index - 1
        next
      } else if(menu_selection == length(menu_displayed) - 2){
        menu_index <- menu_index + 1
        next
      } else {
        pdb_selected <- TRUE
      }

    } else if ((menu_index == length(ff_menu_options)) && (length(ff_menu_options) != 1)){


      if(menu_selection == length(menu_displayed)){
        go_to_more_options <- TRUE
      } else if(menu_selection == (length(menu_displayed) - 1)){
        menu_index <- menu_index - 1
        next
      } else if(menu_selection == length(menu_displayed) - 2){
        #go back to first page of hits
        menu_index <- 1
        next
      } else {
        pdb_selected <- TRUE
      }

    }

    #if user chooses to do more options

    while(go_to_more_options == TRUE){

      pdb_readline_search <- FALSE
      #move to second while loop
      more_options_list <- c("Search the PDB hits by PDB ID",
                             "Do not use any PDB structure -- generate 2D structure from sequence",
                             "Go back to last page")

      #establish menu with more_options_list
      mo_title <- "More Options"
      mo_menu_selection <- menu(more_options_list,title=mo_title)

      if(mo_menu_selection == 1){ #search by PDB ID

        pdb_readline_search <- TRUE

        while(pdb_readline_search == TRUE){
          pdb_readline_input <- readline(prompt = "Enter a PDB ID: ")
          pdb_readline_hit_table <- ff_blast$hit.tbl[grepl(pdb_readline_input,ff_blast$hit.tbl$subjectids),]
          ht_row_text_list <- c()
          for(ht_row_num in 1:nrow(pdb_readline_hit_table)){
            ht_row <- pdb_readline_hit_table[ht_row_num,]
            #get relevant information to output to table

            ht_row_text <- paste(ht_row$subjectids," | ",ht_row$positives,
                                 sep='')

            ht_row_text_list <- c(ht_row_text_list,ht_row_text)
          }

          more_options_lines <- c("Search again","Go back to more options")

          ht_row_text_list <- c(ht_row_text_list,more_options_lines)

          #change title based on if nrow(pdb_readline_hit_table) is 1 or higher
          #if == 1, "There was 1 match"
          ht_menu_title <- paste("There were ",as.character(nrow(pdb_readline_hit_table)),
                                 " matches to your search.\nPDB ID | positives",
                                 sep='')
          ht_menu_selection <- menu(ht_row_text_list,title=ht_menu_title)

          if(ht_menu_selection == length(ht_row_text_list)){ #go back to more options
            #go back to more options
            #change menu_displayed to more_options menu
            #if this menu is within a while loop, exit out of the loop
            pdb_readline_search <- FALSE

          } else if(ht_menu_selection == (length(ht_row_text_list) - 1)){
            #search again
            #send back to beginning of while statement
            next

          } else {

            searched_pdb_ht <- pdb_readline_hit_table[ht_menu_selection,]


            use_searched_pdb <- TRUE

            #exits out of all the while loops
            pdb_readline_search <- FALSE
            go_to_more_options <- FALSE
            pdb_selected <- TRUE

          }

        }

      } else if(mo_menu_selection == 2){ #generate 2D structure

        print('Function not currently available')
        go_to_more_options <- FALSE
        #need to update to include this functionality

      } else if(mo_menu_selection == 3){ #go back to previous page
        #next?
        #make while statement false to exit the while loop
        #if at end of while statement, next is not necessary
        go_to_more_options <- FALSE
      }


    } #end of go_to_more_options while loop


  } #end of pdb_selected while loop


  if(use_searched_pdb == FALSE){
    pdb_hit_index <- ((menu_index-1) * 10) + menu_selection
    ff_pdb_info <- (get_pdb_info(menu_selection = pdb_hit_index,
                                 pdb_hit_table = ff_blast$hit.tbl))
  } else {
    ff_pdb_info <- get_pdb_info(1,searched_pdb_ht)
  }


  return(ff_pdb_info)
}
