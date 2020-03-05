#' Get Binding Sequence from RBDmap Data - Zhang et. al 2019 Method
#'
#' This function gets the binding sequence from RBDmap data
#'
#' @param ms_sequence The input sequence from MS, should be loaded as a string
#' @param protease The protease used for the experiment, either 'ArgC' or 'LysC'
#' @param sequence_for_alignment The sequence to find the binding sequence in. The ms_sequence must be within the sequence_for_alignment or NULL will be returned.
#' @param protein_name Name of the protein to which the sequence_for_alignment belongs. Defaults to NULL
#' @param database_name Name of the database used: 'PDB','UniProt' or 'FASTA' (for the input FASTA sequence). Defaults to NULL
#' @param database_id Identifier for the database. Defaults to NULL
#' @param cleave_offset Cleave offset for the sequence. Defaults to 0
#' @export


rbd.getBindingSeq <- function(ms_sequence,protease,sequence_for_alignment,
                              protein_name = NULL,database_name = NULL,
                              database_id = NULL, cleave_offset=0,last_aa = NULL,
                              include_ambiguous = FALSE, proteolytic_fragments = FALSE){

  if(length(str_locate_all(sequence_for_alignment,ms_sequence)[[1]]) == 0){
    warning("Input sequence not found in the sequence for alignment")
    return(NULL)

  }


  if(is.null(last_aa)){
    last_aa <- substrRight(ms_sequence,1)

  }

  aa_after_peptide_seq <- last_aa

  start_and_end_patterns <- str_locate_all(sequence_for_alignment,ms_sequence)[[1]]

  start_pattern <- start_and_end_patterns[1]
  end_pattern <- start_and_end_patterns[2]

  sequence_for_alignment_split <- strsplit(sequence_for_alignment,'')[[1]]

  go_forwards <- NULL

  #turn this into a function?
  if(grepl('ArgC',protease)){
    end_binding_aa <- 'R'

    if(aa_after_peptide_seq == 'R'){
      #R at the end --> go backwards
      aa_increment <- -1
      current_aa_start <- start_pattern
      go_forwards <- FALSE

    } else if(aa_after_peptide_seq == 'K'){
      #K at the end + ArgC --> go forwards
      aa_increment <- 1
      current_aa_start <- end_pattern - 1
      #end_binding_aa <- 'R'
      go_forwards <- TRUE

    }

  }

  if(grepl('LysC',protease)){
    end_binding_aa <- 'K'

    if(aa_after_peptide_seq == 'R'){
      #R at the end --> go forwards
      aa_increment <- 1
      current_aa_start <- end_pattern - 1
      #end_binding_aa <- 'R'
      go_forwards <- TRUE

    } else if(aa_after_peptide_seq == 'K'){
      #K at the end + LysC --> go backwards
      aa_increment <- -1
      current_aa_start <- start_pattern
      go_forwards <- FALSE

    }


  }

  #need to just add a statement if current_aa == 1 or the last index, need to end the while loop


  #end_binding_aa <- aa_after_peptide_seq
  current_aa_index <- current_aa_start + aa_increment
  binding_sequence <- c()
  current_aa <- sequence_for_alignment_split[current_aa_index]


  while((current_aa != end_binding_aa) && (current_aa_index <= length(sequence_for_alignment_split)) && (current_aa_index > 0)){
    if(aa_increment < 0){
      binding_sequence <- c(current_aa,binding_sequence)
      current_aa_index <- current_aa_index + aa_increment
      current_aa <- sequence_for_alignment_split[current_aa_index]
    } else {
      current_aa_index <- current_aa_index + aa_increment
      current_aa <- sequence_for_alignment_split[current_aa_index]
      binding_sequence <- c(binding_sequence,current_aa)
    }


    #search for another nearby peptide

    #needs to see if there is another

    #if it's 0 --> should return an empty string



    #cleave_offset <- 0

    #may need an if/else statement

    #move current_aa increment down here to get rid of 'R' at end
  }

  #does this need to be within the while() loop so that it will continue to add amino acids
  #to binding sequence if needed

  #or need to make 2 while loops --> then can make whole thing into a function
  #just 2nd as function would make less sense than enclosing both in a function



  if(cleave_offset > 0){

    #set the variable that keeps track of the most recent R/K

    #necesary if the same?
    most_recent_cleave_site <- current_aa_index #current_aa will need to be reset if it goes
    #outside of the loop

    #current_aa_index+1
    #does this need to +/- 1 for the current_aa_index

    end_of_cleave_peptide <- (aa_increment*cleave_offset)+current_aa_index

    cleaved_peptide_indeces <- sort(c(end_of_cleave_peptide,current_aa_index))

    cleave_peptide <- sequence_for_alignment_split[cleaved_peptide_indeces[1]:cleaved_peptide_indeces[2]]

    #check if there is a K or R within the cleave peptide

    if(end_binding_aa %in% cleave_peptide){

      #if the end_binding peptide is in the

      #will need to append it to the current binding sequence and make the new
      #ending the site of the furthest end binding site?
      #will that need to be checked again?

      #can also use the str_locate within the segment before/after the binding sequence
    }

    #fasta_seq <- paste(toupper(fasta_file$`EZH2_Q15910-2`),collapse='')
    #list_of_cleaves <- str_locate_all(fasta_seq,'K')[[1]][,1]
    #find the next one depending on the start/end of the binding sequence
    #can also expand this for originally finding the binding sequence

    #start_pattern <- 40
    #end_pattern <- 80

    #list_of_cleaves_sub <- list_of_cleaves[start_pattern > list_of_cleaves]
    #tail(list_of_cleaves_sub,n=1) #last integer
    #tail(list_of_cleaves_sub,n=2)[1] #second to last integer
    #will have to deal with errors as well --> is.na() probably



    #will have to see if the difference to the next cleave site is
    #more than the cleave_offset


    #would have to account for the beginning/end if it goes out of sequence

    #sort() to make sure they're in the right order


    #get the peptide that's before/after the binding sequence using the indeces

    #will need to keep track of the K or R that's closest to the start/end
    #of the binding site in a variable

    #should put the



    #may need to order



    #check if end_binding_aa is %in% the list of amino acids

    #if TRUE --> need to add the whole peptide to the binding sequence and then go to the very next
    #peptide and continue on with the while loop


    #get a range of the next few
  } #end if(cleave_offset > 0)

  #if 0 then skip
  #will need to be + or - based on which direction it is
  #multiply cleave_offset by the aa_increment and then add

  binding_sequence <- paste(binding_sequence,collapse='')
  #can just do the str_locate_all (?) to get the start and end?

  resno_ms2_peptide <- sequence_for_alignment_split[start_pattern:end_pattern]
  resno_ms2_start_pattern <- start_pattern
  resno_ms2_end_pattern <- end_pattern

  #need to change these for ArgC increments
  #only works for LysC

  if(go_forwards == TRUE){
    binding_site_start_in_pdb <- end_pattern + 1
    binding_site_end_in_pdb <- current_aa_index
    resno_binding_peptide <- sequence_for_alignment_split[binding_site_start_in_pdb:binding_site_end_in_pdb]


  } else if(go_forwards == FALSE){
    binding_site_start_in_pdb <- current_aa_index + 1
    binding_site_end_in_pdb <- start_pattern - 1
    resno_binding_peptide <- sequence_for_alignment_split[binding_site_start_in_pdb:binding_site_end_in_pdb]
  } else {

    binding_site_start_in_pdb <- NA
    binding_site_end_in_pdb <- NA
    resno_binding_peptide <- NA

    #none selected --> have some kind of error --> assign NA to the binding sequence(?)

  }


  #binding_site_start_in_pdb <- min(resno_binding_peptide)
  #binding_site_end_in_pdb <- max(resno_binding_peptide)

  #return(binding_sequence)

  if(include_ambiguous == TRUE){
    if(nchar(as.character(binding_sequence)) == 0){
      binding_sequence <- ms_sequence
      binding_site_start_in_pdb <- resno_ms2_start_pattern
      binding_site_end_in_pdb <- resno_ms2_end_pattern
    }
  }


  if(proteolytic_fragments == TRUE){ #should the include_ambiguous be included in proteloytic fragments anyway??

    if(resno_ms2_start_pattern > binding_site_start_in_pdb){
      #if the tryptic peptide is after the binding site
      binding_sequence <- as.character(paste0(as.character(binding_sequence), as.character(ms_sequence)))
      #need to change the binding start and/or end
      binding_site_end_in_pdb <- resno_ms2_end_pattern

    } else if(resno_ms2_start_pattern < binding_site_start_in_pdb){
      binding_sequence <- as.character(paste0(as.character(ms_sequence), as.character(binding_sequence)))
      binding_site_start_in_pdb <- resno_ms2_start_pattern


    } else if(resno_ms2_start_pattern == binding_site_start_in_pdb){
      #the two are equal
      if(include_ambiguous == TRUE){
        #the two are equal so just ignore it

      } else {
        #error?
        warning('Potential error --> Danger!\n')
      }

      #check if it's include_ambiguous?
      #if not --> there may be an error

    }

    #get the start and end positions
    #get the order of the tryptic and binding peptides
    #exclude binding of the include_ambiguous fragments
    #

  }

  
  if(binding_site_end_in_pdb > nchar(as.character(sequence_for_alignment))){
    binding_site_end_in_pdb <- nchar(as.character(sequence_for_alignment))
    warning('Binding site detected is beyond limits of sequence, switching to last amino acid\n')
  } else if(binding_site_start_in_pdb < 1){
    
    binding_site_end_in_pdb <- 1
    warning('Binding site detected is less than 1, switching to first amino acid in sequence\n')
  }

  binding_site_ouput <- list(binding_site_start=binding_site_start_in_pdb,
                             binding_site_end=binding_site_end_in_pdb,
                             binding_sequence=binding_sequence,
                             ms2_peptide_seq=ms_sequence,
                             ms2_start=resno_ms2_start_pattern,
                             ms2_end=resno_ms2_end_pattern,
                             source_sequence=sequence_for_alignment,
                             pro_name=protein_name,
                             db=database_name,
                             db_id=database_id)


  return(data.frame(binding_site_ouput))
}
