#----BS3 Align to PDB----

#does this need boolean?
#should also accept protein_dict --> but will need to account for the multiple chains
#user should also be able to just put in the PDB ID with no chain --> will have to do a grepl check to see if '_' exists in string

#' PPI Align PDB
#'
#' This functions aligns fasta file to PDB
#'
#' @param fasta_file Name of fasta file or loaded fasta file by seqinr::read.fasta().
#' @param alignIDs A data.frame containing the columns "ProteinName", "UniProtID", and "PDB"
#' @param uniprot2pdb If TRUE, will align to the UniProt sequence before aligning to the PDB. This parameter should be selected if the sequences used are not exactly UniProt (such as a slightly different N-terminal), but are relatively similar.
#' @return A list containing vectors corresponding to the length of the FASTA file sequences. Containments alignments to UniProt and PDB with corresponding chain IDs.
#' @author Emma Gail
#'
#' @export

ppi.alignPDB <- function(fasta_file, alignIDs=NULL, uniprot2pdb=TRUE){

  pdb_vector_match_mega <- list()
  stored_pdbs <- list()
  stored_annos <- list()
  stored_sequences <- list()

  for(protein_name in names(fasta_file)){

    #will need to remove check_download etc if switching to fetch system

    if(is.null(alignIDs)){
      pdb_info <- display_preferred_pdb_structure_menu(protein_name,fasta_file)
      chain <- pdb_info$chain
      pdb_id <- pdb_info$pdb_id
    } else {

      if(!(protein_name %in% as.character(alignIDs$ProteinName))){
        next
      }

      pdb_chain <- strsplit(as.character(alignIDs[alignIDs$ProteinName == protein_name,'PDB']),'_')[[1]]
      chain <- pdb_chain[2]
      pdb_id <- pdb_chain[1]
    }

    if(!(pdb_id %in% names(stored_pdbs))){
      pdb_file <- read.pdb2(pdb_id)
      stored_pdbs[[pdb_id]] <- pdb_file
      pdb_anno <- pdb.annotate(pdb_id)
      stored_annos[[pdb_id]] <- pdb_anno
    } else {
      pdb_anno <- stored_annos[[pdb_id]]
      pdb_file <- stored_pdbs[[pdb_id]]
    }



    #can update this to include option to combine PDB chains or use multiple
    #how to use multiple?


    #pdb_file <- check_download_read_pdb(pdb_info$pdb_id)


    #can update this so that if chain
    if(is.na(chain)){
      #if no chain --> will have to get uniprot_id from the alignIDs
      uniprot_id <- as.character(alignIDs[alignIDs$ProteinName == protein_name,'UniProtID'])
      uniprot_id <- strsplit(uniprot_id,'-')[[1]][1]
    } else {
      uniprot_id <- pdb_anno[pdb_anno$chainId == chain,"db_id"]
    }



    if(uniprot2pdb == TRUE){

      if(uniprot_id %in% names(stored_sequences)){
        uniprot_sequence <- stored_sequences[[uniprot_id]]
      } else {
        uniprot_sequence <- uniprot.fasta(uniprot.id = uniprot_id)
        stored_sequences[[uniprot_id]] <- uniprot_sequence
      }


    } else { #end if(uniprot2pdb == TRUE){
      uniprot_sequence <- toupper(paste0(fasta_file[[protein_name]],collapse=''))

    } #end else to if(uniprot2pdb == TRUE){


    #uniprot_sequence <-uniprot(uniprot_id)$sequence #taking a really long time?

    #if chain is na --> will just use the chain selected as the chain matches
    chain_matches <- pdb_anno[pdb_anno$db_id == uniprot_id,"chainId"]



    #check if seqres == uniprot_id


    #check if the UniProt sequence and the seqres sequence are different?

    #check to see if the seqres of the PDB file matches that of the uniprot sequence --> why would they be different?


    # pwa_results <- quick_pwa_from_pdb(fasta_file[[protein_name]],pdb_file,chain,
    #                                   use_resid_and_resno = TRUE)

    #may need to get the UniProt sequence first before going further down the pipeline
    #and use that sequence for the actual aignment code

    # pwa_results <- quick_pwa_from_pdb(fasta_file[[protein_name]],pdb_file,chain,
    #                                   use_resid_and_resno = FALSE)


    pwa_results <- pairwiseAlignment(uniprot_sequence,toupper(paste0(fasta_file[[protein_name]],collapse='')))

    pwa_strings <- get_pwa_strings(pwa_results)

    #length(get_pwa_strings(pwa_results)$pattern_string)
    pwa_ranges <- get_pwa_ranges(pwa_results) #use the start and end subject to get the range

    #accounting for the dashes that are within the pattern?
    #use the ranges to get the start position
    #make sure that the end position matches the end of the range in the actual
    #should account for gaps in both the subject and pattern
    #start with the position - 1 then +1 in the first iteration

    fasta_seq_alignment_vector <- c()
    subject_start <- pwa_ranges$start_subject - 1
    subject_string <- pwa_strings$subject_string
    for(subject_index in 1:length(subject_string)){

      #go by index or aa??

      #get the index, assign the
      #need to check if there are any dashes within the subject_string
      subject_aa <- subject_string[subject_index]
      if(subject_aa != '-'){
        subject_start <- subject_start + 1
        fasta_seq_alignment_vector <- c(fasta_seq_alignment_vector,subject_start)

      } else {
        #if it does == '-'
        #add a dash
        fasta_seq_alignment_vector <- c(fasta_seq_alignment_vector,'-')
      }

    }

    uniprot_seq_alignment_vector <- c()
    pattern_start <- pwa_ranges$start_pattern -1
    pattern_string <- pwa_strings$pattern_string
    for(pattern_index in 1:length(pattern_string)){

      #go by index or aa??

      #get the index, assign the
      #need to check if there are any dashes within the subject_string
      pattern_aa <- pattern_string[pattern_index]
      if(pattern_aa != '-'){
        pattern_start <- pattern_start + 1
        uniprot_seq_alignment_vector <- c(uniprot_seq_alignment_vector,pattern_start)

      } else {
        #if it does == '-'
        #add a dash
        uniprot_seq_alignment_vector <- c(uniprot_seq_alignment_vector,'-')
      }

    }



    #QC check
    #fasta_seq_alignment_vector <- pwa_ranges$start_subject:pwa_ranges$end_subject #add in dashes? or do that afterwards?
    #does the initial check have to be UniProt to PDB and then that is further condensed using an alignment to the FASTA file
    #uniprot_seq_alignment_vector <- pwa_ranges$start_pattern:pwa_ranges$end_pattern
    #uniprot_seq_alignment_vector <- 1:length(seqres)
    pdb_alignment_vector <- rep('-',length(uniprot_seq_alignment_vector))
    chain_alignment_vector <- rep('-',length(uniprot_seq_alignment_vector))


    #change this to the new function?
    pdb_uniprot_mapping <- GET(paste0('https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=',pdb_id))
    xml_data <- xmlParse(pdb_uniprot_mapping)
    xml_data <- xmlToList(xml_data)

    #alignment_df <- data.frame()

    #go thorugh each alignment?
    alignment_indices <- (1:length(names(xml_data)))[names(xml_data) == 'alignment']
    #(xml_data[[1]][names(xml_data[[1]]) == 'block'])
    #each one named alignment
    for(alignment_index in alignment_indices){

      block_data <- xml_data[[alignment_index]][names(xml_data[[alignment_index]]) == 'block']
      for(block_index in 1:length(block_data)){

        segment_data <- block_data[[block_index]][names(block_data[[block_index]]) == 'segment']
        segment_df <- t(data.frame(segment_data))
        #rownames(segment_df) <- NULL

        #do the block index stuff here
        #make sure there are 2 rows
        if(nrow(segment_df) == 2){

          #check the ID against the vector??
          #vectors will have to be initialized
          #name and then sub vectors --> 3 of them
          #will have to do alignment against the UniProt sequence before doing the
          #alignment
          #

          if(uniprot_id %in% segment_df[,'intObjectId']){

            #if the UniProt ID that corresponds with the input variable or whatever
            #print(segment_df)

            #segment_df <- segment_df2

            uniprot_line <- segment_df[segment_df[,"intObjectId"] == uniprot_id,] #the uniprot line
            pdb_line <- segment_df[segment_df[,"intObjectId"] != uniprot_id,] #the PDB line --> should check to make sure the PDB ID is in there

            #grepl(pdb_info$pdb_id,as.character(pdb_line[names(pdb_line) == "intObjectId"]))

            pdb_line_id <- as.character(pdb_line[names(pdb_line) == "intObjectId"])
            pdb_line_chain <- strsplit(pdb_line_id,'\\.')[[1]][2]
            pdb_line_start <- as.numeric(pdb_line[names(pdb_line) == "start"])
            pdb_line_end <- as.numeric(pdb_line[names(pdb_line) == "end"])

            uniprot_line_start <- as.numeric(uniprot_line[names(pdb_line) == "start"])
            uniprot_line_end <- as.numeric(uniprot_line[names(pdb_line) == "end"])

            #do a QC check here possibly --> check if end is > the



            #put the PDB start and end in the positions where the uniprot_line is

            #should use match() instead of directly the index for better accuracy

            #do a pairwiseAlignment beforehand with the segments from both the UniProt and PDB sequences
            #would use the resno numbering from PDB and the index number

            uniprot_start_index <- match(uniprot_line_start,uniprot_seq_alignment_vector)
            uniprot_end_index <- match(uniprot_line_end,uniprot_seq_alignment_vector)

            #return(c(uniprot_start_index,uniprot_end_index))

            if(uniprot_line_end > max(uniprot_seq_alignment_vector)){
              uniprot_end_index <- max(uniprot_seq_alignment_vector)
            }

            #uniprot_seq_alignment_vector
            pdb_alignment_vector[uniprot_start_index:uniprot_end_index] <- pdb_line_start:pdb_line_end

            pdb_line_fullid <- paste0(pdb_id,'_',pdb_line_chain)
            #can have here that if chain is NA --> will do regardless
            #will need to account for duplicates --> will need guide for how to do this
            #have override option?

            if(is.na(chain)){
              chain_alignment_vector[uniprot_start_index:uniprot_end_index] <- rep(pdb_line_fullid,length(uniprot_line_start:uniprot_line_end))
            } else {
              #if there is a specific chain indicated --> will not do it unless it matches exactly
              if(chain == pdb_line_chain){
                chain_alignment_vector[uniprot_start_index:uniprot_end_index] <- rep(pdb_line_fullid,length(uniprot_line_start:uniprot_line_end))
              }
            }


            #source("https://bioconductor.org/biocLite.R")
            #biocLite('bio3d','Biostrings')


          } #end if(uniprot_id %in% segment_df[,'intObjectId'])

          #check if the UniProt ID is in there?


        } else {
          warning('More than 2 rows found in block.segment in pdb_uniprot_mapping --> check data\n')
        } #end else to if(nrow(segment_df) == 2)
      }

    } #end for(alignment_index in alignment_indices)


    pdb_vector_match_mega[[protein_name]] <- list(fasta=fasta_seq_alignment_vector,
                                                  uniprot=uniprot_seq_alignment_vector,
                                                  pdb=pdb_alignment_vector,
                                                  chain=chain_alignment_vector)

    #once the pdb numbering/chains have been found --> need to use this info for the actual creation
    #of the PyMOL file

    #add all of the info to the mega list


    #should override any extra overlap with the originally chosen chain


    #xml_data[[1]]$block #go through the blocks --> there can be multiple blocks
    #xml_data[[1]]$block[[1]] #go through the segments

    #need to figure out how to make the best sense of the blocks
    #do the same kind of alignment with the UniProt sequence
    #do a kind of frequency vector?
    #have all '-' in it --> then add the PDB numbering to it in the same indices that are in the XML file
    #will also need another vector that contains the chain ID so that it doesn't have to be re-named


    #could then create a vector corresponding to the chain it's on
    #would need 3 vectors --> 1 for the chains, (can have - if no match), 1 for fasta sequence,
    #1 for the numbering of the chain

    #SUZ12 is split between multiple sequences --> option to combine sequences
    #together??
    #quality control --> lots of single amino acids not easily combined


  #   pwa_ranges <- get_pwa_ranges(pwa_results = pwa_results)
  #   pwa_strings <- get_pwa_strings(pwa_results = pwa_results)
  #   #resno_and_resid <- quick_resno_and_resid(pdb_file,chain)
  #
  #   pat_str_list <- paste0(pwa_strings$pattern_string,collapse='')
  #   pat_str_split <- strsplit(pat_str_list,'-')[[1]]
  #   pat_str_split <- pat_str_split[pat_str_split != '']
  #   #measure the length instead
  #   length(pat_str_split[nchar(pat_str_split) == 1])
  #
  #
  #   pwa_strings
  #
  #   fasta_vector <- c()
  #   pdb_vector <- c()
  #
  #   #go through the list
  #   for(pat_str in pat_str_split){
  #
  #     #matching the numbering to the fasta file
  #     fasta_se <- str_locate_all(toupper(paste0(fasta_file[[protein_name]],collapse='')),pat_str)[[1]]
  #     pat_se <- str_locate_all(paste0(resno_and_resid$resid,collapse=''),pat_str)[[1]]
  #     nrow(pat_se) #check if > 1 --> more than 1 match
  #     #can also check the nrow of fasta_se
  #     #can go to the next one numerically?
  #
  #
  #     #what to do if there is more than 1 match?
  #     pat_start <- pat_se[1]
  #     pat_end <- pat_se[2]
  #
  #     #get the numbering from resno from pat_se
  #     resno_start <- pdb_info$resno[pat_start]
  #     resno_end <- pdb_info$resno[pat_end]
  #
  #     fasta_vector <- c(fasta_vector,fasta_se[1]:fasta_se[2])
  #     pdb_vector <- c(pdb_vector,resno_start:resno_end)
  #
  #     #add to vector
  #
  #     #use the original numbering from fasta_se and make into a vector
  #     #both vectors should be the same length --> should do a QC at the end of the loop
  #
  #   } #end for(pat_str in pat_str_split)
  #
  #   length(fasta_vector) == length(pdb_vector)
  #
  #   pdb_vector_match <- list(fasta_vector=fasta_vector,
  #                            pdb_vector=pdb_vector)
  #
  #
  #   pdb_vector_match_mega[[protein_name]] <- pdb_vector_match
  #
  # }  #end for(protein_name in names(fasta_file))
  #


  #output should be the same/similar to the original function
  #need to find a way to handle missing points on the PDB file
  #can then be used for input for AP code
  }

  return(pdb_vector_match_mega)
} #end function ppi.alignPDB



#use the pdb_vector_match to t

# setwd('~/Downloads')
# fasta_file <- seqinr::read.fasta('PRC2_units_5m.fasta')
# pdb_vector_match <- ppi.alignPDB(fasta_file = fasta_file)
#
# #for every line in the loaded data
#
#
#
# xlink_list
#
# for(xlink_index in 1:length(xlink_list)){
#
#
#   #xlink_index <- 94
#   #xlink_list[[xlink_index]]
#
#   seq_xlink <- xlink_list[[xlink_index]]$sequence_xlink
#   pro_xlink <- xlink_list[[xlink_index]]$proteins_xlink
#   pro_xlink <- pro_xlink[1]
#
#   #need to do the check up here first to make sure the protein names are real
#
#   #rename this function?
#   protein_split_list <- extract_protein_name_and_peptide_num(pro_xlink,names(fasta_file))
#
#   #within the main function can make it so that if pdb_numbering == TRUE
#   #do some of the rest of the loop
#   #after the intiial checks
#
#
#   #have boolean to keep track of the proteins?
#
#   for(pro_index in  1:length(protein_split_list)){
#
#     pro_name <- protein_split_list[[pro_index]][1]
#     pro_seq_index <- protein_split_list[[pro_index]][2]
#
#     pdb_vector_pro <- pdb_vector_match[[pro_name]] #need to make sure that this is not NULL?
#     pro_seq_index2 <- match(pro_seq_index,pdb_vector_pro$fasta)
#
#     pdb_index <- as.numeric(pdb_vector_pro$pdb[pro_seq_index2])
#     chain_index <- pdb_vector_pro$chain[pro_seq_index2]
#
#     chain_split <- strsplit(chain_index,'_')[[1]]
#
#     #read.fasta.pdb #look into this further
#     pdb_file <- read.pdb(chain_split[1]) #store the pdb file so that it doesn't need to be re-downloaded every time?
#     pdb_resno_bool <- pdb_file$atom[pdb_file$atom$chain == chain_split[2],'resno'] == pdb_index
#     pdb_resid <- unique(a(firstup(pdb_file$atom[pdb_file$atom$chain == chain_split[2],'resid']))[pdb_resno_bool])
#     #check the index of the seq_xlink?
#     #first check to make sure the length == 1, throw an error if not
#     #will also need to check to make sure that it is not NULL or NA (if the match is not real)
#
#
#     #should double-check the actual index within
#
#
#     #as.numeric('-') #made into NA by coercion
#     #will need to check if is.na()?
#     #if cannot be made into a number --> can make NA for pdb_id
#
#     #does a second check need to happen to make sure the index exists with the
#
#   }
#
#   #get the protein name and match it to the right name within the pdb_vector_match
#   #then get the protein position and then match it within the $fasta list
#   #use match() to get the index
#   #use the index for the $pdb list
#   #use the same index for the $chain list to get the name of the PDB file that will be used for pdb1/2
#
#
# }
#
