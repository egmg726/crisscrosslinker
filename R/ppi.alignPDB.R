#----BS3 Align to PDB----


#does this need boolean?
#should also accept protein_dict --> but will need to account for the multiple chains
#user should also be able to just put in the PDB ID with no chain --> will have to do a grepl check to see if '_' exists in string
ppi.alignPDB <- function(fasta_file){

  pdb_vector_match_mega <- list()

  for(protein_name in names(fasta_file)){

    pdb_info <- display_preferred_pdb_structure_menu(protein_name,fasta_file)
    #can update this to include option to combine PDB chains


    chain <- pdb_info$chain
    pdb_file <- check_download_read_pdb(pdb_info$pdb_id)

    # pwa_results <- quick_pwa_from_pdb(fasta_file[[protein_name]],pdb_file,chain,
    #                                   use_resid_and_resno = TRUE)

    pwa_results <- quick_pwa_from_pdb(fasta_file[[protein_name]],pdb_file,chain,
                                      use_resid_and_resno = FALSE)

    length(get_pwa_strings(pwa_results)$pattern_string)
    pwa_ranges <- get_pwa_ranges(pwa_results) #use the start and end subject to get the range

    #QC check
    fasta_seq_alignment_vector <- pwa_ranges$start_subject:pwa_ranges$end_subject #add in dashes? or do that afterwards?
    uniprot_seq_alignment_vector <- pwa_ranges$start_pattern:pwa_ranges$end_pattern
    pdb_alignment_vector <- rep('-',length(uniprot_seq_alignment_vector))
    chain_alignment_vector <- rep('-',length(uniprot_seq_alignment_vector))
    #account for the dashes that are present within the seqres

    length(get_pwa_strings(pwa_results)$pattern_string) == length(pwa_ranges$start_subject:pwa_ranges$end_subject)



    #after QC check --> get the UniProt ID from the pdb_annotation
    #do another QC check to make sure that the Uniprot sequence is the same as the seqres?


    pdb_anno <- pdb.annotate(pdb_info$pdb_id)
    uniprot_id <- pdb_anno[pdb_anno$chainId == chain,"db_id"]

    chain_matches <- pdb_anno[pdb_anno$db_id == uniprot_id,"chainId"]

    #go through the chain matches? need to have some kind of prompt to see if the user wants to do that
    #if length > 1

    #load the XML data here

    for(chain_name in chain_matches){

      #get the chain name
      #
      chain_name

    }

    #

    #length(fasta_file[[protein_name]])

    #add additional option to use all chains for the menu


    pdb_file$seqres[names(pdb_file$seqres) == chain]
    #get all of the chains that exist in the PDB file


    # chain_mega_list <- c()
    # chains <- unique(names(pdb_file$seqres))
    # for(chain in unique(names(pdb_file$seqres))){
    #
    #   seqres <- paste0(a(firstup(pdb_file$seqres[names(pdb_file$seqres) == chain])),collapse='')
    #   chain_mega_list <- c(chain_mega_list,seqres)
    #   #check to see if any chains == each other
    #   #if they do == each other, add the chain to the identifier?
    #   #add to list of the seqres sequences
    #   #pairwise align the resid to the seqres
    #
    #
    #
    # }



    #pdb_id <- '6C23'
    pdb_anno <- pdb.annotate(pdb_info$pdb_id)

    pdb_anno$compound
    pdb_anno$chainId #will need all

    uniprot_chain_match <- list()


    #this can be used for a menu option
    for(uniprot_id in unique(pdb_anno[pdb_anno$db_name == "UniProt",'db_id'])){
      #get all of the chains that correspond to this uniprot_id
      #get the name of the compound
      #put this into a menu option

      uniprot_chain_match[[uniprot_id]] #or compound ID or both!
      #will have all of the chains


    } #end for(uniprot_id in unique(pdb_anno[pdb_anno$db_name == "UniProt",'db_id']))

    #get the chain IDs
    #trying to measure the overlap??


    #chain_matches <- c()
    #get the uniprot IDs

    # for(useqres in unique(chain_mega_list)){
    #
    #   chains[useqres == chain_mega_list]
    #   #keep chain information somewhere??
    #   #this should be integrated with the chain loop --> use seqres instead the actual sequence
    #
    #   #either do 1 by 1 or can do multiples
    #
    # }


    #combine the chains together --> need to check how that was done before when combining chains together
    paste0(quick_resno_and_resid(pdb_file,'A')$resid,collapse='')

    #xml_data <- xmlParse(paste0('https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=',pdb_id))


    #can make this into its own function
    #input --> pdb_id, optional input: chains and/or UniProt ID?

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

            grepl(pdb_id,as.character(pdb_line[names(pdb_line) == "intObjectId"]))

            pdb_line_id <- as.character(pdb_line[names(pdb_line) == "intObjectId"])
            pdb_line_chain <- strsplit(pdb_line_id,'\\.')[[1]][2]
            pdb_line_start <- as.numeric(pdb_line[names(pdb_line) == "start"])
            pdb_line_end <- as.numeric(pdb_line[names(pdb_line) == "end"])

            uniprot_line_start <- as.numeric(uniprot_line[names(pdb_line) == "start"])
            uniprot_line_end <- as.numeric(uniprot_line[names(pdb_line) == "end"])


            #put the PDB start and end in the positions where the uniprot_line is

            #should use match() instead of directly the index for better accuracy

            #do a pairwiseAlignment beforehand with the segments from both the UniProt and PDB sequences
            #would use the resno numbering from PDB and the index number

            uniprot_start_index <- match(uniprot_line_start,uniprot_seq_alignment_vector)
            uniprot_end_index <- match(uniprot_line_end,uniprot_seq_alignment_vector)


            uniprot_seq_alignment_vector
            pdb_alignment_vector
            chain_alignment_vector

            source("https://bioconductor.org/biocLite.R")
            biocLite('bio3d','Biostrings')


          } #end if(uniprot_id %in% segment_df[,'intObjectId'])

          #check if the UniProt ID is in there?


        } else {
          warning('More than 2 rows found in block.segment in pdb_uniprot_mapping --> check data\n')
        } #end else to if(nrow(segment_df) == 2)
      }

    } #end for(alignment_index in alignment_indices)


    #should override any extra overlap with the originally chosen chain


    xml_data[[1]]$block #go through the blocks --> there can be multiple blocks
    xml_data[[1]]$block[[1]] #go through the segments

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


    pwa_ranges <- get_pwa_ranges(pwa_results = pwa_results)
    pwa_strings <- get_pwa_strings(pwa_results = pwa_results)
    #resno_and_resid <- quick_resno_and_resid(pdb_file,chain)

    pat_str_list <- paste0(pwa_strings$pattern_string,collapse='')
    pat_str_split <- strsplit(pat_str_list,'-')[[1]]
    pat_str_split <- pat_str_split[pat_str_split != '']
    #measure the length instead
    length(pat_str_split[nchar(pat_str_split) == 1])


    pwa_strings

    fasta_vector <- c()
    pdb_vector <- c()

    #go through the list
    for(pat_str in pat_str_split){

      #matching the numbering to the fasta file
      fasta_se <- str_locate_all(toupper(paste0(fasta_file[[protein_name]],collapse='')),pat_str)[[1]]
      pat_se <- str_locate_all(paste0(resno_and_resid$resid,collapse=''),pat_str)[[1]]
      nrow(pat_se) #check if > 1 --> more than 1 match
      #can also check the nrow of fasta_se
      #can go to the next one numerically?


      #what to do if there is more than 1 match?
      pat_start <- pat_se[1]
      pat_end <- pat_se[2]

      #get the numbering from resno from pat_se
      resno_start <- pdb_info$resno[pat_start]
      resno_end <- pdb_info$resno[pat_end]

      fasta_vector <- c(fasta_vector,fasta_se[1]:fasta_se[2])
      pdb_vector <- c(pdb_vector,resno_start:resno_end)

      #add to vector

      #use the original numbering from fasta_se and make into a vector
      #both vectors should be the same length --> should do a QC at the end of the loop

    } #end for(pat_str in pat_str_split)

    length(fasta_vector) == length(pdb_vector)

    pdb_vector_match <- list(fasta_vector=fasta_vector,
                             pdb_vector=pdb_vector)


    pdb_vector_match_mega[[protein_name]] <- pdb_vector_match

  }  #end for(protein_name in names(fasta_file))



  #output should be the same/similar to the original function
  #need to find a way to handle missing points on the PDB file
  #can then be used for input for AP code


  return(pdb_vector_match_mega)
} #end function ppi.alignPDB
