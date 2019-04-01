# library(ggplot2)
# library(bio3d)
# library(Biostrings)
# library(seqinr)
# library(RColorBrewer)
# library(prodlim)
# library(openxlsx)
# library(stringr)
# library(httr)
# library(jsonlite)
# library(xml2)
# library(grDevices)
# library(svglite)
# library(XML)
# library(RCurl)

#need to do importFrom for these packages so that only certain ones are taken out of the
#packages --> no name conflicts


#only import the svg function into the final version



#-----firstup-----

#need a reference here to the stackoverflow page

#' Capitalize first letter
#'
#' This function allows you to capitalize the first letter of a string and returns all other characters as lower case.
#' @param x A string
#' @keywords capitalize
#' @export
#' @examples
#' string <- 'AAA'
#' firstup(string)
#' [1] Aaa


firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}



#-----Get resno and resid list----

#'Get resno and resid from chain in pdb_read
#'
#'This function retrieves the residues and corresponding numbers from a specific chain in a PDB file
#'@param pdb_read PDB file defined by the class in bio3d
#'@param chain Character corresponding to chain ID in PDB file
#'@export

get_resno_and_resid <- function(pdb_read,chain){
  pdb_resno <- pdb_read$atom$resno[pdb_read$atom$chain==chain]
  pdb_resid <- pdb_read$atom$resid[pdb_read$atom$chain==chain]

  resno_and_resid <- list(resno=pdb_resno,
                          resid=pdb_resid)

  return(resno_and_resid)
}



#-----consolidate resno and resid-----

#'Consolidate resno_and_resid
#'
#'This function consolidates the get_resno_and_resid() output to get rid of repeating amino acids
#'@param resno_and_resid Output from get_resno_and_resid()
#'@export

consolidate_resno_and_resid <- function(resno_and_resid){

  pdb_resid <- resno_and_resid$resid
  pdb_resno <- resno_and_resid$resno

  translated_resid <- suppressWarnings(seqinr::a(firstup(pdb_resid)))

  resno_list <- c()
  resid_list <- c()
  resid_count <- 0

  for(resno_num in pdb_resno){
    resid_count <- resid_count + 1
    if(!(resno_num %in% resno_list)){
      resno_list <- c(resno_list,resno_num)
      add_resid <- translated_resid[resid_count]
      resid_list <- c(resid_list,add_resid)
    }
  }

  resid_na <- !is.na(resid_list)
  resno_list <- resno_list[!is.na(resid_list)]
  resid_list <- resid_list[!is.na(resid_list)]

  new_resno_and_resid <- list(resno=resno_list,
                              resid=resid_list,
                              resid_na=resid_na)

  return(new_resno_and_resid)

}



#-----Get fragments------

#'Get missing amino acids/fragments
#'
#'This function finds the missing amino acids between two sequences that have been pairwise aligned
#'@param maal Output from renumber_and_identify_missing_aa()
#'@param pwa_strings Output from Biostrings::pairwiseAlignment()
#'@export

#boolean pattern == TRUE? to distinguish
get_missing_amino_acids <- function(maal,pwa_strings){

  #need to include both pattern and subject
  #missing_animo_acids <- maal$pattern

  if(is.null(maal$pattern) && is.null(maal$subject)){
    aa_fragments <- NULL

    return(aa_fragments)
  }


  missing_animo_acids <- maal$maa_pattern_numbering

  subject_string_vector <- pwa_strings$subject_string
  pattern_string_vector <- pwa_strings$pattern_string

  if(length(subject_string_vector) == 1){
    subject_string_vector <- strsplit(subject_string_vector,split='')[[1]]
    pattern_string_vector <- strsplit(pattern_string_vector,split='')[[1]]
  }

  amino_acid_fragments <- c()
  amino_acid_fragments_main_list <- list()
  list_num <- 0
  start_fragment_list <- list()
  start_fragment2 <- c()


  #supposed start is actually the end of the sequence
  if(length(missing_animo_acids) == 1){
    missing_animo_acids_length <- 1
  } else {
    missing_animo_acids_length <- length(missing_animo_acids)-1
  }

  for(aa_num in 1:(missing_animo_acids_length)){
    amino_acid_fragments <- c(amino_acid_fragments,subject_string_vector[aa_num])

    if(!is.na(missing_animo_acids[aa_num+1])){

      if((missing_animo_acids[aa_num+1] - missing_animo_acids[aa_num]) != 1){
        list_num <- list_num + 1
        amino_acid_fragments_main_list[[list_num]] <- amino_acid_fragments
        start_fragment_list[[list_num]] <- c(missing_animo_acids[aa_num],
                                             length(amino_acid_fragments))

        fragment_start <- missing_animo_acids[aa_num]-length(amino_acid_fragments)+1
        start_fragment2 <- c(start_fragment2,fragment_start)
        amino_acid_fragments <- c()
      }

    } #end if(!is.na(missing_animo_acids[aa_num+1])){


  } #end for(aa_num in 1:(missing_animo_acids_length))

  #if(!is.na(missing_animo_acids[aa_num+1])){
    list_num <- list_num + 1
    amino_acid_fragments_main_list[[list_num]] <- amino_acid_fragments
    start_fragment_list[[list_num]] <- c(missing_animo_acids[aa_num],
                                         length(amino_acid_fragments))

    fragment_start <- missing_animo_acids[aa_num]-length(amino_acid_fragments)+1
    start_fragment2 <- c(start_fragment2,fragment_start)


    aa_fragments <- list(amino_acid_fragments=amino_acid_fragments_main_list,
                         start_fragment=start_fragment_list,
                         start_fragment_num=start_fragment2)

  #}
  return(aa_fragments)
}


#-----Make renumbered vectors-----

#'Make renumbered PDB vectors
#'
#'This function renumbers the PDB file based on a new resno_and_resid
#'@param resno_and_resid Output from get_resno_and_resid() function
#'@param new_resno_and_resid Output from consolidate_resno_and_resid() function
#'@param maal Output from renumber_and_identify_missing_aa() function
#'@export

make_renumbered_pdb_vectors <- function(resno_and_resid, new_resno_and_resid, maal){

  resid_na <- new_resno_and_resid$resid_na

  new_pdb_numbering <- maal$new_pdb_numbering

  pdb_resno <- resno_and_resid$resno
  pdb_resid <- resno_and_resid$resid

  pdb_resno_freq_table <- data.frame(table(pdb_resno))
  rownames(pdb_resno_freq_table) <- pdb_resno_freq_table$pdb_resno
  pdb_resno_freq_table <- pdb_resno_freq_table[resid_na,]
  pdb_resno_list <- as.numeric(as.vector(pdb_resno_freq_table$pdb_resno))

  new_pdb_resnos <- c()
  index_count <- 0
  for(resno_num in pdb_resno_freq_table$pdb_resno){
    index_count <- index_count + 1
    freq <- pdb_resno_freq_table[resno_num,]$Freq
    #new_resno_num <- resno_num
    new_resno_num <- new_pdb_numbering[index_count]
    resno_repeats <- as.numeric(rep(new_resno_num,freq))
    new_pdb_resnos <- c(new_pdb_resnos,resno_repeats)
    if(index_count == 1){
      resno_start <- resno_num
    }
    if(index_count == pdb_resno_list[length(pdb_resno_list)]){
      resno_end <- resno_num
    }
  }

  renumbered_pdb_residues <- list(new_pdb_resnos=new_pdb_resnos,
                                  pdb_resno_list=pdb_resno_list)

  return(renumbered_pdb_residues)
}

#----Get GenBank start lines-----

get_genbank_start_lines <- function(split_genbank_file){

  line_count <- 0
  feature_start_line <- 0
  origin_start_line <- 0
  for(efetch_line in split_genbank_file){
    line_count <- line_count + 1
    if(startsWith(efetch_line,'DEFINITION') == TRUE){
      efetch_definition <- strsplit(efetch_line,'DEFINITION')[[1]][2]
    }
    if(startsWith(efetch_line,'FEATURES') == TRUE){
      feature_start_line <- line_count
    }
    if(startsWith(efetch_line,'ORIGIN') == TRUE){
      origin_start_line <- line_count
    }

  }

  genbank_start_lines <- list(efetch_definition=efetch_definition,
                              feature_start_line=feature_start_line,
                              origin_start_line=origin_start_line)

  return(genbank_start_lines)

}

#-----renumber_feature_subset------

renumber_feature_subset <- function(split_genbank_file,genbank_start_lines,pwa_ranges){

  feature_start_line <- genbank_start_lines$feature_start_line
  origin_start_line <- genbank_start_lines$origin_start_line

  features_subset <- split_genbank_file[(feature_start_line+1):(origin_start_line-1)]

  start_difference_gb <- pwa_ranges$start_difference

  feature_subset_names <- c("source","Protein","Site","Region","CDS")
  five_spaces <- "     "

  new_features_subset <- c()
  feature_line_number <- c()
  for(fs_name in feature_subset_names){
    line_number <- 0
    for(feature_line in features_subset){
      line_number <- line_number + 1
      feature_line_heading <- paste(five_spaces,fs_name,sep = '')
      if(startsWith(feature_line,feature_line_heading)){
        feature_line_number <- c(feature_line_number,line_number)
        for(split_feature in strsplit(feature_line," ")[[1]]){
          if(split_feature != '' && split_feature != fs_name){
            if(grepl('\\.\\.',split_feature)){
              feature_range <- strsplit(split_feature,'\\.\\.')[[1]]
              fs_start <- as.character(as.numeric(feature_range[1]) + start_difference_gb)
              fs_end <- as.character(as.numeric(feature_range[2]) + start_difference_gb)
              new_feature_range <- paste(fs_start,'..',fs_end,sep='')
            } else {
              new_feature_range <- as.character(as.numeric(split_feature) + start_difference_gb)
            }
            num_feature_name_spaces <- 16 - nchar(fs_name)
            feature_name_spaces <- paste(rep(" ",num_feature_name_spaces),collapse = "")
            new_feature_line <- paste(five_spaces,fs_name,feature_name_spaces,new_feature_range,sep = "")
            new_features_subset <- c(new_features_subset,new_feature_line)
          }
        }

      }
    }
  }

  renumbered_feature_lines <- list(new_features_subset=new_features_subset,
                                   feature_line_number=feature_line_number)

  return(renumbered_feature_lines)

}

#-----Create features subset-----

create_features_subset <- function(split_genbank_file,genbank_start_lines){

  feature_start_line <- genbank_start_lines$feature_start_line
  origin_start_line <- genbank_start_lines$origin_start_line

  features_subset <- split_genbank_file[(feature_start_line+1):(origin_start_line-1)]

  return(features_subset)
}

#----Edit features subset----

edit_features_subset <- function(features_subset,renumbered_feature_lines){

  new_features_subset <- renumbered_feature_lines$new_features_subset
  feature_line_number <- renumbered_feature_lines$feature_line_number

  features_subset2 <- features_subset
  for(num in 1:length(new_features_subset)){
    feature_index <- feature_line_number[num]
    nf_line <- new_features_subset[num]
    features_subset2[feature_index] <- nf_line
  }

  return(features_subset2)
}

#------condensed PDB output function-----

#future function


#------split sequences (for xlinking data)-----


#'Split Sequences
#'
#'@export
split_sequences <- function(xlink_string,proteins = TRUE,split_by_parenthesis=FALSE){

  seq_to_split <- xlink_string

  if(proteins == FALSE){
    seq_to_split <- strsplit(xlink_string,':')[[1]][1]
  }

  seq_to_split <- strsplit(seq_to_split,'-')[[1]]
  seq_to_split <- gsub("\\)","",seq_to_split)
  seq_to_split <- strsplit(seq_to_split,'\\(')

  if(split_by_parenthesis == TRUE){
    seq_to_split <- strsplit(seq_to_split,')-')[[1]]
    seq_to_split <- gsub("\\)","",seq_to_split)
    seq_to_split <- strsplit(seq_to_split,'\\(')
  }

  return(seq_to_split)

}


#------split sequences2 (for xlinking data)-----

#'Split sequences from plink output
#'
#'This function splits sequences from the plink output function
#'@param xlink_string String of BS3 XL data made by make_plink2_master_list()
#'@param proteins TRUE/FALSE if string needs to be split by ":"
#'@export

split_sequences2 <- function(xlink_string,proteins = TRUE){

  seq_to_split <- xlink_string

  if(proteins == FALSE){
    seq_to_split <- strsplit(xlink_string,':')[[1]][1]
  }

  seq_to_split <- strsplit(seq_to_split,'\\)-')[[1]]
  #seq_to_split <- strsplit(seq_to_split,'-')[[1]]
  seq_to_split <- gsub("\\)","",seq_to_split)
  seq_to_split <- strsplit(seq_to_split,'\\(')

  return(seq_to_split)

}



#------xyz coordinates-----

#'Get XYZ coordinates from PDB file
#'
#'This function makes a list of the xyz coordinates from a PDB file
#'@param trimmed_pdb PDB file loaded by bio3d::read.pdb2()
#'@export

get_xyz_coordinates_pdb <- function(trimmed_pdb){

  #x
  x_coords <- trimmed_pdb$xyz[c(T,F,F)]
  #y
  y_coords <- trimmed_pdb$xyz[c(F,T,F)]
  #z
  z_coords <- trimmed_pdb$xyz[c(F,F,T)]

  xyz_pdb_coords <- list(x_coords=x_coords,
                         y_coords=y_coords,
                         z_coords=z_coords)

  return(xyz_pdb_coords)

}

#-----get min/max of xyz coordinates-----

#'Get min and max of xyz coordinates
#'
#'This function gets the min and max of get_xyz_coordinates_pdb() function
#'@param xyz_pdb_coords Output of get_xyz_coordinates_pdb()
#'@export

get_xyz_min_max <- function(xyz_pdb_coords){

  min_max_coord_list <- list()
  for(coord_name in names(xyz_pdb_coords)){
    coord_vector <- xyz_pdb_coords[[coord_name]]
    min_coord <- min(coord_vector)
    max_coord <- max(coord_vector)
    coord_list <- c(min_coord,max_coord)
    min_max_coord_list[[coord_name]] <- coord_list
  }

  return(min_max_coord_list)
}




#----------PROTEIN-RNA FUNCTIONS----------


#nick's functions here


#----Combine FASTA files----

#'Combine FASTA files
#'
#'This function combines FASTA files from a list of file names
#'@param fasta_file_names List of file names
#'@param experiment_directory Experiment directory. If NULL, will use current directory using getwd()
#'@param fasta_file_directory If FASTA files are in a different directory. Defaults to NULL.
#'@param fasta_file_in_different_directory If FASTA files are in a different directory. Defaults to FALSE.
#'@export

fasta.combine <- function(fasta_file_names,experiment_directory = NULL,
                                fasta_file_directory = NULL,
                                fasta_file_in_different_directory = FALSE){

  if(is.null(experiment_directory)){
    experiment_directory <- getwd()
  }

  fasta_file <- c()
  for(fasta_file_name in fasta_file_names){
    if(is.null(fasta_file_directory) || fasta_file_in_different_directory == FALSE){
      fasta_file_path <- paste(experiment_directory,'/',fasta_file_name,sep='')
    } else {
      fasta_file_path <- paste(fasta_file_directory,'/',fasta_file_name,sep='')
    }
    fasta_file <- c(fasta_file,seqinr::read.fasta(fasta_file_path))
  }

  return(fasta_file)
}


#----Make Sequence Hit List----

#'Make sequence hit list
#'
#'This function creates a hit list of significant XL in UV XL-MS dataset
#'
#'@param fasta_file FASTA file as loaded by seqinr::read.fasta()
#'@param experiment_directory Directory where experiment files are. If NULL, defaults to current directory using getwd().
#'@param file_format Can be either 'csv' or 'xlsx', the suffix to the experiment files
#'@param cutoff_score Cutoff for Andromeda score. Defaults to 20.
#'@author Emma Gail
#'@export

rbd.makeSeqHitList <- function(fasta_file, experiment_directory = NULL,
                                   file_format='csv',
                                   cutoff_score = 20, RK.cleavage = TRUE){

  if(is.null(experiment_directory)){
    experiment_directory <- getwd()
  }

  protein_prefix <- names(fasta_file)

  file_pattern <- paste('\\.',file_format,'$',sep='')

  list_of_files <- list.files(experiment_directory,pattern=file_pattern)
  list_of_files <- list_of_files[!startsWith(list_of_files,'~')]
  #if startswith '~', need to remove since it messes up the function


  #need to make sure that the makeSeqHitList can also accept already loaded data.frames for
  #the tutorial

  sequence_hit_list <- list()
  for(file_name in list_of_files){

    #should join the file_name and experiment directory
    file_name1 <- paste0(experiment_directory,'/',file_name)

    if(file_format == 'csv'){
      csv_file <- read.csv(file_name1)
    } else if(file_format == 'xlsx'){
      csv_file <- data.frame(read.xlsx(file_name1))
    } else if(file_format == 'txt'){
      csv_file <- read.table(file_name1,sep='\t',header=TRUE)
    }

    csv_file_filtered <- csv_file[grep(paste(protein_prefix, collapse='|'),
                                       csv_file$Proteins, ignore.case=TRUE),]
    if(nrow(csv_file_filtered) == 0){
      sequence_hit_list[[file_name]] <- list(sequence=NULL,
                                             intensity=NULL)
      next
    }
    csv_file_filtered <- csv_file_filtered[csv_file_filtered$Score > cutoff_score,]

    if(nrow(csv_file_filtered) == 0){
      sequence_hit_list[[file_name]] <- list(sequence=NULL,
                                             intensity=NULL)
      next
    }
    list_of_sequences <- c()
    intensity_list <- c()
    protein_list <- c()
    aa_before_list <- c()
    aa_after_list <- c()

    for(row_num in 1:nrow(csv_file_filtered)){
      csv_row <- csv_file_filtered[row_num,]
      amino_acid_before <- as.character(csv_row$Amino.acid.before)
      amino_acid_after <- as.character(csv_row$Last.amino.acid)

      #make sure they're R/K?
      acceptable_aas <- c('R','K')
      add_to_shl <- FALSE
      if((amino_acid_before %in% acceptable_aas) && (amino_acid_after %in% acceptable_aas)){
        if(RK.cleavage == TRUE){
          if((amino_acid_before != amino_acid_after) && (amino_acid_before != '-')){
            add_to_shl <- TRUE
          } else {
            add_to_shl <- FALSE
          }
        } else {
          add_to_shl <- TRUE
        }
      }

      #if((amino_acid_before != amino_acid_after) && (amino_acid_before != '-')){
      if(add_to_shl == TRUE){
        intensity_list <- c(intensity_list,csv_row$Intensity)
        list_of_sequences <- c(list_of_sequences,as.character(csv_row$Sequence))
        protein_list <- c(protein_list,as.character(csv_row$Proteins))
        aa_before_list <- c(aa_before_list,amino_acid_before)
        aa_after_list <- c(aa_after_list,amino_acid_after)

      }

    }

    sequence_hit_list[[file_name]] <- list(sequence=list_of_sequences,
                                           intensity=intensity_list,
                                           protein=protein_list,
                                           amino_acid_before=aa_before_list,
                                           amino_acid_after=aa_after_list)
  }

  if(is.null(unlist(sequence_hit_list))){
    warning('No matches to selected protein names, returning NULL\nTIP: Check to make sure that your FASTA file protein names match the protein names in your MaxQuant file')
    return(NULL)
  } else {
    return(sequence_hit_list)
  }

} #end function rbd.makeHitLists


#----Make Eluate/Input Output Table----

#'Make input and eluate table
#'
#'This function makes an input and eluate table based on the make_sequence_hit_list() function
#'@param prefixes Prefix in file name identifying the experiments
#'@param sequence_hit_list Output from rbd.makeSeqHitList() function
#'@param proteases Proteases used in the experiment. Defaults to c('ArgC','LysC')
#'@param outline_file_name Name of output file. Defaults to 'pipeline_output_file.csv'
#'@author Emma Gail
#'@export

rbd.makeIETable <- function(sequence_hit_list,
                            prefixes,
                            proteases = c('ArgC','LysC'),
                            output_file_name = 'pipeline_output_file.csv',
                            experiment_directory = NULL){
  secondary_prefix <- TRUE

  if(is.null(experiment_directory)){
    experiment_directory <- getwd()
  }

  df_output_list <- list()
  for(prefix in prefixes){
    sub_shl <- sequence_hit_list[grepl(prefix,names(sequence_hit_list))]
    cat(prefix,'\n')

    if(secondary_prefix == TRUE){
      for(sec_prefix in proteases){
        new_sub_shl <- sub_shl[grepl(sec_prefix,names(sub_shl))]
        new_prefix <- paste(prefix,'_',sec_prefix,sep='')
        #cat(paste(sec_prefix,'\n'))
        cat(new_prefix)
        output_df <- make_intensity_output_df(sub_shl = new_sub_shl,prefix = new_prefix)
        df_output_list[[new_prefix]] <- output_df
      }
    } else { #if secondary_prefix == FALSE
      output_df <- make_intensity_output_df(sub_shl = sub_shl,prefix = prefix)
      df_output_list[[prefix]] <- output_df
    }

  }

  #return(df_output_list)

  #another category --> prefix_category for color?
  combined_output_df <- do.call(rbind, df_output_list)
  combined_output_df_nice_colnames <- combined_output_df

  if('input' %in% colnames(combined_output_df)){
    input_col_name <- 'Input'
  } else {
    input_col_name <- NULL
  }
  colnames(combined_output_df_nice_colnames) <- c('Sequence',
                                                  input_col_name,
                                                  'Eluate',
                                                  'Category',
                                                  'Protein Name',
                                                  'Amino Acid Before',
                                                  'Amino Acid After')

  write.csv(combined_output_df_nice_colnames,output_file_name)
  return(combined_output_df)

}


#----Make Input/Eluate Graph----

#'Make input/eluate plot
#'
#'This function takes the data from the input/eluate table and makes it into a plot
#'@param input_eluate_table Output from make_input_eluate_table()
#'@param prefixes List of control and experimental prefixes
#'@param proteases List of secondary prefixes in file names, defaults to proteases used 'ArgC' and 'LysC'
#'@param palette Palette from RColorBrewer to be used for plot. Defaults to 'Dark2'
#'@param fill Same as "fill" from ggplot2 for plot
#'@param colour Same as "colour" from ggplot2 for plot
#'@param size Same as "size" from ggplot2 for plot
#'@param pipeline_generated_plot_title Title for generated plot. Must end in '.png'. Defaults to 'pipeline_generated_plot.png'.
#'@author Emma Gail
#'@export

rbd.makeIEPlot <- function(input_eluate_table,
                           prefixes,
                           proteases = c("LysC",'ArgC'),
                           title = 'Eluate vs. Input',
                           palette = 'Dark2',
                           fill = "#eaf2ff",
                           colour = "#eaf2ff",
                           size = 0.5,
                           pipeline_generated_plot_title='pipeline_generated_plot.png'){


  if(!('input' %in% colnames(input_eluate_table))){
    warning('No input column found, returning NULL\n')
    return(NULL)
  }

  list_of_prefixes <- prefixes
  list_of_secondary_prefixes <- proteases
  # if(length(prefixes) > 1){
  #   secondary_prefix <- TRUE
  # } else {
  #   secondary_prefix <- FALSE
  # }
  secondary_prefix <- TRUE
  combined_output_df <- input_eluate_table
  cols <- colorRampPalette(brewer.pal(8,palette))(length(list_of_secondary_prefixes))
  shapes <- c(16,17,15,18) #add in additional shapes as needed

  #make color list if secondary_prefix == TRUE
  plot_colors <- c()
  plot_shapes <- c()
  for(category_level in levels(factor(combined_output_df$category))){
    if(secondary_prefix == TRUE){

      protease_split <- NULL
      for(protease in proteases){
        if(grepl(protease,category_level)){
          protease_split <- protease
        }
      }

      #split_category_level <- strsplit(category_level,'_')[[1]]

      split_category_level <- strsplit(category_level,paste0('_',protease_split))[[1]]
      #prefix_name <- split_category_level[1]
      #secondary_prefix_name <- split_category_level[2]
      prefix_name <- split_category_level
      #secondary_prefix_name <- protease
      prefix_index <- match(prefix_name,list_of_prefixes)
      #secondary_prefix_index <- match(secondary_prefix_name,list_of_secondary_prefixes)
      secondary_prefix_index <- match(protease_split,list_of_secondary_prefixes)
      plot_colors <- c(plot_colors,cols[secondary_prefix_index])
      plot_shapes <- c(plot_shapes,shapes[prefix_index])
    } else {
      prefix_name <- category_level
      prefix_index <- match(prefix_name,list_of_prefixes)
      plot_shapes <- c(plot_shapes,shapes[prefix_index])

    }

  }

  nudge_variable <- (max(combined_output_df$input)-min(combined_output_df$input))/50


  #change to svg?
  png(pipeline_generated_plot_title, width = 480, height = 480)

  if(secondary_prefix == TRUE){
    p <- ggplot(combined_output_df, aes(x=input,y=eluate)) +
      geom_point(aes(shape = factor(category),colour=factor(category)),size=2) +
      geom_text(aes(label=ifelse(eluate>0,as.character(sequence),'')),hjust=0,vjust=0,nudge_x = nudge_variable) +
      labs(shape = "Category") +
      labs(x = 'Input (Peak Intensity, A.U.)') +
      labs(y = 'Eluate (Peak Intensity, A.U.)') +
      labs(title = title) +
      theme(
        panel.background = element_rect(fill = fill,
                                        colour = colour,
                                        size = 0.5, linetype = "solid"),
        legend.key = element_rect(fill = fill)
      ) +

      scale_colour_manual(name = "",
                          labels = levels(factor(combined_output_df$category)),
                          values = plot_colors) +

      scale_shape_manual(name = "",
                         labels = levels(factor(combined_output_df$category)),
                         values = plot_shapes)

  } else {

    p <- ggplot(combined_output_df, aes(x=input,y=eluate)) +
      geom_point(aes(shape = factor(category)),size=2) +
      geom_text(aes(label=ifelse(eluate>0,as.character(sequence),'')),hjust=0,vjust=0,nudge_x = nudge_variable) +
      labs(shape = "Category") +
      labs(x = 'Input (Peak Intensity, A.U.)') +
      labs(y = 'Eluate (Peak Intensity, A.U.)') +
      labs(title = title) +
      theme(panel.background = element_rect(fill = fill,
                                            colour = colour,
                                            size = size,
                                            linetype = "solid"),
            legend.key = element_rect(fill = fill)
      ) +

      scale_shape_manual(name = "",
                         labels = levels(factor(combined_output_df$category)),
                         values = plot_shapes)

  }

  dev.off()


  return(p)
}





#-----Make intesnity output dataframe-----

#need to include secondary prefix in here?

#'Make intesnity output dataframe
#'
#'This function makes a dataframe out of a subsection of the sequence_hit_list
#'@param sub_shl Subsection of sequence_hit_list, ouput of make_sequence_hit_list()
#'@param prefix Prefix that will be used to categorized the subsection in the output
#'@author Emma Gail
#'@export

make_intensity_output_df <- function(sub_shl,prefix){

  elute_file_name <- names(sub_shl)[grepl('elute',names(sub_shl),ignore.case = TRUE)]
  if(length(elute_file_name) == 0){
    elute_file_name <- names(sub_shl)[grepl('eluate',names(sub_shl),ignore.case = TRUE)]
  }
  input_file_name <- names(sub_shl)[grepl('input',names(sub_shl), ignore.case = TRUE)]



  if(length(input_file_name) == 0){
    input_file_name <- elute_file_name
    cat('No input file detected --> using eluate file\n')
    input <- FALSE
  } else {
    input <- TRUE
  }


  input_sequences <- sub_shl[[input_file_name]]$sequence
  input_intensities <- sub_shl[[input_file_name]]$intensity
  eluate_itensities <- sub_shl[[elute_file_name]]$intensity
  eluate_sequences <- sub_shl[[elute_file_name]]$sequence
  list_of_proteins <- sub_shl[[input_file_name]]$protein
  list_of_eluate_proteins <- sub_shl[[elute_file_name]]$protein
  aa_before_list <- sub_shl[[input_file_name]]$amino_acid_before
  aa_before_list_eluate <- sub_shl[[elute_file_name]]$amino_acid_before
  aa_after_list <- sub_shl[[input_file_name]]$amino_acid_after
  aa_after_list_eluate <- sub_shl[[elute_file_name]]$amino_acid_after

  #prefix_category_list <- rep(prefix,times=length(input_sequences))

  #test_df <- data.frame(list(col1=c(1,1,2),col2=c(2,2,2),col3=c(4,5,6)))
  #duplicated(test_df[c('col1','col2')])
  #test_df[!duplicated(test_df[c('col1','col2')]),]
  #check to make sure last column == 0


  i_e_check_df <- data.frame(list(seqs=c(input_sequences,eluate_sequences),
                                  pros=c(list_of_proteins,list_of_eluate_proteins),
                                  input_ints=c(input_intensities,rep(0,length(eluate_sequences))),
                                  aa_before=c(aa_before_list,aa_before_list_eluate),
                                  aa_after=c(aa_after_list,aa_after_list_eluate)))

  #i_e_check_df <- i_e_check_df[order(i_e_check_df$seqs),]
  i_e_check_df <- i_e_check_df[!duplicated(i_e_check_df[c('seqs','pros')]),]

  #after the check if there are duplicates, can then add the

  #c(input_sequences

  #need to account for the fact that the peptides may show up in eluate but not input
  output_df <- matrix(nrow=nrow(i_e_check_df),ncol=7)
  output_df <- data.frame(output_df)
  rownames(output_df) <- i_e_check_df$seqs
  sequence_label <- 'sequence'
  input_y_axis_label <- 'input'
  elute_x_axis_label <- 'eluate'
  category_label <- 'category'
  protein_label <- 'protein'
  aa_before_label <- 'aa_before'
  aa_after_label <- 'aa_after'
  colnames(output_df) <- c(sequence_label,input_y_axis_label,
                           elute_x_axis_label,category_label,
                           protein_label,aa_before_label,
                           aa_after_label)

  input_sequences2 <- input_sequences
  input_sequences <- as.character(i_e_check_df$seqs)
  input_intensities <- as.numeric(i_e_check_df$input_ints)
  prefix_category_list <- rep(prefix,times=length(input_sequences))

  seq_count <- 0
  for(input_seq in input_sequences){
    seq_count <- seq_count + 1


    if(!(input_seq %in% eluate_sequences) || is.null(eluate_sequences)){
      #if the input sequence is not in the elute file

      output_df[input_seq,elute_x_axis_label] <- 0
      output_df[input_seq,input_y_axis_label] <- input_intensities[seq_count]

    } else {
      intensity_index <- match(input_seq,eluate_sequences)

      output_df[input_seq,elute_x_axis_label] <- eluate_itensities[intensity_index]
      output_df[input_seq,input_y_axis_label] <- input_intensities[seq_count]

    }

  }

  output_df[,category_label] <- prefix_category_list
  output_df[,sequence_label] <- input_sequences
  output_df[,protein_label] <- as.character(i_e_check_df$pros)
  output_df[,aa_after_label] <- as.character(i_e_check_df$aa_after)
  output_df[,aa_before_label] <- as.character(i_e_check_df$aa_before)

  output_df <- output_df[(!((output_df$input) == 0 & (output_df$input==output_df$eluate))),]

  if(input == FALSE){
    output_df[[input_y_axis_label]] <- NULL
  }


  return(output_df)

}

#-----Updated Menu Loop Function------

#can be a useful function on its own so should have an updated name to reflect that
#blast.menu (?)

#'BLAST sequence and go through menu loops
#'
#'This function will BLAST the given sequence and give the user an interactive menu to select the best option.
#'@param fasta_sequence_vector Sequence as a vector to be used to BLAST against
#'@param database Menu selected for blast.pdb(). Options: "pdb","nr", or "swissprot". Defaults to "pdb".
#'@param time_out Integer indicating time in seconds the blast.pdb() will search before time out occurs. Defaults to
#'@export

go_through_menu_loops <- function(fasta_sequence_vector, database = 'pdb', time_out = NULL){

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


#----Quick PWA from PDB (combine 3 functions into one)----

#'Quick PWA from PDB
#'
#'This function is a shortcut to do a BioStrings::pairwiseAlignment() with a selected chain in a PDB file
#'@param query_sequence Sequence to match the PDB sequence against
#'@param pdb_read PDB file loaded by bio3d::read.pdb2()
#'@param chain Character indicating the chain on the PDB file
#'@param use_resid_and_resno Will get the resno_and_resid using get_resno_and_resid() from chain in PDB file rather than the seqres for alignment. Defaults to FALSE.
#'@export

quick_pwa_from_pdb <- function(query_sequence,pdb_read,chain,use_resid_and_resno = FALSE){
  #check if query sequence is a string or a vector of characters
  if(length(query_sequence) != 1){ #not a string
    query_sequence <- paste(query_sequence,collapse='')
  }

  #possibly get rid of the resno and resid portion of this code and replace with seqres if possible

  query_sequence <- toupper(query_sequence)

  seq_on_chain <- a(firstup(as.vector(pdb_read$seqres[names(pdb_read$seqres) == chain])))

  if(length(seq_on_chain) == 0 || use_resid_and_resno == TRUE){
    if(length(seq_on_chain) == 0){
      cat('No SEQRES line in PDB file - using ATOM line instead\n')
    }

    resno_and_resid <- get_resno_and_resid(pdb_read = pdb_read,
                                           chain = chain)
    resno_and_resid2 <- consolidate_resno_and_resid(resno_and_resid = resno_and_resid)
    seq_on_chain <- resno_and_resid2$resid
  }

  sequence_from_pdb <- toupper(paste(seq_on_chain,collapse=''))
  pwa_results <- pairwiseAlignment(sequence_from_pdb,query_sequence)

  #return both the resno_and_resid2 as well as the pwa_results?
  return(pwa_results)
}


#----Quick resno/resid function----

#'resno_and_resid Shortcut Function
#'
#'This function makes the resno_and_resid list directly from a given PDB file and chain
#'@param pdb_read PDB file as loaded by bio3d::read.pdb2()
#'@param chain Character indicating the chain on the PDB file
#'@export

quick_resno_and_resid <- function(pdb_read,chain){

  resno_and_resid <- get_resno_and_resid(pdb_read = pdb_read,
                                         chain = chain)
  resno_and_resid2 <- consolidate_resno_and_resid(resno_and_resid = resno_and_resid)

  #return both the resno_and_resid2 as well as the pwa_results?
  return(resno_and_resid2)
}


#----SeqRes Chain Loop----

#'Match query sequence to chain in PDB file
#'
#'This function will go through all chains in the PDB file to output list of chains that match the sequence
#'@param protein_name Protein name that will be used for the console output
#'@param fasta_seq Sequence from FASTA file
#'@param pdb_read PDB file as loaded by bio3d::seqinr()
#'@param pdb_rl Name of PDB file
#'@param pwa_score_threshold Minimum score to determine a match using Biostrings::pairwiseAlignment()
#'@param peptide_sequence An optional peptide sequence to use for recommendations if trying to match a peptide to a structure. Defaults to NULL
#'@param fasta_file Optional fasta file input
#'@export

#maybe add an option for the fasta file and then
do_sr_chain_loop <- function(protein_name,
                             fasta_seq,
                             pdb_read,
                             pdb_rl,
                             pwa_score_threshold = 500,
                             peptide_sequence = NULL,
                             fasta_file = NULL){

  ff_pdb_read <- pdb_read

  seqres_matching_chains <- list()
  seqres_matching_pdb_info <- list()
  seqres_matching_recommended <- list()

  seqres_chains <- levels(factor(names(ff_pdb_read$seqres)))
  if(length(seqres_chains) == 0){
    seqres_chains <- levels(factor(ff_pdb_read$atom$chain))
  }


  #add in contingency in case seqres doesn't work
  #

  for(sr_chain in seqres_chains){
    cat(paste('Checking chain ',sr_chain,' on sequence ',protein_name,'\n',sep=''))

    sr_pwa_results <- suppressWarnings(quick_pwa_from_pdb(query_sequence = fasta_seq,
                                         pdb_read = ff_pdb_read,
                                         chain = sr_chain))

    if(!is.null(peptide_sequence)){
      #if it's not null --> check if peptide is in the PDB chain
      #get the resid and do the string check to see if it's in there

      pdb_resid <- paste(quick_resno_and_resid(ff_pdb_read,sr_chain)$resid,collapse='')

      if(length(str_locate_all(peptide_sequence,pdb_resid)[[1]]) > 0){
        #if it's greater than 0 --> there is at least 1 suggestion in there
        #add to recommendation list

        seqres_matching_recommended[[protein_name]] <- c(seqres_matching_recommended[[protein_name]],sr_chain)

      }

    }

    seqres_pwa_ranges <- get_pwa_ranges(sr_pwa_results)
    seqres_pdb_info <- list(query_start=seqres_pwa_ranges$start_subject,
                            query_end=seqres_pwa_ranges$end_subject,
                            sequence_start=seqres_pwa_ranges$start_subject,
                            sequence_end=seqres_pwa_ranges$end_subject,
                            pdb_id=pdb_rl,
                            chain=sr_chain,
                            search_method='pwa')

    sr_pwa_score <- score(sr_pwa_results)

    if(sr_pwa_score >= pwa_score_threshold){ #above threshold score

      pdb_and_chain <- paste(seqres_pdb_info$pdb_id,'_',sr_chain,sep='')
      srmc_index_name <- paste(pdb_and_chain,'_',protein_name,sep='')

      seqres_matching_pdb_info[[srmc_index_name]] <- seqres_pdb_info
      seqres_matching_chains[[protein_name]] <- c(seqres_matching_chains[[protein_name]],sr_chain)

    }

  } #end for loop sr_chain

  return(list(pdb_info = seqres_matching_pdb_info,
              chains = seqres_matching_chains,
              recommend = seqres_matching_recommended))

}

#----Make Binding Site Dataframe----

#adding the fasta file sequence as an argument so that it can be compared
#to the binding sequence that was found using the PDB file
#can mark it as a different color for a mutation?

#'Make binding site dataframe
#'
#'This function makes a row of the RBDMap binding_site_df based on a row of the input/eluate table
#'@param pdb_info List of information made from go_through_menu_loops() or do_sr_chain_loop()
#'@param filtered_input_eluate_table_row A row of input/eluate table made from make_in make_input_eluate_table()
#'@param fasta_seq Vector or string of protein sequence to match the input/eluate information against
#'@export

make_binding_site_df <- function(pdb_info,
                                 filtered_input_eluate_table_row,
                                 fasta_seq){

  if(length(fasta_seq) != 1){
    fasta_seq <- toupper(paste(fasta_seq,collapse=''))
  }

  #checking to make sure that pdb_info has a resid string
  #otherwise need to generate it using the id and chain info contained
  #in pdb_info
  #check using 'resid' %in% names(pdb_info)

  #need to reconfigure this to make sure the sequence is the correct one

  #these start/end patterns refer to a pairwise alignment done between the peptide
  #sequence and the chain

  #should probably check if the sequence is in the chain before it is outputted to the
  #user

  #resid_string <- paste(pdb_info$resid,collapse='')
  peptide_sequence <- filtered_input_eluate_table_row$sequence

  #should probably use these results for the rest of the code
  #what to do when pwa_results are way off?


  #should check to make sure that there are no dashes in output

  #pairwiseAlignment('EVSTAPAGTDMPAAK',toupper(paste(fasta_file$EED_O75530,collapse='')))


  # #can always try an additional method of searching for the string within the string
  # #why is


  resid_string <- paste(pdb_info$resid,collapse='')

  if(length(str_locate_all(pattern = peptide_sequence, resid_string)[[1]]) == 0){

    cat('Sequence not found in the chain selected\n')

    #ask user if they want to switch chains?
    #should have this function be earlier in the workflow
    #then can recommend to the user which chain they should pick

    #should then switch to the fasta file to get the binding sequence?
    #the sequence against the fasta file should always be checked --> add to binding sequence output?

    pdb_info$pdb_id <- NA
    pdb_info$chain <- NA
    pdb_info$resid <- strsplit(fasta_seq,'')[[1]]
    pdb_info$sequence_start <- 1
    pdb_info$sequence_end <- nchar(fasta_seq)
    pdb_info$resno <- 1:nchar(fasta_seq)
    pdb_info$search_method <- 'original_fasta'

    #change the pdb_info so that it changes the



  }

  resid_string <- paste(pdb_info$resid,collapse='')
  pwa_results <- pairwiseAlignment(resid_string,peptide_sequence)
  pwa_ranges <- get_pwa_ranges(pwa_results)
  #checking the dashes in the pwa_strings to see if there are any gaps
  #so that it can be noted

  start_pattern <- pwa_ranges$start_pattern
  end_pattern <- pwa_ranges$end_pattern

  aa_before_peptide_seq <- filtered_input_eluate_table_row$aa_before
  aa_after_peptide_seq <- filtered_input_eluate_table_row$aa_after

  #need to check this to see which direction to go to
  #can use && sequences to see which statement

  aa_after_peptide_seq

  #can do it just once and then make the decision based on that
  #(grepl('ArgC',filtered_input_eluate_table_row$category)) && (aa_after_peptide_seq == 'R')
  #if this is true need to go backwards for the
  #if false --> need to see which one is false
  #or



  #need to extract the last amino acid to decide which direction to go instead of the
  #protease used

  go_forwards <- NULL

  #turn this into a function?
  if(grepl('ArgC',filtered_input_eluate_table_row$category)){

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
    # aa_increment <- 1
    # current_aa_start <- end_pattern - 1
    # end_binding_aa <- 'R'
  }
  if(grepl('LysC',filtered_input_eluate_table_row$category)){
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


    # aa_increment <- -1
    # current_aa_start <- start_pattern
    # #end_binding_aa <- 'K'
  }

  end_binding_aa <- aa_after_peptide_seq
  current_aa_index <- current_aa_start + aa_increment
  binding_sequence <- c()
  current_aa <- pdb_info$resid[current_aa_index]

  #maybe can be the same but not use the PDB file in order to get the binding
  #sequence

  if(is.na(current_aa)){
    #not



    # binding_site_start_in_pdb <- NA
    # binding_site_end_in_pdb <- NA
    #
    #
    # pymol_mega_list <- list(protein_name=protein_name,
    #                         binding_site_start_in_pdb=binding_site_start_in_pdb,
    #                         binding_site_end_in_pdb=binding_site_end_in_pdb,
    #                         binding_sequence=binding_sequence,
    #                         pdb_id=pdb_info$pdb_id,
    #                         pdb_chain=pdb_info$chain,
    #                         ms2_peptide_seq=filtered_input_eluate_table_row$sequence,
    #                         ms2_start_in_pdb=resno_ms2_start_pattern,
    #                         ms2_end_in_pdb=resno_ms2_end_pattern,
    #                         category=filtered_input_eluate_table_row$category)
    #
    #

  }


  while((current_aa != end_binding_aa) && (current_aa_index <= length(pdb_info$resid))){
    if(aa_increment < 0){
      binding_sequence <- c(current_aa,binding_sequence)
      current_aa_index <- current_aa_index + aa_increment
      current_aa <- pdb_info$resid[current_aa_index]
    } else {
      current_aa_index <- current_aa_index + aa_increment
      current_aa <- pdb_info$resid[current_aa_index]
      binding_sequence <- c(binding_sequence,current_aa)
    }

    #may need an if/else statement

    #move current_aa increment down here to get rid of 'R' at end
  }

  binding_sequence <- paste(binding_sequence,collapse='')
  #can just do the str_locate_all (?) to get the start and end?

  resno_ms2_peptide <- pdb_info$resno[start_pattern:end_pattern]
  resno_ms2_start_pattern <- min(resno_ms2_peptide)
  resno_ms2_end_pattern <- max(resno_ms2_peptide)

  #need to change these for ArgC increments
  #only works for LysC

  if(is.null(go_forwards)){

    binding_site_start_in_pdb <- NA
    binding_site_end_in_pdb <- NA
    resno_binding_peptide <- NA

    #none selected --> have some kind of error --> assign NA to the binding sequence(?)

  } else if(go_forwards == TRUE){
    binding_site_start_in_pdb <- end_pattern + 1
    binding_site_end_in_pdb <- current_aa_index
    resno_binding_peptide <- pdb_info$resno[binding_site_start_in_pdb:binding_site_end_in_pdb]


  } else if(go_forwards == FALSE){
    binding_site_start_in_pdb <- current_aa_index + 1
    binding_site_end_in_pdb <- start_pattern - 1
    resno_binding_peptide <- pdb_info$resno[binding_site_start_in_pdb:binding_site_end_in_pdb]
  }

  # if(grepl('LysC',filtered_input_eluate_table_row$category)){
  #   #can change these two options to go forward or go backwards
  #   binding_site_start_in_pdb <- current_aa_index + 1
  #   binding_site_end_in_pdb <- start_pattern - 1
  #   resno_binding_peptide <- pdb_info$resno[binding_site_start_in_pdb:binding_site_end_in_pdb]
  # }
  # if(grepl('ArgC',filtered_input_eluate_table_row$category)){
  #   binding_site_start_in_pdb <- end_pattern + 1
  #   binding_site_end_in_pdb <- current_aa_index
  #   resno_binding_peptide <- pdb_info$resno[binding_site_start_in_pdb:binding_site_end_in_pdb]
  # }

  #if there are missing peptides that should not be highlighted,
  #may have to change this code so that it is a vector and not a
  #range of values
  binding_site_start_in_pdb <- min(resno_binding_peptide)
  binding_site_end_in_pdb <- max(resno_binding_peptide)


  #pairwiseAlignment the fasta_seq sequence and the binding sequence?
  #make sure that it is in the right place and there are no gaps
  #can do it with the peptide sequence as well to see if the start
  #and end are only 1 aa apart

  #compare the residues that were in the PDB file to the ones that
  #are in the actual fasta file

  #get actual resnos from this list to be used for the PyMol
  #visualization


  pymol_mega_list <- list(protein_name=protein_name,
                          binding_site_start_in_pdb=binding_site_start_in_pdb,
                          binding_site_end_in_pdb=binding_site_end_in_pdb,
                          binding_sequence=binding_sequence,
                          pdb_id=pdb_info$pdb_id,
                          pdb_chain=pdb_info$chain,
                          ms2_peptide_seq=filtered_input_eluate_table_row$sequence,
                          ms2_start_in_pdb=resno_ms2_start_pattern,
                          ms2_end_in_pdb=resno_ms2_end_pattern,
                          category=filtered_input_eluate_table_row$category)


  return(pymol_mega_list)

} #end making_binding_site_df function

#----Color PyMol-----

#' Color PyMol
#'
#' Generate colors to be used when generating a PyMol script
#'
#' @param vars List of variables to be used
#' @param colors List of either PyMol colors, hexcodes, or standard colors in R. Will also accept color palette names from RColorBrewer and viridis. If left as NULL, will use standard PyMOL tint colors.
#' @param png.name Name of the legend output. Defaults to "pymol_legend.jpg" unless otherwise specified.
#' @param gray0 If TRUE, will shift color palette and replace 1st color as gray (corresponds 0 if using frequency). Defaults to FALSE.
#' @author Emma Gail
#' @export

#zero indices as a parameter?
#have vars as a number? can have boolean that overrides if it is a number
#change png.name to something more general to indicate that it can output different picture file formats
color.pymol <- function(vars, colors = NULL, png.name = 'pymol_legend%03d.svg', gray0 = FALSE, print.legend = TRUE){

  num <- length(vars)

  is_hexcode <- FALSE

  if(is.null(colors)){
    #choose pymol_tints as the chosen colors if no colors are chosen
    #can also automatically choose a spectrum if it is not long enough


    chosen_pymol_colors <- c('wheat','palegreen','lightblue','paleyellow','lightpink',
                             'palecyan','lightorange','bluewhite')

    if(num > length(chosen_pymol_colors)){

      #chosen_pymol_colors <- brewer.pal(num,'Blues')
      chosen_pymol_colors <- colorRampPalette(brewer.pal(9,"Blues"))(num)

      cat('No color chosen --> reverting to RColorBrewer palette "Blues" \n')

    } else {

      chosen_pymol_colors <- chosen_pymol_colors[1:num]

      cat('No color chosen --> reverting to PyMol Tints palette\n')

    } #end else to if(num > length(chosen_pymol_colors))



  } else if((typeof(colors) == 'character') && (length(colors) == 1) && (!startsWith(colors,"#"))){

    #number will have to be the length of the
    #can also have it within a tryCatch loop?

    #length? nrow? --> will depend on intial
    #chosen_pymol_colors <- brewer.pal(num,colors)

    viridis_colors <- c("magma","inferno","plasma","viridis","cividis")

    if(colors %in% rownames(brewer.pal.info)){
      colorbrew_maxcolor <- brewer.pal.info[colors,]$maxcolors

      if(num <= colorbrew_maxcolor){
        chosen_pymol_colors <- brewer.pal(num,colors)
      } else {
        chosen_pymol_colors <- colorRampPalette(brewer.pal(colorbrew_maxcolor,colors))(num)
      }

    } else if(colors %in% viridis_colors){ #end if (colors %in% rownames(brewer.pal.info))

      chosen_pymol_colors <- viridis::viridis_pal(option = colors)(num)

    } else {

      warning("Invalid color name\n")
      return(NULL)

    } #end else if(colors %in% viridis_colors)


    #chosen_pymol_colors <- colorRampPalette(brewer.pal(9,colors))(num)

    is_hexcode <- TRUE

    #brewer.pal(8,'null')

  } else if((typeof(colors) == 'character') && (startsWith(colors,"#"))){

    chosen_pymol_colors <- colors
    is_hexcode <- TRUE

  } else { #end else if(typeof(colors) == 'character')

  chosen_pymol_colors <- colors

  }



  #need to load this as part of the package
  #pymol_color_table <- read.csv('~/Downloads/pymol_colors.csv')
  #devtools::use_data(pymol_color_table,crisscrosslinker)

  pymol_hexcode_list <- c()
  #freq_list <- c()
  #freq_num <- 0
  load_colors <- c()
  all_colors <- c()

  if(is.numeric(vars) & (gray0 == TRUE)){


    all_colors <- c(all_colors, 'gray')
    pymol_hexcode_list <- c(pymol_hexcode_list,'#808080')
  }


  #return(all_colors)
  #return(chosen_pymol_colors)

  #return(chosen_pymol_colors)

  for(pymol_color_name in chosen_pymol_colors){

    #freq_num <- freq_num + 1

    pymol_col_row <- pymol_color_table[pymol_color_table$name == pymol_color_name,]
    if(nrow(pymol_col_row) == 0){
      #no match --> do something else

      #see if it's a hexcode
      if(startsWith(pymol_color_name,'#')){
        #then it's a hexcode
        pymol_hexcode_list <- c(pymol_hexcode_list,pymol_color_name)
        #freq_list <- c(freq_list,freq_num)

        #if it's a hexcode, will have to add the color as a custom color
        #and keep track of the custom colors in a list so they will
        #correspond to the correct


      } else if(pymol_color_name %in% colors()){

        pymol_hexcode_list <- c(pymol_hexcode_list, col2hex(pymol_color_name))

      } else { #end else if(pymol_color_name %in% colors())

        #it's unknown --> return error?

      } #end else to if(startsWith(pymol_color_name,'#'))

    } else { #end if(nrow(pymol_col_row) == 0)
      #there's a match in PyMol colors

      pymol_hexcode <- rgb(as.character(pymol_col_row$red), as.character(pymol_col_row$green), as.character(pymol_col_row$blue), names = NULL, maxColorValue = 1)

      pymol_hexcode_list <- c(pymol_hexcode_list,pymol_hexcode)
      #freq_list <- c(freq_list,freq_num)

      #add hexcode to list of colors
      #update the list for the frequency or whatever the other part of the legend is
      #will either be frequency or the actual sequence depending on if it is a differential analysis
      #for RBDMap -->
      #should have option to make gray == 0, boolean upon entry
      #should the

    } #end to else of if(nrow(pymol_col_row) == 0)

  } #end for(pymol_color_name in chosen_pymol_colors)

  #if the pymol colors are hexcodes will need to translate them into rbg colors
  #for PyMol

  #if the original pymol colors were not actually pymol colors



  if(is_hexcode == TRUE){

    rbg_df <- t(col2rgb(pymol_hexcode_list)/255)

    for(row_num in 1:nrow(rbg_df)){
      rbg_row <- rbg_df[row_num,]

      #paste(rbg_row,collapse=', ')

      #if()
      custom_color_name <- paste0('custom_color_',row_num)
      set_color_line <- paste0("cmd.set_color('",custom_color_name,"', (",paste(rbg_row,collapse=', '),"))")

      load_colors <- c(load_colors,set_color_line)

      #add to list?
      #need to add both

      #need to have list of colors
      #for the

      all_colors <- c(all_colors,custom_color_name)

    }

    #add colors in a list

  } else { #end if(is_hexcode == TRUE)

    #load_colors <- NULL


  } #end else to if(is_hexcode == TRUE)


  #will need to make all of the colors at the beginning of the
  #will be named by palette_number, otherwise customcolor_number

  #add a list of these colors to the very beginning of the PyMol script
  #within the loop
  #cmd.set_color('myfavoritecolor', (1.0, 0.0, 0.0))

  #chosen_pymol_colors

  #return(all_colors)

  # if(!is.null(all_colors)){
  #   if(length(vars) < length(all_colors)){
  #     all_colors <- all_colors[1:length(vars)]
  #     pymol_hexcode_list <- pymol_hexcode_list[1:length(vars)]
  #     load_colors <- load_colors[1:length(vars)]
  #   }
  #
  # }



  #when a color is selected, use the table to get the appropriate RBG numbers to make into a hexcode



  # plot(1, type="n", axes=FALSE, xlab="", ylab="")
  # #can then just add the fill and the legends as lists using the hexcodes
  # legend("topleft", legend = vars, fill=pymol_hexcode_list,
  #        title = '', bty = 'n')


  #if vars is numeric --> should make the legend that corresponds to a more continuous scale
  #create a bar plot that will look like a legend

  #par(mfrow=c(1,1))
  #par(mar=c(7,4,4,6))
  #barplot(pymol_cc$vars,col=pymol_cc$hexcodes)

  #barplot(c(1,2,3),col=c('yellow','red','green'))

  #return(vars)

  if(is.numeric(vars)){ #may want to add a boolean in case they decide to turn this off manually

    #png(png.name, width = 200, height = 500)
    svglite(png.name, width = 4, height = 10)
    par(mar=rep(.5, 4), oma=rep(3, 4), las=1)
    image(1, vars, t(seq_along(vars)), col=pymol_hexcode_list, axes = FALSE)
    axis(4)

  } else {
    #png(png.name, width = 200, height = 500)
    if(typeof(vars) == 'character'){
      max_vars <- max(nchar(vars))
    } else {
      max_vars <- 1
    }


    #should have multple options for output --> svgs, pngs, etc.

    #png(png.name, width = (100 + (max_vars*10)), height = (150+(length(vars)*20)))
    svglite(png.name)
    #gridsvg("gridsvg.svg")
    #par(mfrow=c(1,1))
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    #can then just add the fill and the legends as lists using the hexcodes
    legend("topleft", legend = vars, fill=pymol_hexcode_list,
           title = '', bty = 'n')

    #p1 <- ggplot()
    #p1 +

  }
  dev.off()

  if(print.legend == TRUE){

    if(is.numeric(vars)){
      #par(mar=rep(.5, 4), oma=rep(3, 4), las=1)
      #print(image(1, vars, t(seq_along(vars)), col=pymol_hexcode_list, axes = FALSE))
    } else {
      plot(1, type="n", axes=FALSE, xlab="", ylab="")
      #can then just add the fill and the legends as lists using the hexcodes
      print(legend("topleft", legend = vars, fill=pymol_hexcode_list,
                   title = '', bty = 'n'))

    }

  }


  #library(plotrix)

  #plot(1, type="n", axes=FALSE, xlab="", ylab="")
  #color.legend(0.2,0.2,2,2,legend = c(1,2,3),rect.col = c('blue','green','red'))

  if(is.null(all_colors) || ((all_colors == 'gray') & length(vars > 1))){
    all_colors <- c(all_colors,chosen_pymol_colors)
  }

  if(!is.null(all_colors)){
    if(length(vars) < length(all_colors)){
      all_colors <- all_colors[1:length(vars)]
      pymol_hexcode_list <- pymol_hexcode_list[1:length(vars)]
      load_colors <- load_colors[1:length(vars)]
    }

  }


  color_list <- list(vars=vars,color_names=all_colors,hexcodes=pymol_hexcode_list,set_colors=load_colors)

  return(color_list)

} #end function color.pymol


#----Write PyMol----

#pymol.data should be in either the form of a list
#can check typeof() before input to see if it's a list of data or if it's a single dataframe
#should check original make_diff_analysis function to see what's what

#should have two different types of colors that will be used for the analysis
#ppi.colors --> if null, will use colors from the frequency color
#rbd.colors --> if null, will use a spectrum if more than one file and just default settings if not

write.pymol <- function(pymol.data, colors = NULL, experiment.dir = NULL, write.file = TRUE){

  #boolean to say whether it's ppi or rbd or both
  #accounting for both to combine the two??

  #can have an append statement that would add to the original pymol

  #if there are colnames that means it's a dataframe
  #don't need to go through them like a list

  #if there are no colnames than it is just a list
  is.null(colnames(c(1,2,3)))
  #than go through each of the names
  #check typeof() to see if character --> if character, load the csv or xlsx file
  #than go through it as a data.frame
  #should do some error checking
  #if it is a list than it is a dataframe

  #than check the dataframes to see what the column names are and use that to decide which
  #type of analysis it needs to go through


  #get the
  typeof(pymol.data)
  #if(typeof(pymol.data))

  #should have one format that will be used for the rest of the
  #code


  #extract the vars out of the pymol_mega_list

  #vars <- c(1,2,3)
  #pymol_cols <- color.pymol(vars,colors) #output the hexcolors as well?
  #display the pymol colors as well?
  #should do a par() probably so that both pymol colors can be outputted together in the same
  #plot
  #currently they will override each other as is if done separately

  #can also output this if the user wants
  data.frame(pymol_cols)

  #start of pymol lines
  #if NULL, can still add to it using c()
  pymol_lines <- pymol_cols$set_colors




  #recognize the column names since they are unique to each dataset

  #will have to deal with
  #PPI:
  #xlink_df

  #RBDMap
  #vectors of 0 and 1 (frequency) --> will need information about the PDBs as well
  #or should just have the list of PDB aligned
  #if there are multiple RBDmap dataframes inputted, it may automatically do a frequency analysis
  #bs_output

  #sort all data to start with and then do analysis for the PyMol file
  #can sort into multiple lists

  #will have to track whether or not the PDBS are different as they're being added to their
  #resepctive lists

  #for ppi --> see if '___' exists in the name of the file
  #for rbd --> check if PDB is listed as the db in any of the list dbs (means it wqs aligned)

  ppi_list <- list()
  rbd_list <- list()

  #choose the vars for each of the lists

  #ppi vars --> the frequency
  #will need to input the colors that have been used for the
  #rbd vars --> the binding sequence if there is only 1 RBD file, the frequency if there is more than 1

  #categorize each and then add to each list
  #add numerical identifier to each
  #length(list_name) + 1


  #start with if there is just 1 and then work from there

  #ppi: function already to write

  #if there are multiple ppis in list --> will have to differentially analyze them
  #should automatically make categories for them?

  #if there are multiple rbds in list --> will also have to differentially analyze them


  #should check the PDB IDs because they will probably be different and ask the user if they
  #want to put both on the same structure

  #will need to realign the crosslinking data to the original PDB (probably)
  #should be flexible

  #fasta as input as well?

  #can ignore the fasta if able to use the actual PDB files to align them to each other
  #for instance take

  made_pdb_file_name <- '~/Downloads/ProXL_test/His-TEV-Tub4-yeast___5FM1_C.pdb'

  chain <- substrRight(strsplit(made_pdb_file_name,'.pdb')[[1]],1)
  strsplit(made_pdb_file_name,'.pdb')[[1]]

  made_pdb <- read.pdb2('~/Downloads/ProXL_test/His-TEV-Tub4-yeast___5FM1_C.pdb')
  real_pdb <- read.pdb2('~/Downloads/ProXL_test/5FM1.pdb') #check_download_read?

  made_resno_resid <- quick_resno_and_resid(made_pdb,chain)
  real_resno_resid <- quick_resno_and_resid(real_pdb,chain)

  #check to make sure the sequences are the same between the two
  (paste0(made_resno_resid$resid,collapse='') == paste0(real_resno_resid$resid,collapse=''))

  #make sure the lengths are the same
  #depending on which option has been chosen --> align one to the other

  ###this scenario is more likely --> focus
  #if the real_resno is chosen as the one to align to, only '___' in the xlink_dataframe will be used


  #if made_resno is chosen as the one to align to, will have to use the start_and_end or fasta file

  #with the

  pdb_files <- c()

  #when to do the colors??
  rbd_colors
  ppi_colors


  if(is.null(colnames(pymol.data))){ #can stipulate in documentation that must be a dataframe if 1 and a list if > 1
    #if null, it is a list and should be treated like a list
    #means that there might be more than one file

    #organize each of the data files by if it's a ppi or rbd

    for(index_num in 1:length(pymol.data)){

      data_file <- pymol.data[[index_num]]
      if(typeof(data_file) == 'character'){
        #it's the name of a file name and needs to be loaded into the global environment
        if(endsWith(data_file,'.csv')){
          data_file <- read.csv(data_file)
        } else if(endsWith(data_file,'.xlsx')){
          data_file <- data.frame(read.xlsx(data_file))
        }
      }

      if(('binding_sequence' %in% colnames(data_file)) || ('Binding.Sequence' %in% colnames(data_file))){
        #it's rbd
        #add to list
        rbd_list[[length(rbd_list)+1]] <- data_file

        #add the list of pdb files to a list
        data_file$db_id #remove the '_' from the PDB ID?


      } else {
        #it's ppi
        #add to list
        ppi_list[[length(ppi_list)+1]] <- data_file
        if('pdb1' %in% colnames(data_file)){
          pdb_files <- c(pdb_files,as.character(data_file[['pdb1']]))
          pdb_files <- c(pdb_files,as.character(data_file[['pdb2']]))
        } else {
          pdb_files <- c(pdb_files,as.character(data_file[['PDB1']]))
          pdb_files <- c(pdb_files,as.character(data_file[['PDB2']]))
        }
        pdb_files <- unique(pdb_files)

        #have an or statement in here as well if length of the pdb file name is only 4 chars (real PDB file)
        #what about files that have been generated by the BS3 workflow that do not have the numbering changed?

        pdb_files <- pdb_files[grepl('___',pdb_files)]
        #unlist(strsplit(pdb_files,'.pdb'))
        #need to see what's in there
        #can't just delete the ones that don't have '____'
        #can also see if they're only 4 characters plus .pdb --> can remove .pdb first and then see if length() == 4

      }

      #pdb_files <- c(pdb_files,'6C23')
      real_pdbs <- pdb_files[unlist(lapply(pdb_files,nchar)) == 4]
      made_pdbs <- pdb_files[!unlist(lapply(pdb_files,nchar)) == 4]

      pdb_matches <- list()

      if(length(real_pdbs) > 0){
        for(real_pdb_name in real_pdbs){
          #made_pdbs[grepl(real_pdb_name,made_pdbs)]
          #categorize all of the names into a list
          #what to do about those that don't match
          if(length(made_pdbs[grepl(real_pdb_name,made_pdbs)]) > 0){
            pdb_matches[[real_pdb_name]] <- made_pdbs[grepl(real_pdb_name,made_pdbs)]
          } #end if(length(made_pdbs[grepl(real_pdb_name,made_pdbs)]) > 0){

        } #end for(real_pdb_name in real_pdbs){

      } #end if(length(real_pdbs) > 0){


      if(length(pdb_matches) > 0){
        #if length > 0, it means there is overlap between the two sets of systems
        #alert user and create menu so they can choose where to align
        #for this purpose --> always align to the known structure

        #turn on boolean?

      }

      real_pdbs %in% names(pdb_matches) #can do this because no hits does not get added to the matches
      made_pdbs %in% as.character(unlist(pdb_matches)) #see if there are any excluded matches

      real_pdbs[!(real_pdbs %in% names(pdb_matches))]
      made_pdbs[!(made_pdbs %in% as.character(unlist(pdb_matches)))]



      #this will have to be within the length(pdb_matches) > 0 loop
      # for(pdb_name in names(pdb_matches)){
      #   #will have to replace all of the PDB names with the original PDB ID
      #   #also will have to exclude any created PDB files in the differential analysis
      #   #replace all of the numbering in the BS3 file
      # }
      #

      #extract the chain from the made_pdb


      #when it comes to the


      #grepl(c(real_pdbs,'5WAI'),made_pdbs)

      #need to figure out plan if there are no real PDBs that correspond to one of the made up PDBs




      #after all files are collected in pdb_files
      #need to check if there are any conflicts
      #should there be 2 different lists for this as well?
      #and then it can be a mute point if one list is empty


      #check to see if there's a a 4 character length


      #check if it's rbd or ppi
      #will have to check if it's a character first to see if it's the name of a csv or xlsx file that
      #needs to be loaded
      #make sure that everything is already loaded before going to part 2 for analysis
      #to ensure consistency
      #will do this by the colnames() of the dataframe

      #then will have to load the pdb file names into a list so that the list can check
      #whether or not there are multiples


      #will need 2 different lists
      #sort into the right list


      #for ppi, should only load files that contain '___' in them for the


      #will have to keep a list of the PyMol files to be loaded in the final version before
      #the second part of the function happens



    } #end for(index_num in 1:length(pymol.data))


  } #end if(is.null(colnames(pymol.data)))
  #else statement for if it really is just one dataset


  #once the data has been sorted into the two lists
  #will have two separate procedures for loading onto the PyMol file

  #with the xlink_df --> get all PDB files grepl('___')
  #get protein position
  #get the resno from the PDB file
  #make sure the vectors of both the PDBs are the same length (QC)
  #get the index (index() function??) of the position number
  #match that to the index number of the real PDB --> new position number
  #write the line for the PyMOL file


  #RBDmap
  #filter the results by PDB database
  #get the id
  #get the binding site start and end
  #


  #here is where the pdb_files list will have to be checked
  #do unique()? (maybe should be done earlier, levels in dataframe and then unique)
  #remove all of the .pdbs at the end
  #remove all of the '___'

  if(length(ppi_list) > 0){
    #incorporate BS3 protocol here for loading?

    #use the colors that already exist in the file?
    #may require more user input
    #can do by category

    #so that both can be on the same need par(mfrow=c(2,1)) (is this the right order??)
    #though should check to make sure both

    #ppi_data <- read.csv('/Users/emmagail/Downloads/XL-MS Xlink Visual/BS3 - 6 - PRC2-EPOP/make_diff_analysis Run on 2018-05-11 17-45-04/BS3 - 6 - PRC2-EPOP_xlink_dataframe.csv')

    for(index_num in 1:length(ppi_list)){

      ppi_data <- ppi_list[[index_num]]
      #what vars to use for color.pymol()?
      #frequency?
      ppi.pymol(ppi_data)


    }


    color.pymol
  }



  if(length(rbd_list) > 0){

    #if the length == 1, should use the peptides as vars for color.pymol
    #if length > 1, should use the frequency of the peptides
    #use rbd.freqVector to get the frequency info
    #then get a list corresponding to the unique numbers (in order) that show up in all
    #of the vectors, maybe unlist them unique(unlist(list_here))
    rbd.freqVector(bs_output = bs_output)

    color.pymol #can also change the name of the output legend
  }

  color.pymol #will need 2 different color.pymols for each of the datasets

  for(index_num in 1:length(pymol.data)){
    if(is.null(pymol.data[[index_num]])){
      next
    }

    pdb_id <- pymol_mega_list[[index_num]]$pdb_id
    if(!(pdb_id %in% loaded_pdbs_in_pymol)){
      #pymol_file_line <- paste('fetch ',pdb_id,sep='')
      #pymol_file <- c(pymol_file, pymol_file_line)
      pymol_file_line <- paste('load ',paste(experiment_directory,'/',pdb_id,'.pdb',sep=''),sep='')
      pymol_file <- c(pymol_file, pymol_file_line)
      pymol_file_line <- paste('show surface, ',pdb_id,sep='')
      pymol_file <- c(pymol_file, pymol_file_line)
      pymol_file_line <- paste('color gray, ',pdb_id,sep='')
      pymol_file <- c(pymol_file, pymol_file_line)
      loaded_pdbs_in_pymol <- c(loaded_pdbs_in_pymol,pdb_id)
    }

    #identify where the residues are with the chain
    pdb_chain <- pymol_mega_list[[index_num]]$pdb_chain

    bs_start <- pymol_mega_list[[index_num]]$binding_site_start
    bs_end <- pymol_mega_list[[index_num]]$binding_site_end

    peptide_pymol_selection <- paste('/',pdb_id,'/','/',pdb_chain,
                                     '/',bs_start,'-',bs_end, sep='')

    pymol_file_line <- paste('color ',pymol_colors[index_num],', ',
                             peptide_pymol_selection,
                             sep='')


    pymol_file <- c(pymol_file, pymol_file_line)

    name_of_pymol_selection <- paste(pdb_id,'_',pdb_chain,'_',
                                     bs_start,'-',bs_end, sep='')

    pymol_file_line <- paste('select ',name_of_pymol_selection,', ',
                             peptide_pymol_selection,
                             sep='')

    pymol_file <- c(pymol_file, pymol_file_line)

  }

  pymol_file_name <- 'pymol_output.pml'

  if(write.file == TRUE){
    write(paste(pymol_file,collapse = '\n'),pymol_file_name)
  } else {
    return(pymol_file)
  }




} #end function write.pymol()




#----Make PyMol file----

#writing separate files for each of the PDB files that are being
#loaded into the PyMol space

#will have to rewrite this as write.pymol
#use for both PyMol creations? can have selection 'type' with 'ppi' or 'rbd' as options
write_pymol_file <- function(pymol_mega_list, colors = NULL, experiment_directory = NULL){


  #should have one format that will be used for the rest of the
  #code


  #extract the vars out of the pymol_mega_list

  color.pymol(vars,colors)

  #if null --> can set the colors to pymol_tints
  #brewer.pal(9,"Blues")
  #will have to have a startsWith for the pymol color to see if it's a hexcode

  #save the PyMOL colors to


  #save plot to the directory along with the pymol

  #build sets of colors using the different PyMol colors (much like in RColorBrewer)

  #user should be able to convert any colors they want into colors for PyMol
  #can use standard colors from RColorBrewer as custom colors for PyMol
  #or just use standard PyMol Colors

  # pymol_reds <- c('red','tv_red','raspberry','darksalmon','salmon','deepsalmon',
  #                 'warmpink','firebrick','ruby','chocolate','brown')
  # pymol_greens <- c('green','tv_green','chartreuse','splitpea','smudge',
  #                   'palegreen','limegreen','lime','limon','forest')
  # pymol_blues <- c('blue','tv_blue','marine','slate','lightblue','skyblue',
  #                  'purpleblue','deepblue','density')
  # pymol_yellows <- c('yellow','tv_yellow','paleyellow','yelloworange',
  #                    'limon','wheat','sand')
  # pymol_magentas <- c('magenta','lightmagenta','hotpink','pink','lightpink',
  #                     'dirtyviolet','violet','violetpurple','purple',
  #                     'deeppurple')
  # pymol_cyans <- c('cyan','palecyan','aquamarine','greencyan','teal','deepteal',
  #                  'lightteal')
  # pymol_oranges <- c('orange','tv_orange','brightorange','lightorange','yelloworange',
  #                    'olive','deepolive')
  #
  # pymol_tints <- c('wheat','palegreen','lightblue','paleyellow','lightpink',
  #                  'palecyan','lightorange','bluewhite')
  #
  # #pick a random color from each list
  # pymol_colors <- c('green','blue','red','purple','orange')

  #construct the PyMol file needed to visualize the structure
  #index_num <- 2

  if(is.null(colnames(pymol_mega_list))){
    #if null, it is a list and should be treated like a list
  }

  for(index_num in 1:length(pymol_mega_list)){
    if(is.null(pymol_mega_list[[index_num]])){
      next
    }

    pdb_id <- pymol_mega_list[[index_num]]$pdb_id
    if(!(pdb_id %in% loaded_pdbs_in_pymol)){
      #pymol_file_line <- paste('fetch ',pdb_id,sep='')
      #pymol_file <- c(pymol_file, pymol_file_line)
      pymol_file_line <- paste('load ',paste(experiment_directory,'/',pdb_id,'.pdb',sep=''),sep='')
      pymol_file <- c(pymol_file, pymol_file_line)
      pymol_file_line <- paste('show surface, ',pdb_id,sep='')
      pymol_file <- c(pymol_file, pymol_file_line)
      pymol_file_line <- paste('color gray, ',pdb_id,sep='')
      pymol_file <- c(pymol_file, pymol_file_line)
      loaded_pdbs_in_pymol <- c(loaded_pdbs_in_pymol,pdb_id)
    }

    #identify where the residues are with the chain
    pdb_chain <- pymol_mega_list[[index_num]]$pdb_chain

    bs_start <- pymol_mega_list[[index_num]]$binding_site_start
    bs_end <- pymol_mega_list[[index_num]]$binding_site_end

    peptide_pymol_selection <- paste('/',pdb_id,'/','/',pdb_chain,
                                     '/',bs_start,'-',bs_end, sep='')

    pymol_file_line <- paste('color ',pymol_colors[index_num],', ',
                             peptide_pymol_selection,
                             sep='')


    pymol_file <- c(pymol_file, pymol_file_line)

    name_of_pymol_selection <- paste(pdb_id,'_',pdb_chain,'_',
                                     bs_start,'-',bs_end, sep='')

    pymol_file_line <- paste('select ',name_of_pymol_selection,', ',
                             peptide_pymol_selection,
                             sep='')

    pymol_file <- c(pymol_file, pymol_file_line)

  }

  pymol_file_name <- 'pymol_output.pml'
  write(paste(pymol_file,collapse = '\n'),pymol_file_name)

  return(pymol_file)

}

#----Write Binding Sequence CSV file----

write_binding_sequence_csv <- function(pymol_mega_list){

  pymol_df_column_names <- c('Protein Name',
                             'Binding Site Start (in PDB File)',
                             'Binding Site End (in PDB File)',
                             'Binding Sequence',
                             'PDB ID',
                             'PDB Chain',
                             'MS/MS Peptide Sequence',
                             'MS/MS Site Start (in PDB File)',
                             'MS/MS Site End (in PDB File)',
                             'Category')


  #accounting for nulls in list
  pymol_mega_list_length <- 0
  for(num in 1:length(pymol_mega_list)){
    if(!(is.null(pymol_mega_list[[num]]))){
      pymol_mega_list_length <- pymol_mega_list_length + 1
    }
  }

  pymol_df <- data.frame(matrix(unlist(pymol_mega_list),
                                nrow=pymol_mega_list_length,
                                byrow = TRUE))

  colnames(pymol_df) <- pymol_df_column_names
  write.csv(pymol_df,'binding_site_sequence_data_frame.csv')

}

#----Check, Download, and Read PDB Files----

#'Check Download Read PDB
#'
#'This function is a wrapper for the bio3d PDB functions. It checks if the PDB exists before downloading.
#'
#'@param pdb_id PDB ID
#'@export

check_download_read_pdb <- function(pdb_id){
  if(!file.exists(paste(pdb_id,'.pdb',sep=''))){


    get_pdb_test <- tryCatch(get.pdb(pdb_id),
                         error = function(err){
                           #cat('No results for found for this sequence')
                           return(NULL)
                         })

    if(is.null(get_pdb_test)){
      return(NULL)
    }

  } #end if(!file.exists(paste(pdb_id,'.pdb',sep='')))


  return(read.pdb2(paste(pdb_id,'.pdb',sep='')))
}

#----Find PDB match and Get Info----


find_pdb_match_and_get_info <- function(input_sequence,
                                        protein_name,
                                        full_sequence = NULL){

  if(is.null(full_sequence)){
    full_sequence <- input_sequence
  }

  #aligning the input and full sequences to get the start and end?
  initial_pwa <- pairwiseAlignment(full_sequence,input_sequence)
  initial_ranges <- get_pwa_ranges(initial_pwa)


  input_start <- initial_ranges$start_pattern
  input_end <- initial_ranges$end_pattern


  if(length(input_sequence) != 1){
    input_sequence <- toupper(paste(input_sequence,collapse = ''))
  }

  pdb_info <- go_through_menu_loops(input_sequence)
  if(is.null(pdb_info)){
    return(NULL)
  }
  pdb_read <- check_download_read_pdb(pdb_info$pdb_id)
  pwa_results <- quick_pwa_from_pdb(input_sequence,
                                    pdb_read,
                                    pdb_info$chain)

  pwa_ranges <- get_pwa_ranges(pwa_results)
  pwa_strings <- get_pwa_strings(pwa_results)

  input_range_start <- (input_start:input_end)[pwa_ranges$start_subject]
  input_range_end <- (input_start:input_end)[pwa_ranges$end_subject]

  identifier <- paste(protein_name,'_',pdb_info$pdb_id,
                      '_',pdb_info$chain,'_',input_range_start,
                      '-',input_range_end,sep='')

  updated_pdb_info <- list(query_start=input_range_start,
                           query_end=input_range_end,
                           sequence_start=pwa_ranges$start_pattern,
                           sequence_end=pwa_ranges$end_pattern,
                           pdb_id=pdb_info$pdb_id,
                           protein_name=protein_name,
                           chain=pdb_info$chain)

  updated_pdb_info <- c(updated_pdb_info,pwa_strings)

  mega_list <- list()
  mega_list[[identifier]] <- updated_pdb_info

  return(mega_list)
}

#----Split and Match Peptides----

#accounting for non-matching short sequences
#function will not work if fasta vector is entered --> could account for that
#enclosing in some kind of loop for a generate_2d function
split_and_match_peptides <- function(fasta_sequence,
                                     sequence_name,
                                     non_matching_aa_threshold = 20){

  matching_pdb_master_list <- list()

  #first check if fasta sequence is in proper form
  if(length(fasta_sequence) != 1){
    fasta_sequence <- toupper(paste(fasta_sequence,collapse=''))
  }

  pdb_info <- find_pdb_match_and_get_info(fasta_sequence,sequence_name)
  if(is.null(pdb_info)){
    return(NULL)
  }

  #add pdb_info to mega list
  matching_pdb_master_list <- c(matching_pdb_master_list,pdb_info)

  pdb_info <- pdb_info[[names(pdb_info)]]

  pdb_read <- check_download_read_pdb(pdb_info$pdb_id)

  pwa_results <- quick_pwa_from_pdb(fasta_sequence,pdb_read,pdb_info$chain)
  pwa_ranges <- get_pwa_ranges(pwa_results)


  #start if statement
  if(pwa_ranges$start_subject != 1){
    new_seq_end <- pwa_ranges$start_subject - 1
    new_seq <- strsplit(fasta_sequence,split='')[[1]][1:new_seq_end]
    new_seq_string <- toupper(paste(new_seq,collapse = ''))
    if(length(new_seq) > non_matching_aa_threshold){

      #does it need to be renamed as something else?
      pdb_info2 <- find_pdb_match_and_get_info(new_seq_string,
                                               sequence_name,
                                               fasta_sequence)


      matching_pdb_master_list <- c(matching_pdb_master_list,pdb_info2)
      #add pdb_info to mega list

      #does this need to be a recursive function? --> avoid since computationally expensive

    }

  }

  #end if statement
  if(pwa_ranges$end_subject != nchar(fasta_sequence)){
    new_seq_start <- pwa_ranges$end_subject + 1
    new_seq <- strsplit(fasta_sequence,split='')[[1]][new_seq_start:nchar(fasta_sequence)]
    new_seq_string <- toupper(paste(new_seq,collapse = ''))

    #check if the new_seq is larger than the threshold
    if(length(new_seq) > non_matching_aa_threshold){

      pdb_info2 <- find_pdb_match_and_get_info(new_seq_string,
                                               sequence_name,
                                               fasta_sequence)


      matching_pdb_master_list <- c(matching_pdb_master_list,pdb_info2)

    }

  }

  #contains all of the matches that have been made within the
  #function

  return(matching_pdb_master_list)

} #end of split_and_match_peptides function


#----Load pLink File----

#add xlsx compatibility to this function?
load_plink_file <- function(csv_file_path){

  #if file ends with csv --> read csv
  #if file ends with xlsx --> read xlsx and make into data.frame


  if(typeof(csv_file_path) == 'character'){

    if(endsWith(csv_file_path,'csv')){
      xlink_csv <- read.csv(csv_file_path)
    } else if(endsWith(csv_file_path,'xlsx') || endsWith(csv_file_path,'xls')){
      xlink_csv <- data.frame(read.xlsx(csv_file_path))
    }

  } else if(typeof(csv_file_path) == 'list'){ #end if(typeof(csv_file_path) == 'character'){

      xlink_csv <- csv_file_path

  } else {
    warning('Unsupported file')
  }
  # if(grepl('.csv',csv_file_path)){
  #
  #   xlink_csv <- read.csv(csv_file_path)
  #
  # } else if(grepl('.xlsx',csv_file_path) || grepl('.xls',csv_file_path)){
  #   xlink_csv <- data.frame(read.xlsx(csv_file_path))
  #
  # }

  rows_with_order_nums <- c()
  row_num <- 0
  xlink_list <- list()
  for(order_value in xlink_csv$Order){
    row_num <- row_num + 1
    if(!is.na(order_value)){

      rows_with_order_nums <- c(rows_with_order_nums,row_num)
      next_row_list <- c()
      next_row <- row_num + 1
      #next_row_list <- c(next_row_list,next_row)
      while(is.na(xlink_csv$Order[next_row]) && next_row <= nrow(xlink_csv)){
        #print(c(row_num,next_row))
        next_row_list <- c(next_row_list,next_row)
        next_row <- next_row + 1
      }

      sequence_xlink <- as.character(xlink_csv$Sequence[row_num])

      spectrum_xlink <- c()
      score_xlink <- c()
      calc_m_xlink <- c()
      delta_xlink <- c()
      ppm_xlink <- c()
      modification_xlink <- c()
      sampleid_xlink <- c()
      engine_xlink <- c()
      rank_xlink <- c()
      proteins_xlink <- c()

      for(nr_num in next_row_list){

        spectrum_xlink <- c(spectrum_xlink,as.character(xlink_csv[nr_num,3]))
        score_xlink <- c(score_xlink,as.character(xlink_csv[nr_num,4]))
        calc_m_xlink <- c(calc_m_xlink,as.character(xlink_csv[nr_num,5]))
        delta_xlink <- c(delta_xlink,as.character(xlink_csv[nr_num,6]))
        ppm_xlink <- c(ppm_xlink,as.character(xlink_csv[nr_num,7]))
        modification_xlink <- c(modification_xlink,as.character(xlink_csv[nr_num,8]))
        sampleid_xlink <- c(sampleid_xlink,as.character(xlink_csv[nr_num,9]))
        engine_xlink <- c(engine_xlink,as.character(xlink_csv[nr_num,10]))
        rank_xlink <- c(rank_xlink,as.character(xlink_csv[nr_num,11]))
        proteins_xlink <- c(proteins_xlink,as.character(xlink_csv[nr_num,12]))

      }

      xlink_sub_list <- list(sequence_xlink=sequence_xlink,
                             spectrum_xlink=spectrum_xlink,
                             score_xlink=score_xlink,
                             calc_m_xlink=calc_m_xlink,
                             delta_xlink=delta_xlink,
                             ppm_xlink=ppm_xlink,
                             modification_xlink=modification_xlink,
                             sampleid_xlink=sampleid_xlink,
                             engine_xlink=engine_xlink,
                             rank_xlink=rank_xlink,
                             proteins_xlink=proteins_xlink)

      xlink_list[[order_value]] <- xlink_sub_list
    }
  }

  return(xlink_list)
}

#----Generate 2D structure from missing amino acids----

#how to incorporate this type of function with the splitting function?
find_and_generate_missing_aa_2d_structure <- function(fasta_file,
                                                      fasta_sequence_name,
                                                      pdb_file,pdb_id,
                                                      chain){

  list_of_pdb_file_names <- c()

  pwa_results <- quick_pwa_from_pdb(fasta_file[[fasta_sequence_name]],pdb_file,chain)
  pwa_ranges <- get_pwa_ranges(pwa_results = pwa_results)
  pwa_strings <- get_pwa_strings(pwa_results = pwa_results)
  resno_and_resid <- get_resno_and_resid(pdb_read = pdb_file, chain = chain)
  new_resno_and_resid <- consolidate_resno_and_resid(resno_and_resid)

  maal <- renumber_and_identify_missing_aa(pwa_strings = pwa_strings, start_difference =   pwa_ranges$start_difference)

  aa_fragments <- get_missing_amino_acids(maal = maal, pwa_strings = pwa_strings)

  if(!is.null(aa_fragments)){ #proceed if there exists values in aa_fragments
    aa_count <- 0
    for(amino_acid_fragments in aa_fragments$amino_acid_fragments){
      aa_count <- aa_count + 1
      start_residue <- aa_fragments$start_fragment_num[aa_count]
      new_pdb_file <- generate_2d_pdb_file(amino_acid_fragments = amino_acid_fragments, start_residue = start_residue)
      pdb_file_name2 <- paste(fasta_sequence_name,
                              '_missing_sequence_',
                              as.character(aa_count),
                              '.pdb',sep = '')
      list_of_pdb_file_names <- c(list_of_pdb_file_names,pdb_file_name2)
      write(new_pdb_file,pdb_file_name2,sep='')
    }
  } else {
    cat('There are no missing peptides to generate a 2D file\n')
  }


  start_subject <- pwa_ranges$start_subject
  end_subject <- pwa_ranges$end_subject

  if(start_subject != 1){
    before_alignment_aa <- fasta_file[[fasta_sequence_name]][1:start_subject-1]
    before_alignment_aa <- toupper(as.vector(before_alignment_aa))
    start_residue <- 1
    before_length <- as.character(length(before_alignment_aa))

    before_pdb_file <- generate_2d_pdb_file(before_alignment_aa, start_residue = start_residue)
    pdb_file_name2 <- paste(fasta_sequence_name,
                            '_before_alignment_sequence.pdb',sep='')
    list_of_pdb_file_names <- c(list_of_pdb_file_names,pdb_file_name2)
    write(before_pdb_file,pdb_file_name2,sep='')


  }

  if(end_subject != length(fasta_file[[fasta_sequence_name]])){
    after_alignment_aa <- fasta_file[[fasta_sequence_name]][(end_subject+1):length(fasta_file[[fasta_sequence_name]])]
    after_alignment_aa <- toupper(as.vector(after_alignment_aa))
    start_residue <- end_subject+1
    after_pdb_file <- generate_2d_pdb_file(after_alignment_aa,start_residue = start_residue)
    pdb_file_name2 <- paste(fasta_sequence_name,
                            '_after_alignment_sequence.pdb',sep='')
    list_of_pdb_file_names <- c(list_of_pdb_file_names,pdb_file_name2)
    write(after_pdb_file,pdb_file_name2,sep='')
  }


  renumbered_pdb_residues <- make_renumbered_pdb_vectors(resno_and_resid = resno_and_resid,
                                                         new_resno_and_resid = new_resno_and_resid,
                                                         maal = maal)

  #will need to reorder this in the future to be able to include the xyz
  #function in here that can be applied to the 2d generated PDB files
  pdb_resno_list <- renumbered_pdb_residues$pdb_resno_list
  inds <- atom.select(pdb_file,chain = c(chain),resno = pdb_resno_list[1]:pdb_resno_list[length(pdb_resno_list)])
  trimmed_pdb <- trim.pdb(pdb_file,inds = inds)
  trimmed_pdb$atom$resno <- renumbered_pdb_residues$new_pdb_resnos
  pdb_file_name <- paste(fasta_sequence_name,'___',pdb_id,'_',chain,".pdb",sep='')
  write.pdb(trimmed_pdb,pdb_file_name)
  list_of_pdb_file_names <- c(list_of_pdb_file_names,pdb_file_name)
  #accounting for before and after sequences as well

  #can return a list of the files generated in this function?
  return(list_of_pdb_file_names)
} #end function

#----Move file function----
#originally found from:
#https://stackoverflow.com/questions/10266963/moving-files-between-folders

my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

#create directories
#dir.create(file.path(mainDir, subDir))

#----Save plot as png----
#doesn't work?
save_plot_as_png <- function(plot_object,file_name){
  png(paste(file_name,'.png',sep=''))
  plot_object #print graph to "Plots" tab
  dev.off()
}

#----Get list of duplicated row matches----

get_list_of_duplicated_row_matches <- function(um_row,df_to_match_against){
  match_index <- 1
  row_match_result <- 0
  index_matches <- c()
  while(!is.na(row_match_result)){
    if(match_index <= nrow(df_to_match_against)){
      row_match_result <- row.match(um_row,df_to_match_against[match_index:nrow(df_to_match_against),])
    } else {
      break
    }
    if(!is.na(row_match_result)){
      match_index <- match_index + row_match_result
      index_matches <- c(index_matches,match_index - 1)
    }
  }
  return(index_matches)

}

#----Get list of index matches for dataframe----

make_list_of_index_matches <- function(df_to_match_against){

  if(FALSE %in% (rownames(df_to_match_against) == 1:nrow(df_to_match_against))){
    rownames(df_to_match_against) <- 1:nrow(df_to_match_against)
  }

  list_of_index_matches <- list()
  unique_match_rows <- unique(df_to_match_against[duplicated(df_to_match_against),])
  um_rownames <- rownames(unique_match_rows)
  for(um_row_name in um_rownames){
    um_row <- unique_match_rows[um_row_name,]
    index_matches <- get_list_of_duplicated_row_matches(um_row,df_to_match_against)
    #make this index more sophisticated later on?
    list_of_index_matches[[as.character(index_matches[1])]] <- index_matches
  }
  return(list_of_index_matches)
}

#----Find PDB files from FASTA and directory----

find_pdb_files_from_fasta_and_dir <- function(fasta_file,experiment_directory){
  #advanced iteration may include searching for directories with the name
  #of the fasta_name if the files were not found in the initial directory
  list_pdb_files_dir <- list.files(experiment_directory,pattern='\\.pdb$')
  found_pdb_files <- list()
  for(fasta_name in names(fasta_file)){
    list_matched_pdb <- list_pdb_files_dir[grepl(fasta_name, list_pdb_files_dir)]
    found_pdb_files[[fasta_name]] <- list_matched_pdb
  }
  return(found_pdb_files)
}

#----Make PDB files from FASTA----

make_pdb_files_from_fasta <- function(fasta_file, fasta_names_to_generate_all_2d_structures = NULL){

  list_of_all_pdb_files <- list()
  for(fasta_name in names(fasta_file)){

    #change to fasta_name %in% fasta_names_to_generate_all_2d_structures eventually
    if(fasta_name == 'AEBP2_Q6ZN18'){
      #generate 2D structure for this peptide

      fasta_sequence <- toupper(as.vector(fasta_file[[fasta_name]]))
      generated_pdb_file <- generate_2d_pdb_file(fasta_sequence,1)
      generated_pdb_name <- paste(fasta_name,'_generated_file.pdb',sep='')
      write(generated_pdb_file,generated_pdb_name)

      list_of_all_pdb_files[[fasta_name]] <- generated_pdb_name
      next #don't do rest of loop if it's this sequence
    }

    pdb_info <- go_through_menu_loops(fasta_file[[fasta_name]])
    pdb_read <- check_download_read_pdb(pdb_id = pdb_info$pdb_id)

    list_of_pdb_file_names <- find_and_generate_missing_aa_2d_structure(fasta_file = fasta_file,fasta_sequence_name = fasta_name,
                                                                        pdb_read,pdb_info$pdb_id,pdb_info$chain)

    list_of_all_pdb_files[[fasta_name]] <- list_of_pdb_file_names

  }

  return(list_of_all_pdb_files)

}

#----Make List of all PDBs (with Menu)----

make_list_of_all_pdbs <- function(fasta_file,
                                  experiment_directory = NULL,
                                  make_pdb_files = FALSE){

  if(make_pdb_files == FALSE){
    #put this in a while loop
    valid_directory <- FALSE
    new_directory <- experiment_directory

    while(valid_directory == FALSE){
      list_of_all_pdb_files <- find_pdb_files_from_fasta_and_dir(fasta_file,new_directory)

      if(length(list_of_all_pdb_files) == 0){
        #can also use this if not all fasta files are included

        menu_options <- c('Try a different directory','Make new PDB files')
        menu_title <- paste('Your chosen directory\n',new_directory,
                            '\ndid not turn up any PDB files matching the sequences in your FASTA file. ',
                            '\n\nWhat would you like to do?',
                            sep='')
        menu_selection <- menu(menu_options,title = menu_title)
        if(menu_selection == 1){ #try a different directory

          new_directory <- readline('Enter path: ')

          #input using
          #readline

        } else if(menu_selection == 2){ #make new PDB files
          list_of_all_pdb_files <- make_pdb_files_from_fasta(fasta_file)
          valid_directory <- TRUE
        }

      } else {
        valid_directory <- TRUE
      } #end if/else (length(list_of_all_pdb_files) == 0)

    } #end while


  } else { #end if(make_pdb_files == FALSE)

    list_of_all_pdb_files <- make_pdb_files_from_fasta(fasta_file)

  } #end else to if(make_pdb_files == FALSE)

  return(list_of_all_pdb_files)
}

#----Make Start/End PDB List----

#modifying this function for the single amino acid structures
#can't use the start and end in this way
#can do a boolean statement in this case
#if the file name ends in _single_point or something
#but then how would it get the right intervals?

#would have to get data from the missing amino acid function probably

make_start_end_pdb_list <- function(list_of_all_pdb_files){

  list_of_start_and_end_pdbs <- list()

  for(fasta_name in names(list_of_all_pdb_files)){
    for(pdb_name in list_of_all_pdb_files[[fasta_name]]){

      pdb_read <- read.pdb2(pdb_name)
      resno <- pdb_read$atom$resno
      resno_min <- min(resno)
      resno_max <- max(resno)
      chain_levels <- levels(factor(pdb_read$atom$chain))
      list_of_start_and_end_pdbs[[pdb_name]] <- c(resno_min,resno_max,chain_levels)

    }

  }

  return(list_of_start_and_end_pdbs)
}

#----Get PDB structure from Protein Position----

#slight error in getting the right PDB structure with the code
#accidentally getting the real PDB when it is one of the missing sequences
#start with missing sequences?

get_pdb_structure_from_protein_pos <- function(protein_name,protein_pos,
                                               list_of_all_pdb_files,
                                               list_of_start_and_end_pdbs){

  protein_name_pdb_files <- unique(list_of_all_pdb_files[[protein_name]])
  protein_name_pdb_files <- sort(protein_name_pdb_files,decreasing = TRUE)
  #create variable containing the
  on_this_pdb_structure <- NULL
  #may not need to go through each of the
  for(pdb_name in protein_name_pdb_files){
    start_pos <- list_of_start_and_end_pdbs[[pdb_name]][1]
    end_pos <- list_of_start_and_end_pdbs[[pdb_name]][2]
    chain <- list_of_start_and_end_pdbs[[pdb_name]][3]

    #check if protein_pos is between start and end
    if((as.numeric(protein_pos) >= as.numeric(start_pos)) && (as.numeric(protein_pos) <= as.numeric(end_pos))){
      on_this_pdb_structure <- pdb_name
      break #is this necessary?
    }

  } #end pdb_name in protein_name_pdb_files

  return(on_this_pdb_structure)

}

#----Get last character of string----

#from https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


#----Generate Single Atom PDB File----

generate_2d_pdb_file_single_atom <- function(amino_acid_fragments,start_residue,
                                             x_coord = 82.998, y_coord = -5.760,
                                             z_coord = -20.5){

  res_num <- start_residue
  num <- 1
  new_line <- ''
  new_pdb_file <- c()

  num_atom_spaces <- 7 - nchar(num)
  atom_spaces <- paste(rep(" ",num_atom_spaces),collapse = "")
  two_spaces <- "  "
  atom_number <- num
  residue_name <- toupper(aaa(amino_acid_fragments[num]))
  residue_number <- res_num
  six_spaces <- "      "
  num_residue_spaces <- 4 - nchar(residue_number)

  num_x_y_spaces <- 7 - nchar(y_coord)
  if(num_x_y_spaces > 0){
    x_y_spaces <- paste(rep(" ",num_x_y_spaces),collapse = "")
  } else {
    x_y_spaces <- ""
  }
  num_y_z_spaces <- 7 - nchar(z_coord)
  if(num_y_z_spaces > 0){
    y_z_spaces <- paste(rep(" ",num_y_z_spaces),collapse = "")
  } else {
    y_z_spaces <- ""
  }

  residue_spaces <- paste(rep(" ",num_residue_spaces),collapse = "")
  pdb_line <- paste('ATOM',atom_spaces,as.character(atom_number),two_spaces,
                    "CA",two_spaces,residue_name," A",residue_spaces,
                    as.character(residue_number),six_spaces,
                    as.character(format(x_coord,nsmall=3)),x_y_spaces,
                    as.character(format(y_coord,nsmall = 3)),
                    y_z_spaces, as.character(format(z_coord,nsmall=3)),two_spaces,
                    "1.00109.30           CA",new_line,
                    collapse = "",sep = "")

  new_pdb_file <- c(new_pdb_file,pdb_line)

  new_pdb_file <- paste(new_pdb_file,sep="", collapse="")

  return(new_pdb_file)

}



#----Generate PDB Files and Lists----

generate_pdb_lists_and_files_from_fasta <- function(fasta_file,pdb_suffix,fasta_names_to_generate_all_2d_structures = NULL){

  list_of_all_pdb_files2 <- list()
  list_of_start_and_end_pdbs2 <- list()
  generated_pdb_files <- list()



  for(fasta_sequence_name in names(fasta_file)){

    #in the future, this particular protocol should be included as part of the menu
    #loop as the generate 2D sequence option in the menu loop


    #have two lists that are then joined together


    #put in for loop

    generate_full_2d <- FALSE

    if(!is.null(fasta_names_to_generate_all_2d_structures)){
      for(fn_generate in fasta_names_to_generate_all_2d_structures){

        if(fasta_sequence_name == fn_generate){
          #generate 2D structure for this peptide
          #should there be a list of designated generated fasta sequence names that can be
          #inputted into the function

          fasta_sequence <- toupper(as.vector(fasta_file[[fasta_sequence_name]]))
          generated_pdb_file <- generate_2d_pdb_file_single_atom(fasta_sequence,1)
          generated_pdb_name <- paste(fasta_sequence_name,'_generated_file_',pdb_suffix,'.pdb',sep='')

          #copy the rest of the protocol here
          range_vector <- 1:length(fasta_sequence)

          write(generated_pdb_file,generated_pdb_name)

          #needs to be redone so it matches the rest
          generated_pdb_files[[generated_pdb_name]] <- range_vector

          list_of_all_pdb_files2[[fasta_sequence_name]] <- names(generated_pdb_files)

          for(pdb_file_name in names(generated_pdb_files)){

            resno_min <- min(generated_pdb_files[[pdb_file_name]])
            resno_max <- max(generated_pdb_files[[pdb_file_name]])
            chain_levels <- 'A'

            list_of_start_and_end_pdbs2[[pdb_file_name]] <- c(resno_min,resno_max,chain_levels)

          }

          #change to boolean?

          generate_full_2d <- TRUE

          #next #don't do rest of loop if it's this sequence
        }

      } #end for(fn_generate in fasta_names_to_generate_all_2d_structures)


    } #end if(!is.null)


    if(generate_full_2d == TRUE){
      next #if a full 2d structure has already been created, do not continue with the rest of the loop
    }


    #does the old loop where it
    cat(paste('Now finding PDB structures for:',fasta_sequence_name,'\n'))

    #just need to adjust the go_through_menu_loops so that it will first display
    #the pdb preference menu

    #doesn't work because does not output pdb_info
    pdb_info <- display_preferred_pdb_structure_menu(fasta_sequence_name, fasta_file)
    #pdb_info <- go_through_menu_loops(fasta_file[[fasta_sequence_name]])

    pdb_id <- pdb_info$pdb_id
    chain <- pdb_info$chain
    pdb_file <- check_download_read_pdb(pdb_id = pdb_id)
    pdb_file$atom <- pdb_file$atom[!is.na(a(firstup(pdb_file$atom$resid))),]



    #make the through menu loops output have some kind of indicator if the PDB has been generated by
    #a 2D function

    list_of_pdb_file_names2 <- find_and_generate_missing_aa_2d_structure_all(fasta_file = fasta_file,
                                                                             fasta_sequence_name = fasta_sequence_name,
                                                                             pdb_file = pdb_file, pdb_id = pdb_id,
                                                                             chain = chain, pdb_suffix = pdb_suffix)

    #list_of_pdb_file_names <- c(list_of_pdb_file_names,generated_pdb_files)

    list_of_all_pdb_files2[[fasta_sequence_name]] <- names(list_of_pdb_file_names2)

    for(pdb_file_name in names(list_of_pdb_file_names2)){

      pdb_read <- read.pdb2(pdb_file_name)
      resno_min <- min(list_of_pdb_file_names2[[pdb_file_name]])
      resno_max <- max(list_of_pdb_file_names2[[pdb_file_name]])
      chain_levels <- levels(factor(pdb_read$atom$chain))

      list_of_start_and_end_pdbs2[[pdb_file_name]] <- c(resno_min,resno_max,chain_levels)

    }


  } #end for(fasta_sequence_name in names(fasta_file))


  return(list(all_pdb_files=list_of_all_pdb_files2,
              start_end_chain=list_of_start_and_end_pdbs2))


} #end function generate_pdb_lists_and_files



#---- Find and Generate the Missing AA 2D Structure (All)-----

find_and_generate_missing_aa_2d_structure_all <- function(fasta_file,
                                                          fasta_sequence_name,
                                                          pdb_file,pdb_id,
                                                          chain,pdb_suffix){

  list_of_pdb_file_names <- list() #will contain the names of the files and the start and end

  #need to change this back to the other way of doing PWA to get the better results
  pwa_results <- quick_pwa_from_pdb(fasta_file[[fasta_sequence_name]],pdb_file,chain,
                                    use_resid_and_resno = TRUE)


  pwa_ranges <- get_pwa_ranges(pwa_results = pwa_results)
  pwa_strings <- get_pwa_strings(pwa_results = pwa_results)
  resno_and_resid <- get_resno_and_resid(pdb_read = pdb_file, chain = chain)
  new_resno_and_resid <- consolidate_resno_and_resid(resno_and_resid)

  maal <- renumber_and_identify_missing_aa(pwa_strings = pwa_strings, start_difference =   pwa_ranges$start_difference)

  aa_fragments <- get_missing_amino_acids(maal = maal, pwa_strings = pwa_strings)

  if(!is.null(aa_fragments)){ #proceed if there exists values in aa_fragments
    aa_count <- 0
    for(amino_acid_fragments in aa_fragments$amino_acid_fragments){
      aa_count <- aa_count + 1
      start_residue <- aa_fragments$start_fragment_num[aa_count]
      pdb_file_name2 <- paste(fasta_sequence_name,
                              '_missing_sequence_',
                              as.character(aa_count),
                              "_",pdb_suffix,
                              '.pdb',sep = '')
      #change the xyz coordinates for each of the dots

      if(pdb_suffix == 'single_dot'){
        new_pdb_file <- generate_2d_pdb_file_single_atom(amino_acid_fragments,start_residue)
      } else if(pdb_suffix == 'curved'){

      } else {
        #linear

      }


      range_vector <- start_residue:(start_residue+(aa_fragments$start_fragment[[aa_count]][2] - 1))

      #list_of_pdb_file_names <- c(list_of_pdb_file_names,pdb_file_name2)
      list_of_pdb_file_names[[pdb_file_name2]] <- range_vector
      write(new_pdb_file,pdb_file_name2,sep='')

    }
  } else {
    cat('There are no missing peptides to generate "missing_peptide" 2D files\n\n')
  }


  start_subject <- pwa_ranges$start_subject
  end_subject <- pwa_ranges$end_subject

  if(start_subject != 1){
    before_alignment_aa <- fasta_file[[fasta_sequence_name]][1:start_subject-1]
    before_alignment_aa <- toupper(as.vector(before_alignment_aa))
    start_residue <- 1
    before_length <- as.character(length(before_alignment_aa))


    if(pdb_suffix == 'single_dot'){
      before_pdb_file <- generate_2d_pdb_file_single_atom(before_alignment_aa, start_residue = start_residue)
    } else if(pdb_suffix == 'curved'){ #replace with
      #before_pdb_file <- generate_2d_pdb_file_single_atom(before_alignment_aa, start_residue = start_residue)
    } else {
      #linear
      #before_pdb_file <- generate_2d_pdb_file_single_atom(before_alignment_aa, start_residue = start_residue)
    }

    pdb_file_name2 <- paste(fasta_sequence_name,
                            '_before_alignment_sequence_',pdb_suffix,'.pdb',sep='')

    range_vector <- start_residue:(start_residue+(as.numeric(before_length) - 1))
    list_of_pdb_file_names[[pdb_file_name2]] <- range_vector
    write(before_pdb_file,pdb_file_name2,sep='')


  }


  if(end_subject != length(fasta_file[[fasta_sequence_name]])){
    after_alignment_aa <- fasta_file[[fasta_sequence_name]][(end_subject+1):length(fasta_file[[fasta_sequence_name]])]
    after_alignment_aa <- toupper(as.vector(after_alignment_aa))
    start_residue <- end_subject+1

    if(pdb_suffix == 'single_dot'){
      after_pdb_file <- generate_2d_pdb_file_single_atom(after_alignment_aa,start_residue = start_residue)
    } else if(pdb_suffix == 'curved'){
      #curved function
    } else {
      #linear
    }


    pdb_file_name2 <- paste(fasta_sequence_name,
                            '_after_alignment_sequence_',pdb_suffix,'.pdb',sep='')

    range_vector <- start_residue:length(fasta_file[[fasta_sequence_name]])
    list_of_pdb_file_names[[pdb_file_name2]] <- range_vector
    write(after_pdb_file,pdb_file_name2,sep='')
  }



  renumbered_pdb_residues <- make_renumbered_pdb_vectors(resno_and_resid = resno_and_resid,
                                                         new_resno_and_resid = new_resno_and_resid,
                                                         maal = maal)

  #return(new_resno_and_resid)


  #will need to reorder this in the future to be able to include the xyz
  #function in here that can be applied to the 2d generated PDB files
  pdb_resno_list <- renumbered_pdb_residues$pdb_resno_list
  inds <- atom.select(pdb_file,chain = c(chain),resno = pdb_resno_list[1]:pdb_resno_list[length(pdb_resno_list)])
  trimmed_pdb <- trim.pdb(pdb_file,inds = inds)
  #trimmed_pdb$atom <- trimmed_pdb$atom[!is.na(a(firstup(trimmed_pdb$atom$resid))),]
  trimmed_pdb$atom$resno <- renumbered_pdb_residues$new_pdb_resnos
  #trimmed_pdb$atom <- trimmed_pdb$atom[!is.na(trimmed_pdb$atom$resno),]
  #clean.pdb?
  pdb_file_name <- paste(fasta_sequence_name,'___',pdb_id,'_',chain,".pdb",sep='')
  #need to remove NAs

  #trimmed_pdb <- read.pdb2('H3___1KX5_A.pdb')
  #trimmed_pdb$atom <- trimmed_pdb$atom[!is.na(trimmed_pdb$atom$resno),]
  write.pdb(trimmed_pdb,pdb_file_name)
  resno <- trimmed_pdb$atom$resno[!is.na(trimmed_pdb$atom$resno)]
  resno_min <- min(resno)
  resno_max <- max(resno)
  range_vector <- resno_min:resno_max
  list_of_pdb_file_names[[pdb_file_name]] <- range_vector
  #accounting for before and after sequences as well

  #can return a list of the files generated in this function?
  cat(paste('Generated',length(names(list_of_pdb_file_names)),'PDB files for',fasta_sequence_name,'\n\n'))

  return(list_of_pdb_file_names)
} #end all function






#----Differential Plink Analysis Function----

#'BS3 Differential Analysis Inclusive Function
#'
#'This function combines many of the features of the BS3 XL-MS workflow for ease of use
#'@param list_of_files A list of file names containing pLink data to be used for the analysis
#'@param fasta_file FASTA file as loaded by seqinr::read.fasta()
#'@param file_type_2d The type of file used for PDB creation. "single_dot" creates PDB files with peptides represented by single dots.
#'@param fasta_names_to_generate_all_2d_structures List of protein names corresponding to names in FASTA file to automatically make into 2D PDB structures. Defaults to NULL.
#'@param csv_pdb_input_file csv_pdb_input_file
#'@param category_color_input_file category_color_input_file
#'@param show_only_real_structures show_only_real_structures
#'@param distance_histogram_name distance_histogram_name
#'@param pymol_file_list_file_name pymol_file_list_file_name
#'@param xlink_df_name xlink_df_name
#'@param xlink_view_file_name xlink_view_file_name
#'@param no_pdb_files no_pdb_files
#'@param protein_alternative_names_dict protein_alternative_names_dict
#'@param pdb_directory pdb_directory
#'@param data_input_type data_input_type
#'@export

ppi.analyze <- function(list_of_files,fasta_file,
                        file_type_2d = c('single_dot','linear','curved'),
                        fasta_names_to_generate_all_2d_structures = NULL,
                        csv_pdb_input_file = NULL,
                        category_color_input_file = NULL,
                        show_only_real_structures = NULL,
                        distance_histogram_name = 'distance_histogram.png',
                        pymol_file_list_file_name = 'pymol_ts_output.pml',
                        xlink_df_name = 'xlink_dataframe.csv',
                        xlink_view_file_name = 'xlink_view_df.csv',
                        no_pdb_files = FALSE,
                        protein_alternative_names_dict = NULL,
                        pdb_directory = NULL,
                        data_input_type = 'plink',
                        pdb_numbering = FALSE,
                        pdb_match_vector = NULL){

  if((!is.null(category_color_input_file)) && (typeof(category_color_input_file) == 'character')){
    category_color_input_file <- read.csv(category_color_input_file)
  }


  if(!is.null(pdb_directory)){
    experiment_directory <- pdb_directory
  } else {
    experiment_directory <- getwd()
  }

  xlink_list <- list()

  list_of_protein_names <- names(fasta_file)

  if(!is.null(protein_alternative_names_dict)){

    protein_alternative_names_dict <- read.csv(protein_alternative_names_dict)

    for(pand_name in names(protein_alternative_names_dict)){
      #add entire column to list

      pand_list <- as.character(protein_alternative_names_dict[,pand_name])
      list_of_protein_names <- c(list_of_protein_names,pand_list)

    }

    #add to list of protein names

  }

  list_of_protein_names <- unique(list_of_protein_names)

  current_directory <- getwd()
  Sys.setenv(TZ='Australia/Melbourne')
  #create a system log file that has sys time on it instead?
  #what about the %03d?
  new_directory_title <- paste('make_diff_analysis Run on',Sys.time())
  new_directory_title <- gsub(':','-',new_directory_title)

  new_directory <- paste(current_directory,'/',new_directory_title,sep = '')
  dir.create(new_directory)

  pdb_suffix <- file_type_2d

  # if(file_type_2d == 'single dot'){
  #   pdb_suffix <- 'single_dot'
  # } else {
  #   pdb_suffix <- file_type_2d
  # }



  if(is.null(csv_pdb_input_file)){

    #should have the boolean here that will be able to renumber the plink numbering from the
    #PDB identifier?

    #should bypass this function altogether and create a new one if the user chooses to create a
    #new file

    if(pdb_numbering == TRUE){

      #do PDB alignment?
      #how to do it for each of the proteins?
      #have a mega list for each of the proteins with 2 vectors each

      if(is.null(pdb_match_vector)){
        pdb_match_vector <- ppi.alignPDB(fasta_file = fasta_file)
      }


    } else {
      #pdb_numbering == FALSE --> do what has been done before
      generated_pdb_lists <- generate_pdb_lists_and_files_from_fasta(fasta_file,pdb_suffix,
                                                                     fasta_names_to_generate_all_2d_structures)

      make_start_end_pdb_df_output2(generated_pdb_lists = generated_pdb_lists)

    }



  } else { #if someone has indicated a PDB file to be used

    generated_pdb_lists <- generate_pdb_lists_from_pdb_csv(csv_pdb_input_file)
    #need same input as list of all pdb files
  }

  #return(generated_pdb_lists)


  #if pdb_numbering == TRUE --> will need a substitute for these or otherwise ignore them
  if(pdb_numbering == FALSE){
    list_of_all_pdb_files <- generated_pdb_lists$all_pdb_files
    list_of_start_and_end_pdbs <- generated_pdb_lists$start_end_chain
  } else {
    #end if(pdb_numbering == FALSE)
    list_of_all_pdb_files <- pdb_match_vector
  }

  #return(generated_pdb_lists)

  #generation of this list within the function? or different input variable

  #make the PDBs with the fasta file?
  #function that will produce the two lists

  #
  #

  #make the pdb suffix based on the file type designation
  #remove and replace the space with underscore?

  not_included_xl_sites <- c()

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



  plink_files_already_run <- c()

  # typeof(list_of_files)
  # typeof(list())
  #
  # for(index_num in list(num1=c(1,2,3),num2=c(1,2,3))){
  #   print(index_num)
  # }

  plink_file_index <- 0

  for(plink_file in list_of_files){

    plink_file_index <- plink_file_index + 1

    #need to update this so it will take into account which type of plink file is being uploaded

    #need to check here if plink_file is the name of the file or the actual data

    # if(typeof(plink_file) == 'character'){ #name of file that needs to be loaded
    #
    #   #get files to be loaded earlier in the sequence?
    #
    #
    # } else if(typeof(plink_file) == 'list'){ #dataframe?
    #
    # }

    plink_file_name <- plink_file

    if(typeof(plink_file_name) == 'list'){
      plink_file_name <- names(list_of_files)[plink_file_index]
    }




    if(!(plink_file_name %in% plink_files_already_run)){


      #can put data_input_type == 'proxl'

      if(data_input_type == 'proxl'){

        xlink_list <- load_proxl_data(proxl_data = plink_file)

        plink_files_already_run <- c(plink_files_already_run,plink_file_name)

      }


      #if data_input_type == 'plink'

      if(data_input_type == 'plink'){

        #need to check here if it is a character or a list here since endsWith won't work here
        #way to tell which data structure it is from the lines?

        if(endsWith(plink_file_name,'cross-linked_peptides.csv') || endsWith(plink_file_name,'cross-linked_peptides.xlsx') || endsWith(plink_file_name,'cross-linked_peptides')){

          xl_experiment_name <- strsplit(plink_file_name,'.filtered_cross-linked_peptides.csv')[[1]]
          #loop_linked_name <- paste(xl_experiment_name,'.filtered_loop-linked_peptides.csv',sep='')

          #have a boolean here in case there is onl

          crosslinked_peptides_file_name <- plink_file_name
          #looplinked_peptides_file_name <- loop_linked_name

          crosslinked_peptides <- readLines(crosslinked_peptides_file_name)
          #looplinked_peptides <- readLines(looplinked_peptides_file_name)

          #can probably make the crosslinked master list not have to import the list of protein names
          #since it doesn't use it
          crosslinked_master_list <- make_plink2_master_list(plink2_peptides = crosslinked_peptides,
                                                             list_of_protein_names = list_of_protein_names,
                                                             plink2_type_of_file = 'cross-linked',
                                                             master_list_index = 0)

          # looplinked_master_list <- make_plink2_master_list(plink2_peptides = looplinked_peptides,
          #                                                   list_of_protein_names = list_of_protein_names,
          #                                                   plink2_type_of_file = 'loop-linked',
          #                                                   master_list_index = 0)
          #

          #plink2_master_list <- c(crosslinked_master_list,looplinked_master_list)
          #should do grepl or just loop-linked peptides?


          plink_files_already_run <- c(plink_files_already_run,plink_file_name)

          xlink_list <- crosslinked_master_list



          #cross-link loop
          #see if there is a loop-link file to add
          #add the two together

        } else if(endsWith(plink_file_name,'loop-linked_peptides.csv') || endsWith(plink_file_name,'loop-linked_peptides.xlsx')){

          #plink_file_name <- '/Users/emmagail/Downloads/XL-MS Xlink Visual/BS3 - 2 - PRC2-AEBP2/BS3 - 2 - PRC2-AEBP2 - 2018.03.28_1.filtered_loop-linked_peptides.csv'

          #read.csv(plink_file_name)

          looplinked_peptides <- readLines(plink_file_name)

          looplinked_master_list <- make_plink2_master_list(plink2_peptides = looplinked_peptides,
                                                            list_of_protein_names = list_of_protein_names,
                                                            plink2_type_of_file = 'loop-linked',
                                                            master_list_index = 0)

          plink_files_already_run <- c(plink_files_already_run,plink_file_name)

          xlink_list <- looplinked_master_list

        } else { #assume the file is a plink1 file

          #do the original plink load function

          xlink_list <- load_plink_file(plink_file_name)

          plink_files_already_run <- c(plink_files_already_run,plink_file_name)
        }


      } #end if data_input_type == 'plink


    } else {
      #end if(!(plink_file_name in plink_files_already_run)
      next
    }

    #return(xlink_list)
  #} #end of for(plink_file in list_of_files)


    #return(xlink_list)

    #xlink_list <- load_plink_file(plink_file)

    #return(plink_file_name)

    #return(category_color_input_file)

    for(xlink_index in 1:length(xlink_list)){


      # if(xlink_index == 17){
      #
      #   return(list(proteins_xlink=xlink_list[[xlink_index]]))
      #
      # }


      protein_split_list <- NULL
      seq_xlink <- NULL


      if(xlink_index > length(xlink_list)){
        break
      }




      seq_xlink <- xlink_list[[xlink_index]]$sequence_xlink

      if((length(xlink_list[[xlink_index]]$proteins_xlink) != 0)){
        seq_xlink <- strsplit(seq_xlink,':0')[[1]]
      }


      #pro_xlink <- xlink_list[[xlink_index]]$proteins_xlink[1]
      #this may not necessarily be true if the first is an unusable protein name
      #should switch this so that the loop happens right at the beginning

      #account for the '/' in the original pLink files that were not separated

      #return(xlink_list)

      if((length(xlink_list[[xlink_index]]$proteins_xlink) != 0) && grepl("/",xlink_list[[xlink_index]]$proteins_xlink)){

        xlink_list[[xlink_index]]$proteins_xlink <- strsplit(xlink_list[[xlink_index]]$proteins_xlink,'/')[[1]]

      }

#
#       if((length(xlink_list[[xlink_index]]$proteins_xlink) == 1) && grepl("/",xlink_list[[xlink_index]]$proteins_xlink)){
#
#         xlink_list[[xlink_index]]$proteins_xlink <- strsplit(xlink_list[[xlink_index]]$proteins_xlink,'/')[[1]]
#
#       }

      #xlink_list[[xlink_index]]$proteins_xlink <- strsplit(xlink_list[[xlink_index]]$proteins_xlink,'/')[[1]]




      for(protein_xlink_in_list in xlink_list[[xlink_index]]$proteins_xlink){

        #JAZF1-SUZ12 is extracted but
        protein_split_list <- extract_protein_name_and_peptide_num(split_proteins = protein_xlink_in_list,
                                                                   list_of_protein_names = list_of_protein_names)

        if(!is.null(protein_split_list)){
          pro_output <- protein_split_list
          #pro_xlink <- paste(pro_output[[1]][1],'(',pro_output[[1]][2],')-',pro_output[[2]][1],'(',pro_output[[2]][2],')',sep='')

          #protein_split_list <- extract_protein_name_and_peptide_num('PHF19_Q5T6S3(25)-PHF19_Q5T6S3(32)',list_of_protein_names)

          new_protein_list <- c()

          for(pro_index1 in 1:length(protein_split_list)){

            protein_name <- protein_split_list[[pro_index1]][1]


            protein_name <- check_and_change_protein_name(protein_name,protein_alternative_names_dict,list_of_all_pdb_files,fasta_file)

            protein_split_list[[pro_index1]][1] <- protein_name
            new_protein_list <- c(new_protein_list,protein_name)

          }


          pro_xlink <- paste(new_protein_list[1],'(',pro_output[[1]][2],')-',new_protein_list[2],'(',pro_output[[2]][2],')',sep='')
          #change the names in in protein_split_list


          #go through

          #check_and_change_protein_name(protein_name,protein_alternative_names_dict,list_of_all_pdb_files,fasta_file)


          #protein name hasn't been changed at this point, need to move dictionary check to
          #up here

          #this needs to be changed to the changed version of the pro_xlink instead
          #of the previous version
          break
        }

      } #end for(protein_xlink_in_list in xlink_list[[xlink_index]]$proteins_xlink)

      # protein_split_list <- extract_protein_name_and_peptide_num(split_proteins = xlink_list[[xlink_index]]$proteins_xlink[1],
      #                                                            list_of_protein_names = list_of_protein_names)
      #




      if(is.null(protein_split_list)){

        #write to file

        nixl_site <- paste(xlink_list[[xlink_index]]$proteins_xlink[1],' - ',plink_file)
        not_included_xl_sites <- c(not_included_xl_sites,nixl_site)
        next
      }


      score_xlink <- min(as.numeric(xlink_list[[xlink_index]]$score_xlink))

      #change the seq_xlink to the common name up here?
      #switch the order of some of the following so that it makes more sense

      print(pro_xlink)

      #once the protein name has been changed, change pro_xlink as well so that it matches


      if(is.null(pro_xlink) || is.null(seq_xlink)){
        next
      }

      #print(xlink_index)


      seq_index <- match(seq_xlink,seq_xlink_list) #if == NA, means that it has not been added to the list yet
      pro_index <- match(pro_xlink,pro_xlink_list)

      print(seq_index)
      print(pro_index)

      #NEED TO CHECK if only 1 is
      if(length(unique(is.na(c(seq_index,pro_index)))) == 2){
        #both TRUE and FALSE exist in the list
        #if both are true --> probably should re-add to the list

        seq_index <- NA
        pro_index <- NA

        #should there be a note added?

      }
      #if == 2 --> only 1 shows up so probably should re-add the full thing to
      #


      # if(is.na(seq_index) || is.na(pro_index)){
      #
      #   cat('Error in finding index in sequence list, sending to not included list\n')
      #   not_included_xl_sites <- c(not_included_xl_sites,pro_xlink)
      #
      # }



      if(seq_index != pro_index && !is.na(seq_index) && !is.na(pro_index)){

        if(pro_xlink_list[seq_index] == pro_xlink){

          pro_index <- seq_index

        } else if(seq_xlink_list[pro_index] == seq_xlink){

          seq_index <- pro_index

        } else {

          cat(paste('Index error for',seq_xlink,'and',pro_xlink))
        }

      }


      if((seq_index == pro_index) && !is.na(seq_index)){ #make sure that the indeces match and are real numbers
        #indiciates that they have already been added to the lists respectively

        #check if the file has not already been checked?

        if(!grepl(plink_file_name,files_xlink_list[seq_index])){

          #need to change this so it will change it to the right frequency/color based on
          #plink_file_name

          if(!is.null(category_color_input_file)){
            #if it is not null, will need to change the color based on the
            new_color <- as.character(category_color_input_file[category_color_input_file[['file_name']] == plink_file_name,][['color']])
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

          #will also have to have a contingency here for the plink_file name if a list and
          #not a character



          files_xlink_list[seq_index] <- paste(files_xlink_list[seq_index],'+',plink_file_name,sep='')


        }

      } else { #need to double check to make sure that both are not NA?
        #have not been added to the list
        #seq_xlink_list <- c(seq_xlink_list,seq_xlink)
        #pro_xlink_list <- c(pro_xlink_list,pro_xlink)
        #need to keep track of index or if everything works well, that will be taken care of already?

        seq_xlink2 <- seq_xlink

        sequence_xlink_list <- split_sequences2(xlink_list[[xlink_index]]$sequence_xlink, proteins = FALSE)
        #protein_split_list <- split_sequences2(xlink_list[[xlink_index]]$proteins_xlink[1], proteins = TRUE)

        #should go through the protein split list until there is a valid hit?

        # for(protein_xlink_in_list in xlink_list[[xlink_index]]$proteins_xlink){
        #
        #   protein_split_list <- extract_protein_name_and_peptide_num(split_proteins = protein_xlink_in_list,
        #                                                              list_of_protein_names = list_of_protein_names)
        #
        #   if(!is.null(protein_split_list)){
        #     break
        #   }
        #
        # }
        # # protein_split_list <- extract_protein_name_and_peptide_num(split_proteins = xlink_list[[xlink_index]]$proteins_xlink[1],
        # #                                                            list_of_protein_names = list_of_protein_names)
        # #
        # if(is.null(protein_split_list)){
        #
        #   #write to file
        #   not_included_xl_sites <- c(not_included_xl_sites,pro_xlink)
        #   next
        # }



        otps_list <- c()
        xyz_coord_list <- list()

        #keep these in lists for next part of analysis?
        #if yes --> make into a function so better condense the information
        for(index_num in 1:length(sequence_xlink_list)){ #go through peptides (1 and 2)
          seq_xlink <- sequence_xlink_list[[index_num]]
          pro_split <- protein_split_list[[index_num]]
          protein_name <- pro_split[1]
          protein_pos <- as.numeric(pro_split[2])
          peptide_seq <- seq_xlink[1]
          peptide_pos <- as.numeric(seq_xlink[2])


          if(!(protein_name %in% names(list_of_all_pdb_files))){ #does it have a PDB file?
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

                    cat(paste(protein_name,'changed to',as.character(protein_alternative_names_dict[row_num,'original_name']),
                              'in output\n'))

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


              cat('Is in fasta file but protein name has not been changed yet\n')

              #ask if the user wants to create a separate PDB file for it

              #menu loops function here

            } else {
              #the user has not included it in neither the fasta file nor the dictionary
              #add the whole crosslinking site to a dataframe to be
              #outputted if the xl site does not meet the requirements

              cat('Not in fasta file or in the dictionary\n')


            }

            #should also check protein_pos?

            #can move the is.null if statement within this if statement
            #since it might be helpful to check beforehand if there is a PDB file for
            #the protein

          } #end if(!(protein_name %in% names(list_of_all_pdb_files)))


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
          if(pdb_numbering == TRUE){

            protein_pos_match <- match(protein_pos,pdb_match_vector[[protein_name]]$fasta)

            #how to keep track of the new protein position within the
            pdb_protein_pos <- as.numeric(pdb_match_vector[[protein_name]]$pdb[protein_pos_match])

            #pdb_match_vector[[protein_name]]$chain[protein_pos_match]

            protein_split_list[[index_num]][2] <- pdb_protein_pos

            on_this_pdb_structure <- pdb_match_vector[[protein_name]]$chain[protein_pos_match]


            #if it == '-' --> make it NULL to match the rest of the code


            #need to add to otps_list in order for this to work

            #on_this_pdb_structure
            otps_list <- c(otps_list,on_this_pdb_structure)



          } else {



            on_this_pdb_structure <- get_pdb_structure_from_protein_pos(protein_name,
                                                                        protein_pos,
                                                                        list_of_all_pdb_files,
                                                                        list_of_start_and_end_pdbs)



            if(!is.null(on_this_pdb_structure)){
              otps_list <- c(otps_list,on_this_pdb_structure)
            } #end if(!is.null(on_this_pdb_structure)){


          } #end else to if(pdb_numbering == TRUE){

          #return(list(all_pdb=list_of_all_pdb_files,start_and_end=list_of_start_and_end_pdbs))

          #add protein positions to another list?

          #cat('On this PDB Structure\n')
          #print(on_this_pdb_structure)



        } #end for(index_num in 1:length(sequence_xlink_list))

        if((length(otps_list) != 2) || ('-' %in% otps_list)){

          #not_included_xl_sites <- c(not_included_xl_sites,seq_xlink2)
          nixl_site <- paste(pro_xlink,' - ',plink_file)
          not_included_xl_sites <- c(not_included_xl_sites,nixl_site)
          cat('Skipping this crosslinking site, adding to list of XL sites not listed')
          next
          #cat('Skipping this crosslinking site, adding to list of XL sites not listed')
          #add crosslinking site to list of not included crosslinking sites?
          #list or dataframe?

          #not_included_xl_sites <- c(not_included_xl_sites,seq_xlink2)

        } #end if(length(otps_list) != 2)

        #cat('OTPS List')
        #print(otps_list)

        pdb1_xlink_list <- c(pdb1_xlink_list,otps_list[1])
        pdb2_xlink_list <- c(pdb2_xlink_list,otps_list[2])

        #print(pdb1_xlink_list)
        #print(pdb2_xlink_list)
        #print(grepl('___',otps_list))


        #this only works when pdb_numbering == FALSE
        #can add an or statement
        #pdb_numbering must be true --> is that it?
        if((!(FALSE %in% grepl('___',otps_list))) || pdb_numbering == TRUE){ #make sure that it is all TRUEs to go forward with dist calculation
          #calculate the distance here for the real distances

          #need to get/load the correct PDB file
          #use the protein position numbers to get the right XYZ coordinates

          #use indices instead of going through the list
          #and use something like this instead
          #sequence_xlink_list[[c(2,2)]]

          #return(otps_list)

          for(otps_index in 1:length(otps_list)){ #need to load the protein position as well
            #otps_index <- 1
            otps <- otps_list[otps_index]
            #print(otps)


            if(pdb_numbering == TRUE){

              #can check the temporary directory for the PDB file
              #see if pdb_id.pdb exists within the tempdir
              #if yes, can
              otps_split <- strsplit(otps,'_')[[1]]
              pdb_id <- otps_split[1]
              chain <- otps_split[2]


              if(is.na(pdb_id)){
                break
              }

              #print(pdb_id)
              if(paste0(pdb_id,'.pdb') %in% list.files(tempdir())){

                pdb_read <- read.pdb2(paste0(tempdir(),'/',pdb_id,'.pdb'))

              } else {

                pdb_read <- read.pdb2(pdb_id)

              } #end if(paste0(pdb_id,'.pdb') %in% list.files(tempdir()))


              #once the pdb has been read --> get just the chain?

              pdb_read$atom <- pdb_read$atom[pdb_read$atom$chain == chain,]

              #use the read.pdb function to be able to

            } else { #end if(pdb_numbering == TRUE)

              if(!is.null(pdb_directory)){
                #paste the directory to otps
                if(!endsWith('/',pdb_directory)){
                  pdb_directory <- paste(pdb_directory,'/',sep='')
                }

                pdb_path <- paste(pdb_directory,otps,sep='')
                pdb_read <- read.pdb2(pdb_path)

              } else {
                pdb_read <- read.pdb2(otps)
              }


            } #end  else to if(pdb_numbering == TRUE)



            xyz_coords <- get_xyz_coordinates_pdb(pdb_read)



            #return(protein_split_list)
            protein_pos <- protein_split_list[[c(otps_index,2)]]

            #if pdb_numbering == TRUE, this may need to be changed slightly
            pos_match <- match(protein_pos,pdb_read$atom$resno)

            if(!is.na(pos_match)){ #checking if there is a true match
              #there is a true match

              #loop to get the right atom
              #check to see if at that atom == 'CA' or whatever

              #put in while loop so that it will keep increasing
              #

              atom_match <- pdb_read$atom$elety[pos_match]
              while(atom_match != 'CA'){ #can choose any atom here, can change to N later on
                #make 'CA' as a variable when this is made into a function
                #with 'CA' as potential default in function
                pos_match <- pos_match + 1
                atom_match <- pdb_read$atom$elety[pos_match]
              }

              atom_match
              pos_match

              if(pdb_read$atom$resno[pos_match] == protein_pos){

                #should anything be done in here?


              } else {
                #what to do if it has been overshot
                cat('No match to chosen atom, reverting to first match')
                pos_match <- match(protein_pos,pdb_read$atom$resno)

              }

              #then at end, double-check to make sure that the protein_pos is still the same
              #and didn't overshoot


              #then can use pos_match to get the right xyz coordinates

            } else {
              #if it is NA and there is no match
              #what is the recourse here?

            }

            #can use the following to check what atom it is
            #pdb_read$atom$elety
            #if doing a cycle, need to confirm that it is still the right amino acid number

            #only got the first atom? match() returns the first match
            #check which atom it is before continuing
            #need to select the CA or something

            xyz_matches <- c()
            #make list of the xyz coordinates
            for(xyz_name in names(xyz_coords)){
              xyz1 <- xyz_coords[[xyz_name]][pos_match]
              xyz_matches <- c(xyz_matches,xyz1)
            }

            xyz_coord_list[[otps_index]] <- xyz_matches


          } #end of otps_index in 1:length(opts_list) for loop


          #calculate distance here

          #double-check formula before using in the final version of distance
          if(length(xyz_coord_list) == 2){
            dist_xyz <- dist.xyz(xyz_coord_list[[1]],xyz_coord_list[[2]])[1,1]
            dist_xlink_list <- c(dist_xlink_list,dist_xyz)
          } else {
            dist_xlink_list <- c(dist_xlink_list,NA)
          }



          #making note of numbering?
          #can be used as the row index?
          #cmd.set_name(string old_name, string new_name) --> pymol command that can be used to set new name

          #add to list of distances for the matrix

          #use bio3d functions?


          #get xyz coordinates from protein position
          #use of index or can use residue number?

          #formula for xyz coordinates


        } else { #one or more PDB structures is not "real" so distance should not be calculated


          dist_xlink_list <- c(dist_xlink_list,NA)

          #add NA or some other value to dist_xlink_list

        }

        if(!is.null(category_color_input_file)){
          #need to make the color match the category
          new_color <- as.character(category_color_input_file[category_color_input_file[['file_name']] == plink_file_name,][['color']])
          freq_color_xlink_list <- c(freq_color_xlink_list,new_color)

        } else {

          freq_color_xlink_list <- c(freq_color_xlink_list,frequency_color_list[1])

        }

        freq_xlink_list <- c(freq_xlink_list,1)

        if(typeof(plink_file) == 'list'){
          #if list --> need to have some kind of character or something to replace it instead
          #just have a dummy name for now?
        }
        files_xlink_list <- c(files_xlink_list,plink_file_name)

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


      } #!!!!!ELSE ENDS HERE



    } #end for(xlink_index in 1:length(xlink_list))

  } #end of for(plink_file in list_of_files)

  #theoretically end list_of_files for loop here


  #make data.frame here out of the many lists

  xlink_mega_list <- list(seq=seq_xlink_list,
                          pro=pro_xlink_list,
                          dist=dist_xlink_list,
                          freq=freq_xlink_list,
                          freq_color=freq_color_xlink_list,
                          files=files_xlink_list,
                          pdb1=pdb1_xlink_list,
                          pdb2=pdb2_xlink_list,
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

  for(row_num in 1:nrow(xlink_mega_df)){
    filtered_pdb_xl_df_row <- xlink_mega_df[row_num,]

    #make this into a function
    pdb1 <- as.character(filtered_pdb_xl_df_row$pdb1)
    pdb2 <- as.character(filtered_pdb_xl_df_row$pdb2)

    if((!grepl('___',pdb1)) || (!(grepl('___',pdb2)))){

      next
    }

    pdb1 <- strsplit(pdb1,'___')[[1]][2]
    pdb1 <- strsplit(pdb1,'_')[[1]][1]

    pdb2 <- strsplit(pdb2,'___')[[1]][2]
    pdb2 <- strsplit(pdb2,'_')[[1]][1]

    if(pdb1 != pdb2){
      #reassign distance to NA
      xlink_mega_df[row_num,][['dist']] <- NA
    }

  }


  if(!is.null(category_color_input_file)){
    #read file based on the extension
    #file_input <- read.csv(category_color_input_file)

    #category_colors <- get_category_colors_from_input_file(xlink_mega_df,file_input = file_input)

    #xlink_mega_df$freq_color <- category_colors

    freq_color_pretty_name <- "Category Color"
  } else {
    freq_color_pretty_name <- "Frequency Color"
  }

  #print(freq_color_pretty_name)

  pretty_xlink_colnames <- c('Sequence',
                             'Protein',
                             'Distance',
                             'Frequency',
                             freq_color_pretty_name,
                             'Files',
                             'PDB1',
                             'PDB2',
                             'Protein Position 1',
                             'Protein Position 2',
                             'Protein Name 1',
                             'Protein Name 2',
                             'Peptide Sequence 1',
                             'Peptide Sequence 2',
                             'Peptide Position 1',
                             'Peptide Position 2',
                             'Score')

  #need to write the mega_df to a csv
  #write csv file as output
  #re-order and make with nice headers
  xlink_mega_df_export <- xlink_mega_df

  #bs3_xlink_function_output[order(bs3_xlink_function_output$freq),]
  xlink_mega_df_export <- xlink_mega_df_export[order(xlink_mega_df_export$freq, decreasing = TRUE),]

  colnames(xlink_mega_df_export) <- pretty_xlink_colnames
  #let user choose the name of the dataframe?
  #can have multiple options for this so that people can choose the columns that they want
  #in the export of the csv file
  #as default will export all
  #can have users input a list, default can be null or something
  #only_selected_columns = c('PDB1','PDB2')
  #could have a warning if the user only puts in 1 but not 2 or excessive?
  #or not have the 1 or 2 at all, the user only inputs 'PDB' and it will output both
  #should the user be able to write the file as an xlsx file as well

  setwd(new_directory)
  write.csv(xlink_mega_df_export,file=xlink_df_name)

  if(length(not_included_xl_sites) != 0){

    #should probably update the input file name if this is actually included
    #in the output
    write(not_included_xl_sites,file = 'not_inlcuded_sites.txt')

  }

  #should order by frequency (does not currently do that though it looks like it because of test data sets)

  #pretty_xlink_colnames <- c()
  #duplicate csv file first and then rename columns before exporting


  #Distance Histogram
  #enclose in a tryCatch() or just ensure that ggplot2 is installed?
  qplot(xlink_mega_df$dist,
        main = "Distance Histogram",
        xlab = "Distance",
        ylab = "Frequency")

  ggsave(distance_histogram_name)
  setwd(current_directory)

  ################################
  #### PyMOL part of function ####


  pymol_lines <- c()

  #names(fasta_file)

  #put within
  #is.null(show_only_real_structures)
  #show_only_real_structures will be a list like names(fasta_file)
  #if the name is in the list && !grepl('___')
  #then do not load or include the

  #names(fasta_file)
  #names(fasta_file) %in% "SUZ12_Q15022test"




  #load the files
  pdbs_in_df <- unique(c(levels(xlink_mega_df$pdb1),levels(xlink_mega_df$pdb2)))
  #loop through the pdbs
  #pdb_name <- pdbs_in_df[1]
  for(pdb_name in pdbs_in_df){
    #load pdb files

    if(!is.null(show_only_real_structures) && (pdb_numbering == FALSE)){
      pdb_in_sors <- FALSE
        for(sors in show_only_real_structures){
          #noob solution to current problem
          if(grepl(sors,pdb_name)){
            pdb_in_sors <- TRUE
          }


        } #end first for(sors in show_only_real_structures)


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

      if(pdb_numbering == TRUE){

        #use the fetch command instead
        pdb_name <- strsplit(pdb_name,'_')[[1]][1]

        py_line <- paste0('fetch ',pdb_name,', async=0')


      } else { #end if(pdb_numbering == TRUE)

        py_line <- paste('load ',experiment_directory,'/',pdb_name,
                         sep='')

      } #end else to if(pdb_numbering == TRUE)



      if(!(py_line %in% pymol_lines)){
        pymol_lines <- c(pymol_lines,py_line)
      }




    }

    # py_line <- paste('load ',experiment_directory,'/',pdb_name,
    #                  sep='')
    # pymol_lines <- c(pymol_lines,py_line)


    #if this boolean is turned on,
    #color color_name, protein_name

    #should each of the PDBs be colored based on the protein they come from
    #can make this a user option

  }

  #can just have one option
  #if show_surface, color_gray == TRUE
  #py_line <- 'show surface'
  #pymol_lines <- c(pymol_lines,py_line)
  py_line <- 'color gray'
  pymol_lines <- c(pymol_lines,py_line)


  mega_distance_count <- 0
  for(row_num in 1:nrow(xlink_mega_df)){

    pdb1 <- strsplit(as.character(xlink_mega_df$pdb1[row_num]),'.pdb')[[1]]
    pdb2 <- strsplit(as.character(xlink_mega_df$pdb2[row_num]),'.pdb')[[1]]

    if(pdb_numbering == TRUE){
      chain1 <- strsplit(pdb1,'_')[[1]][2]
      chain2 <- strsplit(pdb2,'_')[[1]][2]
      pdb1 <- strsplit(pdb1,'_')[[1]][1]
      pdb2 <- strsplit(pdb2,'_')[[1]][1]


    }

    #if grepl '___' get last character for the chain
    #otherwise use 'A'

    pro_pos1 <- as.character(xlink_mega_df$pro_pos1[row_num])
    pro_pos2 <- as.character(xlink_mega_df$pro_pos2[row_num])

    #if statement to select the right chain


    #can do the check here in the grepl to do the check if it a single dot PDB structure, linear or curved

    #should make a list of pdbs so that this is not repeated twice?
    #can make output as a list that is then used for the next part of the function

    pdb_names <- c(pdb1,pdb2)
    pro_pos_list <- c(pro_pos1,pro_pos2)
    if(pdb_numbering == TRUE){
      chain_list <- c(chain1,chain2)
    }
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

      if(grepl('___',pdb_num) || pdb_numbering == TRUE){ #if 'real' PDB structure use last character since that is the chain

        if(pdb_numbering == TRUE){
          chain <- chain_list[pdb_num_count]
        } else {
          chain <- substrRight(pdb_num,1)
        }

        selection_name_list[[chain_name]] <- chain

      } else {

        #should it be in this part
        #should have an or statement or excessive?
        #do we need file_type --> unless the file is being created within the function
        #if(file_type_2d == 'single atom' && grepl('single_atom',pdb_num)){
        if(grepl('single_dot',pdb_num)){

          #transform the protein position for the PyMol script so it is the first in the list
          pro_pos <- list_of_start_and_end_pdbs[[paste(pdb_num,'.pdb',sep='')]][1]

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
          freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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

          freq_color <- as.character(xlink_mega_df$freq_color[row_num])
          py_line <- paste('color ',freq_color,', ',dist_name,sep='')
          pymol_lines <- c(pymol_lines,py_line)

        } #end if((grepl(sors,pdb_name) && grepl('___',pdb_name)) || !grepl(sors,pdb_name))

    } else {#end if(!is.null(show_only_real_structures))

      freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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

      freq_color <- as.character(xlink_mega_df$freq_color[row_num])
      py_line <- paste('color ',freq_color,', ',dist_name,sep='')
      pymol_lines <- c(pymol_lines,py_line)

    }
    #color the amino acids
    # freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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
    # freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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

  #write lines to .txt
  #separate by '\n'

  #pymol_file_list_file_name <- 'pymol_ts_output.pml'

  setwd(new_directory)
  write(paste(pymol_lines,collapse = '\n'),pymol_file_list_file_name)
  if(!is.null(xlink_view_file_name)){
    xlink_viewer_csv <- ppi.xinet(xlink_mega_df)
    write.csv(xlink_viewer_csv,xlink_view_file_name)
  }
  setwd(current_directory)

  #rename the freq_color column name if there's a category color instead?
  if(!is.null(category_color_input_file)){

    names(xlink_mega_df)[names(xlink_mega_df) == 'freq_color'] <- 'category_color'

  }


  return(xlink_mega_df)


} #end function ppi.analyze






#----Get new PWA Score-----

get_new_pwa_score <- function(){

  valid_pwa_number <- FALSE
  while(valid_pwa_number == FALSE){
    new_pwa_score_rl <- readline(prompt = "Enter a new cutoff score: ")
    if(is.na(as.numeric(new_pwa_score_rl))){
      cat('Please enter a real number.\n')
    } else {
      valid_pwa_number <- TRUE
    }
  } #end while valid_pwa_number == FALSE

  #turn off warnings for this function?

  return(new_pwa_score_rl)
}

#----Match Sequence to PDB and chain----

match_sequence_to_pdb_and_chain <- function(protein_name,
                                            fasta_seq,
                                            pdb_read,
                                            pdb_rl,
                                            pwa_score_threshold = 500,
                                            peptide_sequence = NULL){

  valid_chain_and_pdb_selected <- FALSE
  while(valid_chain_and_pdb_selected == FALSE){

    seqres_matching_list <- do_sr_chain_loop(protein_name = protein_name,
                                             fasta_seq = fasta_seq,
                                             pdb_read = pdb_read,
                                             pdb_rl = pdb_rl,
                                             pwa_score_threshold = pwa_score_threshold,
                                             peptide_sequence = peptide_sequence)



    #if the list == 0, do the advanced options to get a match
    if(length(seqres_matching_list$chains) == 0){

      #include the score
      no_chains_matched_title <- paste('\nNo chains matched ',protein_name,' against ',pdb_rl,' using your pairiwse alignment score of ',
                                       as.character(pwa_score_threshold),
                                       '\n\nChoose a new option from the list below.',sep='')

      #new list of options that will be displayed to the user
      no_chains_matched_options <- c('Change the pairwise alignment score cutoff',
                                     'Search using different PDB ID',
                                     'BLAST this sequence')


      ncm_menu_selection <- menu(no_chains_matched_options,title = no_chains_matched_title)

      if(ncm_menu_selection == 1){ #change the pwa score

        #have the user enter a new pwa score, go back to the top of the while loop to rerun the function

        pwa_score_threshold <- get_new_pwa_score()

        #should the pwa score be stored in the same variable and then sent back to the top of the while loop


        #enter the new score into the sr_chain_loop?

      } else if(ncm_menu_selection == 2){ #search using different PDB ID


        #should this also be turned into a function for reusability?
        valid_pdb_selected <- FALSE
        while(valid_pdb_selected == FALSE){

          pdb_rl <- readline(prompt = "Enter the PDB ID: ")

          ff_pdb_read <- check_download_read_pdb(pdb_rl)
          if(is.null(ff_pdb_read)){
            cat('Please enter a valid PDB ID\n')
          } else {
            valid_pdb_selected <- TRUE
          }

        } #end while(valid_pdb_selected == FALSE)


        #have the user enter a new pdb id, go back to the top of the while loop

      } else if(ncm_menu_selection == 3){ #blast this sequence

        #just go through menu loops?

        #copy the code for doing the BLAST sequence?
        #should turn into a function for better reproducibility?

        pdb_info <- go_through_menu_loops(fasta_sequence_vector = fasta_seq)

        ff_pdb_read <- check_download_read_pdb(pdb_info$pdb_id)

        resno_and_resid <- quick_resno_and_resid(ff_pdb_read,pdb_info$chain)
        pdb_info <- c(pdb_info,resno_and_resid)

        valid_chain_and_pdb_selected <- TRUE
        #return(pdb_info)
        #boolean to escape from while loop?

      }

    }  else { #end if(length(seqres_matching_list$chains) == 0)


      #if the list actually contains options
      #the user will pick the options here from a list
      pick_chain_title <- paste("\nPick a chain to use for ",protein_name," for PDB ID: ",pdb_rl,sep='')

      seqres_matching_list$recommend[[protein_name]]

      #can just alter the sequence

      pick_chain_menu_selection <- menu(seqres_matching_list$chains[[protein_name]],title=pick_chain_title)

      chain_picked <- seqres_matching_list$chains[[protein_name]][pick_chain_menu_selection]



      pdb_info <- seqres_matching_list$pdb_info[[paste(pdb_rl,'_',chain_picked,'_',protein_name,sep='')]]

      resno_and_resid <- quick_resno_and_resid(pdb_read,chain_picked)
      pdb_info <- c(pdb_info,resno_and_resid)

      #either return statement or boolean to break out of the while loop
      valid_chain_and_pdb_selected <- TRUE
      #return(pdb_info)



    } #end else to if(length(seqres_matching_list$chains) == 0)

  } #end while(valid_chain_and_pdb_selected == FALSE)

  return(pdb_info)

} #end function match_sequence_to_pdb_and_chain






#----Preferred PDB Menu (still under construction)----

#in order to have it so that the PDB structure will be used multiple times,
#the protein name may have to be within the loop
#what then will be the output of the loop?
#will need multiple options for what will happen with the PDB info
#need to create a scenario that will accept pdb_info for plink analysis

#'Display Preferred PDB Structure Menu
#'
#'This function creates a menu interface for preferred PDB structure
#'
#'@param protein_name Protein name that corresponds to name in fasta_file to extract sequence
#'@param fasta_file FASTA file as read by seqinr::read.fasta()
#'@param pwa_score_threshold Minimum score to determine a match using Biostrings::pairwiseAlignment()
#'@param peptide_sequence An optional peptide sequence to use for recommendations if trying to match a peptide to a structure. Defaults to NULL
#'@export

display_preferred_pdb_structure_menu <- function(protein_name, fasta_file, pwa_score_threshold = 500,
                                                 peptide_sequence = NULL){

  fasta_seq <- fasta_file[[protein_name]]
  #can ask user whether or not they want to go through the loops
  #keep in a while loop pdb_selected == FALSE
  run_through_loops_title <- paste('Do you have a preferred PDB structure you would like to match to your sequence ',
                                   protein_name,'?',sep='')
  run_through_loops_options <- c('Yes','No, run BLAST on this sequence','Advanced options')
  rtl_menu_selection <- menu(run_through_loops_options, title=run_through_loops_title)
  if(rtl_menu_selection == 1){ #yes
    #error check: make sure that the PDB ID is valid
    #may require a tryCatch()

    valid_pdb_selected <- FALSE
    while(valid_pdb_selected == FALSE){

      pdb_rl <- readline(prompt = "Enter the PDB ID: ")

      ff_pdb_read <- check_download_read_pdb(pdb_rl)
      if(is.null(ff_pdb_read)){
        cat('Please enter a valid PDB ID\n')
      } else {
        valid_pdb_selected <- TRUE
      }

    } #end while(valid_pdb_selected == FALSE)

    pdb_info <- match_sequence_to_pdb_and_chain(protein_name = protein_name,
                                                fasta_seq = fasta_seq,
                                                pdb_read = ff_pdb_read,
                                                pdb_rl = pdb_rl,
                                                pwa_score_threshold = pwa_score_threshold)

    return(pdb_info)

  } else if(rtl_menu_selection == 2){ #no

    pdb_info <- go_through_menu_loops(fasta_sequence_vector = fasta_seq)

    pdb_read <- check_download_read_pdb(pdb_info$pdb_id)

    resno_and_resid <- quick_resno_and_resid(pdb_read,pdb_info$chain)
    pdb_info <- c(pdb_info,resno_and_resid)
    return(pdb_info)


  } else if(rtl_menu_selection == 3){ #advanced options

    #advanced options menu
    #as part of a while loop?
    #need to revisit the original menu cycle to make sure all of the while loops
    #are appropriate

    advanced_options <- c('Enter a PDB ID to use for all sequences in this fasta file',
                          'Run BLAST for all sequences in this fasta file',
                          'Use previously found PDB file for this sequence',
                          'Go back')
    advanced_options_title <- 'Advanced Options'
    advanced_options_menu_selection <- menu(advanced_options,title = advanced_options_title)

    if(advanced_options_menu_selection == 1){ #enter a PDB ID for all sequences

      #turn on some kind of boolean
      #store the pdb in a variable

    } else if(advanced_options_menu_selection == 2){ #run BLAST for all sequences

      #turn on a different boolean

    } else if(advanced_options_menu_selection == 3){ # use previously found PDB file

      #when a PDB has been searched for, add to a list that will be displayed here

    } else if(advanced_options_menu_selection == 4){#go back


    }

  } #end else if rtl_menu_selection == 3

  return(pdb_info)

} #end function display_preferred_pdb_structure_menu


#----Generate PDB Lists from CSV with PDB Info----

#'Generate PDB Lists from PDB CSV
#'
#'@param csv_file_name String containing path to csv file
#'@export

generate_pdb_lists_from_pdb_csv <- function(csv_file_name){

  master_list_of_pdbs_df <- read.csv(csv_file_name)
  all_pdb_files <- list()
  start_end_chain <- list()

  #separate the rows by the seq_name and then make everything left as a list
  for(seq_name in unique(as.vector(master_list_of_pdbs_df$seq_name))){


    mlop_sub_df <- master_list_of_pdbs_df[master_list_of_pdbs_df$seq_name == seq_name,]
    pdb_names_for_seq <- as.vector(mlop_sub_df$pdb_name)
    all_pdb_files[[seq_name]] <- pdb_names_for_seq
    for(pdb_name in pdb_names_for_seq){

      mlop_sub_df2 <- mlop_sub_df[mlop_sub_df$pdb_name == pdb_name,]

      start_end_chain[[pdb_name]] <- c(as.character(mlop_sub_df2$start_pos),
                                       as.character(mlop_sub_df2$end_pos),
                                       as.character(mlop_sub_df2$chain))

    }

  } #

  return(list(all_pdb_files=all_pdb_files,start_end_chain=start_end_chain))

} #end function generate_pdb_lists_from_pdb_csv


#----Get Category Colors from Input File for Diff. Analysis----

get_category_colors_from_input_file <- function(diff_analysis_df,file_input,multi_cat_color='grey'){

  #can have user change this as needed, but grey will be the defeault color
  #can check this only once at the beginning after file_input has been entered
  if('color' %in% names(file_input)){
    color_spelling <- 'color'
  } else if('colour' %in% names(file_input)){
    color_spelling <- 'colour'
  } else {
    cat('No color specified in the file')
    return(NULL)
  }

  category_colors_list <- c()

  for(files_found in as.vector(diff_analysis_df$files)){
    #use boolean to keep track of if files fall into one or more categories
    #or have list and then see how many unique values there are in the list
    #separate by + in list?
    #can also check if + in list before separation?

    if(grepl('\\+',files_found)){
      files_found <- strsplit(files_found,'\\+')[[1]]
    }


    file_colors <- c()
    for(file_found in files_found){

      file_color <- as.character(file_input[file_input$file_name == file_found,][[color_spelling]])
      file_colors <- c(file_colors,file_color)


    } #end for(file_found in files_found)

    #check file_colors list
    unique_file_colors <- unique(file_colors)
    if(length(unique_file_colors) > 1){
      #not only one color --> assign grey as the color
      category_color <- multi_cat_color
    } else {
      #only one color in the list
      #assign the color in that unique
      category_color <- unique_file_colors
    }

    #add category color to the list of category colors
    category_colors_list <- c(category_colors_list,category_color)

  } #end for(files_found in as.vector(bs3_diff_analysis_with_csv$files))

  return(category_colors_list)

} #end function get_category_colors_from_input_file


#----CSV export from PDB master list-----

make_start_end_pdb_df_output <- function(master_list_of_pdbs,
                                         output_file_name='start_and_end'){

  seq_name_list <- c()
  pdb_name_list <- c()
  start_pos_list <- c()
  end_pos_list <- c()
  chain_list <- c()

  for(fasta_name in names(master_list_of_pdbs)){

    #make the list that will be turned into a data.frame()

    #seq_name_list <- c(seq_name_list,fasta_name)

    sub_mlop <- master_list_of_pdbs[[fasta_name]]

    #go through each of the
    for(pdb_name in names(sub_mlop)){

      seq_name_list <- c(seq_name_list,fasta_name)
      pdb_name_list <- c(pdb_name_list,pdb_name)
      start_pos_list <- c(start_pos_list,min(sub_mlop[[pdb_name]]))
      end_pos_list <- c(end_pos_list,max(sub_mlop[[pdb_name]]))

      if(grepl('___',pdb_name)){
        pdb_name1 <- strsplit(pdb_name,'.pdb')[[1]]
        chain <- substrRight(pdb_name1,1)
      } else {
        chain <- 'A'
      }

      chain_list <- c(chain_list,chain)
    }

  }

  master_list_of_pdbs_df <- list(seq_name=seq_name_list,
                                 pdb_name=pdb_name_list,
                                 start_pos=start_pos_list,
                                 end_pos=end_pos_list,
                                 chain=chain_list)

  master_list_of_pdbs_df <- data.frame(master_list_of_pdbs_df)

  write.csv(master_list_of_pdbs_df,paste(output_file_name,'.csv',sep=''))

  return(master_list_of_pdbs_df)

}

#----Make Start/End Output 2----

#'Make start and end csv file for make_diff_analysis function
#'
#'This function creates the csv file within the make_diff_analysis function for BS3 XL-MS function
#'that can be used as input to direct function to correct PDB files
#'
#'@param generated_pdb_lists list() generated from generate_pdb_lists_and_files_from_fasta() or generate_pdb_lists_and_files_pdb_csv()
#'@param output_file_name String containing name of output file. Defaults to 'start_and_end'
#'@export

make_start_end_pdb_df_output2 <- function(generated_pdb_lists,
                                          output_file_name='start_and_end'){

  seq_name_list <- c()
  pdb_name_list <- c()
  start_pos_list <- c()
  end_pos_list <- c()
  chain_list <- c()



  for(seq_name in names(generated_pdb_lists$all_pdb_files)){

    #make the list that will be turned into a data.frame()

    #seq_name_list <- c(seq_name_list,fasta_name)

    sub_mlop <- generated_pdb_lists$all_pdb_files[[seq_name]]

    #go through each of the
    for(pdb_name in sub_mlop){

      seq_name_list <- c(seq_name_list,seq_name)
      pdb_name_list <- c(pdb_name_list,pdb_name)
      start_pos_list <- c(start_pos_list,generated_pdb_lists$start_end_chain[[pdb_name]][1])
      end_pos_list <- c(end_pos_list,generated_pdb_lists$start_end_chain[[pdb_name]][2])
      chain_list <- c(chain_list,generated_pdb_lists$start_end_chain[[pdb_name]][3])
    }

  }

  master_list_of_pdbs_df <- list(seq_name=seq_name_list,
                                 pdb_name=pdb_name_list,
                                 start_pos=start_pos_list,
                                 end_pos=end_pos_list,
                                 chain=chain_list)

  master_list_of_pdbs_df <- data.frame(master_list_of_pdbs_df)

  write.csv(master_list_of_pdbs_df,paste(output_file_name,'.csv',sep=''))

  return(master_list_of_pdbs_df)

}



#----Make Xlinker Viewer CSV----

#'Makes file for xiNET input
#'
#'This function writes a file suitable for input into the xiNET xlink viewer website
#'
#'@param xlink_df Dataframe created by make_diff_analysis function
#'@param write_file Boolean, TRUE writes file as csv
#'@param xlink_viewer_csv_file_name String indicating file name of csv output if write_file = TRUE
#'@export

ppi.xinet <- function(xlink_df, write_file = TRUE, xlink_viewer_csv_file_name = 'xlink_viewer.csv', add_color = FALSE){
  #should have a conversion tool here that would convert between the

  xlink_viewer_csv <- list()

  #xlink_pep_pos1 <- as.numeric(as.character(xlink_df$pro_pos1)) - as.numeric(as.character(xlink_df$pep_pos1))
  #xlink_pep_pos2 <- as.numeric(as.character(xlink_df$pro_pos2)) - as.numeric(as.character(xlink_df$pep_pos2))

  #xlink.df$pep_start1 <- (xlink.df$pro_pos1 - xlink.df$pep_pos1)

  #xlink_df <- prc2_rna_diff_analysis

  if('pro_pos1' %in% colnames(xlink_df)){
    xlink_pep_pos1 <- as.numeric(as.character(xlink_df$pro_pos1)) - as.numeric(as.character(xlink_df$pep_pos1)) +1
    xlink_pep_pos2 <- as.numeric(as.character(xlink_df$pro_pos2)) - as.numeric(as.character(xlink_df$pep_pos2)) +1
    #xlink_pep_pos1 <- as.numeric(as.character(xlink_df$pro_pos1))
    #xlink_pep_pos2 <- as.numeric(as.character(xlink_df$pro_pos2))


    # xlink_viewer_csv$PepPos1 <- xlink_df$pro_pos1
    # xlink_viewer_csv$PepPos2 <- xlink_df$pro_pos2
    xlink_viewer_csv$PepPos1 <- xlink_pep_pos1
    xlink_viewer_csv$PepPos2 <- xlink_pep_pos2

    # xlink_viewer_csv$Protein1 <- xlink_df$pro_name1
    # xlink_viewer_csv$Protein2 <- xlink_df$pro_name2
    xlink_viewer_csv$PepSeq1 <- xlink_df$pep_seq1
    xlink_viewer_csv$PepSeq2 <- xlink_df$pep_seq2
    xlink_viewer_csv$LinkPos1 <- xlink_df$pep_pos1
    xlink_viewer_csv$LinkPos2 <- xlink_df$pep_pos2
    # xlink_viewer_csv$Score <- xlink_df$score
    xlink_viewer_csv$Protein1 <- xlink_df$pro_name1
    xlink_viewer_csv$Protein2 <- xlink_df$pro_name2
    #xlink_viewer_csv$LinkPos1 <- xlink_df$pro_pos1
    #xlink_viewer_csv$LinkPos2 <- xlink_df$pro_pos2
    xlink_viewer_csv$Score <- xlink_df$score

    if(add_color == TRUE){
      if('category_color' %in% colnames(xlink_df)){
        xlink_viewer_csv$CategoryColor <- xlink_df$category_color
      } else if('freq_color' %in% colnames(xlink_df)){
        xlink_viewer_csv$FrequencyColor <- xlink_df$freq_color
      }

    } #end if(add_color == TRUE)


  } else {
    xlink_pep_pos1 <- as.numeric(as.character(xlink_df[['Protein.Position.1']])) - as.numeric(as.character(xlink_df[['Peptide.Position.1']])) +1
    xlink_pep_pos2 <- as.numeric(as.character(xlink_df[['Protein.Position.2']])) - as.numeric(as.character(xlink_df[['Peptide.Position.2']])) +1
    #xlink_pep_pos1 <- as.numeric(as.character(xlink_df[['Protein.Position.1']]))
    #xlink_pep_pos2 <- as.numeric(as.character(xlink_df[['Protein.Position.2']]))


    # xlink_viewer_csv$PepPos1 <- xlink_df[['Protein.Position.1']]
    # xlink_viewer_csv$PepPos2 <- xlink_df[['Protein.Position.2']]
    xlink_viewer_csv$Protein1 <- xlink_df[['Protein.Name.1']]
    xlink_viewer_csv$Protein2 <- xlink_df[['Protein.Name.2']]
    xlink_viewer_csv$PepPos1 <- xlink_pep_pos1
    xlink_viewer_csv$PepPos2 <- xlink_pep_pos2
    xlink_viewer_csv$PepSeq1 <- xlink_df[['Peptide.Sequence.1']]
    xlink_viewer_csv$PepSeq2 <- xlink_df[['Peptide.Sequence.2']]
    xlink_viewer_csv$LinkPos1 <- xlink_df[['Peptide.Position.1']]
    xlink_viewer_csv$LinkPos2 <- xlink_df[['Peptide.Position.2']]
    #xlink_viewer_csv$LinkPos1 <- xlink_df[['Protein.Position.1']]
    #xlink_viewer_csv$LinkPos2 <- xlink_df[['Protein.Position.2']]
    xlink_viewer_csv$Score <- xlink_df[['Score']]

    if(add_color == TRUE){
      if('Category.Color' %in% colnames(xlink_df)){
        xlink_viewer_csv$CategoryColor <- xlink_df[['Category.Color']]
      } else if('Frequency.Color' %in% colnames(xlink_df)){
        xlink_viewer_csv$FrequencyColor <- xlink_df[['Frequency.Color']]
      }

    } #end if(add_color == TRUE)



  }



  xlink_viewer_csv <- data.frame(xlink_viewer_csv)

  if(write_file == TRUE){

    if(!endsWith(xlink_viewer_csv_file_name,'csv')){
      xlink_viewer_csv_file_name <- paste0(xlink_viewer_csv_file_name,'.csv')
    }
    write.csv(xlink_viewer_csv,xlink_viewer_csv_file_name)
  }

  return(xlink_viewer_csv)
}

#----Make Plink2 Master List-----

#'Make pLink2 master list
#'
#'This function creates the list of
#'@export
make_plink2_master_list <- function(plink2_peptides, list_of_protein_names,
                                    plink2_type_of_file=c('cross-linked','loop-linked'),
                                    master_list_index = 0){

  plink2_master_list <- list()

  first_iteration <- TRUE

  #need to intialize layer1 and layer2 column names
  layer1_column_names <- strsplit(plink2_peptides[1],',')[[1]] #layer 1
  layer2_column_names <- strsplit(plink2_peptides[2],',')[[1]]

  for(index_num in 3:length(plink2_peptides)){

    #index_num <- 56

    #starting with the values make sure the number of columns is the same
    plink2_row <- strsplit(plink2_peptides[index_num],',')[[1]]
    #plink2_row <- strsplit(looplinked_peptides[92],',')[[1]]

    if(length(plink2_row) == length(layer1_column_names)){
      #it is in layer 1

      if(first_iteration == FALSE){

        plink2_master_list[[master_list_index]] <- row_list

      } else {

        first_iteration <- FALSE

      }
      #add previous row_list to the master list before intializing a new list


      master_list_index <- master_list_index + 1 #plink2_row[1]
      #will need to revisit this as if multiple files are being combined into one, might
      #be better to do a count

      row_list <- list()

      #add first column/name of the first column
      order_number <- plink2_row[1] #NEED TO FIX
      peptide <- plink2_row[2]

      if(plink2_type_of_file == 'cross-linked'){

        split_peptide <- split_sequences(peptide)

        proteins <- plink2_row[5] #will need quality checks if doing it this way
        split_proteins <- strsplit(proteins,'/')[[1]]

        row_list$proteins_xlink <- split_proteins

        row_list$sequence_xlink <- peptide

        row_list$modification_xlink <-  plink2_row[layer1_column_names == 'Modifications']

      } else if(plink2_type_of_file == 'loop-linked'){

        #will just need to change the split peptide portion
        #plink2_row2 <- strsplit(looplinked_peptides[3],',')[[1]]

        #need to remove the '2' for the final code

        plink2_row2 <- plink2_row

        proteins <- plink2_row2[layer1_column_names == 'Proteins']
        split_proteins <- strsplit(proteins,'/')[[1]]

        #row_list$proteins_xlink <- split_proteins

        row_list$modification_xlink <-  plink2_row2[layer1_column_names == 'Modifications']

        peptide2 <- plink2_row2[2]



        split_peptide2 <- strsplit(peptide2,'\\)')[[1]]
        split_peptide2 <- strsplit(split_peptide2,'\\(')
        split_peptide2 <- c(split_peptide2[[1]],split_peptide2[[2]])

        split_peptide2 <- split_peptide2[!('' == split_peptide2)]
        pep_sequence2 <- split_peptide2[1]
        pep_pos1_2 <- split_peptide2[2]
        pep_pos2_2 <- split_peptide2[3]

        new_peptide <- paste(pep_sequence2,'(',pep_pos1_2,')-',pep_sequence2,'(',pep_pos2_2,')',sep='')
        row_list$sequence_xlink <- new_peptide

        #need to duplicate the first sequence for the second peptide

        proteins2 <- plink2_row2[layer1_column_names == 'Proteins']
        #proteins2 <- strsplit(looplinked_peptides[118],',')[[1]][layer1_column_names == 'Proteins']


        split_proteins2 <- strsplit(proteins2,'/')[[1]]

        #row_list$sequence_xlink <- split_proteins2

        new_peptide2_list <- c()

        #NEEDS TO BE MODIFIED SINCE EZH2 errors are occuring here


        #proteins2 <- strsplit(looplinked_peptides[118],',')[[1]][layer1_column_names == 'Proteins']

        #split_protein2 <- split_proteins2[1]

        #this should be within a for loop for each of the different proteins
        for(split_protein2 in split_proteins2){

          #can do the same kind of loop here that will remove the actual

          #list_of_protein_names
          #have separate function to check on the protein name that can be included
          #in both other functions

          #only need to check once since there is only 1 protein in the loop-linked file

          loop_pep_numbers <- ''
          loop_protein_name <- ''

          for(protein_name in list_of_protein_names){

            enclosed_protein_name <- gsub('\\)','\\\\)',protein_name)
            enclosed_protein_name <- gsub('\\(','\\\\(',enclosed_protein_name)
            enclosed_protein_name <- gsub('\\|','\\\\|',enclosed_protein_name)

            #see if the protein name is in split_protein2

            if(grepl(enclosed_protein_name,split_protein2)){

              #make a string containing just the peptide numbers?

              loop_pep_numbers <- gsub(enclosed_protein_name,'',split_protein2)
              loop_protein_name <- protein_name
            }

            #if yes, replace the protein name with an empty string '' in the string

            #can break and move on to the next step
            if(!(loop_pep_numbers == '')){
              break
            }

          }

          if(!(loop_pep_numbers == '')){
            #if there is actually something and peptide numbers were extracted


            loop_pep_split <- strsplit(loop_pep_numbers,'\\(')[[1]]
            loop_pep1 <- strsplit(loop_pep_split,'\\)')[[2]]
            loop_pep2 <- strsplit(loop_pep_split,'\\)')[[3]]


            new_peptide2 <- paste(loop_protein_name,'(',loop_pep1,')-',loop_protein_name,'(',loop_pep2,')',sep='')


            loop_qc_check <- paste(loop_protein_name,'(',loop_pep1,')(',loop_pep2,')',sep='')

            #do paste() quality control here to make sure that the previous
            #peptide matches the

            if(loop_qc_check == split_protein2){

              new_peptide2_list <- c(new_peptide2_list,new_peptide2)

            } else {
              loop_qc_check <- paste(loop_protein_name,' (',loop_pep1,')(',loop_pep2,')',sep='')
              if(loop_qc_check == split_protein2){
                new_peptide2_list <- c(new_peptide2_list,new_peptide2)
              } else {
                #cat(paste(loop_qc_check,split_protein2,'\n'))
                new_peptide2_list <- c(new_peptide2_list,paste('failed_quality_check:',split_proteins2))
              }
            }

            #new_peptide2_list <- c(new_peptide2_list,new_peptide2)

          } else {
            #nothing extracted and likely a contaminant
            #do what the original code did and hope for the best
            #(likely won't matter because it's a contaminant and is not critical)

            split_peptide2 <- strsplit(split_protein2,'\\)')[[1]]
            split_peptide2 <- strsplit(split_peptide2,'\\(')
            split_peptide2 <- c(split_peptide2[[1]],split_peptide2[[2]])

            split_peptide2 <- split_peptide2[!('' == split_peptide2)]
            pep_sequence2 <- split_peptide2[1]
            pep_pos1_2 <- split_peptide2[2]
            pep_pos2_2 <- split_peptide2[3]

            new_peptide2 <- paste(pep_sequence2,'(',pep_pos1_2,')-',pep_sequence2,'(',pep_pos2_2,')',sep='')
            #will have to change this to more of a list structure to be able to account for multiples
            #will this need to be minimized to just 1? how to account for this when
            #there are different lengths for peptide or does it not matter in this case?

            new_peptide2_list <- c(new_peptide2_list,new_peptide2)


          } #end else for if(!(loop_pep_numbers == ''))

          #the next step:
          #should account for the fact that the protein may not have showed up in the
          #list of proteins (like in the other set of codes)
          #if it does not show up --> what to do?
          #split it like it was done originally and hope for the best?




          #row_list$sequence_xlink <- new_peptide

        } #end for(split_protein2 in split_proteins2)

        row_list$proteins_xlink <- new_peptide2_list
        #for loop?
        #will then do the same thing as the first one

        #will need to adjust both proteins and sequence to match


      } #end else if(plink2_type_of_file == 'loop-linked')

    } else if((length(plink2_row) == length(layer2_column_names)) || ((length(plink2_row)-1) == length(layer2_column_names))){
      #it is in layer 2
      #add all of this information to the index as defined above
      #plink2_row

      if(((length(plink2_row)-1) == length(layer2_column_names))){
        if(is.na(as.numeric(plink2_row[layer2_column_names == 'Charge']))){
          #if this is true, the Charge column is not actually the charge column
          #need to merge the 2 columns
          plink2_row <- c(plink2_row[1:2],paste0(plink2_row[3:4],collapse=''),plink2_row[5:length(plink2_row)])

        } else { #end if(is.na(as.numeric(plink2_row[layer2_column_names == 'Charge']))){
          cat('Unknown layer detected\n')
          next
          }#end else to if(is.na(as.numeric(plink2_row[layer2_column_names == 'Charge']))){
      } #end if(((length(plink2_row)-1) == length(layer2_column_names))){


      row_list$score_xlink <- c(row_list$score_xlink,plink2_row[layer2_column_names == 'Score'])
      row_list$delta_xlink <- c(row_list$delta_xlink,plink2_row[layer2_column_names == 'Precursor_Mass_Error(Da)'])
      row_list$score_xlink <- c(row_list$score_xlink,plink2_row[layer2_column_names == 'Score'])

      #maybe should have a list of the terms used in this type of file to make sure that
      #there is consistency in the output even if the terms are not in the same order
      #should first make sure that

      #lists for this section:
      #spectrum_xlink (labeled in the second row as "Title")
      #evalue_xlink (not in previous iteration but labled as "Evalue")
      #modification_xlink ("Modification")
      #delta_xlink ('Precursor_Mass_Error(Da)' --> might be slightly different when loaded)
      #score_xlink ('Score')



    } else {
      #it is in an unknown layer
      cat('Unknown layer detected\n')

    } #end else #if(length(plink2_row) == length(layer1_column_names))




  } #end for(index_num in 3:length(crosslinked_peptides))

  plink2_master_list[[master_list_index]] <- row_list

  return(plink2_master_list)

} #end function make_plink2_master_list

#----strcount-----

#from https://aurelienmadouasse.wordpress.com/2012/05/24/r-code-how-the-to-count-the-number-of-occurrences-of-a-substring-within-a-string/

strcount <- function(x, pattern, split){

  unlist(lapply(
    strsplit(x, split),
    function(z) na.omit(length(grep(pattern, z)))
  ))

}

#----Extract proteins and peptide position from string-----

#split_proteins <- 'PHF19_Q5T6S3(25)-PHF19_Q5T6S3(32)'

extract_protein_name_and_peptide_num <- function(split_proteins,list_of_protein_names){

  pro_pos_list <- c()
  matched_pro_list <- c()
  pep_num_to_extract <- split_proteins

  #make sure the list of proteins is in reverse order in case not done before
  list_of_protein_names <- list_of_protein_names[order(nchar(list_of_protein_names),
                                                       list_of_protein_names,decreasing=TRUE)]



  for(protein_name in list_of_protein_names){

    #check if the protein name is within the protein string

    #grepl(protein_name,split_proteins)
    #enclosed doesn't work, will say TRUE for anything
    #grepl("[HMBP_PrS_EZH2(iso2-BsaI)_CD136]",split_proteins)

    #check first if there are parantheses in the string or will it not do anything if there
    #are none

    #can just keep this for all of them since gsub will not throw an error if
    #the substring you want to substitute in not in the string

    #protein_name <- 'tr|Q0MUU8|Q0MUU8_TRINI'
    #protein_name <- 'EED_O75530'

    #protein_name <- 'PHF19_Q5T6S3'

    enclosed_protein_name <- gsub('\\)','\\\\)',protein_name)
    enclosed_protein_name <- gsub('\\(','\\\\(',enclosed_protein_name)

    enclosed_protein_name <- gsub('\\|','\\\\|',enclosed_protein_name)



    #grepl(enclosed_p,'PHF19(333)-tr|Q0MUU8|Q0MUU8_TRINI(190)')

    #gsub('\\(','\\\\','test')

    #enclosed_protein_name

    #grepl(enclosed_protein_name,split_proteins)
    #if true, split the
    #strsplit(split_proteins,enclosed_protein_name)[[1]]

    #split_proteins <- 'PHF19(333)-tr|Q0MUU8|Q0MUU8_TRINI(190)'

    #split_proteins <- 'PHF19(31)-EED_O75530(339)'

    enclosed_split_proteins <- gsub('\\)','\\\\)',split_proteins)
    enclosed_split_proteins <- gsub('\\(','\\\\(',enclosed_split_proteins)

    enclosed_split_proteins <- gsub('\\|','\\\\|',enclosed_split_proteins)

    protein_count <- str_count(split_proteins,enclosed_protein_name)

    #str_count('PHF19(333)-tr\\|Q0MUU8\\|Q0MUU8_TRINI\\(190\\)',enclosed_split_proteins)

    if(protein_count == 0){

      next
      #protein has not showed up at all in this string, move on to the next

    } else if(protein_count == 1){

      #need to enclose the split proteins name as well?

      #startsWith(enclosed_split_proteins,enclosed_protein_name)
      if(startsWith(split_proteins,protein_name)){
        #if it does start with the protein
        #then assign this position 1

        #need to make sure that pro_pos_list does not already have 1 in it

        if(!(1 %in% pro_pos_list)){
          pro_pos_list <- c(pro_pos_list,1)
          do_not_add <- FALSE
        } else {
          do_not_add <- TRUE
        }

        #pro_pos_list <- c(pro_pos_list,1)

      } else {
        #does not start with this protein
        #assign this position 2


        if(!(2 %in% pro_pos_list)){
          pro_pos_list <- c(pro_pos_list,2)
          do_not_add <- FALSE
        } else {
          do_not_add <- TRUE
        }

      }

      if(do_not_add == FALSE){

        matched_pro_list <- c(matched_pro_list,protein_name)
        pep_num_to_extract <- gsub(enclosed_protein_name,'',pep_num_to_extract)

      }


      #use startsWith to check if protein is at beginning of the string
      #then order the proteins based on that
      #question would be when to order it
      #if it's a dataframe then would have to order them after the loop is finished

    } else if(protein_count == 2){

      #add the protein twice to the list and/or dataframe
      #may want to make it a dataframe as well for consistency?
      #then can be ordered the same way without error

      pep_num_to_extract <- gsub(enclosed_protein_name,'',pep_num_to_extract)

      matched_pro_list <- c(protein_name,protein_name)

      pro_pos_list <- c(1,2)


      #break
    }

    if(length(pro_pos_list) == 2){
      break
    }



    #quality control check --> paste everything together and see if it
    #matches the original split_peptides to make sure it's the same and no mistakes were made
    #can cat or do a warning --> 'Possible naming error in proteins' or something to that effect




    #can do a quality control check to make sure the replaced protein names are in fact underscores
    #may be able to check the order by making the order '_1_' and '_2_' to double check that they
    #are in the right place


    #enclosed_protein_name <- (paste('[',protein_name,']',sep=''))

    #replace the parantheses with underscores in the protein name
    #add this protein name to a larger list of protein names
    #or will this always have to be done more than once?

    #will likely have to cycle through the list twice(?)
    #what to do when the protein shows up twice?
    #just remove the protein name and then check it again?
    #
    # gsub('\\(','_',protein_name)
    # gsub('\\)','_',protein_name)

    #something to do with the parantheses is causing it not to recognize the string


    #then use gsub to replace the string
    #make the protein name the official protein name

  } #end for(protein_name in list_of_protein_names)

  #once all of the proteins have been run:

  #make dataframe from the pro_pos_list and the matched_pro_list
  #order the dataframe based on the pro_pos_list

  #pro_sort_df <- data.frame(list(protein_names=c('beta','alpha'),protein_order=c(2,1)))

  #pro_sort_df[order(pro_sort_df$protein_order),]

  pro_sort_df <- data.frame(list(protein_names=matched_pro_list,protein_order=pro_pos_list))


  if(nrow(pro_sort_df) == 0){ #nothing in the dataframe

    cat('Neither sequence in fasta file or dictionary - sending to unused file\n')
    return(NULL)
  } else if (nrow(pro_sort_df) == 1){
    cat('Only 1 protein found in fasta file and/or dictionary - sending to unused file\n')
    return(NULL)
  }

  pro_sort_df <- pro_sort_df[order(pro_sort_df$protein_order),]

  #should also check if only 1 of the proteins in the list

  #pep_num_to_extract
  pep_num_to_extract <- strsplit(pep_num_to_extract,'-')[[1]]
  pep_num_to_extract <- strsplit(pep_num_to_extract,'\\(')
  pep_num_to_extract2 <- c()
  for(pep_num in 1:length(pep_num_to_extract)){
    pnte <- strsplit(pep_num_to_extract[[pep_num]][2],'\\)')[[1]]
    pep_num_to_extract2 <- c(pep_num_to_extract2,pnte)
  }

  #as.character(pro_sort_df$protein_names)
  #pep_num_to_extract2
  #combine data together --> for loop? (so it is split like how it should be)

  #initialize list here?
  pro_output <- list()
  #qc_pro_check <- ''

  for(ind_num in 1:length(pep_num_to_extract2)){

    pro_output[[ind_num]] <- c(as.character(pro_sort_df$protein_names)[ind_num],pep_num_to_extract2[ind_num])

  }

  #QC check

  qc_pro_check <- paste(pro_output[[1]][1],'(',pro_output[[1]][2],')-',pro_output[[2]][1],'(',pro_output[[2]][2],')',sep='')

  #return(c(qc_pro_check,split_proteins))

  if(!(qc_pro_check == split_proteins)){

    qc_pro_check2  <- paste(pro_output[[1]][1],' (',pro_output[[1]][2],')-',pro_output[[2]][1],' (',pro_output[[2]][2],')',sep='')
    if(!(qc_pro_check2 == split_proteins)){
      cat('Protein naming error may have occurred')
      cat(qc_pro_check)
      return(NULL)
    }
  }

  return(pro_output)


  #return protein output variable

} #end function extract_protein_name_and_peptide_num


#----Check protein name-----


check_and_change_protein_name <- function(protein_name,protein_alternative_names_dict,
                                          list_of_all_pdb_files,fasta_file){

  if(protein_name %in% protein_alternative_names_dict$original_name){

    return(protein_name)

  }

  if(!(protein_name %in% names(list_of_all_pdb_files))){


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

            cat(paste(protein_name,'changed to',as.character(protein_alternative_names_dict[row_num,'original_name']),
                      'in output\n'))

            protein_name <- as.character(protein_alternative_names_dict[row_num,'original_name'])
            #activate boolean?
            #what to do if protein name does not show up in this list and does not have a
            #fasta file?

            protein_name_has_been_changed <- TRUE

            #add protein to a list that will be outputted at the end of the loop

            break

          } #end if(protein_name %in% pnames_as_list)


        } #for(row_num in 1:nrow(pnames_alternatives))

        #if it does, rename the protein to the previous

        #if it does not, put up a message to the
      }


    }

    if((protein_name %in% names(fasta_file)) && (protein_name_has_been_changed == FALSE)){
      #if it's not in the protein dictionary, check if it is in the fasta file


      cat('Is in fasta file but protein name has not been changed yet\n')

      #ask if the user wants to create a separate PDB file for it

      #menu loops function here


    } else {
      #the user has not included it in neither the fasta file nor the dictionary
      #add the whole crosslinking site to a dataframe to be
      #outputted if the xl site does not meet the requirements

      cat('Not in fasta file or in the dictionary\n')

    }

    #should also check protein_pos?

    #can move the is.null if statement within this if statement
    #since it might be helpful to check beforehand if there is a PDB file for
    #the protein

  } #if it does not have a PDB file but it does have a fasta file
  #ask if the user does not if they want to run it in the fasta


  return(protein_name)

}


#----Renumber xiNET file-----

renumber_xinet_from_uniprot_fasta <- function(xinet_file,protein_to_uniprot_id,fasta_file){

  for(row_num in 1:nrow(xinet_file)){
    xinet_file_row <- xinet_file[row_num,]

    #can also replace the 1 with 2 in a cycle
    for(pep_num in 1:2){

      peppos_name <- paste('PepPos',as.character(pep_num),sep='')
      peppos1 <- as.numeric(xinet_file_row[[peppos_name]])

      protein_xi_name <- paste('Protein',as.character(pep_num),sep='')
      protein1 <- as.character(xinet_file_row[[protein_xi_name]])

      pepstring_name <- paste('PepSeq',as.character(pep_num),sep='')
      pepstring <- as.character(xinet_file_row[[pepstring_name]])

      linkpos_name <- paste('LinkPos',as.character(pep_num),sep='')
      linkpos <- as.numeric(as.character(xinet_file_row[[linkpos_name]]))

      uniprot_id <- as.character(protein_to_uniprot_id[protein_to_uniprot_id$protein_id == protein1,'uniprot_id'])

      # if(nchar(uniprot_id) == 0){
      #   cat(protein1)
      # }

      #lab_fasta_seq <- toupper(paste(fasta_file[[protein1]],collapse=''))

      uniprot_fasta_file_name <- paste(uniprot_id,'.fasta',sep='')
      uniprot_fasta <- seqinr::read.fasta(uniprot_fasta_file_name)

      uniprot_fasta_seq <- toupper(paste(uniprot_fasta[[names(uniprot_fasta)]],collapse=''))
      #uniprot_fasta_seq <- toupper((uniprot_fasta[[names(uniprot_fasta)]]))

      pepstring_location <- str_locate_all(uniprot_fasta_seq,pepstring)[[1]]

      if(length(pepstring_location) == 0){
        #no match

        xinet_file[row_num,][[peppos_name]] <- NA

      } else {

        pepstring_location_start <- pepstring_location[1]

        new_peppos1 <- pepstring_location_start + linkpos-1

        #can also do a QC check to make sure that it == 'K'
        #uniprot_fasta_seq <- toupper((uniprot_fasta[[names(uniprot_fasta)]]))
        #uniprot_fasta_seq[732]

        levels(xinet_file[[peppos_name]]) <- c(levels(xinet_file[[peppos_name]]),new_peppos1)
        xinet_file[row_num,][[peppos_name]] <- new_peppos1

      }



      # pwa_uniprot_lab <- pairwiseAlignment(uniprot_fasta_seq,lab_fasta_seq)
      # #re-numbering the sequence based on the
      # pwa_ranges_uniprot_lab <- get_pwa_ranges(pwa_results = pwa_uniprot_lab)
      # pwa_strings_uniprot_lab <- get_pwa_strings(pwa_results = pwa_uniprot_lab)
      #
      #
      # maal <- renumber_and_identify_missing_aa(pwa_strings = pwa_strings_uniprot_lab, start_difference =   pwa_ranges_uniprot_lab$start_difference)
      #
      # #use start and end of strings to get the match of the right index
      # #when lengths are different --> uses occur
      #
      # renumbered_uniprot_sequence_vector <- maal$new_pdb_numbering
      #
      # new_peppos1 <- match(peppos1,renumbered_uniprot_sequence_vector)
      # #new_peppos1 <- renumbered_uniprot_sequence_vector[peppos1]
      #
      # xinet_file[row_num,][[peppos_name]] <- new_peppos1
      #
      #

    } #end for(pep_num in 1:2)


  } #end for(row_num in 1:nrow(xinet_file))


  return(xinet_file)

} #end function renumber_xinet_from_uniprot_fasta()


#----get vector of distances by position frequency----

#' Get Vector of Distances by Position Frequency
#'
#'

get_vector_of_distances_by_pos_freq <- function(xl_dataframe){

  xl_freq_df <- xl_dataframe

  xl_freq_df_no_na_dist_list <- c()

  xl_freq_df_no_na <- xl_freq_df[!is.na(xl_freq_df$Distance),]
  #rownames(xl_freq_df_no_na) <- 1:length(xl_freq_df_no_na)
  for(row_num in 1:nrow(xl_freq_df_no_na)){

    xl_freq_df_no_na_row <- xl_freq_df_no_na[row_num,]

    pos_freq <- as.numeric(xl_freq_df_no_na_row$Position_Frequency)
    dist_freq <- as.numeric(as.character(xl_freq_df_no_na_row$Distance))

    xl_freq_df_no_na_dist_list <- c(xl_freq_df_no_na_dist_list,rep(dist_freq,pos_freq))

  }

  return(xl_freq_df_no_na_dist_list)

}


#-----Count Frequencies in XLMS output----

#'Count Frequencies of Sequences in ppi.analyze() output
#'
#'This function counts the frequencies of sequences from ppi.analyze() output in the case of potential repeats of protein positions
#'
#'@param xlink_df Output from ppi.analyze()
#'@export

ppi.freqCount <- function(xlink_df){

  xl_dataframe <- xlink_df

  #can check if xl_dataframe is a character or already a dataframe
  if(typeof(xl_dataframe) == 'character'){
    xl_dataframe <- read.csv(xl_dataframe)
  }

  pos_freq_list <- c()
  seq_freq_list <- c()
  pos_string_freq_list <- c()
  pro_string_freq_list <- c()
  score_freq_list <- c()
  distance_freq_list <- c()
  pdb_freq_list <- c()

  bs3_colname <- ppi.colnamesConvert(xl_dataframe)
  pro <- as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'Protein', bs3_colname])
  freq <- as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'Frequency', bs3_colname])
  seq <- as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'Sequence', bs3_colname])
  score <- as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'Score', bs3_colname])
  dist <-as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'Distance', bs3_colname])
  pdb <- as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'PDB1', bs3_colname])

  #should be a way of getting the most recently made folder?
  #organize by date? choose the one that starts with make_diff_analysis and is at the top or bottom of the list

  for(protein_string in as.character(xl_dataframe[[pro]])){

    xl_sub_df <- xl_dataframe[xl_dataframe[[pro]] == protein_string,]
    #add up the frequencies
    #"Position Frequency"
    pos_freq <- sum(xl_sub_df[[freq]]) #number of times protein position shows up in all experiments
    pos_freq_list <- c(pos_freq_list,pos_freq)
    #"Sequence Frequency"
    seq_freq <- nrow(xl_sub_df) #number of sequence "types"
    seq_freq_list <- c(seq_freq_list,seq_freq)

    pos_string <- paste(xl_sub_df[[seq]],collapse='+')
    pos_string_freq_list <- c(pos_string_freq_list,pos_string)

    pro_string_freq_list <- c(pro_string_freq_list,protein_string)

    #lowest score?
    score_f <- min(xl_sub_df[[score]])
    score_freq_list <- c(score_freq_list,score_f)

    #distance_freq_list <- c(distance_freq_list,as.character(xl_dataframe[[dist]]))
    distance_freq_list <- c(distance_freq_list,unique(as.character(xl_sub_df[[dist]])))

    #have the PDB here
    if(is.na(unique(xl_sub_df[[dist]]))){
      #if is
      pdb_freq_list <- c(pdb_freq_list,NA)

    } else {#end if(is.na(xl_sub_df[[dist]])){

      pdb_freq_list <- c(pdb_freq_list,strsplit(unique(xl_sub_df[[pdb]]),'_')[[1]][1])
    }
    #list of the sequence types?
    #columns: score,

  }

  #create list() here of all of the lists
  xl_freq_list <- list(Protein_String=pro_string_freq_list,
                       Sequence_Frequency=seq_freq_list,
                       Position_Frequency=pos_freq_list,
                       Sequence_Strings=pos_string_freq_list,
                       Lowest_Score=score_freq_list,
                       Distance=distance_freq_list,
                       PDB=pdb_freq_list)

  #return(xl_freq_list)

  xl_freq_df <- data.frame(xl_freq_list)

  #need to remove repeating rows in the dataframe

  xl_freq_df <- xl_freq_df[!duplicated(xl_freq_df),]


  return(xl_freq_df)
}


#-----Filter XLink DF by Protein Names----

#need to have a dataframe or something that contains each of the lists to be able to easily change
#between the two

colnames1 <- c("seq","pro","dist","freq","freq_color","files","pdb1","pdb2",
               "pro_pos1","pro_pos2","pro_name1","pro_name2","pep_seq1","pep_seq2",
               "pep_pos1","pep_pos2","score")
colnames2 <-  c('Sequence', 'Protein', 'Distance', 'Frequency', "Frequency.Color",
                'Files', 'PDB1', 'PDB2','Protein.Position.1','Protein.Position.2',
                'Protein.Name.1','Protein.Name.2','Peptide.Sequence.1','Peptide.Sequence.2',
                'Peptide.Position.1','Peptide.Position.2','Score')
colnames3 <-  c('Sequence', 'Protein', 'Distance', 'Frequency', "Category.Color",
                'Files', 'PDB1', 'PDB2','Protein.Position.1','Protein.Position.2',
                'Protein.Name.1','Protein.Name.2','Peptide.Sequence.1','Peptide.Sequence.2',
                'Peptide.Position.1','Peptide.Position.2','Score')

bs3_colnames_converter <- data.frame(colnames1=colnames1,colnames2=colnames2,colnames3=colnames3)

#save data to be able to loaded by the function

ppi.filter <- function(xlink_df,filter_vars,filter_by=c('frequency','score','proteins'),condition='and'){

  #get the filter_vars
  #filter vars should be in the format of a list() corresponding to filter_by variable or a regular
  #c() list that is equal to the length of the filter_by list

  if(typeof(filter_vars) == 'list'){
    #it's a list -->

  } else {
    #it's a character, double --> c()
    #check to make sure the length
    length(filter_vars) == length(filter_by)

  }



} #end function ppi.filter


ppi.colnamesConvert <- function(xl_dataframe){

  bs3_colname <- NULL

  for(col_name in colnames(bs3_colnames_converter)){
    tf_unique <- unique(colnames(xl_dataframe) %in% bs3_colnames_converter[[col_name]])
    if(tf_unique == TRUE){
      bs3_colname <- col_name
      break
    }
  }

  return(bs3_colname)

} #end function ppi.colnamesConvert




filter_xlink_df_by_protein_names <- function(xl_dataframe,list_of_protein_names){

  #turn this into a function?
  bs3_colname <- NULL

  for(col_name in colnames(bs3_colnames_converter)){
    tf_unique <- unique(colnames(xl_dataframe) %in% bs3_colnames_converter[[col_name]])
    if(tf_unique == TRUE){
      bs3_colname <- col_name
      break
    }
  }

  #get the row associated with the actual name --> convert to the appropriate colname column

  pro_name1 <- as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'Protein.Name.1',bs3_colname])
  pro_name2 <- as.character(bs3_colnames_converter[bs3_colnames_converter$colnames2 == 'Protein.Name.2',bs3_colname])


  core_subunit_names <- list_of_protein_names

  boolean_selection_list <- c()

  #go through each row
  for(row_num in 1:nrow(xl_dataframe)){

    xl_sub_df <- xl_dataframe[row_num,]
    if((as.character(xl_sub_df[[pro_name1]]) %in% core_subunit_names) && (as.character(xl_sub_df[[pro_name2]]) %in% core_subunit_names)){

      #if it's in both --> add to list?
      #can also make a boolean list that will extract all of the correct lists

      boolean_selection_list <- c(boolean_selection_list,TRUE)

    } else {

      boolean_selection_list <- c(boolean_selection_list,FALSE)

    }

  }

  xl_df_core_pro <- xl_dataframe[boolean_selection_list,]

  return(xl_df_core_pro)

}

#-----Generate Random Lysine Distances in PDB----

#what is the frequency vector?? need this info for documentation
generate_random_lysine_distances_in_pdb <- function(pdb_id,frequency_vector,chains=NULL){

  #pdb_6c23 <- check_download_read_pdb(pdb_id = pdb_id)
  pdb_6c23 <- read.pdb2(pdb_id)

  pdb_6c23_atom_filtered <- pdb_6c23$atom[pdb_6c23$atom$elety == 'CA',]
  pdb_6c23_atom_filtered <- pdb_6c23_atom_filtered[pdb_6c23_atom_filtered$resid == 'LYS',]

  if(!is.null(chains)){

    pdb_6c23_atom_filtered <- pdb_6c23_atom_filtered[pdb_6c23_atom_filtered$chain %in% chains,]

  }

  length_of_freq_vector <- length(frequency_vector)

  random_row_vector_list <- list()

  for(vector_num in 1:2){

    random_row_vector_list[[vector_num]] <- sample(1:nrow(pdb_6c23_atom_filtered),length_of_freq_vector,replace = TRUE)

  }

  lysine_random_xyz_distances <- c()

  for(freq_vec_num in 1:length_of_freq_vector){

    #get info from both random_row_vectors

    xyz_coordinates <- list()

    for(vector_num in 1:length(random_row_vector_list)){

      random_row_vector <- random_row_vector_list[[vector_num]]

      row_num <- random_row_vector[freq_vec_num]

      pdb_6c23_atom_filtered_row <- pdb_6c23_atom_filtered[row_num,]

      xyz_coordinates[[vector_num]] <- c(pdb_6c23_atom_filtered_row$x,
                                         pdb_6c23_atom_filtered_row$y,
                                         pdb_6c23_atom_filtered_row$z)


    } #end for(vector_name in names(random_row_vector_list))

    #calculate the distance here

    #calculate the xyz coordinates
    calculated_dist <- dist.xyz(xyz_coordinates[[1]],xyz_coordinates[[2]])[[1]]
    lysine_random_xyz_distances <- c(lysine_random_xyz_distances,calculated_dist)

  } #end for(row_num in 1:length_of_freq_vector)

  #return histogram or list?
  return(lysine_random_xyz_distances)

} #end function generate_random_lysine


#----Renumber Binding Sequence----


#'Renumber binding site df from uniprot fasta
#'
#'Renumber binding site df from uniprot fasta
#'@param binding_site_df binding_site_df
#'@param protein_to_uniprot_id protein_to_uniprot_id
#'@param fasta_file fasta_file
#'@param fasta_file_directory fasta_file_directory
#'@param use_pipeline_output_df use_pipeline_output_df
#'@export

renumber_binding_site_df_from_uniprot_fasta <- function(binding_site_df,protein_to_uniprot_id,fasta_file,
                                                        fasta_file_directory = NULL,
                                                        use_pipeline_output_df = FALSE){

  #rename the columns?

  if(use_pipeline_output_df == TRUE){
    binding_site_df[c('Start','End')] <- NA
  }


  for(row_num in 1:nrow(binding_site_df)){
    binding_site_df_row <- binding_site_df[row_num,]

    protein1 <- as.character(binding_site_df_row[['Protein.Name']])

    if(use_pipeline_output_df == TRUE){
      pepstring <- as.character(binding_site_df_row[['Sequence']])
      #add a new column for the start and end positions


    } else {
      pepstring <- as.character(binding_site_df_row[['Binding.Sequence']])
    }


    uniprot_id <- as.character(protein_to_uniprot_id[protein_to_uniprot_id$protein_id == protein1,'uniprot_id'])

      # if(nchar(uniprot_id) == 0){
      #   cat(protein1)
      # }

      #lab_fasta_seq <- toupper(paste(fasta_file[[protein1]],collapse=''))

    if(is.null(fasta_file_directory)){
      fasta_file_directory <- getwd()
    }

    uniprot_fasta_file_name <- paste(fasta_file_directory,'/',uniprot_id,'.fasta',sep='')
    uniprot_fasta <- seqinr::read.fasta(uniprot_fasta_file_name)

    uniprot_fasta_seq <- toupper(paste(uniprot_fasta[[names(uniprot_fasta)]],collapse=''))

    pepstring_location <- str_locate_all(uniprot_fasta_seq,pepstring)[[1]]

      if(length(pepstring_location) == 0){
        #no match

        if(use_pipeline_output_df == FALSE){
          binding_site_df[row_num,]$Binding.Site.Start..in.PDB.File. <- NA
          binding_site_df[row_num,]$Binding.Site.End..in.PDB.File. <- NA
        }

      } else {

        pepstring_location_start <- pepstring_location[1]

        pepstring_location_end <- pepstring_location[2]

        if(use_pipeline_output_df == TRUE){
          binding_site_df[row_num,]$Start <- pepstring_location_start
          binding_site_df[row_num,]$End <- pepstring_location_end


        } else {
          binding_site_df[row_num,]$Binding.Site.Start..in.PDB.File. <- pepstring_location_start
          binding_site_df[row_num,]$Binding.Site.End..in.PDB.File. <- pepstring_location_end

        }


      }


  } #end for(row_num in 1:nrow(binding_site_df))

  return(binding_site_df)

} #end function renumber_xinet_from_uniprot_fasta()

#----Get binding Sequence-----

#ms_sequence <- 'EVSTAPAGTDMPAAK'
#sequence_for_alignment <- toupper(paste(fasta_file$EED_O75530,collapse=''))
#last_aa <- 'K'


#need to make it so it will also accept the input eluate table as well
#can make a wrapper for the function so that people can use the input eluate table for it
#for RBDmap data --> should be part of a list

#' Get Binding Sequence from RBDmap Data
#'
#' This funtion gets the binding sequence from RBDmap data
#'
#' @param ms_sequence The input sequence from MS, should be loaded as a string
#' @param protease The protease used for the experiment, either 'ArgC' or 'LysC'
#' @param sequence_for_alignment The sequence to find the binding sequence in. The ms_sequence must be within the sequence_for_alignment or NULL will be returned.
#' @param protein_name Name of the protein to which the sequence_for_alignment belongs. Defaults to NA
#' @param database_name Name of the database used: 'PDB','UniProt' or 'FASTA' (for the input FASTA sequence). Defaults to NA
#' @param database_id Identifier for the database. Defaults to NA
#' @param cleave_offset Cleave offset for the sequence. Defaults to 0
#' @export


rbd.getBindingSeq <- function(ms_sequence,protease,sequence_for_alignment,
                                 protein_name = NA,database_name = NA,
                                 database_id = NA, cleave_offset=0,last_aa = NULL,
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


  #will need to do this again if the n_terminal boolean is activated

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

  #set the
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
  } #end while((current_aa != end_binding_aa) && (current_aa_index <= length(sequence_for_alignment_split)) && (current_aa_index > 0)){

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

    fasta_seq <- paste(toupper(fasta_file$`EZH2_Q15910-2`),collapse='')
    list_of_cleaves <- str_locate_all(fasta_seq,'K')[[1]][,1]
    #find the next one depending on the start/end of the binding sequence
    #can also expand this for originally finding the binding sequence

    start_pattern <- 40
    end_pattern <- 80

    list_of_cleaves_sub <- list_of_cleaves[start_pattern > list_of_cleaves]

    cleave_site_found <- FALSE

    last_cleave <- tail(list_of_cleaves_sub,n=1) #last integer
    last_cleave2 <- tail(list_of_cleaves_sub,n=2)[1] #second to last integer

    #can have boolean out here that is false
    #while the boolean is false --> the function will keep searching for the nearby
    #cleave sites

    #turn this into a function?


    #put inside a while loop?
    if((last_cleave-last_cleave2) < cleave_offset){
      #if TRUE --> the difference between the two is less than

    }

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


  #return(binding_site_ouput)
  return(data.frame(binding_site_ouput))
}

#----Count UV/NoXL experiments in input/eluate dataframes----

count_input_eluate_in_supp_df <- function(input_dataframe){

  #input_dataframe_og <- input_dataframe

  input_dataframe2 <- input_dataframe

  empty_uv_df <- data.frame(UV_count=rep(NA,nrow(input_dataframe)),NoXL_count=rep(NA,nrow(input_dataframe)))

  input_dataframe <- cbind(input_dataframe,empty_uv_df)

  # input_dataframe2[c('UV_count','NoXL_count')] <- rep(NA,nrow(input_dataframe))
  #
  # if(('UV_count' %in% colnames(input_dataframe)) && ('NoXL_count' %in% colnames(input_dataframe))){
  #   input_dataframe[c('UV_count','NoXL_count')] <- NA
  # }


  input_dataframe_og <- input_dataframe[,!(colnames(input_dataframe) %in% c('UV_count','NoXL_count'))]

  for(row_num in 1:nrow(input_dataframe_og)){
    #filter which ones are

    noxl_cols <- grepl('-',colnames(input_dataframe_og))
    uv_cols <- !noxl_cols
    #colnames(input_dataframe)[uv_cols]

    values_in_row <- (t(input_dataframe_og[row_num,]))
    uv_values <- values_in_row[uv_cols]

    uv_count <- length(uv_values[uv_values > 0])

    input_dataframe[row_num,]$UV_count <- uv_count

    noxl_values <- values_in_row[noxl_cols]

    noxl_count <- length(noxl_values[noxl_values > 0])

    input_dataframe[row_num,]$NoXL_count <- noxl_count


    #add uv and noxl count

  }

  return(input_dataframe)

}

#----Make Proteins and Intensities List-----

#'Make Proteins and Intensities List
#'
#'Make Proteins and Intensities List
#'@param experiment_prefixes experiment_prefixes
#'@param protein_to_uniprot_id protein_to_uniprot_id
#'@param fasta_file fasta_file
#'@param fasta_file_directory fasta_file_directory
#'@param identifier_is_binding_seq identifier_is_binding_seq
#'@export

make_proteins_and_intensity_list_from_pipeline_output <- function(experiment_prefixes,
                                                                  protein_to_uniprot_id,
                                                                  fasta_file,
                                                                  fasta_file_directory,
                                                                  identifier_is_binding_seq = FALSE){

  bs_output_df <- data.frame(stringsAsFactors = FALSE)

  experiment_number <- 0
  experiment_number_list <- c()
  uv_or_noxl_list <- c()
  protease_list <- c()
  proteins_and_intensity_list <- list()

  for(e_prefix in experiment_prefixes){

    #binding_site_df <- read.csv(paste(e_prefix,'_binding_site_sequence_data_frame.csv',sep=''))
    pipeline_output_df <- read.csv(paste(e_prefix,'_pipeline_output_file.csv',sep=''))

    #need to change this so that it uses the pipeline output file instead of the binding site
    #output file


    #rownames(pipeline_output_df) <- pipeline_output_df$X
    #pipeline_output_df$X <- NULL
    #colnames(pipeline_output_df) <- c("sequence","input","eluate","category","protein","aa_before","aa_after")

    binding_site_df <- renumber_binding_site_df_from_uniprot_fasta(binding_site_df = pipeline_output_df,
                                                                   protein_to_uniprot_id = protein_to_uniprot_id,
                                                                   fasta_file = fasta_file,
                                                                   fasta_file_directory = fasta_file_directory,
                                                                   use_pipeline_output_df = TRUE)

    rownames(binding_site_df) <- binding_site_df$X
    binding_site_df$X <- NULL
    #colnames(binding_site_df) <- c("sequence","input","eluate","category","protein","aa_before","aa_after",'start','end')

    experiment_number <- experiment_number + 1

    #add in columns to the pipeline_output_df to get the binding sites --> will have to put
    #the ArgC and LysC detection in its own function

    for(row_num in 1:nrow(binding_site_df)){
      input_list <- c()
      eluate_list <- c()
      #noxl_count <- 0
      #uv_count <- 0

      bs_df_row <- binding_site_df[row_num,]

      bs_pro_name <- as.character(bs_df_row$Protein.Name)
      bs_seq <- as.character(bs_df_row$Sequence)
      #bs_ms_pep_seq <- as.character(bs_df_row$MS.MS.Peptide.Sequence)
      bs_category <- as.character(bs_df_row$Category)
      bs_start <- as.character(bs_df_row$Start)
      bs_end <- as.character(bs_df_row$End)

      #add each of these to the list
      bs_identifier <- paste(bs_pro_name,'_',bs_seq,'_',bs_start,'-',bs_end,sep='')

      filtered_po_df <- bs_df_row

      #filtered_po_df <- pipeline_output_df[pipeline_output_df$Protein.Name == bs_pro_name,]
      #filtered_po_df <- filtered_po_df[filtered_po_df$Sequence == bs_ms_pep_seq,]
      #filtered_po_df <- filtered_po_df[filtered_po_df$Category == bs_category,]

      if(nrow(filtered_po_df) == 0){
        cat("no matches found")
        return(NULL)
      } else if(nrow(filtered_po_df) == 1){


      } else {
        #multiple matches found
        filtered_po_df <- filtered_po_df[1,]
      }

      #put each one of these in their own tables
      input_value <- filtered_po_df$Input
      eluate_value <- filtered_po_df$Eluate

      category_filtered <- as.character(filtered_po_df$Category)
      category_filtered_split <- strsplit(category_filtered,'_')[[1]]
      f_uv <- category_filtered_split[1]
      f_protease <- category_filtered_split[2]

      experiment_number_list <- c(experiment_number_list,as.character(experiment_number))
      protease_list <- c(protease_list,f_protease)

      #need to first check to see if the grouping already exists
      #or can just delete any duplicates afterwards

      if(f_uv == 'NoXL'){

        #add '-' to some list
        uv_or_noxl <- '-'
        uv_or_noxl_list <- c(uv_or_noxl_list,'-')
        #add to some count for that protein (in list() form probably)

      } else if(f_uv == 'UV' || f_uv == 'UVCD' || f_uv == 'UVMA'){

        #add '+' to some list
        uv_or_noxl <- '+'
        uv_or_noxl_list <- c(uv_or_noxl_list,'+')
        #add to some count for that protein (in list() form probably)

      } else {
        cat(paste('Unknown UV category',f_uv))
        uv_or_noxl <- '?'
        uv_or_noxl_list <- c(uv_or_noxl_list,'?')
        #could add a '?' to this category just so the code will continue
      } #end control flow that starts with if(f_uv == 'NoXL')

      #need to go through the other experiments and find the right input/eluate?
      #check if that shows up in another experiment?

      #add all of the relevant data to the correct identifier

      #can have two similanteous lists: one that keeps track of the pattern for the
      #header and one that keeps track of the intensities that actually show up in the dataframe
      #then create a function that will combine all of the information together in a dataframe


      #bs_df_row <- make_binding_site_df(pdb_info,bs_df_row,fasta_seq)

      uniprot_id <- as.character(protein_to_uniprot_id[protein_to_uniprot_id$protein_id == bs_pro_name,'uniprot_id'])

      if(is.null(fasta_file_directory)){
        fasta_file_directory <- getwd()
      }

      uniprot_fasta_file_name <- paste(fasta_file_directory,'/',uniprot_id,'.fasta',sep='')
      uniprot_fasta <- seqinr::read.fasta(uniprot_fasta_file_name)

      uniprot_fasta_seq <- toupper(paste(uniprot_fasta[[names(uniprot_fasta)]],collapse=''))

      protease <- strsplit(as.character(bs_df_row$Category),'_')[[1]][2]

      bs_output <- get_binding_sequence(ms_sequence = as.character(bs_df_row$Sequence),
                                        protease = protease,
                                        last_aa = as.character(bs_df_row$Amino.Acid.After),
                                        sequence_for_alignment = uniprot_fasta_seq)


      if(identifier_is_binding_seq == TRUE){

        # as.character(bs_output$binding_site_start)
        # as.character(bs_output$binding_site_end)
        # as.character(bs_output$binding_sequence)

        bs_identifier <- paste(bs_pro_name,'_',as.character(bs_output$binding_sequence),
                               '_',as.character(bs_output$binding_site_start),'-',
                               as.character(bs_output$binding_site_end),sep='')

      }

      #can make a true/false statement here to make the identifier the binding sequence instead
      #of the

      bs_output$experiment_number <- as.character(experiment_number)
      bs_output$uv_or_noxl <- uv_or_noxl
      bs_output$protease <- f_protease
      bs_output$input <- input_value
      bs_output$eluate <- eluate_value
      bs_output$pro_name <- bs_pro_name
      bs_output$experiment_prefix <- experiment_prefixes[experiment_number]

      bs_output_df <- rbind(bs_output_df,data.frame(bs_output))

      if(bs_identifier %in% names(proteins_and_intensity_list)){
        #if the identifier alrady exists in the list --> need to add to the list instead of
        #making a new one

        proteins_and_intensity_list[[bs_identifier]]$experiment_number <- c(proteins_and_intensity_list[[bs_identifier]]$experiment_number,
                                                                            as.character(experiment_number))
        proteins_and_intensity_list[[bs_identifier]]$uv_or_noxl <- c(proteins_and_intensity_list[[bs_identifier]]$uv_or_noxl,
                                                                     uv_or_noxl)
        proteins_and_intensity_list[[bs_identifier]]$protease <- c(proteins_and_intensity_list[[bs_identifier]]$protease,
                                                                   f_protease)
        proteins_and_intensity_list[[bs_identifier]]$input <- c(proteins_and_intensity_list[[bs_identifier]]$input,
                                                                input_value)
        proteins_and_intensity_list[[bs_identifier]]$eluate <- c(proteins_and_intensity_list[[bs_identifier]]$eluate,
                                                                 eluate_value)

        proteins_and_intensity_list[[bs_identifier]]$pro_name <- c(proteins_and_intensity_list[[bs_identifier]]$pro_name,
                                                                   bs_pro_name)

        proteins_and_intensity_list[[bs_identifier]]$pep_seq <- c(proteins_and_intensity_list[[bs_identifier]]$pep_seq,
                                                                  bs_seq)

        proteins_and_intensity_list[[bs_identifier]]$pep_start <- c(proteins_and_intensity_list[[bs_identifier]]$pep_start,
                                                                    bs_start)

        proteins_and_intensity_list[[bs_identifier]]$pep_end <- c(proteins_and_intensity_list[[bs_identifier]]$pep_end,
                                                                  bs_end)

        proteins_and_intensity_list[[bs_identifier]]$bs_seq <- c(proteins_and_intensity_list[[bs_identifier]]$bs_seq,
                                                                 bs_output$binding_sequence)

        proteins_and_intensity_list[[bs_identifier]]$bs_start <- c(proteins_and_intensity_list[[bs_identifier]]$bs_start,
                                                                   bs_output$binding_site_start)

        proteins_and_intensity_list[[bs_identifier]]$bs_end<- c(proteins_and_intensity_list[[bs_identifier]]$bs_end,
                                                                bs_output$binding_site_end)
        #paste(bs_pro_name,'_',bs_seq,'_',bs_start,'-',bs_end,sep='')

      } else {
        #initialize the identifier in the list

        proteins_and_intensity_list[[bs_identifier]] <- list(experiment_number=as.character(experiment_number),
                                                             uv_or_noxl=uv_or_noxl,
                                                             protease=f_protease,
                                                             input=input_value,
                                                             eluate=eluate_value,
                                                             pro_name=bs_pro_name,
                                                             pep_seq=bs_seq,
                                                             pep_start=bs_start,
                                                             pep_end=bs_end,
                                                             bs_seq=bs_output$binding_sequence,
                                                             bs_start=bs_output$binding_site_start,
                                                             bs_end=bs_output$binding_site_end)

      } #end if(bs_identifier %in% names(proteins_and_intensity_list))

      #proteins_and_intensity_list[bs_identifier]


      #list names: "experiment_number", "uv", "protease", "key_for_protein_start_end"


      #make all of these into lists and then transpose the lists to make them
      #into the format that CD wants


      #keep track of the experiment number and the UV status


      #split the protein name for the first part of the

      #get info from

      #match sequence, protein name, and category

      #need to use ms sequence to match with pipeline output file
      #output the binding sequence to

    } #end for(row_num in nrow(binding_site_df))



  } #end for(e_prefix in experiment_prefixes)


  supp_header_df <- list(experiment_number=experiment_number_list,
                         uv_or_noxl=uv_or_noxl_list,
                         protease=protease_list)


  supp_header_df <- data.frame(supp_header_df)
  supp_header_df <- supp_header_df[!duplicated(supp_header_df),]


  return(list(proteins_and_intensity_list=proteins_and_intensity_list,
              supp_header_df=supp_header_df,
              bs_output_df=bs_output_df))

} #end make_proteins_and_intensity_list_from_pipeline_output function

#----Make input and eluate dataframes from proteins_and_intensities----

#proteins_and_intensity_list_with_supp_header <- proteins_and_intensity_list_with_supp_header2

make_input_and_eluate_dfs_from_pi_list <- function(proteins_and_intensity_list_with_supp_header){

  supp_header_df <- proteins_and_intensity_list_with_supp_header$supp_header_df
  proteins_and_intensity_list <- proteins_and_intensity_list_with_supp_header$proteins_and_intensity_list

  input_dataframe <- list()
  eluate_dataframe <- list()

  supp_header_identifier_list <- c()
  for(supp_row_num in 1:nrow(supp_header_df)){
    supp_header_df_row <- supp_header_df[supp_row_num,]
    #can do a similar filtering step
    supp_header_identifier <- as.character(apply(supp_header_df_row, 1, paste, collapse="_"))
    supp_header_identifier_list <- c(supp_header_identifier_list,supp_header_identifier)
  }


  #maybe should go by protein instead and then transpose the dataframe to get the right orientation
  for(protein_binding_identifier in names(proteins_and_intensity_list)){

    pi_list_row  <- proteins_and_intensity_list[[protein_binding_identifier]]
    pi_list_row_df <- data.frame(pi_list_row)
    #levels(pi_list_row_df$pep_seq)
    #could have this be different depending on a boolean (such as binding sequence instead)

    #get the levels of the difference peptide sequences to get the different rows
    for(peptide_level in levels(pi_list_row_df$pep_seq)){

      #should add peptide level to identifier

      peptide_level_identifier <- paste(protein_binding_identifier,'_',peptide_level,sep='')
      input_list_for_df <- c()
      eluate_list_for_df <- c()
      #input_dataframe[[peptide_level_identifier]] <- input_list_for_df

      pi_list_row_df_filtered <- pi_list_row_df[pi_list_row_df$pep_seq == peptide_level,]
      #go through each row in the supp_header_df to see if the values are in input or eluate
      #then add to each of the input and eluate df lists

      for(supp_row_num in 1:nrow(supp_header_df)){
        supp_header_df_row <- supp_header_df[supp_row_num,]
        #can do a similar filtering step
        #supp_header_identifier <- as.character(apply(supp_header_df_row, 1, paste, collapse="_"))
        #supp_header_identifier_list <- c(supp_header_identifier_list,supp_header_identifier)

        #should assign to a new variable so that it does not get destroyed in for loop
        pi_list_row_df_filtered2 <- pi_list_row_df_filtered
        for(supp_header_name in names(supp_header_df_row)){

          pi_list_row_df_filtered2 <- pi_list_row_df_filtered2[pi_list_row_df_filtered2[[supp_header_name]] == as.character(supp_header_df_row[[supp_header_name]]),]

        }

        if(nrow(pi_list_row_df_filtered2) == 0){
          #add zeros to the input and eluate lists
          input_list_for_df <- c(input_list_for_df,as.numeric(0))
          eluate_list_for_df <- c(eluate_list_for_df,as.numeric(0))


        } else {
          #should check to make sure that it really only equals == 1
          input_list_for_df <- c(input_list_for_df,as.numeric(pi_list_row_df_filtered2$input))
          eluate_list_for_df <- c(eluate_list_for_df,as.numeric(pi_list_row_df_filtered2$eluate))


        }


      } #end for supp_row_num in 1:nrow(supp_header_df)



      input_dataframe[[peptide_level_identifier]] <- input_list_for_df
      eluate_dataframe[[peptide_level_identifier]] <- eluate_list_for_df


    } #end for(peptide_level in levels(pi_list_row_df$pep_seq))




  } #end for(protein_identifier in names(proteins_and_intensity_list))

  #should also make a list of the peptide identifiers for renaming the columns?
  input_dataframe <- data.frame(input_dataframe)
  eluate_dataframe <- data.frame(eluate_dataframe)

  input_dataframe <- t(input_dataframe)
  eluate_dataframe <- t(eluate_dataframe)

  colnames(input_dataframe) <- supp_header_identifier_list
  colnames(eluate_dataframe) <- supp_header_identifier_list

  input_dataframe <- count_input_eluate_in_supp_df(input_dataframe)
  eluate_dataframe <- count_input_eluate_in_supp_df(eluate_dataframe)

  input_and_eluate_dfs <- list(input_df=input_dataframe,
                               eluate_df=eluate_dataframe)

  return(input_and_eluate_dfs)


  # for(row_num in 1:nrow(supp_header_df)){
  #   supp_header_df_row <- supp_header_df[row_num,]
  #
  #   supp_header_identifier <- as.character(apply(supp_header_df_row, 1, paste, collapse="_"))
  #   supp_header_identifier_list <- c(supp_header_identifier_list,supp_header_identifier)
  #
  #   input_list_for_df <- c()
  #   eluate_list_for_df <- c()
  #
  #   #go through rows of pro_int_list here
  #
  #   for(protein_binding_identifier in names(proteins_and_intensity_list)){
  #
  #     #SUZ12_Q15022_PLATR_376-380
  #     pro_int_list_row <- proteins_and_intensity_list[[protein_binding_identifier]]
  #     pro_int_list_row_df <- data.frame(pro_int_list_row)
  #     #can make this into a dataframe structure and then check the number of
  #     #matching rows
  #     #can do a subsequent for loop to go through each of the
  #     for(supp_header_name in names(supp_header_df_row)){
  #
  #       #pro_int_list_row[[supp_header_name]]
  #       #need to account for the fact that there may be more than one
  #
  #       #starts off as TRUE and then will change to FALSE if any of the statements are FALSE
  #       #narrow down the list
  #       pro_int_list_row_df <- pro_int_list_row_df[pro_int_list_row_df[[supp_header_name]] == as.character(supp_header_df_row[[supp_header_name]]),]
  #
  #       if(nrow(pro_int_list_row_df) == 0){ #break the cycle if it has been narrowed to 0
  #         break
  #       }
  #     }
  #
  #     if(nrow(pro_int_list_row_df) == 0){
  #       #make the values 0
  #     } else {
  #       #do a for loop and go through each of the rows
  #       #need to add a row??
  #
  #
  #     }
  #
  #     #not getting both of the numbers here
  #     #maybe should break if it's T for one of them
  #     #would it be too long if it was happening more than once?
  #
  #     input_list_for_df_mini <- c()
  #     eluate_list_for_df_mini <- c()
  #
  #     for(index_num in 1:length(pro_int_list_row$experiment_number)){
  #
  #       all_supp_headers_match <- TRUE
  #
  #       for(supp_header_name in names(supp_header_df_row)){
  #
  #         #pro_int_list_row[[supp_header_name]]
  #         #need to account for the fact that there may be more than one
  #
  #         #starts off as TRUE and then will change to FALSE if any of the statements are FALSE
  #         if(!(as.character(supp_header_df_row[[supp_header_name]]) == pro_int_list_row[[supp_header_name]][index_num])){
  #           all_supp_headers_match <- FALSE
  #
  #         }
  #
  #
  #
  #         #can break if no match
  #         #if no match --> assign value of zero
  #         #will need to put into two different tables for input and eluate
  #       }
  #
  #       if(all_supp_headers_match == TRUE){
  #         #if all headers match --> add the input and eluate values to the appropriate
  #         #dataframes
  #
  #         #can add boolean here so if input_or_eluate == 'eluate' --> get that column
  #         #or just have the input_or_eluate as the the selected column in the rows to add
  #
  #         input_list_for_df_mini <- c(input_list_for_df_mini,pro_int_list_row$input[index_num])
  #         eluate_list_for_df_mini <- c(eluate_list_for_df_mini,pro_int_list_row$eluate[index_num])
  #
  #         #input_list_for_df <- c(input_list_for_df,pro_int_list_row$input[index_num])
  #         #eluate_list_for_df <- c(eluate_list_for_df,pro_int_list_row$eluate[index_num])
  #
  #       } else {
  #         #if one of the headers does not match --> put 0 for the input and eluate value
  #
  #         input_list_for_df_mini <- c(input_list_for_df_mini,0)
  #         eluate_list_for_df_mini <- c(eluate_list_for_df_mini,0)
  #
  #         #input_list_for_df <- c(input_list_for_df,0)
  #         #eluate_list_for_df <- c(eluate_list_for_df,0)
  #
  #       }
  #
  #
  #
  #       #Add a 'break' statement if it's TRUE down here still --> don't need to check again
  #
  #
  #     } #end for(index_num in length(pro_int_list_row$experiment_number))
  #
  #     if(TRUE %in% (input_list_for_df_mini > 0)){
  #       #at least one of the values is not zero
  #
  #       #could it be that 2 values are being added here sometimes --> should try to reflect this
  #       #better in the
  #       #need to account for the different pep_seqs
  #       #maybe do a round for one of the pep_seqs then do the next one --> BOOLEAN
  #
  #       input_list_for_df <- c(input_list_for_df,input_list_for_df_mini[(input_list_for_df_mini > 0)])
  #
  #       if(length(input_list_for_df_mini[(input_list_for_df_mini > 0)]) >= 2){
  #         print(protein_binding_identifier)
  #       }
  #
  #     } else {
  #       #add zero to list
  #       input_list_for_df <- c(input_list_for_df,0)
  #
  #     }
  #
  #
  #     if(TRUE %in% (eluate_list_for_df_mini > 0)){
  #       #at least one of the values is not zero
  #
  #       eluate_list_for_df <- c(eluate_list_for_df,eluate_list_for_df_mini[(eluate_list_for_df_mini > 0)])
  #
  #     } else {
  #       #add zero to list
  #       eluate_list_for_df <- c(eluate_list_for_df,0)
  #
  #     }
  #
  #     #can check here if anything was added to the list and if not --> add a 0
  #
  #
  #   } #end for(protein_binding_identifier in names(proteins_and_intensity_list))
  #
  #   #add the input and eluate lists
  #   input_dataframe[[supp_header_identifier]] <- input_list_for_df
  #   eluate_dataframe[[supp_header_identifier]] <- eluate_list_for_df
  #
  # } #end for row_num in 1:nrow(supp_header_df)
  #
  # #can make this into a function as well that can be used to generate both the input and eluate
  #
  # # for(row_name in names(input_dataframe)){
  # #   print(length(input_dataframe[[row_name]]))
  # #
  # # }
  #
  # #differing row lengths between the
  #
  # #can make the actual dataframe here
  # #use the names of the proteins and intensity list for the row names in the input and eluate
  # #dataframes
  # #names(proteins_and_intensity_list)
  #
  # input_dataframe <- data.frame(input_dataframe)
  # rownames(input_dataframe) <- names(proteins_and_intensity_list)
  # colnames(input_dataframe) <- supp_header_identifier_list
  #
  # eluate_dataframe <- data.frame(eluate_dataframe)
  # rownames(eluate_dataframe) <- names(proteins_and_intensity_list)
  # colnames(eluate_dataframe) <- supp_header_identifier_list
  #
  # input_dataframe <- count_input_eluate_in_supp_df(input_dataframe)
  # eluate_dataframe <- count_input_eluate_in_supp_df(eluate_dataframe)
  #
  # input_and_eluate_dfs <- list(input_df=input_dataframe,
  #                              eluate_df=eluate_dataframe)
  #
  # return(input_and_eluate_dfs)
  #
} #end function make_input_and_eluate_dfs_from_pi_list


#----Make binding site rows for supp df table-----

#'make_binding_site_rows_for_supp_df
#'
#'make_binding_site_rows_for_supp_df
#'@param proteins_and_intensity_list proteins_and_intensity_list
#'@param identifier_is_binding_seq identifier_is_binding_seq
#'@export

make_binding_site_rows_for_supp_df <- function(proteins_and_intensity_list,identifier_is_binding_seq = FALSE){

  protein_names_list <- c()
  binding_sequences_list <- c()
  binding_seq_start_list <- c()
  binding_seq_end_list <- c()
  uv_count_list <- c()
  noxl_count_list <- c()
  binding_sequences_list_real <- c()
  binding_seq_start_list_real <- c()
  binding_seq_end_list_real <- c()



  #proteins_and_intensity_list[['AEBP2_Q6ZN18_HMLTHSGDKPFK_318-329']]

  for(protein_binding_identifier in names(proteins_and_intensity_list)){

    #protein_binding_identifier <- names(proteins_and_intensity_list)[10]
    pro_int_list_row <- proteins_and_intensity_list[[protein_binding_identifier]]

    #may need to split this up further

    if(identifier_is_binding_seq == TRUE){

      #establish variable to go through levels
      #may not want to do it yet for FALSE since it has not been set up for that yet
      pro_int_list_row_df1 <- data.frame(pro_int_list_row)

      for(peptide_level in levels(pro_int_list_row_df1$pep_seq)){
        #do basically the same thing as the else but

        pro_int_list_row_df <- pro_int_list_row_df1
        #pro_int_list_row_df <- data.frame(pro_int_list_row)
        pro_int_list_row_df <- pro_int_list_row_df[pro_int_list_row_df$pep_seq == peptide_level,]

        protein_names_list <- c(protein_names_list,as.character(pro_int_list_row_df$pro_name[1]))
        binding_sequences_list <- c(binding_sequences_list,peptide_level)
        binding_seq_start_list <- c(binding_seq_start_list,as.character(pro_int_list_row_df$pep_start[1]))
        binding_seq_end_list <- c(binding_seq_end_list,as.character(pro_int_list_row_df$pep_end[1]))
        binding_sequences_list_real <- c(binding_sequences_list_real,as.character(pro_int_list_row_df$bs_seq[1]))
        binding_seq_start_list_real <- c(binding_seq_start_list_real,as.character(pro_int_list_row_df$bs_start[1]))
        binding_seq_end_list_real <- c(binding_seq_end_list_real,as.character(pro_int_list_row_df$bs_end[1]))

      }

    } else {
      protein_names_list <- c(protein_names_list,pro_int_list_row$pro_name[1])
      binding_sequences_list <- c(binding_sequences_list,pro_int_list_row$pep_seq[1])
      binding_seq_start_list <- c(binding_seq_start_list,pro_int_list_row$pep_start[1])
      binding_seq_end_list <- c(binding_seq_end_list,pro_int_list_row$pep_end[1])
      binding_sequences_list_real <- c(binding_sequences_list_real,paste(unique(pro_int_list_row$bs_seq),collapse='_'))
      binding_seq_start_list_real <- c(binding_seq_start_list_real,paste(unique(pro_int_list_row$bs_start),collapse='_'))
      binding_seq_end_list_real <- c(binding_seq_end_list_real,paste(unique(pro_int_list_row$bs_end),collapse='_'))

    }
    #
    # protein_names_list <- c(protein_names_list,pro_int_list_row$pro_name[1])
    # binding_sequences_list <- c(binding_sequences_list,pro_int_list_row$pep_seq[1])
    # binding_seq_start_list <- c(binding_seq_start_list,pro_int_list_row$pep_start[1])
    # binding_seq_end_list <- c(binding_seq_end_list,pro_int_list_row$pep_end[1])
    # binding_sequences_list_real <- c(binding_sequences_list_real,paste(unique(pro_int_list_row$bs_seq),collapse='_'))
    # binding_seq_start_list_real <- c(binding_seq_start_list_real,paste(unique(pro_int_list_row$bs_start),collapse='_'))
    # binding_seq_end_list_real <- c(binding_seq_end_list_real,paste(unique(pro_int_list_row$bs_end),collapse='_'))
    #

    #initially set both values to 0

    uv_count <- as.character(0)
    noxl_count <- as.character(0)

    for(uv_or_noxl_value in names((table(pro_int_list_row$uv_or_noxl)))){
      if(uv_or_noxl_value == '+'){

        uv_count <- as.character((table(pro_int_list_row$uv_or_noxl)[[uv_or_noxl_value]]))

      } else if(uv_or_noxl_value == '-')

        noxl_count <- as.character((table(pro_int_list_row$uv_or_noxl)[[uv_or_noxl_value]]))
    }

    uv_count_list <- c(uv_count_list,uv_count)
    noxl_count_list <- c(noxl_count_list,noxl_count)

    #need to add the UV count to this

  }

  #length(binding_site_for_supp_df$Binding_End)

  binding_site_for_supp_df <- list(Protein_Name=protein_names_list,
                                   MS_Start=binding_seq_start_list,
                                   MS_End=binding_seq_end_list,
                                   MS_Sequence=binding_sequences_list,
                                   Binding_Sequence=binding_sequences_list_real,
                                   Binding_Start=binding_seq_start_list_real,
                                   Binding_End=binding_seq_end_list_real)


  binding_site_for_supp_df <- data.frame(binding_site_for_supp_df, stringsAsFactors = FALSE)

  return(binding_site_for_supp_df)

} #end function make_binding_site_rows_for_supp_df


#----Get Domains/Variations Info from Proteins API (Still needs work!)----

#need to split by '-'
#put in a for loop, combine the xinet_domain_outputs


#'get uniprot info from proteins api
#'
#'get uniprot info from proteins api
#'@param uniprot_id uniprot_id
#'@param type_of_info type_of_info
#'@param select_categories select_categories
#'@param protein_name protein_name
#'@param fasta_file fasta_file
#'@param start_and_end_pos start_and_end_pos
#'@param output_xinet_domain_df output_xinet_domain_df
#'@export



get_uniprot_info_from_proteins_api <- function(uniprot_id,
                                               type_of_info = c('variation','feature'),
                                               select_categories = NULL,
                                               protein_name = NULL,
                                               fasta_file = NULL,
                                               start_and_end_pos = NULL,
                                               output_xinet_domain_df = FALSE){



  #still need to add custom colors
  #should separate so that the domain df is different?

  #start_and_end_pos
  #will be a list of two positions
  #can have it also be a single number
  #check by seeing if length == 1 || length == 2


  if(type_of_info == 'variation'){


    #still need to include the variation workflow here
    requestURL <- paste("https://www.ebi.ac.uk/proteins/api/variation/",uniprot_id,sep='')

    r <- GET(requestURL, accept("application/json"))

    stop_for_status(r)

    json <- toJSON(content(r))
    #head(fromJSON(json))

    json_variation <- fromJSON(json)




    #json_variation$features


  } else if(type_of_info == 'feature'){ #end if(type_of_info == 'variation')

    #protein_to_uniprot_id$uniprot_id
    #uniprot_id <- 'Q15910'

    requestURL <- paste("https://www.ebi.ac.uk/proteins/api/features/",uniprot_id,sep='')

    r <- GET(requestURL, accept("application/json"))

    stop_for_status(r)

    json <- toJSON(content(r))

    #fromJSON(toJSON(content(r)))

    #head(fromJSON(json))



    json_features <- jsonlite::fromJSON(json)

    #return(json_features)

    if(is.null(json_features$features)){
      return(NULL)
      #can also try to split by '-' and try again(?)
    }

    #json_features$features$category == c('STRUCTURAL','DOMAINS_AND_SITES')
    #

    #unlisted_features <- (unlist(json_features$features))
    if(is.null(select_categories)){
      #if select categories is NULL, make it so that it is equal to all of the categories
      #unlisted_features <- (unlist(json_features$features))
      #select_categories <- (unique(unlisted_features[names(unlisted_features) == 'category']))
      select_categories <- unique(unlist(json_features$features)[names(unlist(json_features$features)) == 'category'])
      #select_categories <- unlist(unique(json_features$features$category))

    }


    #add to blank dataframe or establish the
    for(category in select_categories){

      #json_features$features[json_features$features$category == category,]


    }

    #subset()


    #json_features2_df <- data.frame(json_features$features)
    #json_features2_df$category == select_categories

    #go through the entire list of select categories and add to a larger list
    #domain_features <- json_features$features[unlist(json_features$features$category) == select_categories,]


    #unlisted_features
    domain_features <- json_features$features[unlist(json_features$features$category) == 'DOMAINS_AND_SITES',]


    #make fasta sequence based on the fasta input
    #if is.null() --> don't do a pairwise alignment, just use Uniprot numbering instead
    #fasta_sequence <- toupper(paste(fasta_file$SUZ12_Q15022,collapse=''))

    #protein_name <- 'EZH2_Q15910-2'

    if(!is.null(fasta_file)){

      fasta_sequence <- toupper(paste(fasta_file[[protein_name]],collapse=''))
      pwa_results_feature <- pairwiseAlignment(json_features$sequence,fasta_sequence)
      pwa_ranges_feature <- get_pwa_ranges(pwa_results = pwa_results_feature)
      pwa_strings_feature <- get_pwa_strings(pwa_results = pwa_results_feature)
      maal_feature <- renumber_and_identify_missing_aa(pwa_strings = pwa_strings_feature, start_difference =   pwa_ranges_feature$start_difference)

      aa_fragments_feature <- get_missing_amino_acids(maal = maal_feature, pwa_strings = pwa_strings_feature)


      feature_index_in_sequence <- 0
      new_pdb_numbering <- maal_feature$new_pdb_numbering
      #new_pdb_numbering corresponds to the Uniprot sequence
      #length of new_pdb_numbering should == json_features$sequence
      #use the number in the domains as the index for the new_pdb_numbering
      #re-number the index for the output

      #renumber in json_features or establish new ?


    } else {
      #no fasta file provided --> use regular numbering


      new_pdb_numbering <- 1:nchar(json_features$sequence)

    }




    #use the new



    if(output_xinet_domain_df == TRUE){

      proteinid_list <- c()
      annotname_list <- c()
      startres_list <- c()
      endres_list <- c()
      color_list <- c()

      #proteinid <- 'SUZ12_Q15022'
      hexcodes <- c('#B0E0E6','#66CDAA','#F9D4F2','#FFCC99','#FE6F5E')

      for(row_num in 1:nrow(domain_features)){

        domain_row <- domain_features[row_num,]
        begin_index <- as.numeric(domain_row$begin)
        end_index <- as.numeric(domain_row$end)
        domain_description <- domain_row$description[[1]]

        new_domain_begin <- new_pdb_numbering[begin_index]
        new_domain_end <- new_pdb_numbering[end_index]

        if(!is.null(protein_name)){
          proteinid_list <- c(proteinid_list,protein_name)
        } else {
          proteinid_list <- c(proteinid_list,json_features$entryName)
        }

        annotname_list <- c(annotname_list,domain_description)
        startres_list <- c(startres_list,new_domain_begin)
        endres_list <- c(endres_list,new_domain_end)

        #can change the color based on the domain
        color_list <- c(color_list,'#B0E0E6')

      }

      #c('ProteinId','AnnotName','StartRes','EndRes','Color')
      xinet_annot_df <- list(ProteinId=proteinid_list,
                             AnnotName=annotname_list,
                             StartRes=startres_list,
                             EndRes=endres_list,
                             Color=color_list)

      xinet_annot_df <- data.frame(xinet_annot_df)

      return(xinet_annot_df)

    } else {
      #if FALSE -->

      return(domain_features)

    }

    #c('ProteinId','AnnotName','StartRes','EndRes','Color')




  } else { #end else if(type_of_info == 'feature'){

    cat('Unknown type_of_info selected')

  }#end else to if(type_of_info == 'variation')



  #if type of info == variation
  #do variation loop
  #if select_type_in_column is null --> do all of the columns?
  #select_type_in_column should be a list of all of the columns that you want to
  #select from the output of the inquiry


} #end get uniprot info function

#----Make xiNET/PyMOL visualization for combined UV XL-MS data(pymol still needs work)----

#'Make UV XLMS Visualization
#'
#'make_uv_xlms_visualization
#'@param uniprot_fasta_file_numeric uniprot_fasta_file_numeric
#'@param output_type output_type
#'@export

make_uv_xlms_visualization <- function(uniprot_fasta_file_numeric,output_type=c('xinet','pymol')){

  list_of_colors <- c('gray','blue','purple','red')
  #should also have option that will do the PyMol heatmap colors for large numbers of uniprot_fasta_file_numeric

  #step 1: if pymol is selected --> need to align the sequence with the PDB file
  #if pymol is not selected --> step is not necessary as just need the uniprot sequence

  xinet_mega_df <- data.frame()

  for(protein_name in names(uniprot_fasta_file_numeric)){

    #will need to potentially change the numbering (use resid/resno and indeces potentially)
    ff_numeric_pdb_select <- uniprot_fasta_file_numeric[[protein_name]]#add start and end here for pymol

    pdb_color_list <- list()
    start_num <- 1
    index_num <- 1

    while(index_num <= length(ff_numeric_pdb_select)){

      while(ff_numeric_pdb_select[start_num] == ff_numeric_pdb_select[index_num]){
        #add +1 to index num
        index_num <- index_num + 1
        #get the right resid from the resid vector
        #need to make it
        if(index_num == length(ff_numeric_pdb_select)){
          break
        }
      }

      pdb_color_list_index <- length(pdb_color_list) + 1

      #get
      #get color

      if(index_num == length(ff_numeric_pdb_select)){
        end_num <- index_num
      } else {
        end_num <- (index_num-1)
      }

      if(output_type == 'pymol'){
        start_resno <- resid_and_resno$resno[start_num]
        end_resno <- resid_and_resno$resno[end_num]
      } else {

        start_resno <- start_num
        end_resno <- end_num
      }


      # color_for_pdb <- list_of_colors[as.numeric(ff_numeric_pdb_select[end_num]) + 1]
      # if(output_type == 'xinet'){
      #   color_for_pdb <- col2hex(color_for_pdb)
      # }

      pdb_color_list[[pdb_color_list_index]] <- c(start_resno,end_resno,(as.numeric(ff_numeric_pdb_select[end_num])))
      start_num <- index_num
      #once the while loop is done --> add the start and end to list
      index_num <- index_num + 1
      #add one to index_num

    }

    proteinid_list <- c()
    annotname_list <- c()
    startres_list <- c()
    endres_list <- c()
    color_list <- c()

    for(pdbcl_index in 1:length(pdb_color_list)){
      pdb_color_row <- pdb_color_list[[pdbcl_index]]

      pdb_color_row[1] #start
      pdb_color_row[2] #end
      pdb_color_row[3] #color

      #color code to translate so both colors are the same
      #need both a PyMOL and xiNET output for UV XL-MS

      #protein_name will be part of function input
      proteinid_list <- c(proteinid_list,protein_name)
      annotname_list <- c(annotname_list,paste('Frequency:',as.character(pdb_color_row[3]))) #list frequency number here paste0('Frequency',number)(?)
      startres_list <- c(startres_list,pdb_color_row[1])
      endres_list <- c(endres_list,pdb_color_row[2])
      color_for_pdb <- list_of_colors[(pdb_color_row[3] + 1)]
      if(output_type == 'xinet'){
        color_for_pdb <- col2hex(color_for_pdb)
      }
      color_list <- c(color_list,color_for_pdb)


    }

    #note: need to add

    xinet_annot_df <- list(ProteinId=proteinid_list,
                           AnnotName=annotname_list,
                           StartRes=startres_list,
                           EndRes=endres_list,
                           Color=color_list)



    #at end of loop, add the lists to the main output file
    #rbind()?

    xinet_mega_df <- rbind(xinet_mega_df,data.frame(xinet_annot_df))

  } #end for for(protein_name in names(uniprot_fasta_file_numeric))

  return(xinet_mega_df)

} #end function make uv_xlms_visualization

#----Load ProXL Data----
#' Load ProXL Data
#'
#' This function takes data as downloaded from ProXL's public demo
#'
#' @param proxl_data File name or data.frame to ProXL data
#' @export

load_proxl_data <- function(proxl_data){

  if(typeof(proxl_data) == 'character'){
    proxl_data <- read.table(proxl_data,sep='\t',header = TRUE)
  }

  xlink_list <- list()

  for(row_num in 1:nrow(proxl_data)){
    proxl_data_row <- proxl_data[row_num,]

    peptide_1 <- as.character(proxl_data_row$PEPTIDE.1)
    peptide_2 <- as.character(proxl_data_row$PEPTIDE.2)
    position_1 <- as.character(proxl_data_row$POSITION)
    position_2 <- as.character(proxl_data_row$POSITION.1)

    sequence_xlink <- paste0(peptide_1,'(',position_1,')-',peptide_2,'(',position_2,')')

    protein_1 <- as.character(proxl_data_row$PROTEIN.1)
    protein_2 <- as.character(proxl_data_row$PROTEIN.2)

    proteins_xlink <- paste0(protein_1,'-',protein_2)

    score <- as.character(proxl_data_row$Best.PSM.Value.E.value.SEARCH.ID..20.) #score

    xlink_list[[row_num]] <- list(sequence_xlink=sequence_xlink,
                                  proteins_xlink=proteins_xlink,
                                  score_xlink=score)

  }


  return(xlink_list)

} #end function load_proxl_data


#---RBD Frequency and Heatmap----

#rbdfreqvec <- rbd.freq_vector(data.frame(bs_output))

#should there be a saved version with all of the sequences already saved with the full sequence?
#the full sequence can be in bs_output and then removed if needed or too long
#update the color scheme

#'Get Frequency Vector for RBDmap
#'
#'This function gets the frequency vector
#'@param bs_output bs_output
#'@param name_by name_by
#'@param heatmap heatmap
#'@param db_selection db_selection
#'@param colors colors
#'@param save_plot save_plot
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



ppi.pymol2 <- function(xlink_mega_df,list_of_start_and_end_pdbs = NULL,show_only_real_structures = NULL, write_file = FALSE,
                       custom.color = FALSE, colors = NULL, color_by = 'freq', write.df = FALSE, experiment_directory = NULL,
                       pdb_numbering = FALSE, pymol_file_list_file_name = 'pymol_output.pml'){

  if(is.null(experiment_directory)){
    experiment_directory <- getwd()
  }

  #if there is no start and end PDBs can just get the info from the PDB file itself?
  #sd_pdb <- read.pdb2('Spc98-yeast_missing_sequence_3_single_dot.pdb')
  #if(endsWith(pdb_name, 'single_dot.pdb'))
  #sd_pdb$atom$resno #this will be the position for the making of the pymol file

  #colnames(xlink_mega_df)

  if('X' %in% colnames(xlink_mega_df)){
    xlink_mega_df$X <- NULL
  }

  og_xlink_colnames <- c("seq","pro","dist","freq","freq_color","files","pdb1","pdb2",
                         "pro_pos1","pro_pos2","pro_name1","pro_name2","pep_seq1","pep_seq2",
                         "pep_pos1","pep_pos2","score")

  if('Protein.Position.1' %in% colnames(xlink_mega_df)){
    colnames(xlink_mega_df) <- og_xlink_colnames
  }

  if('category_color' %in% colnames(xlink_mega_df)){
    colnames(xlink_mega_df)[colnames(xlink_mega_df) == 'category_color'] <- 'freq_color'

  }

  pymol_lines <- c()

  #if custom color is TRUE
  #would vars be frequency?
  #xlink_mega_df <- xl_df2

  if(custom.color == TRUE){

    pymol_cc <- color.pymol(sort(unique(xlink_mega_df[[color_by]])), colors = colors)
    #pymol_cc$vars
    #pymol_cc$color_names

    pymol_lines <- c(pymol_lines,pymol_cc$set_colors)

    levels(xlink_mega_df$freq_color) <- c(levels(xlink_mega_df$freq_color),pymol_cc$color_names)

    #xlink_mega_df$freq

    for(index_num in 1:length(pymol_cc$vars)){

      var_name <- pymol_cc$vars[index_num]

      #will rename the freq_color by whatever the var name is


      xlink_mega_df[xlink_mega_df[[color_by]] == var_name, 'freq_color'] <- pymol_cc$color_names[index_num]

      #within(xlink_mega_df, freq_color[color_by == var_name] <- pymol_cc$color_names[index_num])

      #need to select and change the variable name based on the value


    } #end for(index_num in 1:length(pymol_cc$vars))

  } else {#end if custom color == TRUE
    #custom_color == FALSE

    #make the PyMOL legend still even if they're using the original colors
    #just input those colors into the colors/vars

  }
  #should still produce a legend even if they have not opted for a custom color


  #would have to replace the values in freq_color



  #load the files
  pdbs_in_df <- unique(c(levels(xlink_mega_df$pdb1),levels(xlink_mega_df$pdb2)))
  #loop through the pdbs
  #pdb_name <- pdbs_in_df[1]
  for(pdb_name in pdbs_in_df){
    #load pdb files

    if(!is.null(show_only_real_structures) && (pdb_numbering == FALSE)){
      pdb_in_sors <- FALSE
      for(sors in show_only_real_structures){
        #noob solution to current problem
        if(grepl(sors,pdb_name)){
          pdb_in_sors <- TRUE
        }


      } #end first for(sors in show_only_real_structures)


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

      if(pdb_numbering == TRUE){

        #use the fetch command instead
        pdb_name <- strsplit(pdb_name,'_')[[1]][1]

        py_line <- paste0('fetch ',pdb_name,', async=0')


      } else { #end if(pdb_numbering == TRUE)

        py_line <- paste('load ',experiment_directory,'/',pdb_name,
                         sep='')

      } #end else to if(pdb_numbering == TRUE)



      if(!(py_line %in% pymol_lines)){
        pymol_lines <- c(pymol_lines,py_line)
      }




    } #end else to if(!is.null(show_only_real_structures))

    # py_line <- paste('load ',experiment_directory,'/',pdb_name,
    #                  sep='')
    # pymol_lines <- c(pymol_lines,py_line)


    #if this boolean is turned on,
    #color color_name, protein_name

    #should each of the PDBs be colored based on the protein they come from
    #can make this a user option

  }

  #can just have one option
  #if show_surface, color_gray == TRUE
  #py_line <- 'show surface'
  #pymol_lines <- c(pymol_lines,py_line)
  py_line <- 'color gray'
  pymol_lines <- c(pymol_lines,py_line)


  mega_distance_count <- 0
  for(row_num in 1:nrow(xlink_mega_df)){

    pdb1 <- strsplit(as.character(xlink_mega_df$pdb1[row_num]),'.pdb')[[1]]
    pdb2 <- strsplit(as.character(xlink_mega_df$pdb2[row_num]),'.pdb')[[1]]

    if(pdb_numbering == TRUE){
      chain1 <- strsplit(pdb1,'_')[[1]][2]
      chain2 <- strsplit(pdb2,'_')[[1]][2]
      pdb1 <- strsplit(pdb1,'_')[[1]][1]
      pdb2 <- strsplit(pdb2,'_')[[1]][1]


    }

    #if grepl '___' get last character for the chain
    #otherwise use 'A'

    pro_pos1 <- as.character(xlink_mega_df$pro_pos1[row_num])
    pro_pos2 <- as.character(xlink_mega_df$pro_pos2[row_num])

    #if statement to select the right chain


    #can do the check here in the grepl to do the check if it a single dot PDB structure, linear or curved

    #should make a list of pdbs so that this is not repeated twice?
    #can make output as a list that is then used for the next part of the function

    pdb_names <- c(pdb1,pdb2)
    pro_pos_list <- c(pro_pos1,pro_pos2)
    if(pdb_numbering == TRUE){
      chain_list <- c(chain1,chain2)
    }
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

      if(grepl('___',pdb_num) || pdb_numbering == TRUE){ #if 'real' PDB structure use last character since that is the chain

        if(pdb_numbering == TRUE){
          chain <- chain_list[pdb_num_count]
        } else {
          chain <- substrRight(pdb_num,1)
        }

        selection_name_list[[chain_name]] <- chain

      } else {

        #should it be in this part
        #should have an or statement or excessive?
        #do we need file_type --> unless the file is being created within the function
        #if(file_type_2d == 'single atom' && grepl('single_atom',pdb_num)){
        if(grepl('single_dot',pdb_num)){

          #transform the protein position for the PyMol script so it is the first in the list
          pro_pos <- list_of_start_and_end_pdbs[[paste(pdb_num,'.pdb',sep='')]][1]

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
        freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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

        freq_color <- as.character(xlink_mega_df$freq_color[row_num])
        py_line <- paste('color ',freq_color,', ',dist_name,sep='')
        pymol_lines <- c(pymol_lines,py_line)

      } #end if((grepl(sors,pdb_name) && grepl('___',pdb_name)) || !grepl(sors,pdb_name))

    } else {#end if(!is.null(show_only_real_structures))

      freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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

      freq_color <- as.character(xlink_mega_df$freq_color[row_num])
      py_line <- paste('color ',freq_color,', ',dist_name,sep='')
      pymol_lines <- c(pymol_lines,py_line)

    }
    #color the amino acids
    # freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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
    # freq_color <- as.character(xlink_mega_df$freq_color[row_num])
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

  #write lines to .txt
  #separate by '\n'

  #pymol_file_list_file_name <- 'pymol_ts_output.pml'

  #setwd(new_directory)
  write(paste(pymol_lines,collapse = '\n'),pymol_file_list_file_name)




} #end function ppi.pymol2




#----Write PPI PyMOL----



#----Write RBD PyMol-----



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
    pymol_cc <- color.pymol(bsfv_vars, colors = colors, png.name = 'freqvector_legend_test.svg', gray0 = gray0)



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


#----Get Binding Sequence from Input/Eluate Table----

#change protein_to_uniprot_id to a protein dictionary?
#can use it for both the uniprot and pdb ids

#'Get Binding Sequences from Input/Eluate Table
#'
#'This function gets the binding sequences from the input eluate table
#'
#'@param filtered_input_eluate_table filtered_input_eluate_table
#'@param fasta_file fasta_file
#'@param align_to align_to
#'@param protein_to_uniprot_id protein_to_uniprot_id
#'@param protein_dict protein_dict
#'@param include_ambiguous include_ambiguous
#'@export

rbd.getBSfromIET <- function(filtered_input_eluate_table,
                             fasta_file,
                             align_to = 'PDB',
                             protein_to_uniprot_id = NULL,
                             protein_dict = NULL,
                             include_ambiguous = FALSE){

  #filter the table beforehand as well for QC??
  filtered_input_eluate_table <- filtered_input_eluate_table[filtered_input_eluate_table$eluate > 0,]

  #assign protein dict to protein to uniprot id before
  #do a search and replace in the future to get rid of it later??



  #need
  bs_output_mega <- data.frame()
  protein_names_already_run <- list()



  #num <- 2
  for(num in 1:nrow(filtered_input_eluate_table)){

    protein_name <- filtered_input_eluate_table$protein[num]
    # if(protein_name == 'SUZ12_Q15022'){ #need to make changes later to make sure SUZ12 can be accounted for later
    #   next
    # }

    filtered_input_eluate_table_row <- filtered_input_eluate_table[num,]
    input_sequence <- filtered_input_eluate_table_row$sequence
    protein_name <- filtered_input_eluate_table_row$protein
    fasta_seq <- fasta_file[[protein_name]]

    #check if this has a '_' in it before splitting it?
    #need to get last
    protease_split <- strsplit(filtered_input_eluate_table_row$category,'_')[[1]]
    protease <- protease_split[length(protease_split)]

    bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                   protease = protease,
                                   sequence_for_alignment=toupper(paste0(fasta_seq,collapse='')),
                                   protein_name = protein_name,
                                   database_name = 'FASTA',
                                   database_id = protein_name,
                                   include_ambiguous = include_ambiguous)


    binding_seq <- as.character(bs_output$binding_sequence)

    #make a list of sequences that have already been run through the
    #loop so that it doesn't occur multiple times




    #can have an align_to variable


    #database_name <- align_to
    #if database_name is not changed before rbd.

    if(align_to == 'PDB'){

      #see if the protein has already been run (same as before)
      #add in protein_dict functionality

      if(!is.null(protein_dict)){
        #protein_name
        pdb_id <- as.character(protein_dict[protein_dict$protein_id == protein_name,'pdb_id'])
        #need to split pdb_id into pdb_id and chain before going into pdb_info
        pdb_split <- strsplit(pdb_id,'_')[[1]]
        pdb_id <- pdb_split[1]
        chain <- pdb_split[2]
        pdb_read <- check_download_read_pdb(pdb_id = pdb_id)
        resno_and_resid <- quick_resno_and_resid(pdb_read = pdb_read, chain = chain)

        pdb_info <- list(pdb_id=pdb_id,
                         chain=chain,
                         resno=resno_and_resid$resno,
                         resid=resno_and_resid$resid)

        #if not NULL --> skip the protein_names already run?
        #will need to get pdb_info from the data
        #may just need the id, resid, and resno
        #do pairwise alignment anyway with the PDB file?
        #check_download with strsplit(pdb_id,'_)[[1]][1]?

      } else { #end if(!is.null(protein_dict)){

        #no protein dictionary entered --> display the menus to get pdb_info that way

        if(!(protein_name %in% names(protein_names_already_run))){
          pdb_info <- display_preferred_pdb_structure_menu(protein_name = protein_name, fasta_file)
          protein_names_already_run[[protein_name]] <- pdb_info
        } else { #end if(!(protein_name %in% names(protein_names_already_run)))
          pdb_info <- protein_names_already_run[[protein_name]]

        } #end else to if(!(protein_name %in% names(protein_names_already_run)))


      } #end else to if(!is.null(protein_dict)){

      #if alignment can't be done --> need to change database_name


      #if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
      if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),binding_seq)[[1]]) == 0){

        #try a pairwise alignment
        # if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) > 0){
        #   #the input sequence is found in the
        #
        # }


        #pwa_results <- pairwiseAlignment(paste0(pdb_info$resid,collapse=''),binding_seq)
        #score(pwa_results) --> do by the number of dashes?



        cat('Binding sequence found in FASTA not found in PDB file\n')
        #bs_output_db <- rbd.menuDBSearch(input_sequence,fasta_file,protein_name,pdb_info = pdb_info)

        #bs_output <- rbd.menuDBSearch('NEFISEYCGEIISQDEADRR',fasta_file,protein_name)
        #next
      } else { #end if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
        #the input sequence really has shown up in the pdb_info

        bs_s_e <- str_locate_all(paste0(pdb_info$resid,collapse=''),binding_seq)[[1]]
        bs_output$binding_site_start <- pdb_info$resno[bs_s_e[1]]
        bs_output$binding_site_end <- pdb_info$resno[bs_s_e[2]]
        ms_s_e <- str_locate_all(paste0(pdb_info$resid,collapse=''),as.character(bs_output$ms2_peptide_seq))[[1]]
        if(length(ms_s_e) == 0){
          bs_output$ms2_start <- NA
          bs_output$ms2_end <- NA
        } else {
          bs_output$ms2_start <- pdb_info$resno[ms_s_e[1]]
          bs_output$ms2_end <- pdb_info$resno[ms_s_e[2]]
        }

        bs_output$source_sequence <- paste0(pdb_info$resid,collapse='')

        bs_output$db <- 'PDB'
        bs_output$db_id <- paste0(pdb_info$pdb_id,'_',pdb_info$chain)

        # bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
        #                                protease = protease,
        #                                sequence_for_alignment=paste0(pdb_info$resid,collapse=''),
        #                                protein_name = protein_name,
        #                                database_name = 'PDB',
        #                                database_id = paste0(pdb_info$pdb_id,'_',pdb_info$chain))
      } #end else to if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0)

    } else if(align_to == 'UniProt'){ #end if(align_to == 'PDB')

      #how to get the UniProt ID --> blast.pdb with the swissprot database?
      #will need to update the menu function
      #check if there's a

      if(is.null(protein_to_uniprot_id)){
        #if it is NULL --> need to search the database
        #should compile a list for future use? or it will be part of the bs_output so no need?

        if(!(protein_name %in% names(protein_names_already_run))){
          pdb_info <- go_through_menu_loops(fasta_sequence_vector = fasta_seq, database='swissprot')
          protein_names_already_run[[protein_name]] <- pdb_info
        } else { #end if(!(protein_name %in% names(protein_names_already_run)))
          pdb_info <- protein_names_already_run[[protein_name]]

        } #end else to if(!(protein_name %in% names(protein_names_already_run)))


        isoform_num <- substrRight(pdb_info$pdb_id,2)
        if((startsWith(isoform_num,'.')) & (!is.na(as.numeric(substrRight(pdb_info$pdb_id,1))))){
          if(as.numeric(substrRight(pdb_info$pdb_id,1)) > 1){
            cat('Different isoform than 1 detected but using isoform 1 from Uniprot. Recommend downloading FASTA file for analysis.')
          }
        }

        uniprot_id <- strsplit(pdb_info$pdb_id,'\\.')[[1]][1]

        uniprot_info <- uniprot(pdb_info$pdb_id)
        uniprot_seq <- uniprot_info$sequence

      } else { #end if(is.null(protein_to_uniprot_id))
        #something is here and should match the FASTA files to those in the file

        #will search for the files ending in .fa or .fasta

        uniprot_id <- as.character(protein_to_uniprot_id[protein_to_uniprot_id$protein_id == protein_name,'uniprot_id'])

        if(paste0(uniprot_id,'.fa') %in% list.files()){
          uniprot_seq <- paste0(toupper(seqinr::read.fasta(paste0(uniprot_id,'.fa'))),collapse='')
        } else if(paste0(uniprot_id,'.fasta') %in% list.files()){
          uniprot_seq <- paste0(toupper(seqinr::read.fasta(paste0(uniprot_id,'.fasta'))),collapse='')
        } else {
          #does not exist in current directory --> fetching from UniProt
          cat("\nCan't find FASTA file in current directory --> fetching from UniProt\n")
          uniprot_seq <- uniprot(uniprot_id)$sequence
        } #end else to if(paste0(uniprot_id,'.fa') %in% list.files()){




      } #end else to if(is.null(protein_to_uniprot_id))

      #return(uniprot_seq)

      #grepl(uniprot_info$accession,pdb_info$pdb_id)

      #return(length(str_locate_all(uniprot_seq,binding_seq)[[1]]))

      #if(length(str_locate_all(uniprot_seq,input_sequence)[[1]]) == 0){
      if(length(str_locate_all(uniprot_seq,binding_seq)[[1]]) == 0){
        #bs_output <- rbd.menuDBSearch(input_sequence,fasta_file,protein_name,pdb_info = pdb_info)
        cat('Binding sequence not found in chosen UniProt sequence')

      } else { #end if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
        # bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
        #                                protease = protease,
        #                                sequence_for_alignment=uniprot_seq,
        #                                protein_name = protein_name,
        #                                database_name = 'UniProt',
        #                                database_id = uniprot_id)

        bs_s_e <- str_locate_all(uniprot_seq,binding_seq)[[1]]
        #return(bs_s_e)

        #return(bs_output)
        bs_output$binding_site_start <- bs_s_e[1]
        bs_output$binding_site_end <- bs_s_e[2]
        ms_s_e <- str_locate_all(uniprot_seq,as.character(bs_output$ms2_peptide_seq))[[1]]
        #return(ms_s_e)
        bs_output$ms2_start <- ms_s_e[1]
        bs_output$ms2_end <- ms_s_e[2]
        bs_output$source_sequence <- uniprot_seq

        bs_output$db <- 'UniProt'
        bs_output$db_id <- uniprot_id


        #return(bs_output)

      } #end else to if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){


    } else if(align_to == 'FASTA') { #end else if(align_to == 'UniProt')

      #if chosen --> get the FASTA sequence and set as the sequence for rbd.getBindingSeq
      bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                     protease = protease,
                                     sequence_for_alignment=toupper(paste0(fasta_file[[protein_name]],collapse='')),
                                     protein_name = protein_name,
                                     database_name = 'FASTA',
                                     database_id = protein_name,
                                     include_ambiguous = include_ambiguous)


      # if(include_ambiguous == TRUE){
      #   if(nchar(as.character(bs_output$binding_sequence)) == 0){
      #
      #     #levels(bs_output$binding_sequence) <- c(levels(bs_output$binding_sequence),as.character(bs_output$ms2_peptide_seq))
      #     bs_output$binding_sequence <- as.character(bs_output$ms2_peptide_seq)
      #     bs_output$binding_site_start <- as.character(bs_output$ms2_start)
      #     bs_output$binding_site_end <- as.character(bs_output$ms2_end)
      #
      #   } #end if(nchar(as.character(bs_output$binding_sequence))){
      # } #end if(include_ambigious == TRUE){


    } #end else if(align_to == 'FASTA')

    #can then run rbd.getBindingSeq out here

    bs_output_mega <- rbind(bs_output_mega,bs_output)

    #return(bs_output_mega)
  } #end for(num in 1:nrow(filtered_input_eluate_table)){

  return(bs_output_mega)

} #end function rbd.getBSfromIET


#----RBDmap database search----

#'Menu DB Search
#'
#'Menu DB Search
#'
#'@param input_sequence input_sequence
#'@param fasta_file fasta_file
#'@param protein_name protein_name
#'@param pdb_info pdb_info
#'@export

rbd.menuDBSearch <- function(input_sequence,fasta_file,protein_name,pdb_info = NULL){

  bs_output <- NULL

  while(is.null(bs_output)){

    if(is.null(pdb_info)){

      menu_choice <- menu(choices = c('PDB','Uniprot','Original FASTA sequence'),
                          title = paste0('Choose a database to align your sequence against:\n'))


      if(menu_choice == 1){ #PDB

        pdb_info <- display_preferred_pdb_structure_menu(protein_name = protein_name, fasta_file)
        if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
          next
        } else {
          bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                         protease = protease,
                                         sequence_for_alignment=paste0(pdb_info$resid,collapse=''),
                                         protein_name = protein_name,
                                         database_name = 'PDB',
                                         database_id = paste0(pdb_info$pdb_id,'_',pdb_info$chain))
        }

      } else if(menu_choice == 2){ #UniProt

        pdb_info <- go_through_menu_loops(fasta_file[[protein_name]], database = 'swissprot')

        uniprot_info <- uniprot(pdb_info$pdb_id)

        if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
          next
        } else {
          bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                         protease = protease,
                                         sequence_for_alignment=uniprot_info$sequence,
                                         protein_name = protein_name,
                                         database_name = 'UniProt',
                                         database_id = pdb_info$pdb_id)
        }

      } else if(menu_choice == 3){ #Original FASTA

        bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                       protease = protease,
                                       sequence_for_alignment=toupper(paste0(fasta_file[[protein_name]],collapse='')),
                                       protein_name = protein_name,
                                       database_name = 'FASTA',
                                       database_id = protein_name)

      } #end menu_choices


    } else {
      #not null --> can use the original menu


    menu_choice <- menu(choices = c('Yes, choose a different chain','Yes, choose a different PDB','No, match against Uniprot ID corresponding to current PDB and chain','No'),
                        title = paste0('This sequence was not found on your selection for\n',input_sequence,'\non protein ',protein_name,'\nWould you like to choose a different chain or PDB?'))

    if(menu_choice == 1) { #yes, chosoe a different chain

      #if yes --> choose a different a chain


      pdb_info <- match_sequence_to_pdb_and_chain(protein_name = protein_name,
                                                  fasta_seq = fasta_file[[protein_name]],
                                                  pdb_read = check_download_read_pdb(pdb_info$pdb_id),
                                                  pdb_rl = pdb_info$pdb_id)

      #automatically go through all of the chains again?
      #readline(prompt = "Enter a new chain: ")



    } else if(menu_choice == 2){ #yes, choose a different PDB

      #go to the display menu again to choose a PDB

      #can change this so that you can just put in the sequence and it will just count the sequences
      #and have that be the display name? (or it will just output the sequence and truncate it to like 10 or so chars)
      pdb_info <- display_preferred_pdb_structure_menu(protein_name = protein_name, fasta_file)
      if(length(str_locate_all(paste0(pdb_info$resid,collapse=''),input_sequence)[[1]]) == 0){
        next
      } else {
        bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                       protease = protease,
                                       sequence_for_alignment=paste0(pdb_info$resid,collapse=''),
                                       protein_name = protein_name,
                                       database_name = 'PDB',
                                       database_id = paste0(pdb_info$pdb_id,'_',pdb_info$chain))
      }


      #should have some kind of boolean in here so that they can't escape
      #the while loop (before the menu choice)

    } else if(menu_choice == 3){ #match against current UniProt ID

      #pdb_anno1 <- pdb.annotate(pdb_info$pdb_id)

      pdb_anno <- pdb.annotate(pdb_info$pdb_id)[paste0(pdb_info$pdb_id,'_',pdb_info$chain),]

      if(length(str_locate_all(pdb_anno$sequence,input_sequence)[[1]]) > 0){
        #if length > 0, the input sequence exists within the sequence obtained from pdb.annotate()
        #if true, will need to use that for get_binding_sequence
        #use db_name and db_id for those fields in the bs_output

        #add db, db_id to input variable

        #just get the bs_output here?

        bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                       protease = protease,
                                       sequence_for_alignment=pdb_anno$sequence,
                                       protein_name = protein_name,
                                       database_name = pdb_anno$db_name,
                                       database_id = pdb_anno$db_id)


      } else {
        #if length == 0, the input sequence does not exist within the sequence obtained from pdb.annotate()
        #what to do here? should be sent to the top of the loop?

        cat('Input sequence not found in the PDB annotation --> using FASTA sequence')
        bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                       protease = protease,
                                       sequence_for_alignment=toupper(paste0(fasta_file[[protein_name]],collapse='')),
                                       protein_name = protein_name,
                                       database_name = 'FASTA',
                                       database_id = protein_name)

      }


    } else if(menu_choice == 4){ #no

      mc2 <- menu(c('Uniprot -- not working atm','Input FASTA'),title = 'Choose a source to compare your sequence to: \n')

      if(mc2 == 1){ #Uniprot

        #two options here?
        if(is.null(protein_to_uniprot_id)){
          #ask the user what they want
          #for now, should just blast it and they'll just have to deal with it
          #BLAST the sequence
          #pick the first one for right now? (unless functionality for menu function is working)
          #get the Uniprot ID


        } else {
          #get the Uniprot ID from the dictionary, see if there's a
        }

        #put Uniprot code here
        #can have 2 options for Uniprot code: either through the uniprot() function from bio3d
        #or can have their own Uniprot FASTA file

        #for isoforms, recommend putting in your own Uniprot FASTA file

        #Uniprot code: can check if there is a file that contains the

        #do a display_preferred_pdb menu?
        #blast.pdb(toupper(paste0(fasta_file[[protein_name]],collapse='')),database = 'swissprot')

        #should turn uniprot thing into a function to be used here

      } else if(mc2 == 2){ #Input FASTA

        #just use the sequence that came with the fasta file
        bs_output <- rbd.getBindingSeq(ms_sequence = input_sequence,
                                       protease = protease,
                                       sequence_for_alignment=toupper(paste0(fasta_file[[protein_name]],collapse='')),
                                       protein_name = protein_name,
                                       database_name = 'FASTA',
                                       database_id = protein_name)

        #break out of the while loop

      }

    } #end if(menu_choice)

  } #end else to if(is.null(pdb_info))


  } #end while(is.null(bs_output))



  return(bs_output)

} #end function rbd.menuDBSearch

#'PyMol Colors
#'
#'Data on PyMol colors
#'
#'@docType data
#'
#'@usage data(pymol_color_table)
#'
"pymol_color_table"




#add in pymol wiki info, etc. for the sake of documentation

#first save the data to
#need to also include rbdmap data and protein-protein interaction data
#will the code be able to handle pre-loaded data? --> should do this for the getSeqHitList function as well
#may need a different code


#save(pymol_color_table,file='data/pymol_color_table.RData')

#setwd('/Users/emmagail/Documents/monash/project.functions')
#setwd('/Users/emmagail/Documents/crisscrosslinker')

#roxygen2::roxygenise()
#requireNamespace("crisscrosslinker")

#build manual
#system("R CMD Rd2pdf ~/Documents/monash/project.functions")
