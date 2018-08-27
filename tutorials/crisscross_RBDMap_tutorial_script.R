###########################
#----IMPORTANT REMINDER----

#Don't forget to import the functions from 'pdb_functions.R' if not
#included as a library in "Import Libraries"

#Always run this standalone script as "Source" not "Run" otherwise
#the menus will not work properly

###########################

#----Import Libraries----
library(ggplot2)
library(bio3d)
library(Biostrings)
library(seqinr)
library(RColorBrewer)

#----Set Variables----
experiment_directory <- '/Users/emmagail/Downloads/20171122_PRC2_XLMS'
control_file_prefix <- 'UVCD'
experimental_file_prefix <- 'UVMA'
secondary_prefix <- TRUE
preferred_pdb <- NULL #'5HYN'
preferred_chain <- NULL
#list_of_secondary_prefixes <- c('LysC')
list_of_secondary_prefixes <- c('LysC','ArgC')
cutoff_score <- 20
name_of_experiment <- "UVCD and UVMA Peptide Intensity: Eluate vs. Input"
fasta_file_in_different_directory <- TRUE
fasta_file_directory <- '/Users/emmagail/Downloads/prc1and2fastafiles'
fasta_file_names <- c('PRC2_5m.fasta')

#----Set Working Directory----
setwd(experiment_directory)


#----Combine FASTA File----
fasta_file <- combine_fasta_files(fasta_file_names = fasta_file_names,
                                  experiment_directory = experiment_directory,
                                  fasta_file_directory = fasta_file_directory,
                                  fasta_file_in_different_directory = TRUE)



#----PART ONE: Make Graphs Comparing Input and Eluate----

#----Make Sequence Hit List----

sequence_hit_list <- make_sequence_hit_list(experiment_directory = experiment_directory,
                                            fasta_file = fasta_file)

#---Make Input/Eluate Output Table----

#still need to double check amino acid before and amino acid after columns
#edit make_intensity_output_df to include secondary prefix
#or at least make sure new_shl only contains the secondary prefix
input_eluate_table <- make_input_eluate_table(control_file_prefix,
                                              experimental_file_prefix,
                                              sequence_hit_list,
                                              secondary_prefix,
                                              list_of_secondary_prefixes)


list_of_prefixes <- c(control_file_prefix,experimental_file_prefix)


#----Make Input/Elute Graph----

input_eluate_graph <- make_input_eluate_graph(input_eluate_table,
                                              list_of_prefixes,
                                              list_of_secondary_prefixes,
                                              secondary_prefix)

#need to edit this to make sure the plot saves to the directory
#save_plot_as_png(input_eluate_graph,name_of_experiment)
#function for this part of the experiment?

png(paste(name_of_experiment,'.png',sep=''))
input_eluate_graph #print graph to "Plots" tab
dev.off()


#----PART TWO: GET THE RNA-BINDING SEQUENCES----

#get only the rows where the eluate variable is greater than 0
filtered_input_eluate_table <- input_eluate_table[input_eluate_table$eluate > 0,]

pymol_mega_list <- list()
protein_names_already_run <- list()
pymol_mega_list <- list()

#input for function:
#filtered_input_eluate_table
#fasta_file



num <- 1
for(num in 1:nrow(filtered_input_eluate_table)){
  protein_name <- filtered_input_eluate_table$protein[num]
  if(protein_name == 'SUZ12_Q15022'){ #need to make changes later to make sure SUZ12 can be accounted for later
    next
  }
  
  #make a list of sequences that have already been run through the 
  #loop so that it doesn't occur multiple times
  if(!(protein_name %in% protein_names_already_run)){
    
    fasta_seq <- fasta_file[[protein_name]]
    #can ask user whether or not they want to go through the loops
    run_through_loops_title <- paste('Do you have a preferred PDB structure you would like to match to your sequence ',
                                     protein_name,'?',sep='')
    run_through_loops_options <- c('Yes','No, run BLAST on this sequence','Advanced options')
    rtl_menu_selection <- menu(run_through_loops_options, title=run_through_loops_title)
    if(rtl_menu_selection == 1){ #yes
      #error check: make sure that the PDB ID is valid
      #may require a tryCatch()
      pdb_rl <- readline(prompt = "Enter the PDB ID: ")
      
      ff_pdb_read <- check_download_read_pdb(pdb_rl)
      
      
      #if checking multiple sequences, enclose in a loop
      #what if there are no matches in the chain loop?
      seqres_matching_list <- do_sr_chain_loop(protein_name = protein_name,
                                               fasta_seq = fasta_seq,
                                               pdb_read = ff_pdb_read)
      
      #have output be null if it does not match
      #if is.null() --> go to menu selection 2
      #turn into a function?
      
      #ask user to pick one for their analysis with a function
      pick_chain_title <- paste("Pick a chain to use for ",protein_name," for PDB ID: ",pdb_rl,sep='')
      pick_chain_menu_selection <- menu(seqres_matching_list$chains[[protein_name]],title=pick_chain_title)
      
      chain_picked <- seqres_matching_list$chains[[protein_name]][pick_chain_menu_selection]
      
      pdb_info <- seqres_matching_list$pdb_info[[paste(pdb_rl,'_',chain_picked,'_',protein_name,sep='')]]
      
      resno_and_resid <- quick_resno_and_resid(ff_pdb_read,chain_picked)
      pdb_info <- c(pdb_info,resno_and_resid)
      
      
    } else if(rtl_menu_selection == 2){ #no
      
      pdb_info <- go_through_menu_loops(fasta_sequence_vector = fasta_seq)
      
      ff_pdb_read <- check_download_read_pdb(pdb_info$pdb_id)
      
      resno_and_resid <- quick_resno_and_resid(ff_pdb_read,pdb_info$chain)
      pdb_info <- c(pdb_info,resno_and_resid)
      
      
    } else if(rtl_menu_selection == 3){ #advanced options
      
      #advanced options menu
      #as part of a while loop?
      
      advanced_options <- c('Enter a PDB ID to use for all sequences in this fasta file',
                            'Run BLAST for all sequences in this fasta file',
                            'Use previously found PDB file for this sequence',
                            'Go back')
      advanced_options_title <- 'Advanced Options'
      advanced_options_menu_selection <- menu(advanced_options,title = advanced_options_title)
      
      if(advanced_options_menu_selection == 1){ #enter a PDB ID for all sequences
        
        
      } else if(advanced_options_menu_selection == 2){ #run BLAST for all sequences
        
        
      } else if(advanced_options_menu_selection == 3){ # use previously found PDB file
        
        
      } else if(advanced_options_menu_selection == 4){#go back 
        
        
      }
      
    } #end else if rtl_menu_selection == 3
    
    protein_names_already_run[[protein_name]] <- pdb_info
    
    
  } else {#end if !protein_name %in% protein_names_already_run
    
    pdb_info <- protein_names_already_run[[protein_name]]
    
  }
  
  filtered_input_eluate_table_row <- filtered_input_eluate_table[num,]
  
  binding_sequence_row <- make_binding_site_df(pdb_info,
                                               filtered_input_eluate_table_row,
                                               fasta_seq = fasta_seq)
  pymol_mega_list[[num]] <- binding_sequence_row
  
  
} #end of for loop num in 1:nrow(filtered_input_eluate_table)

#----Write Output Files----

write_binding_sequence_csv(pymol_mega_list)
write_pymol_file(pymol_mega_list)



#----The Table------

experiment_directory <- '~/Downloads/filenamesforxlmspapersupplementarydata/Making The Table'
setwd(experiment_directory)

protein_to_uniprot_id <- read.csv('~/Downloads/XL-MS Xlink Visual/protein_to_uniprot_id2.csv')
fasta_file_directory <- '~/Downloads/XL-MS Xlink Visual/xiNET Uniprot Correction'
fasta_file <- seqinr::read.fasta('~/Downloads/PRC2_5m.fasta.txt')

#experiment_number <- 0
#experiment_number_list <- c()
uv_or_noxl_list <- c()
protease_list <- c()
proteins_and_intensity_list <- list()

#need to change the order so they are in chronological order like the ones
#in Nick's email

#experiment_prefixes <- c('20170817','UVCD','UVMA','20171026')
#experiment_prefixes <- c('20170817','UVCD','UVMA','20170828')
experiment_prefixes <- c('20170817','20170828','UVMA','UVCD')

#e_prefix <- 'UVCD'


#should turn this into its own function
#add everything else as functions as well

#make mega function?

proteins_and_intensity_list_with_supp_header2 <- make_proteins_and_intensity_list_from_pipeline_output(experiment_prefixes,
                                                                                                       protein_to_uniprot_id,
                                                                                                       fasta_file,
                                                                                                       fasta_file_directory,
                                                                                                       identifier_is_binding_seq = TRUE)


input_and_eluate_dfs2 <- make_input_and_eluate_dfs_from_pi_list(proteins_and_intensity_list_with_supp_header2)

# proteins_and_intensity_list_with_supp_header <- make_proteins_and_intensity_list_from_pipeline_output(experiment_prefixes,
#                                                                   protein_to_uniprot_id,
#                                                                   fasta_file,
#                                                                   fasta_file_directory)


#proteins_and_intensity_list_with_supp_header$bs_output_df
#write.csv(proteins_and_intensity_list_with_supp_header$bs_output_df,'bs_output_df.csv')

bs_output_df <- proteins_and_intensity_list_with_supp_header2$bs_output_df

bs_output_df_e_filtered <- bs_output_df[bs_output_df$eluate > 0,]

write.csv(bs_output_df_e_filtered,'bs_output_filtered.csv')

proteins_and_intensity_list <- proteins_and_intensity_list_with_supp_header2$proteins_and_intensity_list
supp_header_df <- proteins_and_intensity_list_with_supp_header2$supp_header_df

input_and_eluate_dfs <- make_input_and_eluate_dfs_from_pi_list(proteins_and_intensity_list_with_supp_header2)

input_dataframe <- input_and_eluate_dfs$input_df
eluate_dataframe <- input_and_eluate_dfs$eluate_df



#this will need to be altered to comply with the new input/eluate output
#maybe make another boolean as an input (similar to the boolean when making the proteins and intensity list)
#binding_site_for_supp_df <- make_binding_site_rows_for_supp_df(proteins_and_intensity_list = proteins_and_intensity_list_with_supp_header$proteins_and_intensity_list)

binding_site_for_supp_df2 <- make_binding_site_rows_for_supp_df(proteins_and_intensity_list = proteins_and_intensity_list_with_supp_header2$proteins_and_intensity_list,
                                                                identifier_is_binding_seq = TRUE)

write.csv(binding_site_for_supp_df2,'binding_site_for_supp.csv')

#write.csv(binding_site_for_supp_df,'binding_site_for_supp_final.csv')
write.csv(t(supp_header_df),'supp_header_df_final.csv')
write.csv(input_dataframe,'input_supp_df_final.csv')
write.csv(eluate_dataframe,'eluate_supp_final.csv')


#----The Heatmap----


list_of_aa_vectors <- list()

aa_mega_df <- data.frame()

for(protein_level in levels(bs_output_df_e_filtered$pro_name)){
  
  bs_o_df_rows <- bs_output_df_e_filtered[bs_output_df_e_filtered$pro_name == protein_level,]
  rownames(bs_o_df_rows) <- 1:nrow(bs_o_df_rows)
  
  uniprot_id <- protein_to_uniprot_id[protein_to_uniprot_id$protein_id == protein_level,'uniprot_id']
  uniprot_fasta <- seqinr::read.fasta(paste(fasta_file_directory,'/',uniprot_id,'.fasta',sep=''))
  
  uniprot_fasta_seq <- toupper(uniprot_fasta[[names(uniprot_fasta)]])
  #make vector of 0s?
  #use the indices as a count
  aa_count_vector <- rep(0,length(uniprot_fasta_seq))
  
  #aa_count_vector[1:6] <- aa_count_vector[1:6] + 1
  
  for(row_num in 1:nrow(bs_o_df_rows)){
    
    bs_o_df_row <- bs_o_df_rows[row_num,]
    
    bs_start <- as.numeric(bs_o_df_row$binding_site_start)
    bs_end <- as.numeric(bs_o_df_row$binding_site_end)
    bs_seq <- as.character(bs_o_df_row$binding_sequence)
    
    
    aa_count_vector[bs_start:bs_end] <- aa_count_vector[bs_start:bs_end] + 1
    
  }
  
  #once the vector is done --> need to add to list
  list_of_aa_vectors[[protein_level]] <- aa_count_vector
  
  aa_df <- data.frame(xx=1:length(list_of_aa_vectors[[protein_level]]),
                      Count=list_of_aa_vectors[[protein_level]],
                      zz=rep(protein_level,length(list_of_aa_vectors[[protein_level]])))
  
  aa_mega_df <- rbind(aa_mega_df,data.frame(aa_df))
  
  # ggplot(aa_df, aes(xx,zz)) + geom_tile(aes(fill=Count)) + xlab('') + ylab('') + coord_fixed(ratio = 20) +
  #   theme(
  #     panel.background = element_rect(fill = "#ffffff",
  #                                     colour = "#000000",
  #                                     size = 0.5, linetype = "solid"),
  #     legend.key = element_rect(fill = "#eaf2ff")
  #   )
  # 
  # 
  # #heatmap(matrix(1:length(list_of_aa_vectors[[protein_level]]),list_of_aa_vectors[[protein_level]]))
  # aa_heatmap_file_name <- paste(protein_level,'_aa_count_heatmap.png',sep='')
  # ggsave(aa_heatmap_file_name)
  
  
}

list_of_aa_vectors

#make this heatmap for all of the input and eluate as well
#4 heatmaps --> + and - UV  for both input and eluate

color_palette1 <- colorRampPalette(c('#d0d0d0','#2f5ac6','#e50000'),3)
#color_palette1 <- colorRampPalette(c('#d0d0d0','#2f5ac6','#e50000'))

#hm.palette <- colorRampPalette(rev(brewer.pal(4, 'Spectral')), space='Lab')
hm.palette <- colorRampPalette(c('#d0d0d0','#2f5ac6','#e50000'), space='Lab')

unique(aa_mega_df$Count)

ggplot(aa_mega_df, aes(xx,zz)) + geom_tile(aes(fill=Count)) + xlab('AA Position in Protein Sequence') + ylab('Protein Name') + coord_fixed(ratio = 50) +
  theme(
    panel.background = element_rect(fill = "#ffffff",
                                    colour = "#000000",
                                    size = 0.5, linetype = "solid"),
    legend.key = element_rect(fill = "#eaf2ff")
  ) + scale_fill_gradientn(colours = hm.palette(100))

ggsave('uvxlms_heatmap_all_proteins_gray_blue_red.png')
