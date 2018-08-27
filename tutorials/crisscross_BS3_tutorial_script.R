setwd()

# library(devtools)
# library(roxygen2)

#should just load crisscross package here

#----Load Libraries/Files-----

#Load all FASTA files
fasta_prc2_accessories <- seqinr::read.fasta('~/Downloads/PRC2_accessories2.txt')
fasta_file <- seqinr::read.fasta('~/Downloads/PRC2_5m.fasta.txt')
nucleosome_fasta <- seqinr::read.fasta("~/Downloads/Nucleosome_seq.fasta.txt")

#can combine them into one FASTA file if needed
combined_fasta_file <- c(fasta_file,fasta_prc2_accessories,nucleosome_fasta)

protein_to_uniprot_id <- read.csv('~/Downloads/XL-MS Xlink Visual/protein_to_uniprot_id2.csv')
head(protein_to_uniprot_id)

#Establish the main experiment directory
main_experiment_directory <- '~/Downloads/XL-MS Xlink Visual/'
#Set the working directory
setwd(main_experiment_directory)

#If experiments are categorized into different folders
experiment_folder_names <- c("BS3 - 2 - PRC2-AEBP2",
                             "BS3 - 5 - PRC2-PHF19",
                             "BS3 - 7 - PRC2-MTF2-EPOP")

core_subunits <- c("EZH2_Q15910-2","SUZ12_Q15022","EED_O75530","RBBP4_Q09028")

freq_dist_vectors <- c()
experiment_name_list <- c()
random_lysine_vectors <- c()

chains <- c('A','M','Q','K','C','L','N')

freq_dist_vectors_list <- list()
random_lysine_vectors_list <- list()


for(experiment_folder in experiment_folder_names){
  
  specific_experiment_directory <- paste(main_experiment_directory,experiment_folder,sep='')
  setwd(specific_experiment_directory)
  
  list_of_files_in_directory <- list.files(specific_experiment_directory)
  
  enclosed_experiment_folder <- gsub('\\+','\\\\+',experiment_folder)
  enclosed_experiment_folder <- gsub('\\(','\\\\(',enclosed_experiment_folder)
  enclosed_experiment_folder <- gsub('\\)','\\\\)',enclosed_experiment_folder)
  
  #grepl(paste(plink_file_name,'_xinet_input_file.csv',sep=''),list_of_files_in_directory)
  
  list_of_files_in_directory <- list_of_files_in_directory[grepl(enclosed_experiment_folder,list_of_files_in_directory)]
  directory_boolean <- !(grepl(paste(enclosed_experiment_folder,'_xinet_input_file.csv',sep=''),list_of_files_in_directory))
  list_of_files_in_directory <- list_of_files_in_directory[directory_boolean]
  
  list_of_files_in_directory <- list_of_files_in_directory[!startsWith(list_of_files_in_directory,'~$')]
  
  #should make option to opt out of making a new directory during each run
  xlink_file <- make_diff_analysis_plink_csv_and_pymol(list_of_files_in_directory,fasta_file = combined_fasta_file,
                                                       file_type_2d = 'single_dot',
                                                       fasta_names_to_generate_all_2d_structures = NULL,
                                                       csv_pdb_input_file = '~/Downloads/XL-MS Xlink Visual/master_start_and_end_pdb.csv',
                                                       show_only_real_structures = names(fasta_file),
                                                       distance_histogram_name = paste(experiment_folder,'_distance_histogram.png',sep=''),
                                                       pymol_file_list_file_name = paste(experiment_folder,'_pymol_ts_output.pml',sep=''),
                                                       xlink_df_name = paste(experiment_folder,'_xlink_dataframe.csv',sep=''),
                                                       xlink_view_file_name = paste(experiment_folder,'_xinet_input.csv',sep=''),
                                                       protein_alternative_names_dict = '~/Downloads/protein_fasta_names_and_alternatives.csv',
                                                       pdb_directory = '~/Downloads/XL-MS Xlink Visual/')
  
  #can also have an option where if they name the different "names" NULL that they will not be outputted
  #should include xlink_viewer_csv as part of the make_diff_analysis function and then just add
  #a boolean that will turn it on or off --> boolean could be NULL like the other options
  
  
  #Remove the distances between real structures that are not in the same complex
  #Should be integrated into make_diff_analysis
  
  filtered_core_subunit_df <- filter_xlink_df_by_protein_names(xlink_file,core_subunits)
  
  
  #need to fix the xinet_file output so that it's only the 4 columns since it doesn't work with all of them
  #should have option to include/exclude score in the output since it can't be used in the current state
  xinet_file <- make_xlink_viewer_csv(xlink_file)
  #will have to also alter the renumber_xinet function if the make_xlink_viewer function output is altered
  #change the column names that need to be renumbered
  
  #will need preview of protwein_to_uniprot_id
  
  renumbered_xinet_file <- renumber_xinet_from_uniprot_fasta(xinet_file,protein_to_uniprot_id,combined_fasta_file)
  
  #add in function to make domain file as well
  
  
  #can have checks here to show the before and after of the protein names once they have been filtered
  
  filtered_xl_freq_df <- count_frequencies_in_xlms_output(xl_dataframe = filtered_core_subunit_df)
  
  xl_freq_df <- filtered_xl_freq_df
  xl_freq_df_filtered <- xl_freq_df[((xl_freq_df$Lowest_Score <= 1e-05) | (xl_freq_df$Position_Frequency >= 2)),]
  xl_freq_df_distances <- get_vector_of_distances_by_pos_freq(xl_freq_df_filtered)
  #add to list() using [[experiment_folder]]
  #freq_dist_vectors_for_hist[[experiment_folder]] <- xl_freq_df_distances
  
  experiment_folder_split <- strsplit(experiment_folder,'- ')[[1]][3]
  
  
  freq_dist_vectors <- c(freq_dist_vectors,xl_freq_df_distances)
  ef_repeats <- rep(experiment_folder_split,length(xl_freq_df_distances))
  experiment_name_list <- c(experiment_name_list,ef_repeats)
  
  freq_dist_vectors_list[[experiment_folder]] <- xl_freq_df_distances
  
  
  #repeat 
  
  random_lysines_in_list <- generate_random_lysine_distances_in_pdb('6C23',frequency_vector=xl_freq_df_distances,chains=chains)
  random_lysine_vectors <- c(random_lysine_vectors,random_lysines_in_list)
  
  random_lysine_vectors_list[[experiment_folder]] <- random_lysines_in_list
  
  
  #choose filtration criteria here --> should be within the initial code but can do it outside for the
  #purposes of this tutorial
  
  #filter the categories here for the filtered_xl_freq_df(?)
  
  
  
  #filtered_xl_freq_df <- count_frequencies_in_xlms_output(xlink_file = filtered_core_subunit_df)
  
  
  
  
  #switch to experiment directory
  #diff_analysis
  
  #uniprot correction (can be combined in diff_analysis)
  #freq_count
  #xinet_creation
  

  #combine data into a dataframe for the histogram?
  
}

freq_dist_vectors_for_hist <- list(xx=freq_dist_vectors,
                                   yy=experiment_name_list)

random_lysine_vectors_for_hist <- list(xx=random_lysine_vectors,
                                       yy=experiment_name_list)


freq_dist_vectors_for_hist <- data.frame(freq_dist_vectors_for_hist)
random_lysine_vectors_for_hist <- data.frame(random_lysine_vectors_for_hist)

#Style for Supplementary Files

ggplot(freq_dist_vectors_for_hist,aes(x=xx,fill=yy)) +
  geom_histogram(position='dodge') +
  labs(shape = "Category") +
  labs(x = 'Distance') +
  labs(y = 'Frequency') +
  ylim(0,30) +
  xlim(0,135) +
  labs(title = 'Distance Histogram') +
  theme(
    panel.background = element_rect(fill = "#ffffff",
                                    colour = "#000000",
                                    size = 0.5, linetype = "solid"),
    legend.key = element_rect(fill = "#eaf2ff")
  ) +
  
  scale_fill_manual(name = "",
                    labels = levels(factor(freq_dist_vectors_for_hist$yy)),
                    values = c('red','blue','green'))

ggsave('distance_histogram_dodge_freq_2_p_e05_core.png')




ggplot(random_lysine_vectors_for_hist,aes(x=xx,fill=yy)) +
  geom_histogram(position='dodge') +
  labs(shape = "Category") +
  labs(x = 'Distance') +
  labs(y = 'Frequency') +
  labs(title = 'Distance Histogram') +
  ylim(0,30) +
  xlim(0,135) +
  theme(
    panel.background = element_rect(fill = "#ffffff",
                                    colour = "#000000",
                                    size = 0.5, linetype = "solid"),
    legend.key = element_rect(fill = "#eaf2ff")
  ) +
  
  scale_fill_manual(name = "",
                    labels = levels(factor(freq_dist_vectors_for_hist$yy)),
                    values = c('red','blue','green'))


ggsave('distance_histogram_control_dodge_freq_2_p_e05_core.png')





#needs to be a clean file with defined sections

#load files
#define variables at the top (probably)
#add in any selection criteria

#should have the experimental_files loop?


#have initial setup from plink1_and_2 analysis script

#should potentially have option to change numbering to match PDB structure --> would need to add boolean to function
#add in elements from both freq and xinet_renumbering scripts
#essential --> renumbering based on Uniprot

#histogram creation
#control lysines as well as original histogram


#will also need to update the UV XL-MS workflow as well 
#can have final version contain both
#the final version should be 1 cohesive document 

#need to add the functions that made the heatmap

#should include examples of the files --> load the data to the tutorial and use head() to show example
#use fake data for the tutorial for the plink data(?)
#add the heatmap of the different proteins used in the experiment used for the UV XL-MS data
#should have both versions of the heatmap in the tutorial

#how should the package be installed?
#can just upload to github and can be installed using devtools
#will need to stress that this a work in progress and the public should direct inquiries to my email address
#for specific conderns or tutorials that they would like to see
#documentation is still in its infancy and will need to be updated 



