#what does this function do:
#produce the plot of the distance
#if plot is turned off, will produce the data to make your own plot


#input xlink.df.pdb

#if plot == TRUE, will produce the plots (have grid option as well?)
#if plot == FALSE, it will return the dataframe that was used to create the plots
#should include both options within the tutorial to be able to showcase how to make
#a custom plot

#add in wes anderson plots to tutorial
#wes_palettes



#ppi.distAnalysis(xlink.df.pdb = xlink.df.pdb)

#'PPI Distance Analysis
#'
#'This function analyzes the distribution of calculated distances of your PDB matches by comparing it to a random distribution of lysines.
#'
#'@param xlink.df.pdb Dataframe created by ppi.matchPDB2()
#'@param plot If TRUE, will produce the plots. If FALSE, will return a set of dataframes that were used to generate the plots. Defaults to TRUE.
#'
#'@author Emma Gail
#'
#'@export

ppi.distAnalysis <- function(xlink.df.pdb,plot=TRUE){
  #xlink.df.pdb$pdb1
  #xlink.df.pdb$pdb2

  #xlink.df.pdb$dist


  unlisted_pdbs <- unlist(strsplit(unique(c(xlink.df.pdb$pdb1,xlink.df.pdb$pdb2)),'_'))

  #need to get all of the ones in here that are not from the same pdb_id

  pdbs <- unlisted_pdbs[unlisted_pdbs != '-'][c(T,F)]
  chains <- unlisted_pdbs[unlisted_pdbs != '-'][c(F,T)]

  filtered_xl_freq_df <- ppi.freqCount(xlink_df = xlink.df.pdb)
  filtered_xl_freq_df <- na.omit(filtered_xl_freq_df)

  #for loop here of the pdbs

  #if there are multiple PDBS (make multiple plots?)
  #or should have option to have all of it on one plot?
  #will still have to add them together in the end because they will have different chains

  #nrow(na.omit(xlink.df.pdb))


  dist.hist.df.mega <- list()

  for(pdb_id in unique(pdbs)){

    filtered_xl_freq_df_pdb <- filtered_xl_freq_df[filtered_xl_freq_df$PDB == pdb_id,]
    xl_freq_df_distances <- get_vector_of_distances_by_pos_freq(filtered_xl_freq_df_pdb)


    #get the chains for this pdb_id
    chains_pdb <- chains[pdbs == pdb_id]


    random_lysines_in_list <- generate_random_lysine_distances_in_pdb(pdb_id = pdb_id,
                                                                      frequency_vector=xl_freq_df_distances,
                                                                      chains=chains)

    #can add this to
    dist.hist.df <- data.frame(Experimental=xl_freq_df_distances,
                               Control=random_lysines_in_list)

    dist.hist.df.mega[[pdb_id]] <- dist.hist.df

    if(plot == TRUE){

      p <- ggplot(melt(dist.hist.df), aes(value, fill = variable)) +
        geom_histogram(position = "dodge") +
        labs(x=expression(paste('Distance (',ring(A),')')),y='Frequency',title=pdb_id) +
        theme(legend.title = element_blank())

    } #end if(plot == TRUE)


    #add both of the vectors to the

  } #end for(pdb_id in unique(pdbs))

  #need to include the PDB ID in here as well in order to filter by
  #filtered_xl_freq_df <- ppi.freqCount(xlink_df = xlink.df.pdb)
  #xl_freq_df_distances <- get_vector_of_distances_by_pos_freq(filtered_xl_freq_df)



  if(plot == FALSE){
    return(dist.hist.df.mega)
  } else {
    return(p)
  }

} #end function ppi.distAnalysis
