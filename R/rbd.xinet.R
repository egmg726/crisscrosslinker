#'RBDmap: xiNET Output
#'
#'This function outputs an annotation file for import into the web server xiNET. Designed to be used in conjunction with output from ppi.xinet().
#'
#'@param bs.output Output from function rbd.getBSfromIET()
#'@param colors Colors to be used for the annotation file. Can be an RColorBrewer palette name, virdis palette name, or list of hexcodes.
#'@param write.file If TRUE, will write to a .csv file
#'@param file.name File name for the written out file.
#'@export

rbd.xinet <- function(bs.output,colors=c('#d0d0d0','#2f5ac6','#e50000'),write.file=TRUE,file.name='rbd_xinet.csv'){
  
  bs_freqVector <- rbd.freqVector(bs.output,heatmap = FALSE)
  
  viridis_colors <- c("magma","inferno","plasma","viridis","cividis")
  
  
  if((colors %in% rownames(brewer.pal.info)) || (colors %in% viridis_colors)){
    freq_vars <- sort(unique(unlist(bs_freqVector)))
    freq_colors <- color.pymol(freq_vars,colors=colors)
    hm.palette <- colorRampPalette(colors = freq_colors$hexcodes)
  } else {
    hm.palette <- colorRampPalette(colors = colors)
  }
  
  color.list <- hm.palette(length(unique(unlist(bs_freqVector))))
  
  fvl_si_df <- data.frame()
  for(protein_name in names(bs_output_repeats_uniprot_freqvector)){
    freq_vector_list <- bs_output_repeats_uniprot_freqvector[[protein_name]]
    for(n_repeats in sort(unique(freq_vector_list))){
      
      fvl_sub_ind <-(1:length(freq_vector_list))[freq_vector_list == n_repeats]
      fvl_si <- c()
      #fvl_si_df <- data.frame()
      for(i in 1:length(fvl_sub_ind)){
        fvl_si <- c(fvl_si,fvl_sub_ind[i])
        if(is.sequential(fvl_si) == TRUE){
          next #if still sequentil, no break yet
        } else {
          fvl_si <- head(fvl_si,-1)
          fvl_si_min <- min(fvl_si) #get start
          fvl_si_max <- max(fvl_si) #get end
          #fvl_si_df <- rbind(fvl_si_df,data.frame(min=fvl_si_min,max=fvl_si_max,num=1))
          fvl_si_df <- rbind(fvl_si_df,data.frame(ProteinId=protein_name,AnnotName=paste0('Num of Repeats: ',as.character(n_repeats)),StartRes=fvl_si_min,EndRes=fvl_si_max,Color=color.list[n_repeats+1]))
          fvl_si <- c(fvl_sub_ind[i]) #reset the list 
        } #end else to if(is.sequential(fvl_si) == TRUE){
      } #end for(i in 1:length(fvl_sub_ind)){
      
      #repeat to get the last list 
      fvl_si_min <- min(fvl_si)
      fvl_si_max <- max(fvl_si)
      
      #what to do about the annotation name? (just the number of repeats)   
      fvl_si_df <- rbind(fvl_si_df,data.frame(ProteinId=protein_name,AnnotName=paste0('Num of Repeats: ',as.character(n_repeats)),StartRes=fvl_si_min,EndRes=fvl_si_max,Color=color.list[n_repeats+1]))
      
    } #end for(num_repeats in sort(unique(freq_vector_list))){
  } #end for(protein_name in names(bs_output_repeats_uniprot_freqvector)){
  
  if(write.file == TRUE){
    write.csv(fvl_si_df, file = file.name)
  }
  
  return(fvl_si_df)
  
  
} #end function xinet


#'Is Sequential
#'
#'Shows TRUE/FALSE if list contains sequential values
#'
#'@param x List of values to be checked.
#'@export


is.sequential <- function(x){
  all(abs(diff(x)) == 1)
} 
