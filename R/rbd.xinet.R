#'RBDmap: xiNET Output
#'
#'This function outputs an annotation file for import into the web server xiNET. Designed to be used in conjunction with output from ppi.xinet().
#'
#'@param bs.output Output from function rbd.getBSfromIET()
#'@param colors Colors to be used for the annotation file. Can be an RColorBrewer palette name, virdis palette name, or list of hexcodes.
#'@param database Database to be used for the visualization. Can choose betwen "PDB", "UniProt", or "FASTA". Defaults to UniProt.
#'@param write.file If TRUE, will write to a .csv file
#'@param file.name File name for the written out file.
#'@export


rbd.xinet <- function(bs.output,colors=c('#d0d0d0','#2f5ac6','#e50000'),
                      database='UniProt',
                      write.file=TRUE,file.name='rbd_xinet.csv'){
  
  #needs to be all the same source, choose a resource
  #still need to write in resource
  #will send out a warning message if not 
  
  bs.output2 <- bs.output[bs.output$db == database,]
  bs_freqVector <- rbd.freqVector(bs.output2,heatmap = FALSE)
  viridis_colors <- c("magma","inferno","plasma","viridis","cividis")
  
  if((colors %in% rownames(brewer.pal.info)) || (colors %in% viridis_colors)){
    freq_vars <- sort(unique(unlist(bs_freqVector)))
    freq_colors <- color.pymol(freq_vars,colors=colors)
    hm.palette <- colorRampPalette(colors = freq_colors$hexcodes)
  } else {
    hm.palette <- colorRampPalette(colors = colors)
  }
  
  color.list <- hm.palette(length(unique(unlist(bs_freqVector))))
  
  anno.df <- data.frame()
  
  for(protein_name in names(bs_freqVector)){
    
    freq_vector <- bs_freqVector[[protein_name]]
    unique_nums <- as.numeric(unique(freq_vector))
    
    for(unum in unique_nums){
      
      freq_vector2 <- (1:length(freq_vector))[freq_vector == unum]
      list_ranges <- tapply(freq_vector2, cumsum(c(TRUE, diff(freq_vector2) > 2)), range)
      
      for(lr_num in 1:length(list_ranges)){
        start_end <- list_ranges[[lr_num]]
        lr_start <- start_end[1]
        lr_end <- start_end[2]
        
        anno_list <- list(ProteinId=protein_name,
                          AnnotName=paste0('Number of Repeats: ',unum), #string here establishing the number of repeats
                          StartRes=lr_start,
                          EndRes=lr_end,
                          Color=color.list[unum+1]) 
        
        anno.df <- rbind(anno.df,(data.frame(anno_list)))
        
      } #end for(lr_num in 1:length(list_ranges)){
    } #end for(unum in unique_nums){
  } #end for(protein_name in names(bs_freqVector)){
  
  
  if(write.file == TRUE){
    write.csv(anno.df,file.name,row.names=FALSE) #need to make sure that the number column is removed?
  } else {
    return(anno.df)
  }
  
  #make sure the write.csv function works
  
  
} #end function rbd.xinet


