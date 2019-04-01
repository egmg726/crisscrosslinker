#' Match PDB 3
#'
#' @param xlms.df xlms.df
#' @param pdb_numbering pdb_numbering
#' @export

ppi.matchPDB3 <- function(xlms.df,pdb_numbering){

  stored_pdbs <- list()

  pdb_match_vector <- pdb_numbering

  xlink_df <- xlms.df
  xlink_df$dist <- rep(NA,nrow(xlink_df))
  xlink_df$pdb1 <- rep(NA,nrow(xlink_df))
  xlink_df$pdb2 <- rep(NA,nrow(xlink_df))
  xlink_df$pro_pos1 <- rep(NA,nrow(xlink_df))
  xlink_df$pro_pos2 <- rep(NA,nrow(xlink_df))

  for(row_num in 1:nrow(xlms.df)){
    xlink_row <- xlms.df[row_num,]
    otps_list <- c()
    xyz_coord_list <- list()

    for(pro_num in 1:2){
      protein_name <- as.character(xlink_row[[paste0('pro_name',pro_num)]])
      protein_pos <- as.character(xlink_row[[paste0('pro_pos',pro_num)]])

      protein_pos_match <- match(protein_pos,pdb_match_vector[[protein_name]]$fasta)
      pdb_protein_pos <- as.numeric(as.character(pdb_match_vector[[protein_name]]$pdb[protein_pos_match]))
      #xlink_row[[paste0('pro_pos',pro_num)]] <- pdb_protein_pos
      xlink_df[row_num,paste0('pro_pos',pro_num)] <- pdb_protein_pos

      on_this_pdb_structure <- pdb_match_vector[[protein_name]]$chain[protein_pos_match]
      xlink_df[row_num,paste0('pdb',pro_num)] <- on_this_pdb_structure

      #get the xyz coordinates from here and then add them to list
      pdb_chain <- strsplit(on_this_pdb_structure,'_')[[1]]
      pdb <- pdb_chain[1]
      chain <- pdb_chain[2]

      if(!is.na(pdb) && (!is.na(chain))){

        if(pdb %in% names(stored_pdbs)){
          pdb_read <- stored_pdbs[[pdb]]
        } else {
          pdb_read <- read.pdb2(pdb)
          stored_pdbs[[pdb]] <- pdb_read
        }

        pdb_read$atom <- pdb_read$atom[pdb_read$atom$chain == chain,]
        pdb_read_chain <- pdb_read$atom[pdb_read$atom$resno == pdb_protein_pos,]


        if(unique(pdb_read_chain$resid) != 'LYS'){
          cat(paste0('Possible mismatch,',unique(pdb_read_chain$resid),'\n'))
          #return(xlink_row)
          xyz_coord_list[[pro_num]] <- NA

        } else {
          xyz_matches <- as.numeric(pdb_read_chain[pdb_read_chain$elety == 'CA',][c('x','y','z')])
          xyz_coord_list[[pro_num]] <- xyz_matches
        }
        #cat(paste0(unique(pdb_read_chain$resid),'\n'))


      } else {#end if(!is.na(pdb) && (!is.na(chain)))

        xyz_coord_list[[pro_num]] <- NA

      } #end else to if(!is.na(pdb) && (!is.na(chain)))

    } #end for(pro_num in 1:2){

    if(!(NA %in% xyz_coord_list)){
      dist_xyz <- dist.xyz(xyz_coord_list[[1]],xyz_coord_list[[2]])[1,1]
      xlink_df[row_num,'dist'] <- dist_xyz

    }#end if(!(NA %in% xyz_coord_list)){


  } #end for(row_num in 1:nrow(xlms.df)){

  return(xlink_df)

} #end function ppi.matchPDB3

#
# xlms.df.mpdb3 <- ppi.matchPDB3(xlms.df,pdb_numbering)
#
# xlms.df.mpdb3 <- na.omit(xlms.df.mpdb3)
#
# xlms.df.mpdb3[xlms.df.mpdb3$dist > 80,]

# na.omit(xlms.df.mpdb3)
#
# ppi.distAnalysis(na.omit(xlms.df.mpdb3))
#
# ppi.pymol(na.omit(xlms.df.mpdb3))




