

#-----Get PDB info based on menu selection-----

#'Get PDB Info
#'
#'Get information about PDB selection from menu selection
#'
#'@param menu_selection Numerical selection that has been selected corresponding to the row number of pdb_hit_table
#'@param pdb_hit_table Hit table from bio3d::blast.pdb() function
#'@export

get_pdb_info <- function(menu_selection,pdb_hit_table){

  query_start <- pdb_hit_table$q.start[menu_selection]
  query_end <- pdb_hit_table$q.end[menu_selection]
  sequence_start <- pdb_hit_table$s.start[menu_selection]
  sequence_end <- pdb_hit_table$s.end[menu_selection]

  split_hit <- strsplit(pdb_hit_table$subjectids[menu_selection],'_')[[1]]
  pdb_id <- split_hit[1]
  chain <- split_hit[2]

  pdb_info <- list(query_start=query_start,
                   query_end=query_end,
                   sequence_start=sequence_start,
                   sequence_end=sequence_end,
                   pdb_id=pdb_id,
                   chain=chain)

  return(pdb_info)

}
