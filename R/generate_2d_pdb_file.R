#' Generate 2D PDB File
#'
#' This function allows you to turn a list of amino acids into a 2D PDB file that can be downloaded.
#' @param amino_acid_fragments A vector containing an amino acid sequence.
#' @param start_residue Start residue
#' @param x_coord X coordinate
#' @param y_coord Y coordinate
#' @param Z coordinate
#' @keywords pdb
#' @author Emma Gail
#' @export

generate_2d_pdb_file <- function(amino_acid_fragments,start_residue,
                                 x_coord = 82.998, y_coord = -5.760,
                                 z_coord = -20.345){
  #x_coord <- 82.998
  #y_coord <- -42.760
  #y_coord <- -5.760
  #z_coord <- -10.345
  #z_coord <- -20.345
  downward <- TRUE
  z_negative <- FALSE
  res_num <- start_residue

  new_pdb_file <- c()

  for(num in 1:length(amino_acid_fragments)){
    if(num != length(amino_acid_fragments)){
      new_line <- '\n'
    } else {
      new_line <- ''
    }

    num_atom_spaces <- 7 - nchar(num)
    atom_spaces <- paste(rep(" ",num_atom_spaces),collapse = "")
    two_spaces <- "  "
    atom_number <- num
    residue_name <- toupper(aaa(amino_acid_fragments[num]))
    residue_number <- res_num
    six_spaces <- "      "
    num_residue_spaces <- 4 - nchar(residue_number)

    if(downward == TRUE){
      x_coord <- x_coord - 0.5
      y_coord <- y_coord - 0.6
      z_coord <- z_coord - 1.3
      downward <- FALSE

    } else {
      x_coord <- x_coord + 0.2
      y_coord <- y_coord - 1.5
      if(z_negative == TRUE){
        z_coord <- z_coord - 0.05
        z_negative <- FALSE
      } else {
        z_coord <- z_coord + 0.05
        z_negative <- TRUE
      }
      downward <- TRUE
    }


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
    res_num <- res_num + 1
  }

  new_pdb_file <- paste(new_pdb_file,sep="", collapse="")

  return(new_pdb_file)

}
