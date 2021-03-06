% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdb.pymol.R
\name{rbd.pymol}
\alias{rbd.pymol}
\title{Write PyMol file for RBDmap Data}
\usage{
rbd.pymol(
  bs_output,
  color_by = "binding_sequence",
  colors = NULL,
  file.name = "rbd_pymol.pml",
  write.file = TRUE,
  experiment.dir = NULL,
  gray0 = FALSE,
  heatmap = TRUE,
  assembly = 0,
  fetch = TRUE
)
}
\arguments{
\item{bs_output}{Output from rbd.getBindingSeq()}

\item{color_by}{The variable that will be used to determine color. Must a column in bs_output or use "freq" to color by frequency of amino acid position in binding sequences. Defaults to "binding_sequence".}

\item{colors}{Colors to be used for color.pymol() function.}

\item{file.name}{If write.file = T, name of the output file. Defaults to 'rbd_pymol.pml'}

\item{write.file}{Boolean. If TRUE, will output file with name of file.name. If FALSE, will return file as list.}

\item{experiment.dir}{Directory of experiment, where PDB files are. If left as NULL, will use current directory.}

\item{gray0}{Boolean. If TRUE, will shift color scheme for frequency analysis and label 0s as "gray" in PyMOL.}

\item{assembly}{Assembly number for PyMOL output that corresponds to RCSB database. Works for PyMOL 1.8 and above.}
}
\description{
This function makes a PyMOL file (.pml) that can be used to visualize data obtained from rdb.getBindingSeq()
}
