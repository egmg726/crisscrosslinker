% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdb_functions2.R
\name{rbd.makeIEPlot}
\alias{rbd.makeIEPlot}
\title{Make input/eluate plot}
\usage{
rbd.makeIEPlot(
  input_eluate_table,
  prefixes,
  proteases = c("LysC", "ArgC"),
  title = "Eluate vs. Input",
  palette = "Dark2",
  fill = "#eaf2ff",
  colour = "#eaf2ff",
  size = 0.5,
  pipeline_generated_plot_title = "pipeline_generated_plot.png"
)
}
\arguments{
\item{input_eluate_table}{Output from make_input_eluate_table()}

\item{prefixes}{List of control and experimental prefixes}

\item{proteases}{List of secondary prefixes in file names, defaults to proteases used 'ArgC' and 'LysC'}

\item{palette}{Palette from RColorBrewer to be used for plot. Defaults to 'Dark2'}

\item{fill}{Same as "fill" from ggplot2 for plot}

\item{colour}{Same as "colour" from ggplot2 for plot}

\item{size}{Same as "size" from ggplot2 for plot}

\item{pipeline_generated_plot_title}{Title for generated plot. Must end in '.png'. Defaults to 'pipeline_generated_plot.png'.}
}
\description{
This function takes the data from the input/eluate table and makes it into a plot
}
\author{
Emma Gail
}
