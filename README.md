**NSMB Tutorials are now available under the 'tutorials' folder as RMarkdown files. The HTML workflows for [protein-protein interactions](https://egmg726.github.io/crisscrosslinker/nsmb_pipeline_ppi.html) and [protein-RNA interactions](https://egmg726.github.io/crisscrosslinker/nsmb_pipeline_rbdmap.html) are available via the links provided.**

![logo](https://raw.githubusercontent.com/egmg726/crisscrosslinker/master/images/crisscrosslinker_logo.png)

This package is designed to analyze and visualize crosslinking (XL) data from multiple sources (mainly XL-MS and RBDmap).

**Analysis**

-   Differential analysis

-   Crosslinking frequency

-   Creation of PDB files for missing peptides

-   Mutation and Uniprot Integration (*coming soon!*)

**XL-MS**

-   Control distances for PDB structures

**RBDMap**

-   Frequency of RNA-likened amino acids

**Visualization**

-   Distance histograms

-   Amino acid frequency heatmaps

-   2D: xiNET output

-   3D: PyMOL output

Installation
------------


To download the functions with documentation:
``` r
library('devtools')
install_github("egmg726/crisscrosslinker")
```

**Note**

This is still the beta release of the package/documentation so please contact me if you come across any issues!
