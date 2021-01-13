# Madidi-project IN PROGRESS not working currently
Scripts used in analysis of phylogenetic turnover across environmental gradients
## intended use 
These scripts represent those used for our recent manucript:   
"The evolutionary assembly of forest communities along environmental gradients: recent diversification or sorting of pre-adapted clades?"   
https://www.biorxiv.org/content/10.1101/2020.12.22.424032v1   
Once you have the appropriate packages (see list below) download the entire 02_RegionalAnalyses directory and you should be able to run the scripts within /02_RegionalAnalyses/03_Rcode/ (after specifying the appropriate working directory for you machine). Run the scripts sequentially (01 and 02, 03 ... etc).     
A results folder will be created with figures, stats, etc. 
Note: We have provided a toy dataset (made up of fictional data) to run the scripts. Original data used in our study can be found at https://doi.org/10.5281/zenodo.4276558. 

## what each script does
'01_and_02_make_phylogeny.R' uses the tree dataset to generate a phylogeny (via V.phylomaker)   
INPUTS: 'ToyData_trees.txt'   
OUTPUTS: 
## what R-packages you will need to run these scripts
- library(dplyr)
- library(data.table)
- library(V.PhyloMaker)
- library(picante)
- library(phytools)
- library(vegan)
- library(latticeExtra)
- library(plyr)
- library(gtools)
- library(geosphere)
- library(grr)
