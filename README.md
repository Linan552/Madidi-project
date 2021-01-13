# Madidi-project IN PROGRESS not working currently
Scripts used in analysis of phylogenetic turnover across environmental gradients (elevation in this case)
## intended use 
These scripts represent those used for our recent manucript:   
"The evolutionary assembly of forest communities along environmental gradients: recent diversification or sorting of pre-adapted clades?"   
https://www.biorxiv.org/content/10.1101/2020.12.22.424032v1   
Once you have the appropriate packages (see list below) download the entire 02_RegionalAnalyses directory and you should be able to run the scripts within /02_RegionalAnalyses/03_Rcode/ (after specifying the appropriate working directory for you machine). Run the scripts sequentially (01 and 02, 03 ... etc).     
A results folder will be created with figures, stats, etc. 
Note: We have provided a toy dataset (made up of fictional data) to run the scripts. Original data used in our study can be found at https://doi.org/10.5281/zenodo.4276558. 

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
---
## what each script does
### '01_and_02_make_phylogeny.R' 
uses the tree dataset to generate a time calibrated phylogeny by placing taxa into Smith and Brown’s (2018) global mega-phylogeny of seed plants  (via V.phylomaker).   
**INPUTS:** 'ToyData_trees.txt' *A dataset representing individuals trees surveyed within a set of plots*   
**OUTPUTS:** 
- 'diferrentGeneraList.txt' *A list of genera that are not found in the Smith et al's mega phylogeny*
- 'tip_info.txt' *A table of infomation for each of the "tips" of the phylogeny (species names, genera, family, etc)*   
- 'species_data.txt' *A list of unique species and morphospecies found within 'ToyData_trees.txt'*
- 'phylomaker_phylo_all_S3.tre' *phylogeny using scenario 3 (see V.phylomaker manual)*
- 'phylomaker_phylo_all_S1.tre' *phylogeny using scenario 1 (see V.phylomaker manual)*
---
### '03_turnover_decomposition_2021-01-11.R' 
divides phylogeny into clades based on specified age (30 and 60 my old) and decomposes turnover in two ways (Bray-Curtis and Sorenson), into among- and within-clade components.  In addition, the script creates a null model of turnover.    
NOTE: this script will do all these calculatations of turnover on a variety of datasets (with or without morphospecies, varying clade ages) and will output results of all these in the /04_Results/ dierctory.    
**INPUTS:** 
- 'ToyData_trees.txt' *A dataset representing individuals trees surveyed within a set of plots*  
- 'ToyData_plots.txt' *A dataset containing info about the plots from which trees were surveyed* 
- 'phylomaker_phylo_all_S3.tre' or 'phylomaker_phylo_all_S1.tre' *Depending on options used*   
**OUTPUTS:** 
- 'BC.total.emp.rds' *Bray-Curtis turnover between plots (your empirical data)*
- 'S.total.emp.rds' *Sorenson turnover between plots (your empirical data)*
- 'BC.ag.null.rds'  *null model Bray-Curtis among-clade component (proportion) of turnover*
- 'BC.wg.null.rds'  *null model Bray-Curtis within-clade component (proportion) of turnover*
- 'BC.BC.ws.null.rds'  *null model Bray-Curtis within-species component (proportion) of turnover*
- 'S.ag.null.rds'  *null model Sorenson among-clade component (proportion) of turnover*
- 'S.wg.null.rds'  *null model Sorenson within-clade component (proportion) of turnover*
- 'SpeciesData.csv'  
- 'CladeData.csv'  
- 'clade.data.null.rds' 
---
### '04_visualize_EVERYTHING_2021-01-12.R' 
generates a variety of figures, visualizing components of turnover across a variety of datasets generated in above scripts.   
NOTE: plot comparisons are sub-selected to create datasets that compare turnover components across elevational distance while minimizing effects of geographic distance and vice-versa.    
**INPUTS:** 
- 'BC.total.emp.rds' **or** 'S.total.emp.rds' *(must specify Sorenson or Bray-Curtis)*
- 'BC.ag.null.rds', BC.wg.null.rds', and BC.ws.null.rds'; **or** 'S.ag.null.rds' and 'S.wg.null.rds' 
- 'SpeciesData.csv'  
- 'CladeData.csv'  
- 'clade.data.null.rds'  
- 'S.wg.null.rds'  

**OUTPUTS:** 
- 'TurnoverGradients.tiff'  *graph of total turnover of plots across elevational distance*
- 'ProportionOfTurnoverGradients.tiff'  *graphs comparing components of turnover to null models. Examining across elevational gradient (top) and geographcic distance (bottom)*
- 'DifferenceBetweenTurnoverComponentsGradients.tiff'  *graphs comparing just the among-component of turnover to null models (similar to 'ProportionOfTurnoverGradients.tiff')*
- 'LogitProportionOfTurnoverGradients.tiff' *same graphs as 'ProportionOfTurnoverGradients.tiff' but with logit transformed regressions.*
- 'GeoVsElevDists.tiff' *graph showing the selection of plot comparisions used across geographic and elevational distance*
- 'SppNByClade.tiff'  *graphs depicting the distribution of species across clades.*
- 'CladeRangeSizes.tiff'  *elevational ranges of clades* 
- 'Results_summary.txt' 
