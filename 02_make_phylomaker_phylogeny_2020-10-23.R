
rm(list = objects())
## the plan is to combine all the relevant info for my my project from madidi and ABERG
## small and large plots plots. 


################################################################################
#### 1. PRELIMINARIES ##########################################################
################################################################################

## Set the seed for randomizations - so that no elev.rangess change
set.seed(1981)

# setwd("C:/Users/alx55/Dropbox/ESCvsRAD_CentralAndes/02_RegionalAnalyses/")
# setwd("C:/Users/stello/Dropbox/Work/01_MyProjects/ESCvsRAD_CentralAndes/02_RegionalAnalyses/")
# setwd("C:/Users/alinan/Documents/02_RegionalAnalyses/")
setwd("/mnt/drive_2_500gb/ESCvsRAD_CentralAndes/02_RegionalAnalyses/")


## Packages to open
library(dplyr)
library(data.table)
library(V.PhyloMaker)
library(picante)
library(phytools)


################################################################################
#### 2. OPEN DATA ##############################################################
################################################################################

## Merged dataset (all plots) 
trees.merged <- read.csv("./02_Data/MergedCleanData/trees_merged.txt", 
  header=TRUE, as.is=TRUE, sep = "\t", fileEncoding="UTF-8")


## Species data for phylomaker #note: this version uses new version of input file excluding species found in excluded plots.
## note also: for some reason cyathea lechleri was in the input twice (because in one the genera Cyathea was capitolized and the other was not)
phylomaker.input <- read.csv("02_Data/MergedCleanData/input_phylomaker_list_OCT2020.csv") 



################################################################################
#### 3. CREATE PHYLOMAKER PHYLOGENY ############################################
################################################################################

plot.genera <- unique(trees.merged$Genus)
smith.genera <- unique(tips.info$genus)

# get list of genera that differ 
diff.genera <- setdiff(plot.genera, smith.genera)

write.table(diff.genera, "02_Data/MergedCleanData/diferrentGeneraList.txt", 
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")


# get list of families that differ between smith phylogeny and our dataset
plot.fam <- unique(trees.merged$Family)
smith.fam <- unique(tips.info$family)

setdiff(plot.fam, smith.fam) #should be two familys that don't match

write.table(tips.info, "02_Data/MergedCleanData/tip_info.txt", 
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t", 
  fileEncoding="UTF-8")


#### RUN PHYLOMAKER ############################################################
## note: also remeMber to rename to column headings with lower case family, 
##       species etc (see below). 
## note: input_phylomaker.list.csv was edited such that some species not  
##       included in the phylomaker phylogeny were assigned the "genus.relative" 
##       value was used to place the species at the base of the genus in the 
##       tree based on literature searches of more recent phylogenetic work.


# rename columns to be appropriate for phylmaker
names(phylomaker.input)[1]<-"species"
names(phylomaker.input)[2]<-"genus"
names(phylomaker.input)[3]<-"family"


# bind the "genus.relative" to related genera in the tree
rel <- bind.relative(sp.list = phylomaker.input, tree = GBOTB.extended, nodes = nodes.info.1)


# generate the tree
phylomaker.phylo.all.S1 <- phylo.maker(sp.list = rel$species.list, tree = rel$phylo, 
  nodes = rel$nodes.info, scenarios="S1") 
 
write.tree(phylomaker.phylo.all.S1$scenario.1, 
  "04_Results/phylomaker_phylo_all_S1.tre")


phylomaker.phylo.all.S3 <- phylo.maker(sp.list = rel$species.list, tree = rel$phylo, 
  nodes = rel$nodes.info, scenarios="S3") 

write.tree(phylomaker.phylo.all.S3$scenario.3, 
  "04_Results/phylomaker_phylo_all_S3.tre")

#these both should not be different. (return "charecter(0)")
setdiff(phylomaker.phylo.all.S3$scenario.3$tip.label, 
  unique(paste(trees.merged$Genus, trees.merged$Species, sep="_")))

setdiff(unique(paste(trees.merged$Genus, trees.merged$Species, sep="_")),
  phylomaker.phylo.all.S1$scenario.1$tip.label)

