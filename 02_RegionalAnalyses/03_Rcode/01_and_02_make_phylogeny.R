# this script is intended for toydata only 
rm(list = objects())
## the plan is to combine all the relevant info for my my project from madidi and ABERG
## small and large plots plots. 


################################################################################
#### 1. PRELIMINARIES ##########################################################
################################################################################

## Set the seed for randomizations - so that no elev.rangess change
set.seed(1981)


#Set your working directory to 02_RegionalAnalyses folder.
setwd("~/Madidi_ToyData/02_RegionalAnalyses")
# setwd("C:/Users/stello/Dropbox/Work/01_MyProjects/ESCvsRAD_CentralAndes/02_RegionalAnalyses/")


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
trees.merged <- read.csv("./02_Data/MergedCleanData/ToyData_trees_merged.txt", 
                         header=TRUE, as.is=TRUE, sep = "\t", fileEncoding="UTF-8")


# uniqueify my combined plot tree data based on SpeciesName
unique.species <- distinct(trees.merged, SpeciesName, .keep_all= TRUE)

species.data <- unique.species


################################################################################
#### 3. CREATE PHYLOMAKER PHYLOGENY ############################################
################################################################################
# take only data for phylomaker input into new dataframe
phylomaker.input <- data.frame("Species" = species.data$SpeciesName, 
                               "Genus" = species.data$Genus, "Family" = species.data$Family, 
                               "species.relative" = NA, "genus.relative" = NA)

# write csv for phylomaker input for merged.tree data from all plots
write.csv(phylomaker.input, 
          file = "02_Data/input_phylomaker_list.csv", 
          row.names=FALSE, na="")

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


# generate the tree using scenario 1 settings (see V.phylomaker manual for details)
phylomaker.phylo.all.S1 <- phylo.maker(sp.list = rel$species.list, tree = rel$phylo, 
                                       nodes = rel$nodes.info, scenarios="S1") 

write.tree(phylomaker.phylo.all.S1$scenario.1, 
           "04_Results/phylomaker_phylo_all_S1.tre")

# generate the tree using scenario 3 settings (see V.phylomaker manual for details)
phylomaker.phylo.all.S3 <- phylo.maker(sp.list = rel$species.list, tree = rel$phylo, 
                                       nodes = rel$nodes.info, scenarios="S3") 
# write output phylogeny
write.tree(phylomaker.phylo.all.S3$scenario.3, 
           "04_Results/phylomaker_phylo_all_S3.tre")

# write tsv for your whole dataset of unique species
write.table(species.data, "02_Data/MergedCleanData/species_data.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")

#these both should not be different. (return "charecter(0)")
setdiff(phylomaker.phylo.all.S3$scenario.3$tip.label, 
        unique(paste(trees.merged$Genus, trees.merged$Species, sep="_")))

setdiff(unique(paste(trees.merged$Genus, trees.merged$Species, sep="_")),
        phylomaker.phylo.all.S1$scenario.1$tip.label)



