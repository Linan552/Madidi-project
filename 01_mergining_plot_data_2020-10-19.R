
rm(list = objects())
## the plan is to combine all the relevant info for my my project from madidi PT, PP, and ABERG plots. 


################################################################################
#### 1. PRELIMINARIES ##########################################################
################################################################################

## Set the seed for randomizations - so that no elev.rangess change
set.seed(1981)


# setwd("C:/Users/alx55/Dropbox/ESCvsRAD_CentralAndes/02_RegionalAnalyses/")
# setwd("C:/Users/stello/Dropbox/Work/01_MyProjects/ESCvsRAD_CentralAndes/02_RegionalAnalyses/")
setwd("/mnt/drive_2_500gb/ESCvsRAD_CentralAndes/02_RegionalAnalyses/")

## Packages to open
library(dplyr)
library(data.table)
library(stringr)



################################################################################
#### 2. OPEN DATA ##############################################################
################################################################################


### MADIDI DATA ################################################################

tree.data.PP <- 
  read.csv("02_Data/MadidiRawAndCleanData_v4.1/03_CleanMadidiData/04_TreeData_PP_Clean_v4.1_2020-09-07.txt",
    header=TRUE, as.is=TRUE, sep = "\t")

tree.data.PT <- 
  read.csv("02_Data/MadidiRawAndCleanData_v4.1/03_CleanMadidiData/05_TreeData_PT_Clean_v4.1_2020-09-07.txt",
    header=TRUE, as.is=TRUE, sep = "\t", fileEncoding="UTF-8")

plot.data <- read.csv("02_Data/MadidiRawAndCleanData_v4.1/03_CleanMadidiData/01_PlotData_Clean_v4.1_2020-09-07.txt", 
  header=TRUE, as.is=TRUE, sep = "\t", fileEncoding="UTF-8")


### ABERG DATA #################################################################

# aberg data was modified slightly from dropbox download. some headings were set 
# to match madidi plotdata file 
aberg.data.PP <- read.csv("./02_Data/ABERG/ABERG_1-haPlots_WFR_2019-Sebas.csv", 
  header=TRUE, as.is=TRUE, sep = ",", fileEncoding="UTF-8")

aberg.data.PT <- read.csv("./02_Data/ABERG/0.1ha-plots-WFR-2020-edited-headers.csv", 
                          header=TRUE, as.is=TRUE, sep = ",", fileEncoding="UTF-8")

################################################################################
#### 3. MODIFY DATA ############################################################
################################################################################

### REPLACE SPECIAL CHARACTERS #################################################

# this can't display tilde n so I replace the unknown charecter with n 
# (this is definately not UTF-8)
tree.data.PP$PlotTreeCode<- str_replace(tree.data.PP$PlotTreeCode, "\\??", "n")
tree.data.PP$PlotName<- str_replace(tree.data.PP$PlotName, "\\??", "n")

tree.data.PP$PlotTreeCode<- str_replace(tree.data.PP$PlotTreeCode, "?", "n")
tree.data.PP$PlotName<- str_replace(tree.data.PP$PlotName, "?", "n")


tree.data.PT$PlotTreeCode<- str_replace(tree.data.PT$PlotTreeCode, "\\??", "n")
tree.data.PT$PlotName<- str_replace(tree.data.PT$PlotName, "\\??", "n")

tree.data.PT$PlotTreeCode<- str_replace(tree.data.PT$PlotTreeCode, "?", "n")
tree.data.PT$PlotName<- str_replace(tree.data.PT$PlotName, "?", "n")


# replace the tilde n in the the plotfile (is actually coded in UTF-8) so 
# that we can match to tree data
plot.data$PlotName<- str_replace(plot.data$PlotName, "\\??", "n")
plot.data$PlotName<- str_replace(plot.data$PlotName, "?", "n")


### SUBSET RELEVANT COLUMNS AND FILTER OUT UNIDENTIFIED INDIVS. ################
tre.pp <- tree.data.PP[tree.data.PP$TypeOfName!="InvalidSpeciesName", ]
tre.pp <- tre.pp[,c("PlotName","TypeOfName","Family","Genus","Species")]
# tre.pp <- filter(tree.data.PP %>% 
#     select(PlotName, TypeOfName, Family, Genus, Species), 
#   TypeOfName != "InvalidSpeciesName")

tre.pp <- filter(tre.pp, Species != ".NoName")
tre.pp <- filter(tre.pp, Genus != ".NoName")
tre.pp$Genus <- gsub("Monteverdia", "Maytenus", tre.pp$Genus)

tre.pp$PlotType <- "PP" #add PlotType

tre.pt <- tree.data.PT[tree.data.PT$TypeOfName!="InvalidSpeciesName", ]
tre.pt <- tre.pt[,c("PlotName","TypeOfName","Family","Genus","Species")]
# tre.pt <- filter(tree.data.PT %>% 
#     select(PlotName, TypeOfName, Family, Genus, Species), 
#   TypeOfName != "InvalidSpeciesName")

tre.pt <- filter(tre.pt, Species != ".NoName")
tre.pt <- filter(tre.pt, Genus != ".NoName")
tre.pt$Genus <- gsub("Monteverdia", "Maytenus", tre.pt$Genus)
tre.pt$PlotType <- "PT" # add plottype


################################################################################
#### 4. MERGE  MADIDI DATASETS ########################################################
################################################################################

### MERGE THE TWO MADIDI TREE DATASETS AND ADD ELEVATION #######################

tree.madidi.merged <- bind_rows(tre.pp, tre.pt)

#add new column with elevation, searching plot.data elevation column matching the plotname in merged madidi
# this gives the elevation for every individual tree in the dataset 
tree.madidi.merged$MinimumElevationMeters <- 
  plot.data$MinimumElevationMeters[match(tree.madidi.merged$PlotName, plot.data$PlotName)]


### PREPARE ABERG TREE DATA TO MERGE ###########################################
#remove herbaceous plants in PT plots (keep, tree ferns, palms, and lianas)
nrow(aberg.data.PT)
aberg.data.PT <- aberg.data.PT[aberg.data.PT$life.form != "herbaceous", ]

# aberg.data.PP <- aberg.data.PP[aberg.data.PP$life.form != "PT", ]
# aberg.data.PP <- aberg.data.PP[aberg.data.PP$life.form != "TF", ]
# aberg.data.PP <- aberg.data.PP[aberg.data.PP$life.form != "LI", ]
# aberg.data.PT <- aberg.data.PT[aberg.data.PT$life.form != "fern", ]
# aberg.data.PT <- aberg.data.PT[aberg.data.PT$life.form != "palm", ]
# aberg.data.PT <- aberg.data.PT[aberg.data.PT$life.form != "liana", ]

# filter data to select columns of interest and remove indet samples
tree.aberg.PP <- aberg.data.PP[aberg.data.PP$validname!="FALSE", ]

tree.aberg.PP <- tree.aberg.PP[,c("plot","elev","family","genus","specie")]
tree.aberg.PT <- aberg.data.PT[,c("plot","elev","family","genus","specie")]
# tree.aberg <- filter(aberg.data %>% 
#     select(plot, elev, family, genus, specie, validname), 
#   validname != "FALSE")

colnames(tree.aberg.PP) <- c("PlotName", "MinimumElevationMeters", 
  "Family", "Genus", "Species")

colnames(tree.aberg.PT) <- c("PlotName", "MinimumElevationMeters", 
                             "Family", "Genus", "Species")
tree.aberg.PP$PlotType <- "abergPP" #add PlotType
tree.aberg.PT$PlotType <- "abergPT" #add PlotType

tree.aberg.PP <- filter(tree.aberg.PP, Species != "indet")

tree.aberg.PT <- filter(tree.aberg.PT, Species != "indet")
## merge ABERG plots 
tree.aberg.merged <- bind_rows(tree.aberg.PP, tree.aberg.PT)

# create new column (TypeOfName) and if there is a number in the species name column, 
#  then write "morphospecies" and if not write in "speciesname"
tree.aberg.merged$TypeOfName <- ifelse(grepl("\\d", tree.aberg.merged$Species), 
  "MorphoSpeciesName", "SpeciesName")

# since the aberg PT plots data does not have a "Type of Name" section and genera are occationally
# given a "morpho genus" number, I will get rid of this as they are essentially indets
# create new column (TypeOfName) and if there is a number in the species name column, 
#  then write "morphoGenus" and if not write in "RealGenus"
tree.aberg.merged$Tempcol <- ifelse(grepl("\\d", tree.aberg.merged$Genus), 
                                    "MorphoGenusName", "RealGenus")
tree.aberg.merged <- filter(tree.aberg.merged, Tempcol != "MorphoGenusName")
tree.aberg.merged <- select(tree.aberg.merged, -Tempcol) #remove our temporary column

#reorganize columns to match mergered.madi 
tree.aberg.org <- tree.aberg.merged %>% 
  dplyr::select(PlotName, TypeOfName, Family, Genus, Species, MinimumElevationMeters, PlotType)


### MERGE MADIDI AND ABERG TREE DATA ###########################################

trees.merged <- bind_rows(tree.madidi.merged, tree.aberg.org)

## NOTE: This removes all the specimen numbers from the ABERG morphospecies
trees.merged$Species <- gsub("\\(.*", "", trees.merged$Species)

## NOTE: This removes all empty spaces in species names
trees.merged$Family <- gsub(" ", "", trees.merged$Family)
trees.merged$Genus <- gsub(" ", "", trees.merged$Genus)
trees.merged$Species <- gsub(" ", "", trees.merged$Species)

trees.merged$SpeciesName <- 
  paste(trees.merged$Genus, trees.merged$Species, sep=" ")

## NOTE: This standardizes some family names between our data and phylomaker
trees.merged$Family[trees.merged$Family == "Ximeniaceae"] <- "Olacaceae"
trees.merged$Family[trees.merged$Family == "Viburnaceae"] <- "Adoxaceae"
trees.merged$Family[trees.merged$Genus == "Chaetocarpus"] <- "Euphorbiaceae"

## ELIMINATE SPECIES OF CACTI, AND, BAMBU, from merged dataset #########################
if(sum(trees.merged$Family == "Cactaceae")>0)
  trees.merged <- trees.merged[-which(trees.merged$Family == "Cactaceae"),]

if(sum(trees.merged$Family == "Poaceae")>0)
  trees.merged <- trees.merged[-which(trees.merged$Family == "Poaceae"),]
### MERGE MADIDI AND ABERG PLOT DATA ###########################################

# take the information of interest from aberg dataset
aberg.plots.PT <- unique(aberg.data.PT[,c("plot", "elev", "lat", "long")])
colnames(aberg.plots.PT) <- c("PlotName", "MinimumElevationMeters", "lat", "long")
aberg.plots.PT$PlotType <- "abergPT"

aberg.plots.PP <- unique(aberg.data.PP[,c("plot", "elev", "lat", "long")])
colnames(aberg.plots.PP) <- c("PlotName", "MinimumElevationMeters", "lat", "long")
aberg.plots.PP$PlotType <- "abergPP"

# take the information of interest from madidi plots dataset
madidi.plots <- plot.data[,c("PlotName", "MinimumElevationMeters", "LatitudeDecDeg", "LongitudeDecDeg", "PlotType")]
names(madidi.plots)[3]<-"lat"
names(madidi.plots)[4]<-"long"

# merge the plot data 
plots.merged <- bind_rows(aberg.plots.PP, aberg.plots.PT, madidi.plots)


### REMOVE PLOTS THAT WONT BE USED IN ANALYSES #################################

spp.before <- unique(trees.merged$SpeciesName)

compo.temp <- table(trees.merged$PlotName, trees.merged$SpeciesName)
rich.temp <- rowSums(compo.temp>0)

plot(plots.merged$MinimumElevationMeters ~ plots.merged$lat)
abline(h=3800)

# remove plots too far from study areas, too high above treeline, and with low species diversiy
remove.far <- plots.merged$PlotName[plots.merged$long > -67]
remove.high <- plots.merged$PlotName[plots.merged$MinimumElevationMeters > 3800]
remove.lowdiv <- names(rich.temp)[rich.temp<=3]
remove.all <- unique(c(remove.far, remove.high, remove.lowdiv))

plots.merged <- plots.merged[-which(plots.merged$PlotName %in% remove.all),]
trees.merged <- trees.merged[!(trees.merged$PlotName %in% remove.all),]

spp.after <- unique(trees.merged$SpeciesName)

#see the species that have been removed 
setdiff(spp.before, spp.after)

# DONT WORRY ABOUT THIS# write a small list of species that need to be removed from phyomaker input file. 
# write.table(setdiff(spp.before, spp.after), "02_Data/MergedCleanData/spp_to_delete.txt", 
#   row.names = FALSE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")


### WRITE FILES ################################################################

#some family names in our dataset have to be corrected (subfamily names used)
trees.merged[] <- lapply(trees.merged, gsub, pattern = "Cordiaceae", replacement = "Boraginaceae", fixed = TRUE)
trees.merged[] <- lapply(trees.merged, gsub, pattern = "Actinidaceae", replacement = "Actinidiaceae", fixed = TRUE)
trees.merged[] <- lapply(trees.merged, gsub, pattern = "Staphylleaceae", replacement = "Staphyleaceae", fixed = TRUE)
trees.merged[] <- lapply(trees.merged, gsub, pattern = "Lissocarpaceae", replacement = "Ebenaceae", fixed = TRUE)


# write csv for your whole dataset of merged plots
write.table(trees.merged, "02_Data/MergedCleanData/trees_merged.txt", 
  row.names = FALSE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")


# write csv for your whole dataset of merged plots
write.table(plots.merged, "02_Data/MergedCleanData/plots_merged.txt", 
  row.names = FALSE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")



################################################################################
#### 5. CALCULATE ELEVATIONAL RANGES AND CREATE PHYLOMAKER INPUT ###############
################################################################################

### ELEVATIONAL RANGES ###
## This caluculats the range in elevation of each species in plots by using
## the MinElevationMeters info for each indivdual tree

spp.list <- unique(trees.merged$SpeciesName)

elev.ranges <- data.frame(matrix(NA, ncol=2, nrow=length(spp.list)))
rownames(elev.ranges) <- spp.list
colnames(elev.ranges) <- c("minElevation", "maxElevation")

for(i in 1:length(spp.list))
{
  sp.i <- spp.list[i]
  data.i <- base::subset(trees.merged, SpeciesName == sp.i, 
    select = "MinimumElevationMeters") 
  elev.ranges[i,] <- range(as.numeric(data.i$MinimumElevationMeters), na.rm = TRUE)
}

###################################
#### 9. combine range dataframe with unique species dataframe #######
################################################################################

# uniqueify my combined plot tree data based on SpeciesName
unique.species <- distinct(trees.merged, SpeciesName, .keep_all= TRUE)
identical(unique.species$SpeciesName, rownames(elev.ranges))

species.data <- bind_cols(unique.species, elev.ranges)

# write tsv for your whole dataset of unique species with elevational ranges
write.table(species.data, "02_Data/MergedCleanData/species_data.txt", 
  row.names = FALSE, quote = FALSE, sep = "\t", fileEncoding="UTF-8")


################################################################################
#### 10. write to file and write phylomaker input ##############################
################################################################################

# take only data for phylomaker input into new dataframe
phylomaker.input <- data.frame("Species" = species.data$SpeciesName, 
  "Genus" = species.data$Genus, "Family" = species.data$Family, 
  "species.relative" = NA, "genus.relative" = NA)

# write csv for phylomaker input for merged.tree data from all plots
write.csv(phylomaker.input, 
  file = "02_Data/MergedCleanData/input_phylomaker_list.csv", 
  row.names=FALSE, na="")

