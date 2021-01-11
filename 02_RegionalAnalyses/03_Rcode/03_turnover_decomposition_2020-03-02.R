
rm(list = objects())
## thia script generates turnover metrics for more combined Madidi and 
## ABERG plot-data. 
## NOTE: the Null model randomizations has been reduced in this toy script to 100 from 999 to facilitate speed.



################################################################################
#### 1. PRELIMINARIES ##########################################################
################################################################################

## Set the seed for randomizations - so that no results change
set.seed(1981)



# project.directory <- "C:/Users/stello/Dropbox/Work/01_MyProjects/ESCvsRAD_CentralAndes/02_RegionalAnalyses/"
# project.directory <- "C:/Users/alx55/Dropbox/ESCvsRAD_CentralAndes/02_RegionalAnalyses/"
# project.directory <- "C:/Users/alinan/Documents/02_RegionalAnalyses/"
project.directory <- "/mnt/drive_2_500gb/ESCvsRAD_CentralAndes/02_RegionalAnalyses/"

setwd(project.directory)


## Packages to open
library(dplyr)
library(data.table)
library(phytools)
library(vegan)
library(latticeExtra)

source("03_Rcode/decompose_my_beta_2019-09-30.R")


datasets <- c("LargePlots", "SmallPlots")
morpho.opts <- c("WithMorphos", "NoMorphos")
cutoffs <- c(60, 30)

#phylo.opt <- "S1"
phylo.opt <- "S3"



################################################################################
#### 2. OPEN DATA ##############################################################
################################################################################

#open raw tree data
trees.data <- read.csv("02_Data/MergedCleanData/trees_merged.txt", 
  header=TRUE, as.is=TRUE, sep = "\t", fileEncoding="UTF-8")

# edit species names in the raw tree data to include and underscore instead of a space
trees.data$SpeciesName <- gsub('\\s+', '_', trees.data$SpeciesName)


#open merged plot info
plots.data <- read.csv("02_Data/MergedCleanData/plots_merged.txt", 
  header=TRUE, as.is=TRUE, sep = "\t", fileEncoding="UTF-8")

#open merged plot info
if(phylo.opt == "S1")
  phylo.all <- read.tree("04_Results/phylomaker_phylo_all_S1.tre")

if(phylo.opt == "S3")
  phylo.all <- read.tree("04_Results/phylomaker_phylo_all_S3.tre")


setdiff(phylo.all$tip.label, unique(trees.data$SpeciesName))
sort(setdiff(unique(trees.data$SpeciesName), phylo.all$tip.label))



################################################################################
#### 3. RUN ANALYSES ###########################################################
################################################################################

for(i in 1:length(datasets))
{
  
  dataset.i <- datasets[i]
    
  
  #### SUBSET TREE DATA TO CUTOFF i AND DATASET j ############################
    
  if(dataset.i == "LargePlots")
    trees.data.i <- trees.data[trees.data$PlotType %in% c("PP", "abergPP"),]
    
  if(dataset.i == "SmallPlots")
    trees.data.i <- trees.data[trees.data$PlotType %in% c("PT", "abergPT"),]
    

  for(m in 1:length(morpho.opts))
  {   
    morpho.opt.m <- morpho.opts[m]

    if(morpho.opt.m == "WithMorphos")
      trees.data.im <- trees.data.i
      
    if(morpho.opt.m == "NoMorphos")
      trees.data.im <- trees.data.i[trees.data.i$TypeOfName!="MorphoSpeciesName",]
  
   
    for(j in 1:length(cutoffs))
    {
      cutoff.j <- cutoffs[j]
  
      print(c(dataset.i, morpho.opt.m, cutoff.j))
          
  
      #### MAKE CUTS IN THE PHYLOGENY ##############################################
      
      phylo.cuts.j <- treeSlice(phylo.all, trivial=TRUE, 390.703-cutoff.j)
      
      # make my lists of lists just a list of of the tip.label vectors in taken 
      # from each of the other lists. 
      clades.j <- lapply(phylo.cuts.j, "[[", "tip.label")
      names(clades.j) <- seq_along(clades.j)
      names(clades.j) <- paste0("clade_", names(clades.j))
      
      
      #merge the list of vectors into a dataframe
      groups.spp.j <- stack(clades.j)
      groups.spp.j <- groups.spp.j[c(2,1)]
      names(groups.spp.j) <- c("groups", "spp")
      groups.spp.j <- groups.spp.j[order(groups.spp.j$spp),]
          

      #### SUBSET THE DATA TO options i, m and j ###############################
      
      trees.data.imj <- 
        trees.data.im[as.character(trees.data.im$SpeciesName) %in% groups.spp.j$spp,]
      
      groups.spp.imj <- groups.spp.j[groups.spp.j$spp %in% 
          as.character(trees.data.imj$SpeciesName),]
  
      plots.data.imj <- plots.data[plots.data$PlotName %in% 
          as.character(trees.data.imj$PlotName),]  
      plots.data.imj <- plots.data.imj[order(plots.data.imj$PlotName),]
      
      ## Add the clade data to the tree data
      trees.data.imj <- cbind(trees.data.imj, 
        groups.spp.imj[match(trees.data.imj$SpeciesName, groups.spp.imj$spp),])       
  

      #### EMPIRICAL SPECIES COMPOSITION TABLE #################################
      
      compo.imj <- table(trees.data.imj$PlotName, 
        trees.data.imj$SpeciesName)
      
      compo.imj <- compo.imj[order(rownames(compo.imj)),]
      compo.imj <- compo.imj[,order(colnames(compo.imj))]
      
      # identical(plots.data.imj$PlotName, rownames(compo.imj))
  
      #### CREATES A SPECIES DATA TABLE ##########################################
  
      species.limits.imj <- aggregate(trees.data.imj$MinimumElevationMeters, 
        by=list(trees.data.imj$SpeciesName), FUN=range, na.rm=TRUE)
  
      species.w.mean.imj <- aggregate(trees.data.imj$MinimumElevationMeters, 
        by=list(trees.data.imj$SpeciesName), FUN=mean, na.rm=TRUE)    
      
      species.midpoint.imj <- apply(species.limits.imj$x, 1, mean) 
      
      species.range.imj <- apply(species.limits.imj$x, 1, diff)     
  
      species.indN.imj <- aggregate(trees.data.imj$MinimumElevationMeters, 
        by=list(trees.data.imj$SpeciesName), FUN=function(x) sum(!is.na(x)))     
      
      
      species.data.imj <- groups.spp.imj[,c(2,1)]
      
      #check to make sure species align with data. specifically, that names match
      
      if(!identical(as.character(species.limits.imj$Group.1), species.data.imj$spp))
        stop("Error, species names don't coincide among objects")
  
      species.data.imj <- data.frame(species.data.imj, species.indN.imj$x, 
        species.limits.imj$x, species.range.imj, species.midpoint.imj, species.w.mean.imj$x)
        
      colnames(species.data.imj) <- c("species", "clade", "indN", 
        "min.elev", "max.elev", "range.elev", "midpoint.elev", "mean.w.elev")
      
  
      #### EMPIRICAL CLADE DATA ##################################################
  
      #clade.limits.imj <- aggregate(trees.data.imj$MinimumElevationMeters, 
      #  by=list(trees.data.imj$groups), FUN=range, na.rm=TRUE)    
  
      clade.limits.imj <- 
        aggregate(species.data.imj$midpoint.elev, 
        by=list(species.data.imj$clade), FUN=range, na.rm=TRUE)   
      
      clade.w.mean.imj <- aggregate(species.data.imj$midpoint.elev, 
        by=list(species.data.imj$clade), FUN=mean, na.rm=TRUE)
      
      clade.midpoint.imj <- apply(clade.limits.imj$x, 1, mean) 
      
      clade.range.imj <- apply(clade.limits.imj$x, 1, diff)     
      
      clade.sppN.imj <- aggregate(species.data.imj$midpoint.elev, 
        by=list(species.data.imj$clade), FUN=function(x) sum(!is.na(x)))     
      
      
      if(!identical(as.character(clade.limits.imj$Group.1), 
        as.character(clade.w.mean.imj$Group.1)))
        stop("Error 2, species names don't coincide among objects")
  
      clade.data.imj <- data.frame(clade.sppN.imj$x)
      rownames(clade.data.imj) <- clade.limits.imj$Group.1
      
      clade.data.imj <- data.frame(clade.data.imj, clade.limits.imj$x, 
        clade.range.imj, clade.midpoint.imj, clade.w.mean.imj$x)  
        
      colnames(clade.data.imj) <- c("sppN", "min.elev", "max.elev", 
        "range.elev", "midpoint.elev", "mean.w.elev")      
      
      
      #### EMPIRICAL SPECIES DECOMPOSITION #######################################
              
      decomp.res.imj <- decompose.my.beta(compo=compo.imj, sp.by.group=groups.spp.imj)
      
      BC.empirical <- as.matrix(decomp.res.imj$bray.curtis$BC)
      S.empirical <- as.matrix(decomp.res.imj$sorensen$S)
      
      
      #### RUN NULL MODEL ANALYSES ###############################################
      
      # Null model randomizations 
      rand.n <- 100
      #rand.n <- 3
  
  
      # Creates empty arrays for elev. ranges null values
      
      clade.data.null <- array(data=NA, 
        dim= c(nrow(clade.data.imj), ncol(clade.data.imj), rand.n+1))
      
      dimnames(clade.data.null) <- list(rownames(clade.data.imj), 
        colnames(clade.data.imj), c("empirical", paste("rand", 1:rand.n)))    
    
          
      # Creates empty arrays for decomposition null values
      
      BC.ag.null <- array(data=NA, 
        dim= c(nrow(BC.empirical), ncol(BC.empirical), rand.n+1))
        
      dimnames(BC.ag.null) <- 
        list(rownames(BC.empirical), rownames(BC.empirical), 
        c("empirical", paste("rand", 1:rand.n)))
      
    
      BC.wg.null <- BC.ag.null
      BC.ws.null <- BC.ag.null
      
      S.ag.null <- BC.ag.null
      S.wg.null <- BC.ag.null
      
      
      # Enters empirical values in the first position of the null arrays
  
      clade.data.null[,,1] <- as.matrix(clade.data.imj)
  
      BC.ag.null[,,1] <- as.matrix(decomp.res.imj$bray.curtis$BC.ag)
      BC.wg.null[,,1] <- as.matrix(decomp.res.imj$bray.curtis$BC.wg)
      BC.ws.null[,,1] <- as.matrix(decomp.res.imj$bray.curtis$BC.ws)
      
      S.ag.null[,,1] <- as.matrix(decomp.res.imj$sorensen$S.ag)
      S.wg.null[,,1] <- as.matrix(decomp.res.imj$sorensen$S.wg)
      
      
      ## Loop through iterations to use in for loop
      rand.groups.spp.imj <- groups.spp.imj
      rand.species.data.imj <- species.data.imj
      
      #this for loop does it all. randomizes species clade assignments, decomposes turnover, and makes stats for each randomization 
      for(n in 1:rand.n)
      {
        
        # Randomize species across groups
        rand.groups.spp.imj$groups <- sample(rand.groups.spp.imj$groups)
        rand.species.data.imj$clade <- rand.groups.spp.imj$groups
        
  
        # calculate clade stats for Null clade data
        null.clade.limits.imj <- 
          aggregate(rand.species.data.imj$midpoint.elev, 
          by=list(rand.species.data.imj$clade), FUN=range, na.rm=TRUE)   
        
        null.clade.w.mean.imj <- aggregate(rand.species.data.imj$midpoint.elev, 
          by=list(rand.species.data.imj$clade), FUN=mean, na.rm=TRUE)
        
        null.clade.midpoint.imj <- apply(null.clade.limits.imj$x, 1, mean) 
        
        null.clade.range.imj <- apply(null.clade.limits.imj$x, 1, diff)     
        
        null.clade.sppN.imj <- aggregate(rand.species.data.imj$midpoint.elev, 
          by=list(rand.species.data.imj$clade), FUN=function(x) sum(!is.na(x)))     
    
        null.clade.data.imj <- data.frame(null.clade.sppN.imj$x)
        rownames(null.clade.data.imj) <- null.clade.limits.imj$Group.1      
          
        null.clade.data.imj <- data.frame(null.clade.data.imj, null.clade.limits.imj$x, 
          null.clade.range.imj, null.clade.midpoint.imj, null.clade.w.mean.imj$x)  
               
        clade.data.null[,,n+1] <- as.matrix(null.clade.data.imj)      
        
  
        # Null decomposition of turnover
        
        null.decomp.res.imj <- decompose.my.beta(compo=compo.imj, 
          sp.by.group=rand.groups.spp.imj)
      
        BC.ag.null[,,n+1] <- as.matrix(null.decomp.res.imj$bray.curtis$BC.ag)
        BC.wg.null[,,n+1] <- as.matrix(null.decomp.res.imj$bray.curtis$BC.wg)
        BC.ws.null[,,n+1] <- as.matrix(null.decomp.res.imj$bray.curtis$BC.ws)
      
        S.ag.null[,,n+1] <- as.matrix(null.decomp.res.imj$sorensen$S.ag)
        S.wg.null[,,n+1] <- as.matrix(null.decomp.res.imj$sorensen$S.wg)
        
      }
  
      ## write all output filts 
      write.csv(x=plots.data.imj, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/PlotData.csv"))
      
      write.csv(x=compo.imj, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/CompoData.csv"))  
      
      ##also save compo as R object because its easier to load
      save(x=compo.imj, 
           file = paste0(project.directory,"04_Results/", 
                         dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/CompoData"))
      
      
      write.csv(x=species.data.imj, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/SpeciesData.csv"))
  
      write.csv(x=clade.data.imj, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/CladeData.csv"))
  
    
      write.csv(x=BC.empirical, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC_empirical.csv"))
  
      write.csv(x=S.empirical, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Sor_empirical.csv"))
  
  
      save(x=clade.data.null, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/clade.data.null"))    
      
          
      save(x=BC.ag.null, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC.ag.null"))
  
      save(x=BC.wg.null, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC.wg.null"))
      
      save(x=BC.ws.null, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC.ws.null"))
      
      
      save(x=S.ag.null, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/S.ag.null"))
      
      save(x=S.wg.null, 
        file = paste0(project.directory,"04_Results/", 
          dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/S.wg.null"))    
      
    }
  }  
}

