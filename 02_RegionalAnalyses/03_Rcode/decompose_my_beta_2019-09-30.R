
#rm(list=objects())


### FUNCTION DEFINITION ########################################################

## This function calculates components of the Bray-Curtis and Sorensen 
## measuremens of distance between communities (i.e. tunover). For Sorensen, 
## there are two components: (a) within and (b) among groups. For Bray-Curtis  
## there are three components: (a) within species, (b) among species but 
## within groups and (c) among species among groups. 

## Authors: J. Sebastian Tello & Ivan Jimenez


## ARGUMENTS ##

## compo - this is a composition table with sites in rows and species in 
##         columns

## sp.by.group - this is a data.frame with two columns: 'spp' and "groups'. The  
##               column 'spp' contains names that must match the column names  
##               in 'compo'. The column 'groups' indicates to what group 
##               (i.e clade) each species belongs.


## OUTPUT ##

## The output is a list containing another two lists, each of which has the 
## components of either Bray-Curtis or Sorensen distance metrics. 

## 'bray-curtis' - the first list, contains four elements:

##  -'BC': is a 'dist' object with the total Bray-Curtis distace. The other 
##          three elements are additive fractions of 'BC'. 
##  -'BC.ws': is the component of 'BC' owing to within-species variation in 
##             abundance. 
##  -'BC.wg': is the component of 'BC' owing to turnover of species WITHIN groups.   
##  -'BC.ag': is the component of 'BC' owing to turnover of species AMONG groups.


## 'sorensen' - the second list, contains three elements:

##  -'S': is a 'dist' object with the total Sorensen distace. 
##  -'S.wg': is the component of "S" owing to turnover of species WITHIN groups.
##  -'S.ag': is the component of "S" owing to turnover of species AMONG groups. 


decompose.my.beta <- function(compo, sp.by.group)
{
  require(vegan)
  require(R.utils)
  
  
  ## Makes sure this is a dataframe
  sp.by.group <- data.frame(sp.by.group, stringsAsFactors = FALSE)
  

  ## Checks that both input objects match
  if(!identical(sp.by.group$spp, colnames(compo)))
    stop("Your column names and species in 'sp.by.group' do not match")
  
  
  ## Calculates components W, A and B - See Legendre and Legendre 2012
  W <- as.matrix(designdist(compo, method = "J", terms = c("minimum"), abcd = F))
  A <- as.matrix(designdist(compo, method = "A", terms = c("minimum"), abcd = F))
  B <- as.matrix(designdist(compo, method = "B", terms = c("minimum"), abcd = F))

  
  ## Calculates components a, b and c - See Legendre and Legendre 2012
  ## a: number of shared species; b and c: unique species in sites i and k
  a <- as.matrix(designdist(compo, method = "a", terms = c("binary"), abcd = T))
  b <- as.matrix(designdist(compo, method = "b", terms = c("binary"), abcd = T))
  c <- as.matrix(designdist(compo, method = "c", terms = c("binary"), abcd = T) ) 

  
  ## Creates a composition table for groups
  colnames. <- colnames(compo)
  rownames. <- rownames(compo)
  compo <- apply(compo,2,as.numeric)
  colnames(compo) <- colnames.
  rownames(compo) <- rownames.
  
  compo.t <- t(compo)
  compo.groups.t <- aggregate(compo.t, by=list(sp.by.group$groups), FUN=sum)
  rownames(compo.groups.t) <- compo.groups.t$Group.1
  compo.groups.t <- compo.groups.t[,-1]
  compo.groups <- t(compo.groups.t)


  ## Makes lists of the names of species and groups
  spp.names <- colnames(compo)
  group.names <- colnames(compo.groups)
  
   
  ## Creates empty matrices that will be filled by the loop
  bc.dist <- vegdist(compo, method="bray")
  
  BC.ws <- as.matrix(bc.dist)
  BC.ws[] <- NA
  BC <- BC.ag <- BC.wg <- BC.ws
  S <- S.ag <- S.wg <- BC.ws
  
  
  pb <- txtProgressBar(min=0, max=(nrow(compo)-1), style=3)  
  for(i in 1:(nrow(compo)-1))
  {
    setTxtProgressBar(pb, i)
    
    
    ## Finds species and groups in site i
    spp.site.i <- spp.names[compo[i,] > 0]
    groups.site.i <- group.names[compo.groups[i,] > 0]
    
    for(k in (i+1):nrow(compo))
    {

      ## Finds species and groups in site k
      spp.site.k <- spp.names[compo[k,] > 0]
      groups.site.k <- group.names[compo.groups[k,] > 0]
      
      ## Finds shared and unique species in sites i and k
      shared.spp <- intersect(spp.site.i, spp.site.k)
      unique.spp.i <- setdiff(spp.site.i, spp.site.k)
      unique.spp.k <- setdiff(spp.site.k, spp.site.i)
      
      ## Finds shared and unique groups in sites i and k
      shared.groups <- intersect(groups.site.i, groups.site.k)
      unique.groups.i <- setdiff(groups.site.i, groups.site.k)
      unique.groups.k <- setdiff(groups.site.k, groups.site.i)
      
      ## Finds all species of groups shared between sites i and k
      spp.shared.groups <- 
        sp.by.group$spp[sp.by.group$groups %in% shared.groups]
      
      ## Finds unique species to sites i and k of shared groups
      unique.spp.shared.groups.i <- intersect(unique.spp.i, spp.shared.groups)
      unique.spp.shared.groups.k <- intersect(unique.spp.k, spp.shared.groups)
      
      ## Finds unique species of unique groups in site i
      unique.spp.unique.groups.i <- 
        sp.by.group$spp[sp.by.group$groups %in% unique.groups.i]
      
      ## Finds unique species of unique groups in site k
      unique.spp.unique.groups.k <- 
        sp.by.group$spp[sp.by.group$groups %in% unique.groups.k]
      

      ### SORENSEN DISTANCE CALCULATIONS ###
      
      ## Calculares the number of unique species in shared groups
      b.wg <- length(unique.spp.shared.groups.i)
      c.wg <- length(unique.spp.shared.groups.k)
      
      ## Calculares the number of unique species in unique groups
      b.ag <- b[k,i] - b.wg
      c.ag <- c[k,i] - c.wg
      
      ## Sorensen distance decomposition
      S.wg[k,i] <- (b.wg + c.wg)/(2 * a[k,i] + (b[k,i] + c[k,i]))
      S.ag[k,i] <- (b.ag + c.ag)/(2 * a[k,i] + (b[k,i] + c[k,i]))
      
      ## Sums components to the full Sorensen distance
      S[k,i] <- S.wg[k,i] + S.ag[k,i]
            

      ### BRAY-CURTIS DISTANCE CALCULATIONS ###
            
      ## Calculares the site-sums of shared species
      si <- sum(compo[i,shared.spp])
      sk <- sum(compo[k,shared.spp])

      ## Calculares the site-sums of unique species in shared groups
      ei.wg <- sum(compo[i,unique.spp.shared.groups.i])
      ek.wg <- sum(compo[k,unique.spp.shared.groups.k])

      ## Calculares the site-sums of unique species in unique groups
      ei.ag <- sum(compo[i,unique.spp.unique.groups.i])
      ek.ag <- sum(compo[k,unique.spp.unique.groups.k])

      ## Calculates the components of Bray-Curtis
      BC.ws[k,i] <- (si + sk - 2*W[k,i]) / (A[k,i] + B[k,i])
      BC.wg[k,i] <- (ei.wg + ek.wg) / (A[k,i] + B[k,i])
      BC.ag[k,i] <- (ei.ag + ek.ag) / (A[k,i] + B[k,i]) 

      ## Sums components to the full Bray-Curtis distance
      BC[k,i] <- BC.ws[k,i] + BC.wg[k,i] + BC.ag[k,i]
    }
  }
  close(pb)

  output <- list(
    bray.curtis=list(BC=as.dist(BC), BC.ws=as.dist(BC.ws), 
    BC.wg=as.dist(BC.wg), BC.ag=as.dist(BC.ag)), 
    sorensen=list(S=as.dist(S), S.wg=as.dist(S.wg), S.ag=as.dist(S.ag)))
    
  output
}



### EXAMPLE ####################################################################

if(FALSE)
{
  library(vegan)
  
  n.col <- 100
  n.row <- 50
 
  
  compo <- matrix(rpois(n=n.col*n.row, 1.25), ncol=n.col)
  colnames(compo) <- paste0("sp", 1:ncol(compo))  
  rownames(compo) <- paste0("site", 1:nrow(compo))  
  
  dim(compo)
   
  
  groups <- sample(paste0("clade_", LETTERS), size=ncol(compo), replace=TRUE)
  sp.by.group <- data.frame(groups, spp=colnames(compo), stringsAsFactors = FALSE)
  
  
  res2 <- decompose.my.beta(compo=compo, sp.by.group=sp.by.group)
  
  
  BrayCurtis <- vegdist(compo, method="bray")
  Sorensen <- vegdist(compo, method="bray", binary=TRUE)
  
  
  plot(res$sorensen$S ~ Sorensen, ylim=c(0,max(Sorensen)))
  points(res$sorensen$S.wg ~ Sorensen, col="blue")
  points(res$sorensen$S.ag ~ Sorensen, col="green")
  abline(c(0,1))
  
  
  plot(res$bray.curtis$BC ~ BrayCurtis, ylim=c(0, max(BrayCurtis)))
  points(res$bray.curtis$BC.ws ~ BrayCurtis, col="red")
  points(res$bray.curtis$BC.wg ~ BrayCurtis, col="blue")
  points(res$bray.curtis$BC.ag ~ BrayCurtis, col="green")
  abline(c(0,1))

}

  
 