
# the goal of this script is to produce publication ready figures of results 
# for 30 and 60 MY old clades as well as for Large and Small plot data.

rm(list = objects())

#you will need a lot of memory to run this
# memory.size(max=16235.29)
# memory.limit()
# ?memory.limit
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
library(plyr)
library(gtools)
library(geosphere)
library(grr)
source("03_Rcode/PlotMyBetaProps_2020-03-20.R")


datasets <- c("LargePlots", "SmallPlots")
morpho.opts <- c("WithMorphos", "NoMorphos")
cutoffs <- c(60, 30)

metric.to.use <- "S"
#metric.to.use <- "BC"




## make empty summary table that will be filled as each dataset loops (see end of loop for how the table is filled)
## note: expand.grid gives combinations of a set of objects.
summary.table <- expand.grid(datasets, morpho.opts, cutoffs)
summary.table <- data.frame(summary.table, Nspecies=NA, totalNplots=NA, Nplots.used=NA, Nclades=NA, N1sppclades=NA,
                            Meanclades=NA, p.val.ag.elev=NA, p.val.wg.elev=NA, p.val.ag.geo=NA, p.val.wg.geo=NA)


for(i in 1:length(datasets))
{
  dataset.i <- datasets[i]

  
  for(m in 1:length(morpho.opts))
  {   
    morpho.opt.m <- morpho.opts[m]  

    
    for(j in 1:length(cutoffs))
    {
      cutoff.j <- cutoffs[j]
      
      ##########################################################################
      #### 2. OPEN DATA AND SET VARIABLES ######################################
      ##########################################################################
      
      plots.data <- read.csv(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/PlotData.csv"), 
        header=TRUE, as.is=TRUE, sep = ",", fileEncoding="UTF-8")
      
      plots.data <- plots.data %>% arrange(PlotName)
      
      
      species.data <- read.csv(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/SpeciesData.csv"), 
        header=TRUE, as.is=TRUE, sep = ",", fileEncoding="UTF-8")
      
      
      clade.data <- read.csv(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/CladeData.csv"), 
        header=TRUE, as.is=TRUE, sep = ",", fileEncoding="UTF-8")
      
      load(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/clade.data.null"))
      
      
      BC.total.emp <- as.matrix(read.csv(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC_empirical.csv"), 
        header=TRUE, row.names=1, as.is=TRUE, sep = ",", fileEncoding="UTF-8"))
      
      load(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC.ag.null"))
      
      load(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC.wg.null"))
      
      #load(file=paste0(project.directory, "04_Results/", 
      #  dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/BC.ws.null"))
      
      
      S.total.emp <- as.matrix(read.csv(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Sor_empirical.csv"), 
        header=TRUE, row.names=1, as.is=TRUE, sep = ",", fileEncoding="UTF-8"))
      
      load(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/S.ag.null"))
      
      load(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/S.wg.null"))
      
      
      if(metric.to.use == "BC")
      {
        tovr.total.emp <- BC.total.emp
        tovr.ag.null <- BC.ag.null
        tovr.wg.null <- BC.wg.null
        #tovr.ws.null <- BC.ws.null
        
      }
      
      if(metric.to.use == "S")
      {
        tovr.total.emp <- S.total.emp
        tovr.ag.null <- S.ag.null
        tovr.wg.null <- S.wg.null
        #tovr.ws.null <- NA
        
      }
      
      rm(BC.total.emp, BC.ag.null, BC.wg.null, S.total.emp, S.ag.null, S.wg.null)
      
      
      ##########################################################################
      #### 3. CALCULATIONS OF TURNOVER #####################################################
      ##########################################################################
      	
      ### ELEVATIONAL AND GEOGRAPHICAL DISTANCES ###############################
      	
      if(dataset.i == "LargePlots")
      {  
      	elevdistmin <- 0 #in meters (used to minimize effect of elevational 
      	                 #distance when analyzing turnover across geographic distance)
      	
      	elevdistmax <- 200 #in meters (used to minimize effect of elevational 
      	                   #distance when analyzing turnover across geographic distance)
      	
      	geodistmin <- 50 #in km (used to minimize effect of geographic distance 
      	                 #when analyzing turnover across elevations)
      	
      	geodistmax <- 90 #in km (used to minimize effect of geographic distance 
      	                 #when analyzing turnover across elevations)
      }
      
      if(dataset.i == "SmallPlots")
      { 		
      	elevdistmin <- 0 #in meters
      	elevdistmax <- 200 #in meters
      	
      	geodistmin <- 110 #in km
      	geodistmax <- 160 #in km		
      }
      
      
      #extract elevations from plot data
      elevs <- plots.data$MinimumElevationMeters
      names(elevs) <- plots.data$PlotName
      	
      
      #create a distance matrix of elevations 
      elevdist <- dist(elevs)
      
      
      #create a matrix of geographic distances 
      coords <- plots.data[, c("long", "lat")] 
      colnames(coords) <- c("long", "lat")
      	
      geo.dists <- distm(x=coords, fun=distGeo)
      geo.dists <- geo.dists/1000
      colnames(geo.dists) <- plots.data$PlotName 
      rownames(geo.dists) <- plots.data$PlotName  
      geo.dists <- as.dist(geo.dists)    
      
   
      #create subsets of elevational and geographic distances that are
      #appropriate for our analyses
         	
      elevdists.to.use <- elevdist  # the elevations to use for selected subset of plots 
      elevdists.to.use[geo.dists < geodistmin] <- NA
      elevdists.to.use[geo.dists > geodistmax] <- NA
      
      identical(elevdist, elevdists.to.use)
      sum(is.na(elevdist))
      sum(is.na(elevdists.to.use))
      	
      geodists.to.use <- geo.dists
      geodists.to.use[elevdist < elevdistmin] <-NA
      geodists.to.use[elevdist > elevdistmax] <-NA
      
      
      ### TURNOVER COMPONENTS INTO PROPORTIONS OF TURNOVER #####################
      
      tovr.ag.null.prop <- tovr.ag.null
      for(k in 1:dim(tovr.ag.null)[3])
        tovr.ag.null.prop[, , k] <- tovr.ag.null[, , k]/tovr.total.emp
      
      tovr.wg.null.prop <- tovr.wg.null
      for(k in 1:dim(tovr.wg.null)[3])
        tovr.wg.null.prop[, , k] <- tovr.wg.null[, , k]/tovr.total.emp
  
        
      tovr.ag.null.prop[tovr.ag.null.prop < 0.00000001] <- 0.00000001
      tovr.ag.null.prop[tovr.ag.null.prop > 0.99999999] <- 0.99999999 
      	
      tovr.wg.null.prop[tovr.wg.null.prop < 0.00000001] <- 0.00000001
      tovr.wg.null.prop[tovr.wg.null.prop > 0.99999999] <- 0.99999999 
      
      ### WITHIN AND AMONG TURNOVER COMPONENTS INTO A DIFFERENCE ###############
      
      tovr.null.diff <- tovr.ag.null - tovr.wg.null
      dim(tovr.null.diff)
      
      tovr.emp.diff <- tovr.null.diff[, , 1]
      dim(tovr.emp.diff)
      
      #reduce memory usage 
      rm(tovr.ag.null, tovr.wg.null)
  
      ##########################################################################
      #### 4. LOGIT TRANSFORM AND lm SLOPE CALC AGAINST ELEVATION AND GEOGRAPHIC DISTANCES ################
      ##########################################################################
      
      ## slope comparison null vs. empirical--################
      #first transform my data into proportions for all matrices in array 
      #(i could not figure out how to do this as an array, so I made for loop 
      #into a list)
      
      #so the for loop is (for each matrix) first getting our data into 
      #proportions of bray curtis/sorenson second, tranforming the data to be more 
      #linear, then fitting a linear regression slope and saving them to a 
      #vector
      
      #so logit doesnt like 0s in the data and because these are randomizations 
      #a few (less than 5) matrices have some value with a zero turnover
      #to correct this I will replace them with 0.1 so that the logit 
      #tranformation does not produce inf. and regressions can be made 
      
      #make empty vectors
      slopes.ag.elev <- rep(NA, 1000)
      slopes.wg.elev <- slopes.ag.elev
      slopes.wg.geo <- slopes.ag.geo <- slopes.ag.elev
      slopes.diff.geo <- slopes.diff.elev <- slopes.ag.elev
      #make your logit props arrays
      logit.ag.props <- tovr.ag.null.prop
      logit.wg.props <- tovr.wg.null.prop
      
      for (k in 1:dim(tovr.ag.null.prop)[3])
      {
        #print(k)
        
        logit.ag.prop.k <- tovr.ag.null.prop[, , k]
        logit.wg.prop.k <- tovr.wg.null.prop[, , k]
        diff.k <- tovr.null.diff[, , k]
        #transform data using logit
        logit.ag.prop.k <- logit(as.dist(logit.ag.prop.k))
        logit.wg.prop.k <- logit(as.dist(logit.wg.prop.k))
        diff.k <- as.dist(diff.k)
        #put those logits into a growing array! 
        logit.ag.props[, , k] <- as.matrix(logit.ag.prop.k)
        logit.wg.props[, , k] <- as.matrix(logit.wg.prop.k)
        #linear model of null values, elevational distances
        ag.elev.lm.k <- lm(logit.ag.prop.k ~ scale(as.numeric(elevdists.to.use)))
        wg.elev.lm.k <- lm(logit.wg.prop.k ~ scale(as.numeric(elevdists.to.use)))
        diff.elev.lm.k <- lm(diff.k ~ scale(as.numeric(elevdists.to.use)))
        #select coefficients linear model for elevation
        slopes.ag.elev[k] <- ag.elev.lm.k$coefficients[2]
        slopes.wg.elev[k] <- wg.elev.lm.k$coefficients[2]
        slopes.diff.elev[k] <- diff.elev.lm.k$coefficients[2]
        #linear model of null values, geographic distances
        ag.geo.lm.k <- lm(logit.ag.prop.k ~ scale(as.numeric(geodists.to.use)))
        wg.geo.lm.k <- lm(logit.wg.prop.k ~ scale(as.numeric(geodists.to.use)))
        diff.geo.lm.k <- lm(diff.k ~ scale(as.numeric(geodists.to.use)))	  
        #select coefficients linear model 
        slopes.ag.geo[k] <- ag.geo.lm.k$coefficients[2]
        slopes.wg.geo[k] <- wg.geo.lm.k$coefficients[2]
        slopes.diff.geo[k] <- diff.geo.lm.k$coefficients[2]
      }
      
      #pick out the empirical logit prop
      tovr.logit.ag.emp.prop <- as.dist(logit.ag.props[, , 1])
      tovr.logit.wg.emp.prop <- as.dist(logit.wg.props[, , 1])
      
      
      ### EMPIRICAL raw PROPORTIONS ###
            
      tovr.ag.emp.prop <- as.dist(tovr.ag.null.prop[, , 1])
      tovr.wg.emp.prop <- as.dist(tovr.wg.null.prop[, , 1])
      #tovr.ws.emp.prop <- as.dist(tovr.ws.null.prop[, , 1])
      
     
     
      
      #############################################################################
      ### STANDARDIZED EFFECT SIZES ############################################
      ###############################################################################
      ##  find mean and SD then calculate SES (dont include the first matrix 
      ##  in array as it is empirical)
      mean.tovr.ag.null.prop <- 
        apply(tovr.ag.null.prop[, , -1], c(1, 2), mean, na.rm=TRUE)
      
      mean.tovr.wg.null.prop <- 
        apply(tovr.wg.null.prop[, , -1], c(1, 2), mean, na.rm=TRUE)
      
      mean.tovr.null.diff <- 
        apply(tovr.null.diff[, , -1], c(1, 2), mean, na.rm=TRUE)
      
      	
      sd.tovr.ag.null.prop <- 
        apply(tovr.ag.null.prop[, , -1], c(1, 2), sd, na.rm=TRUE)
      
      sd.tovr.wg.null.prop <- 
        apply(tovr.wg.null.prop[, , -1], c(1, 2), sd, na.rm=TRUE)
      	
      sd.tovr.null.diff <- 
        apply(tovr.null.diff[, , -1], c(1, 2), sd, na.rm=TRUE)
      
      ### do the same for LOGIT transformed data
      
      ###### calculate mean logit props for null model 
      mean.logit.tovr.ag.null.prop <- 
        apply(logit.ag.props[, , -1], c(1, 2), mean, na.rm=TRUE)
      
      mean.logit.tovr.wg.null.prop <- 
        apply(logit.wg.props[, , -1], c(1, 2), mean, na.rm=TRUE)
      
      sd.logit.tovr.ag.null.prop <- 
        apply(logit.ag.props[, , -1], c(1, 2), sd, na.rm=TRUE)
      
      sd.logit.tovr.wg.null.prop <- 
        apply(logit.wg.props[, , -1], c(1, 2), sd, na.rm=TRUE)
      
    
      
      #####################################
      #calculate standard effect size (ses)
      ses.ag.prop <- 
        (tovr.ag.null.prop[, , "empirical"] - mean.tovr.ag.null.prop) / sd.tovr.ag.null.prop
      
      ses.wg.prop <- 
        (tovr.wg.null.prop[, , "empirical"] - mean.tovr.wg.null.prop) / sd.tovr.wg.null.prop
      	
      ses.diff <- (tovr.emp.diff - mean.tovr.null.diff) / sd.tovr.null.diff
      
      ########################################
      ## do the same for LOGIT transformed data. 
      ses.logit.ag.prop <- 
        (logit.ag.props[, , "empirical"] - mean.logit.tovr.ag.null.prop) / sd.logit.tovr.ag.null.prop
      
      ses.logit.wg.prop <- 
        (logit.wg.props[, , "empirical"] - mean.logit.tovr.wg.null.prop) / sd.logit.tovr.wg.null.prop
      
      #for any pair of sites that share the same species composition we 
      #repalced NAs with zeros so that SES has no NA's 
      ses.ag.prop[sd.tovr.ag.null.prop==0] <- 0
      ses.wg.prop[sd.tovr.wg.null.prop==0] <- 0 
      ses.logit.ag.prop[sd.logit.tovr.ag.null.prop==0] <- 0
      ses.logit.wg.prop[sd.logit.tovr.wg.null.prop==0] <- 0 
      
      ### manage memory 
      rm(tovr.ag.null.prop, tovr.wg.null.prop, tovr.null.diff, logit.ag.props, logit.wg.props)

      ###NEW stuff###################################################################
      ###########################################################################


      #compare AIC of lm  (linear) vs gam (non-linear) models
      # library("MASS")
      # sresid <- studres(lin.model)
      # hist(sresid, freq=FALSE,
      #      main="Distribution of Studentized Residuals for ag.emp.prop")
      # xfit<-seq(min(sresid),max(sresid),length=40)
      # yfit<-dnorm(xfit)
      # lines(xfit, yfit)
      # now test AIC of linear regression vs general additive model (GAM; AKA loess)
      library(mgcv)
      
      ## raw proportions data 
      #null data 
      lm.model <- lm(as.dist(mean.tovr.ag.null.prop) ~ elevdists.to.use)
      gam.model <- gam(as.dist(mean.tovr.ag.null.prop) ~ s(scale(as.numeric(elevdists.to.use))))
      AIC(lm.model)
      AIC(gam.model)
      #empircal data 
      lm.model <- lm(as.dist(tovr.ag.emp.prop) ~ elevdists.to.use)
      gam.model<- gam(as.dist(tovr.ag.emp.prop) ~ s(scale(as.numeric(elevdists.to.use))))
      AIC(lm.model)
      AIC(gam.model)
      
      ##logit data
      #null data logit transformed
      lm.model.logit <- lm(as.dist(mean.logit.tovr.ag.null.prop) ~ elevdists.to.use)
      gam.model.logit <- gam(as.dist(mean.logit.tovr.ag.null.prop) ~ s(scale(as.numeric(elevdists.to.use))))
      AIC(lm.model.logit)
      AIC(gam.model.logit)
      #empircal data logit transformed
      lm.model.logit <- lm(as.dist(tovr.logit.ag.emp.prop) ~ elevdists.to.use)
      gam.model.logit <- gam(as.dist(tovr.logit.ag.emp.prop) ~ s(scale(as.numeric(elevdists.to.use))))
      AIC(lm.model.logit)
      AIC(gam.model.logit)

 
      
      ##########################################################################
      #### 5. plot lm of Logit transformed data  ################################################
      ##########################################################################

      ## Logit Transformed PROPORTION OF VARIATION AGAINS ELEVATIONAL DIFFERENCE ## 
      tiff(file=paste0(project.directory, "04_Results/", 
                       dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Figures/",
                       "LogitProportionOfTurnoverGradients.tiff"),  
           width = 45, height = 19, units = 'cm', res = 300, pointsize = 9)
      
      par(mfrow=c(2, 4), las=1, cex=1, 
          mar=c(4.5, 4.5, 4.1, 2.5), oma= c(1, 1, 1, 12.5))
      
 
      
      #plot empirical
      plot.my.beta.props(ag.dist=as.dist(tovr.logit.ag.emp.prop), 
                         wg.dist=as.dist(tovr.logit.wg.emp.prop), smooth.option="lm", x.dist=elevdists.to.use, 
                         ylab="Logit Proportion of Turnover", xlab="Difference in Elevation (m)", 
                         lwd=3, ylim=c(-5, 5))
      mtext(text="A", side=3, adj=0, font=2, cex=1.5, line=0.2)
      mtext(text="Empirical", side=3, font=2, cex=1.75, line=1.5)
      
      
      #plot null
      plot.my.beta.props(ag.dist=as.dist(mean.logit.tovr.ag.null.prop), 
                         wg.dist=as.dist(mean.logit.tovr.wg.null.prop), smooth.option="lm", x.dist=elevdists.to.use, 
                         ylab="Logit Proportion of Turnover", xlab="Difference in Elevation (m)", 
                         lwd=3, ylim=c(-5, 5))
      mtext(text="B", side=3, adj=0, font=2, cex=1.5, line=0.2)
      mtext(text="Null", side=3, font=2, cex=1.75, line=1.5)
      
      #plot SES
      plot.my.beta.props(ag.dist=as.dist(ses.logit.ag.prop), 
                         wg.dist=as.dist(ses.logit.wg.prop), smooth.option="lm", x.dist=elevdists.to.use, 
                         ylab="Standardized Effect Size", xlab="Difference in Elevation (m)", 
                         lwd=3)
      abline(h=0, col="gray45", lty=2.5)
      mtext(text="C", side=3, adj=0, font=2, cex=1.5, line=0.2)
      mtext(text="Deviations from Null", side=3, font=2, cex=1.75, line=1.5)
      
      
      # slope comparison
      par(cex.axis=1.5)
      
      slope.range <- range(c(slopes.ag.elev, slopes.wg.elev, 0), na.rm = TRUE)
      slope.range <- slope.range + max(slope.range, na.rm = TRUE)*c(-0.05, 0.05)
      
      boxplot(xaxs="i", slopes.ag.elev[-1], slopes.wg.elev[-1], 
              col=c("orange", "royalblue4"), names=c("", ""), 
              ylim= slope.range, ylab= "Slope", cex.lab=1.5, cex.axis=1.25, 
              cex.main=2, cex.sub=2)
      axis(1, at=1:2, c("Among\nClades", "Within\nClades"), line=2, lwd=0)
      
      points(1, slopes.ag.elev[1], pch=24, cex=2.5, 
             col="black", bg="orange", lwd=1)
      points(2, slopes.wg.elev[1], pch=24, cex=2.5, 
             col="black", bg="royalblue4", lwd=1)
      abline(h=0, col="gray45", lty=2)
      
      mtext(text="D", side=3, adj=0, font=2, cex=1.5, line=0.2)
      mtext(text="Slope Comparison", side=3, font=2, cex=1.75, line=1.5)
      
      text(x=1.25, y=slopes.ag.elev[1], 
           labels=expression(italic("empirical slope")), 
           cex=1.1, font=2, pos=4, offset = 1)
      arrows(x0=1.115, x1=1.3, y0=slopes.ag.elev[1], y1=slopes.ag.elev[1],
             length=0, lwd=1.5)
      
      text(x=1, y=mean(slopes.ag.elev[-1])*1.5, 
           labels=expression(italic("null distribution")), 
           cex=1.1, font=2, pos=3, offset = 1)
      
      
      # legend
      
      legend("topright", xpd=NA, inset=c(-0.7, 0), 
             legend=c("Among Clades", "Within Clades"), 
             col=c("orange", "royalblue4"), 
             lty=1, cex=1.5, lwd = 3, box.lty=0)
      
      
      ## LOGIT PROPORTION OF VARIATION AGAINST GEOGRAPHIC DISTANCE ##  
      
      #plot empirical
      plot.my.beta.props(ag.dist=as.dist(tovr.logit.ag.emp.prop), 
                         wg.dist=as.dist(tovr.logit.wg.emp.prop), smooth.option="lm", x.dist=geodists.to.use, 
                         ylab="Logit Proportion of Turnover", xlab="Geographic Distance (km)", 
                         lwd=3, ylim=c(-5, 5))
      mtext(text="E", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
      
      #plot null
      plot.my.beta.props(ag.dist=as.dist(mean.logit.tovr.ag.null.prop), 
                         wg.dist=as.dist(mean.logit.tovr.wg.null.prop), smooth.option="lm", x.dist=geodists.to.use, 
                         ylab="Logit Proportion of Turnover", xlab="Geographic Distance (km)", 
                         lwd=3, ylim=c(-5, 5))
      mtext(text="F", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
      #plot SES
      plot.my.beta.props(ag.dist=as.dist(ses.logit.ag.prop), 
                         wg.dist=as.dist(ses.logit.wg.prop), smooth.option="lm", x.dist=geodists.to.use, 
                         ylab="Standardized Effect Size", xlab="Geographic Distance (km)", lwd=3)
      abline(h=0, col="gray45", lty=2.5)
      mtext(text="G", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
      
      # slope comparison
      par(cex.axis=1.5)
      
      slope.range <- range(c(slopes.ag.geo, slopes.wg.geo, 0), na.rm = TRUE)
      slope.range <- slope.range + max(slope.range, na.rm = TRUE)*c(-0.05, 0.05)
      
      boxplot(xaxs="i", slopes.ag.geo[-1], slopes.wg.geo[-1], 
              col=c("orange", "royalblue4"), names=c("", ""), 
              ylim= slope.range, ylab= "Slope", cex.lab=1.5, cex.axis=1.25, 
              cex.main=2, cex.sub=2)
      axis(1, at=1:2, c("Among\nClades", "Within\nClades"), line=2, lwd=0)
      
      points(1, slopes.ag.geo[1], pch=24, cex=2.5, 
             col="black", bg="orange", lwd=1)
      points(2, slopes.wg.geo[1], pch=24, cex=2.5, 
             col="black", bg="royalblue4", lwd=1)
      abline(h=0, col="gray45", lty=2)
      
      mtext(text="H", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
      dev.off()
      
      
      
      
      ######################################################################
      ###calculate p-values, if empirical is sigificantly different from null (based on linear regressions)

      slopes.ag.elev.center <- 
        abs(slopes.ag.elev - mean(slopes.ag.elev, na.rm = TRUE))
      
      slopes.wg.elev.center <- 
        abs(slopes.wg.elev - mean(slopes.wg.elev, na.rm = TRUE))
      	
      slopes.ag.geo.center <- 
        abs(slopes.ag.geo - mean(slopes.ag.geo, na.rm = TRUE))
      
      slopes.wg.geo.center <- 
        abs(slopes.wg.geo- mean(slopes.wg.geo, na.rm = TRUE))
      
      ### these objects will be used in the summary table generation at end of loop      	
      p.val.ag.elev <- 
        sum(slopes.ag.elev.center >= slopes.ag.elev.center[1]) / 
        length(na.omit(slopes.ag.elev.center))
      
      p.val.wg.elev <- 
        sum(slopes.wg.elev.center >= slopes.wg.elev.center[1]) / 
        length(na.omit(slopes.wg.elev.center))
      
      p.val.ag.geo <- 
        sum(slopes.ag.geo.center >= slopes.ag.geo.center[1]) / 
        length(na.omit(slopes.ag.geo.center))
      
      p.val.wg.geo <- 
        sum(slopes.wg.geo.center >= slopes.wg.geo.center[1]) / 
        length(na.omit(slopes.wg.geo.center))
      	
      	
  
      ##########################################################################
      #### 5. DESCRIPTIVE PLOTS ################################################
      ##########################################################################
      	
      basic.col <- "grey30" 	
      	
      	
      ### GEOGRAPHIC DISTANCES VS. ELEVATIONAL DIFFERENCES #####################
      ## output D : PLOT 1x3 panel of plots used to show the subsets of plots  
      ## used for analyses
      	
      tiff(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Figures/", 
        "GeoVsElevDists.tiff"),  
        width = 35, height = 12, units = 'cm', res = 300, pointsize = 10)
      	
        par(mfrow=c(1, 3), mgp=c(3.5, 1, 0), 
          mar=c(4.5, 5, 4.1, 2.5), oma=c(2, 2, 1, 1))
        par(cex.lab=1.5, cex.axis=1.25, las=1, cex=1)
        
          plot(elevdist ~ geo.dists, pch=21, cex=1.25, 
            bg=basic.col, col="white", lwd=0.75, 
            xlab="Geographic Distance (km)", ylab="Difference in Elevation (m)")
          polygon(x=c(0, 0, 500, 500), 
            y=c(elevdistmin, elevdistmax, elevdistmax, elevdistmin), 
            col=adjustcolor("darkolivegreen4", alpha.f = 0.5), border = NA)
          polygon(x=c(geodistmin, geodistmin, geodistmax, geodistmax), 
            y=c(0, 3500, 3500, 0), col=adjustcolor("darkorange3", alpha.f = 0.5), 
            border = NA)
          mtext(text="A", side=3, adj=0, font=2, cex=1.75, line=0.2)
          
          
          hist(elevdist, breaks = seq(0, 4250, 250), 
            xlab="Difference in Elevation (m)", ylab="", main=NA)
          hist(elevdist[geo.dists>geodistmin & geo.dists<geodistmax], add=TRUE, 
            col="darkorange3", breaks = seq(0, 4250, 250))
          box()
          mtext(text="B", side=3, adj=0, font=2, cex=1.75, line=0.2)
          mtext(text="Frequency", side=2, cex=1.75, las=3, line=4.75)
          ?mtext
          hist(geo.dists, breaks = seq(0, 500, 25), 
            xlab="Geographic Distance (km)", ylab="", main=NA)
          hist(geo.dists[elevdist>elevdistmin & elevdist<elevdistmax], add=TRUE, 
            col="darkolivegreen4", breaks = seq(0, 500, 25))
          box()
          mtext(text="C", side=3, adj=0, font=2, cex=1.75, line=0.2)
          mtext(text="Frequency", side=2, cex=1.75, las=3, line=4.75)
          
      dev.off()
          	
      	
      ### NUMBER OF SPECIES PER CLADE ##########################################
      
      tiff(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Figures/", 
        "SppNByClade.tiff"),  
        width = 13, height = 12, units = 'cm', res = 300, pointsize = 10)
      	
        par(mgp=c(3.25, 1, 0), mar=c(4.5, 5, 4.1, 2.5), oma=c(2, 2, 1, 1))
        par(cex.lab=1.5, cex.axis=1.25, las=1, cex=1)
        
        hist(clade.data$sppN, main=NA, xlab="Number of Species", 
          xlim=c(0,220), col=basic.col)
      
        t1 <- paste0("Number of clades = ", length(clade.data$sppN))
        t2 <- paste0("Clades with a single species = ", sum(clade.data$sppN==1))
        t3 <- paste0("Mean species per clade = ", round(mean(clade.data$sppN), 2))

        text(y=50, x=50, pos=4, cex=1, 
          labels = paste0(t1, "\n", t2, "\n", t3))      
                  
        box()
      
      dev.off()
      
      ### ELEVATIONAL DISTRIBUTION OF CLADES ###################################
      
      null.clade.range <- apply(clade.data.null, c(1,2), mean)[,"range.elev"]
      
      low.ci.null.clade.range <- apply(clade.data.null, c(1,2), quantile, 
        probs = 0.025)[,"range.elev"]
      high.ci.null.clade.range <- apply(clade.data.null, c(1,2), quantile, 
        probs = 0.975)[,"range.elev"]
      
      sd.null.clade.range <- apply(clade.data.null, c(1,2), sd)[,"range.elev"]
      
      ses.clade.range <- 
        (clade.data.null[,"range.elev",1] - null.clade.range) / 
        sd.null.clade.range
      
      
      tiff(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Figures/",
        "CladeRangeSizes.tiff"),  
        width = 22, height = 11, units = 'cm', res = 300, pointsize = 9)
      	
        par(mfrow=c(1,2), mgp=c(2.5, 1, 0), mar=c(4.5, 5, 4.1, 2.5), 
          oma=c(2, 2, 1, 1))
        par(cex.lab=1.5, cex.axis=1.25, las=1, cex=1)
        
          plot((null.clade.range+1) ~ log(clade.data$sppN), pch=20,
            ylim=(range(c(low.ci.null.clade.range+1, high.ci.null.clade.range+1), na.rm = TRUE)),
            xlim=c(log(2), log(max(clade.data$sppN, na.rm = TRUE))),
            xlab="Number of Species (log.)", ylab="", cex=1.25)
          title(ylab="Clade Range Size (m)", line=3.85)
          
          for(k in 1:length(clade.data$sppN))
            arrows(x0=log(clade.data$sppN[k]), x1=log(clade.data$sppN[k]), 
              y0=(low.ci.null.clade.range[k]+1), y1=(high.ci.null.clade.range[k]+1),
              code=3, angle=90, length=0, col="grey30")
          
          points((clade.data$range.elev+1) ~ log(clade.data$sppN), 
            pch=20, cex=1.75, col="orange", lwd=0.75)
          
          mtext(text="A", side=3, adj=0, font=2, cex=1.75, line=0.2)
          
          
          plot(ses.clade.range[clade.data$sppN>1] ~ 
            log(clade.data$sppN[clade.data$sppN>1]), 
            ylim=range(c(max(abs(ses.clade.range[clade.data$sppN>1]), na.rm=TRUE)*c(-1,1), -2, 2)),
            pch=20, cex=1.75, col="orange", lwd=0.75,
            xlab="Number of Species (log.)", ylab="SES of Range Size")
          
          abline(h=0, lwd=2)
          abline(h=-1.96, lwd=1.5, lty=2)
          abline(h=1.96, lwd=1.5, lty=2)
      
          mtext(text="B", side=3, adj=0, font=2, cex=1.75, line=0.2)
          
      dev.off()
      
      
  
      ##########################################################################
      #### 6. MAIN PLOTS #######################################################
      ##########################################################################
              
      ### TOTAL TURNOVER ALONG ELEVATIONAL AND GEOGRAPHIC GRADIENTS ############
      
      tiff(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Figures/",
        "TurnoverGradients.tiff"),  
        width = 11, height = 19, units = 'cm', res = 300, pointsize = 9)
      	
        par(mfrow=c(2,1), mgp=c(3, 1, 0), mar=c(4.5, 5, 4.1, 2.5), 
          oma=c(2, 2, 1, 1))
        par(cex.lab=1.5, cex.axis=1.25, las=1, cex=1)
          
          plot(as.dist(tovr.total.emp) ~ elevdist, type="n",
            ylab="Species Turnover", xlab="Difference in Elevation (m)", 
            ylim=c(0, 1), xlim=c(0, max(elevdist, na.rm = TRUE)))
          
          which.na <- which(is.na(elevdists.to.use))
          points(as.dist(tovr.total.emp)[-which.na] ~ elevdists.to.use[-which.na], 
            pch=21, cex=1.75, bg=basic.col, col="white", lwd=0.75)
          
          ordered.elev <- c(order(elevdists.to.use[-which.na]))
          points(y=predict(loess(as.dist(tovr.total.emp)[-which.na][ordered.elev] ~ 
            as.numeric(elevdists.to.use)[-which.na][ordered.elev])), 
              x=as.numeric(elevdists.to.use)[-which.na][ordered.elev], 
            type="l", col = "darkorange3", lwd=2.5)
      
          mtext(text="A", side=3, adj=0, font=2, cex=1.75, line=0.2)
              
          
          plot(as.dist(tovr.total.emp) ~ geo.dists, type="n",
            ylab="Species Turnover", xlab="Geographic Distance (km)", 
            ylim=c(0, 1), xlim=c(0, max(geo.dists, na.rm = TRUE)))
          
          which.na <- which(is.na(geodists.to.use))
          points(as.dist(tovr.total.emp)[-which.na] ~ geodists.to.use[-which.na], 
            pch=21, cex=1.75, bg=basic.col, col="white", lwd=0.75)
          
          ordered.geo <- c(order(geodists.to.use[-which.na]))
          points(y=predict(loess(as.dist(tovr.total.emp)[-which.na][ordered.geo] ~ 
            as.numeric(geodists.to.use)[-which.na][ordered.geo])), 
              x=as.numeric(geodists.to.use)[-which.na][ordered.geo], 
            type="l", col = "darkorange3", lwd=2.5)
      
          mtext(text="B", side=3, adj=0, font=2, cex=1.75, line=0.2)
              
      dev.off()    
      
      
      tiff(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Figures/",
        "ProportionOfTurnoverGradients.tiff"),  
        width = 45, height = 19, units = 'cm', res = 300, pointsize = 9)
      	
        par(mfrow=c(2, 4), las=1, cex=1, 
          mar=c(4.5, 4.5, 4.1, 2.5), oma= c(1, 1, 1, 12.5))
          
       
        ## PROPORTION OF VARIATION AGAINS ELEVATIONAL DIFFERENCE ##  
         
          #plot empirical
          plot.my.beta.props(ag.dist=as.dist(tovr.ag.emp.prop), 
            wg.dist=as.dist(tovr.wg.emp.prop), x.dist=elevdists.to.use, 
            ylab="Proportion of Turnover", xlab="Difference in Elevation (m)", 
            lwd=3, ylim=c(0, 1))
          mtext(text="A", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Empirical", side=3, font=2, cex=1.75, line=1.5)
          
          
          #plot null
          plot.my.beta.props(ag.dist=as.dist(mean.tovr.ag.null.prop), 
            wg.dist=as.dist(mean.tovr.wg.null.prop), x.dist=elevdists.to.use, 
            ylab="Proportion of Turnover", xlab="Difference in Elevation (m)", 
            lwd=3, ylim=c(0, 1))
          mtext(text="B", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Null", side=3, font=2, cex=1.75, line=1.5)
          
          #plot SES
          plot.my.beta.props(ag.dist=as.dist(ses.ag.prop), 
            wg.dist=as.dist(ses.wg.prop), x.dist=elevdists.to.use, 
            ylab="Standardized Effect Size", xlab="Difference in Elevation (m)", 
            lwd=3)
          abline(h=0, col="gray45", lty=2.5)
          mtext(text="C", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Deviations from Null", side=3, font=2, cex=1.75, line=1.5)
          
          
          # slope comparison
          par(cex.axis=1.5)
          
          slope.range <- range(c(slopes.ag.elev, slopes.wg.elev, 0), na.rm = TRUE)
          slope.range <- slope.range + max(slope.range, na.rm = TRUE)*c(-0.05, 0.05)
          
          boxplot(xaxs="i", slopes.ag.elev[-1], slopes.wg.elev[-1], 
            col=c("orange", "royalblue4"), names=c("", ""), 
            ylim= slope.range, ylab= "Slope", cex.lab=1.5, cex.axis=1.25, 
            cex.main=2, cex.sub=2)
          axis(1, at=1:2, c("Among\nClades", "Within\nClades"), line=2, lwd=0)
      
          points(1, slopes.ag.elev[1], pch=24, cex=2.5, 
            col="black", bg="orange", lwd=1)
          points(2, slopes.wg.elev[1], pch=24, cex=2.5, 
            col="black", bg="royalblue4", lwd=1)
          abline(h=0, col="gray45", lty=2)
          
          mtext(text="D", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Slope Comparison", side=3, font=2, cex=1.75, line=1.5)
          
          text(x=1.25, y=slopes.ag.elev[1], 
            labels=expression(italic("empirical slope")), 
            cex=1.1, font=2, pos=4, offset = 1)
          arrows(x0=1.115, x1=1.3, y0=slopes.ag.elev[1], y1=slopes.ag.elev[1],
            length=0, lwd=1.5)
          
          text(x=1, y=mean(slopes.ag.elev[-1])*1.5, 
            labels=expression(italic("null distribution")), 
            cex=1.1, font=2, pos=3, offset = 1)

          ###### legend
              
          legend("topright", xpd=NA, inset=c(-0.7, 0), 
            legend=c("Among Clades", "Within Clades"), 
            col=c("orange", "royalblue4"), 
            lty=1, cex=1.5, lwd = 3, box.lty=0)
      
          
        ## PROPORTION OF VARIATION AGAINST GEOGRAPHIC DISTANCE ##  
         
          #plot empirical
          plot.my.beta.props(ag.dist=as.dist(tovr.ag.emp.prop), 
            wg.dist=as.dist(tovr.wg.emp.prop), x.dist=geodists.to.use, 
            ylab="Proportion of Turnover", xlab="Geographic Distance (km)", 
            lwd=3, ylim=c(0, 1))
          mtext(text="E", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
          
          #plot null
          plot.my.beta.props(ag.dist=as.dist(mean.tovr.ag.null.prop), 
            wg.dist=as.dist(mean.tovr.wg.null.prop), x.dist=geodists.to.use, 
            ylab="Proportion of Turnover", xlab="Geographic Distance (km)", 
            lwd=3, ylim=c(0, 1))
          mtext(text="F", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
          #plot SES
          plot.my.beta.props(ag.dist=as.dist(ses.ag.prop), 
            wg.dist=as.dist(ses.wg.prop), x.dist=geodists.to.use, 
            ylab="Standardized Effect Size", xlab="Geographic Distance (km)", lwd=3)
          abline(h=0, col="gray45", lty=2.5)
          mtext(text="G", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
          
          # slope comparison
          par(cex.axis=1.5)
          
          slope.range <- range(c(slopes.ag.geo, slopes.wg.geo, 0), na.rm = TRUE)
          slope.range <- slope.range + max(slope.range, na.rm = TRUE)*c(-0.05, 0.05)
          
          boxplot(xaxs="i", slopes.ag.geo[-1], slopes.wg.geo[-1], 
            col=c("orange", "royalblue4"), names=c("", ""), 
            ylim= slope.range, ylab= "Slope", cex.lab=1.5, cex.axis=1.25, 
            cex.main=2, cex.sub=2)
          axis(1, at=1:2, c("Among\nClades", "Within\nClades"), line=2, lwd=0)
      
          points(1, slopes.ag.geo[1], pch=24, cex=2.5, 
            col="black", bg="orange", lwd=1)
          points(2, slopes.wg.geo[1], pch=24, cex=2.5, 
            col="black", bg="royalblue4", lwd=1)
          abline(h=0, col="gray45", lty=2)
          
          mtext(text="H", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
      dev.off()
      
      
      #graph the amoung-within group difference in proportions
      tiff(file=paste0(project.directory, "04_Results/", 
        dataset.i, "/", morpho.opt.m, "/", cutoff.j, "my/Figures/",
        "DifferenceBetweenTurnoverComponentsGradients.tiff"),  
        width = 45, height = 19, units = 'cm', res = 300, pointsize = 9)
      	
        par(mfrow=c(2, 4), las=1, cex=1, 
          mar=c(4.5, 4.5, 4.1, 2.5), oma= c(1, 1, 1, 12.5))
          
       
        ## PROPORTION OF VARIATION AGAINS ELEVATIONAL DIFFERENCE ##  
         
          #plot empirical
          plot.my.beta.props(ag.dist=as.dist(tovr.emp.diff), 
            x.dist=elevdists.to.use,
            ylab="Among - Within Difference", xlab="Difference in Elevation (m)", 
            lwd=3, ylim=c(-1, 1), ag.col="grey50", ag.line.col="grey20")
          mtext(text="A", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Empirical", side=3, font=2, cex=1.75, line=1.5)
          
          
          #plot null
          plot.my.beta.props(ag.dist=as.dist(mean.tovr.null.diff), 
            x.dist=elevdists.to.use,
            ylab="Among - Within Difference", xlab="Difference in Elevation (m)", 
            lwd=3, ylim=c(-1, 1), ag.col="grey50", ag.line.col="grey20")
          mtext(text="B", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Null", side=3, font=2, cex=1.75, line=1.5)
          
          #plot SES
          plot.my.beta.props(ag.dist=as.dist(ses.diff), 
            x.dist=elevdists.to.use,
            ylab="Standardized Effect Size", xlab="Difference in Elevation (m)", 
            lwd=3, ag.col="grey50", ag.line.col="grey20")
          abline(h=0, col="gray45", lty=2.5)
          mtext(text="C", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Deviations from Null", side=3, font=2, cex=1.75, line=1.5)
          
          
          # slope comparison
          par(cex.axis=1.5)
          
          slope.range <- range(c(slopes.diff.elev, 0), na.rm = TRUE)
          slope.range <- slope.range + max(slope.range, na.rm = TRUE)*c(-0.05, 0.05)
          
          hist(slopes.diff.elev, col="grey50", breaks=50, 
            xlim= slope.range, xlab= "Slope", main="", 
            cex.lab=1.5, cex.axis=1.25)
          abline(v=slopes.diff.elev[1], lwd=3)
          box()
      
          mtext(text="D", side=3, adj=0, font=2, cex=1.5, line=0.2)
          mtext(text="Slope Comparison", side=3, font=2, cex=1.75, line=1.5)
          
          text(x=1.25, y=slopes.ag.elev[1], 
            labels=expression(italic("empirical slope")), cex=1.1, font=2, 
            pos=4, offset = 1)
          arrows(x0=1.115, x1=1.3, y0=slopes.ag.elev[1], y1=slopes.ag.elev[1],
            length=0, lwd=1.5)
          
          text(x=1, y=mean(slopes.ag.elev[-1])*1.5, 
            labels=expression(italic("null distribution")), 
            cex=1.1, font=2, pos=3, offset = 1)
          #arrows(x0=1, x1=1, 
          #  y0=mean(slopes.ag.elev[-1]), y1=mean(slopes.ag.elev[-1])*-1,
          #  length=0, lwd=1.5)
      
      
        ## PROPORTION OF VARIATION AGAINST GEOGRAPHIC DISTANCE ##  
         
          #plot empirical
          plot.my.beta.props(ag.dist=as.dist(tovr.emp.diff), 
            x.dist=geodists.to.use, 
            ylab="Among - Within Difference", xlab="Geographic Distance (km)", 
            lwd=3, ylim=c(-1, 1), ag.col="grey50", ag.line.col="grey20")
          mtext(text="A", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
          
          #plot null
          plot.my.beta.props(ag.dist=as.dist(mean.tovr.null.diff), 
            x.dist=geodists.to.use, 
            ylab="Among - Within Difference", xlab="Geographic Distance (km)", 
            lwd=3, ylim=c(-1, 1), ag.col="grey50", ag.line.col="grey20")
          mtext(text="B", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
          #plot SES
          plot.my.beta.props(ag.dist=as.dist(ses.diff), 
            x.dist=geodists.to.use, 
            ylab="Standardized Effect Size", xlab="Geographic Distance (km)", 
            lwd=3, ag.col="grey50", ag.line.col="grey20")
          abline(h=0, col="gray45", lty=2.5)
          mtext(text="C", side=3, adj=0, font=2, cex=1.5, line=0.2)
      
          
          # slope comparison
          a <- as.dist(tovr.emp.diff)[!is.na(geodists.to.use)]
          b <- as.dist(mean.tovr.null.diff)[!is.na(geodists.to.use)]
          
          boxplot(list(a, b), 
            col=c("grey50", "grey50"), names=c("", ""), 
            ylim= range(c(a,b), na.rm = TRUE), ylab= "Among - Within Difference", 
            cex.lab=1.5, cex.axis=1.25, cex.main=2, cex.sub=2)
          axis(1, at=1:2, c("Empirical", "Mean of Null"), line=2, lwd=0)
      
          mtext(text="H", side=3, adj=0, font=2, cex=1.5, line=0.2)    
          
      dev.off()
      
      #########################################
      ####SUMMARIZE RESULTS ACROSS DATASETS #####
      #########################################
      
      # make basic object to use in summary table at end of loop.
      # note: totalNplots and Nplots.used are the number of plot X plot comparisons that were used to calculate turnover etc.
      Nspecies <- length(species.data$species)
      totalNplots <- length(plots.data$PlotName)
      Nplots.used <- length(na.omit(elevdists.to.use))
      Nclades <- length(clade.data$sppN)
      N1sppclades <- sum(clade.data$sppN==1)
      Meanclades <- round(mean(clade.data$sppN), 2)
      ## build a table of summary metrics, 
      which.imj <- which(summary.table$Var1==dataset.i & 
          summary.table$Var2==morpho.opt.m & summary.table$Var3==cutoff.j)
      
      summary.table[which.imj,] <- c(dataset.i, morpho.opt.m, cutoff.j, Nspecies, totalNplots, Nplots.used, Nclades, N1sppclades, Meanclades, p.val.ag.elev, p.val.wg.elev, p.val.ag.geo, p.val.wg.geo)
  
    }
  }  
}



## write table to file 

write.table(summary.table, "04_Results/Results_summary.txt", row.names = FALSE, sep = "\t", fileEncoding="UTF-8")
