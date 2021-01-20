
plot.my.beta.props <- function(ag.dist, wg.dist, ws.dist, x.dist, 
  ylab="Proportion of Turnover", xlab="Difference in Elevation (m)", 
  ylim=c(0,1), xlim, cex.lab=2, cex.axis=1.25, 
  pch=20, cex=1.25, cex.2=1.5, lwd=1.5, 
  smooth.option="lowess", smooth.factor=0.75, degree=2, x.n.breaks=15, 
  ag.col="lightgoldenrod", ag.line.col="orange", 
  wg.col="lightblue", wg.line.col="royalblue4", 
  ws.col="indianred1", ws.line.col="red3")
{


  use.wg <- TRUE
  if(missing(wg.dist))
  {  
    use.wg <- FALSE
    wg.dist <- NA
  }  
    
  use.ws <- TRUE
  if(missing(ws.dist))
  {  
    use.ws <- FALSE
    ws.dist <- NA
  }
  
  if(missing(ylim))
    ylim <- range(c(ag.dist, wg.dist, ws.dist), na.rm=TRUE)

  if(missing(xlim))
    xlim <- range(x.dist, na.rm=TRUE)
  
  if(any(is.na(x.dist)))
  {
    #which.na <- which(is.na(x.dist))
    
    which.na <- unique(c(which(is.na(x.dist)), which(is.na(ag.dist))))
    
    if(use.wg)
      which.na <- unique(c(which.na, which(is.na(wg.dist))))
      
    if(use.ws)
      which.na <- unique(c(which.na, which(is.na(ws.dist))))
    
    
    x.dist <- x.dist[-which.na]
    ag.dist <- ag.dist[-which.na]

    if(use.wg)
      wg.dist <- wg.dist[-which.na]
  
    if(use.ws)
      ws.dist <- ws.dist[-which.na]
  }
	x.dist.order <- order(x.dist)

	x.dist <- x.dist[x.dist.order]
	ag.dist <- ag.dist[x.dist.order]
	
  if(use.wg)
	  wg.dist <- wg.dist[x.dist.order]
	
  if(use.ws)
	  ws.dist <- ws.dist[x.dist.order]


	plot(ag.dist ~ x.dist, type="n", ylab=ylab, xlab=xlab, 
	  ylim=ylim, xlim=xlim, cex.lab=1.5, cex.axis=1.25, las=1)

	points(ag.dist ~ x.dist, pch=pch, col=ag.col, cex=cex, lwd=0.75)
	
  if(use.wg)
	  points(wg.dist ~ x.dist, pch=pch, col=wg.col, cex=cex, lwd=0.75)
	
  if(use.ws)
	  points(ws.dist ~ x.dist, pch=pch, col=ws.col, cex=cex, lwd=0.75)


	if(smooth.option=="lowess")
	{

  	points(y=predict(loess(ag.dist ~ x.dist, 
  	  span=smooth.factor, degree=degree)),
  	  x=as.numeric(x.dist), type="l",
  	  col = ag.line.col, lwd=lwd)

    if(use.wg)
    	points(y=predict(loess(as.numeric(wg.dist) ~ as.numeric(x.dist), 
    	  span=smooth.factor, degree=degree)),
    	  x=na.omit(as.numeric(x.dist)), type="l",
    	  col = wg.line.col, lwd=lwd)

    if(use.ws)
     	points(y=predict(loess(as.numeric(ws.dist) ~ as.numeric(x.dist), 
     	  span=smooth.factor, degree=degree)),
     	  x=na.omit(as.numeric(x.dist)), type="l",
    	  col = ws.line.col, lwd=lwd)
  }

	if(smooth.option=="lm")
	{
	  lm.ag.model <- lm(ag.dist ~ x.dist)
	  
	  abline(lm.ag.model, lwd=lwd, col=ag.line.col)
	  
	  
	  if(use.wg)
	    lm.wg.model <- lm(wg.dist ~ x.dist)
	  
	    abline(lm.wg.model, lwd=lwd, col=wg.line.col)
	  
	  # if(use.ws)
	  #   lm.ws.model <- lm(as.dist(ws.dist) ~ x.dist)
	  #   
	  #   abline(lm.ws.model, lwd=lwd, col=ws.line.col)
	}
	

	if(smooth.option=="group.means")
	{
		range.in.x <- range(x.dist, na.rm=TRUE)

    x.breaks <- seq(range.in.x[1], range.in.x[2], length.out=x.n.breaks)

	  x.groups <- cut(x.dist, x.breaks)

  	x.means <- tapply(x.dist, x.groups, mean)
  	ag.means <- tapply(ag.dist, x.groups, mean)
    if(use.wg)
  	  wg.means <- tapply(wg.dist, x.groups, mean)
    if(use.ws)
  	  ws.means <- tapply(ws.dist, x.groups, mean)
 
  	points(ag.means ~ x.means, type="l", col=ag.line.col, lwd=lwd)
    if(use.wg)
    	points(wg.means ~ x.means, type="l", col=wg.line.col, lwd=lwd)
    if(use.ws)
    	points(ws.means ~ x.means, type="l", col=ws.line.col, lwd=lwd)  	
  	 
  	points(ag.means ~ x.means, pch=pch, cex=cex.2, col=ag.line.col, lwd=lwd)
    if(use.wg)
  		points(wg.means ~ x.means, pch=pch, cex=cex.2, col=wg.line.col, lwd=lwd)
  	if(use.ws)
      points(ws.means ~ x.means, pch=pch, cex=cex.2, col=ws.line.col, lwd=lwd)
	}

}



## EXAMPLE ##
if(FALSE)
{
  library(scales)
  
  size <- 50
  
  x.example <- matrix(runif(n=size^2), nrow=size, ncol=size)
  
  ag.example <- matrix(1.5*x.example+rnorm(size^2, mean=0, sd=0.5), nrow=size, ncol=size)
  ag.example <- rescale(ag.example)
  wg.example <- 1-ag.example
  
  x.example <- as.dist(x.example)
  ag.example <- as.dist(ag.example)
  wg.example <- as.dist(wg.example)
  
  plot(ag.example ~ x.example)
  plot(wg.example ~ x.example)
  
  
  plot.my.beta.props(ag.dist=ag.example, wg.dist=wg.example, x.dist=x.example, 
  	ylab="Proportion of Turnover", xlab="GeoDists", ylim=c(0,1))
  

  x.example.nas <- as.matrix(x.example)
  x.example.nas[sample(1:length(x.example.nas), 500)] <- NA  
  x.example.nas <- as.dist(x.example.nas)
  sum(is.na(x.example.nas))

  plot.my.beta.props(ag.dist=ag.example, wg.dist=wg.example, x.dist=x.example.nas, 
  	ylab="Proportion of Turnover", xlab="GeoDists", ylim=c(0,1))
    
    
}

