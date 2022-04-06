# DO NOT USE THIS FILE FOR THE LAB!!!
# Keep in your working directory. 

if (!require("RColorBrewer")){install.packages("RColorBrewer");library(RColorBrewer)}
if (!require("forecast")){install.packages("forecast");library(forecast)}


cmav <- function(y,ma=NULL,fill=c(TRUE,FALSE),outplot=c(FALSE,TRUE),fast=c(TRUE,FALSE)){
  # Calculate centred moving average
  #
  # Inputs:
  #   y             Time series vector (can be ts object)
  #   ma            Length of centred moving average. If y is a ts object 
  #                 then the default is its frequency
  #   fill          If TRUE then fill first and last ma/2 observations using ETS
  #   outplot       If TRUE then produce plot
  #   fast          If TRUE then only a limited set of models are evaluated for CMA extrapolation
  #
  # Outputs:
  #   cma           Centred moving average. If y is a ts object, then cma has the
  #                 same properties
  #
  # Example:
  #   cmav(wineind,outplot=TRUE)
  #
  # Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Set defaults
  fill <- fill[1]
  outplot <- outplot[1]
  fast <- fast[1]
  
  # Get MA length
  if (is.null(ma)){
    if (any(class(y) == "ts")){
      ma <- frequency(y)
    } else {
      stop("MA length not defined (y not ts object).")
    }
  }
  
  # Get bounds for MA and correct length
  n <- length(y)
  mlbounds <- c(floor(ma/2)+1, n-floor(ma/2))
  isodd = ma %% 2 
  if (isodd == 0){
    ml <- ma+1
  } else {
    ml <- ma
  }
  
  # Calculate MA
  # Loop across MA-order to speed up things
  mamat <- matrix(NA,nrow=ml,ncol=(n-ml+1))
  for (i in 1:ml){
    mamat[i,] <- y[i:(n-ml+i)]
  }
  if (isodd == 0){
    mamat[c(1,ml),] <- mamat[c(1,ml),]/2
  }
  mamat <- colSums(mamat)/ma
  
  cma <- y
  cma[] <- NA
  cma[mlbounds[1]:mlbounds[2]] <- mamat
  
  # Fill MA is requested
  if (fill == TRUE){
    if (fast == FALSE){
      if ((n-mlbounds[2]) >= 1){   
        cma[(mlbounds[2]+1):n] <- as.vector(forecast::forecast(forecast::ets(cma[(mlbounds[1]:mlbounds[2])], 
                                                                             model="ZZN"),h=(n-mlbounds[2]))$mean)
      }
      if ((mlbounds[1]-1) >= 1){   
        cma[1:(mlbounds[1]-1)] <- rev(as.vector(forecast::forecast(forecast::ets(rev(cma[(mlbounds[1]:mlbounds[2])]),
                                                                                 model="ZZN"),h=(mlbounds[1]-1))$mean))
      }
    } else {
      fit <- forecast::ets(cma[(mlbounds[1]:mlbounds[2])],model="AZN")
      if ((n-mlbounds[2]) >= 1){   
        cma[(mlbounds[2]+1):n] <- as.vector(forecast::forecast(fit,h=(n-mlbounds[2]))$mean)
      }
      if ((mlbounds[1]-1) >= 1){   
        cma[1:(mlbounds[1]-1)] <- rev(as.vector(forecast::forecast(forecast::ets(rev(cma[mlbounds[1]:mlbounds[2]]),
                                                                                 fit,use.initial.values=FALSE),h=(mlbounds[1]-1))$mean))
      }
    }
  }
  
  if (outplot == TRUE){
    plot(y)
    lines(cma,col="red")
  }
  
  return(cma)
  
}

coxstuart <- function(y,type=c("trend","deviation","dispersion"),alpha=0.05){
  # Perform Cox - Stuart test for location and dispersion  
  #
  # H0 --> There is no trend present in location/dispersion
  # H1 --> There is trend (upwards or downwards)  
  #
  # Inputs
  #   y         Vector of data.
  #   type      Type of test:
  #               trend - test for changes in trend [default];
  #               deviation - test for changes in deviation;
  #               dispersion - test for changes in dispersion (range).
  #   alpha     Significance level.
  #
  # Outputs
  #   H         Hypothesis (H0/H1).
  #   p.value   P-value.
  #   Htxt      Description of the result.
  #
  # Example
  #   coxstuart(referrals)
  #
  # Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  type <- match.arg(type,c("trend","deviation","dispersion"))
  
  switch(type,
         "trend" = {
           data <- y
         },
         "deviation"={
           data <- splitdata(y)
           n <- dim(data)[2]
           v.sd <- Vectorize(function(i){temp <- sd(data[,i])})
           data <- v.sd(1:n)
         },
         "dispersion"={
           data <- splitdata(y)
           n <- dim(data)[2]
           v.range <- Vectorize(function(i){temp <- range(data[,i])
           temp <- max(temp)-min(temp)})
           data <- v.range(1:n)
         })
  
  # Find number of pairs for comparison
  n <- length(data)
  C <- ceiling(n/2)
  
  # Create pairs
  idx1 <- 1:(n-C)
  idx2 <- (1+C):n
  pair <- data[idx1] - data[idx2]
  
  # Calculate statistic
  Nplus <- sum(pair>0)
  Nminus <- sum(pair<0)
  stat <- min(c(Nplus,Nminus))
  
  # P-value
  if (sum(pair)!=0){
    p <- pbinom(stat,C,0.5)
  } else {
    p <- 1
  }
  
  if (p <= alpha/2){
    H <- 1
    txt <- "H1: There is trend (upwards or downwards)"
  } else {
    H <- 0
    txt <- "H0: There is no trend present in location/dispersion"
  }
  
  return(list(H=H,p.value=p,Htxt=txt))
  
}

# ------------------------------------------
splitdata <- function(y){
  # Helper function
  
  n <- length(y)
  
  # Find K according to Cox Stuart guidelines
  if (n < 48){
    k <- 2
  } else if(n < 64){
    k <- 3
  } else if(n < 90){
    k <- 4
  } else {
    k <- 5
  }
  
  # Split time series to subsamples
  srem <- n %% k
  ktimes <- floor(n/k)
  idx <- array(1:(n-srem),c(k,ktimes))
  idx[,(round(ktimes/2)+1):ktimes] <- idx[,(round(ktimes/2)+1):ktimes] + srem   
  data <- array(y[idx], c(k,ktimes))
  
}

decomp <- function(y,m=NULL,s=NULL,trend=NULL,outplot=c(FALSE,TRUE),
                   decomposition=c("multiplicative","additive"),
                   h=0,type=c("mean","median","pure.seasonal"),w=NULL)
{
  # Decomposition of series
  #
  # Inputs:
  #   y             Time series vector (can be ts object).
  #   m             Seasonal period. If y is a ts object then the default is its frequency.
  #   s             Starting period in the season. If y is a ts object then default is read.
  #   trend         Vector of level/trend of the time series.
  #                 If NULL then the level/trend is calculated using CMA.
  #   outplot       If TRUE provide a plot of the decomposed components.
  #   decomposition Type of seasonal decomposition: "multiplicative" or "additive".
  #   h             Forecast horizon for seasonal component. 
  #   type          Type of calculation for seasonal component:
  #                   "mean"          - The mean of each seasonal period
  #                   "median"        - The median of each seasonal period
  #                   "pure.seasonal" - Model using a pure seasonal model
  #   w             Percentage or number of observations to winsorise in the calculation
  #                 of mean seasonal indices. If w>1 then it is the number of observations, 
  #                 otherwise it is a percentage. If type != "mean" then this is ignored.
  #
  # Outputs:
  #   List with the following elements:
  #     trend         Trend component
  #     season        Seasonal component
  #     irregular     Irregular component
  #     f.season      Forecasted seasonal component if h>0
  #     g             Purse seasonal model parameters
  #
  # Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Defaults
  outplot <- outplot[1]
  decomposition <- decomposition[1]
  type <- type[1]
  
  # Get m (seasonality)
  if (is.null(m)){
    if (class(y)[1] == "ts"){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }
  
  # Get starting period in seasonality if available
  if (is.null(s)){
    if (class(y)[1] == "ts"){
      s <- start(y)
      s <- s[2]
    } else {
      s <- 1
    }
  } 
  
  n <- length(y)
  
  if ((decomposition == "multiplicative") && (min(y)<=0)){
    decomposition <- "additive"
  }
  
  # If trend is not given then calculate CMA
  if (is.null(trend)){
    trend <- cmav(y=y,ma=m,fill=TRUE,outplot=FALSE,fast=TRUE)
  } else {
    if (n != length(trend)){
      stop("Length of series and trend input do not match.")
    }
  }
  
  if (decomposition == "multiplicative"){
    ynt <- y/trend  
  } else {
    ynt <- y - trend
  }
  
  if (outplot != FALSE){
    ymin.s <- min(ynt)
    ymax.s <- max(ynt)
    yminmax.s <- c(ymin.s-0.1*(ymax.s-ymin.s),ymax.s+0.1*(ymax.s-ymin.s))
    ymin.y <- min(y)
    ymax.y <- max(y)
    yminmax.y <- c(ymin.y-0.1*(ymax.y-ymin.y),ymax.y+0.1*(ymax.y-ymin.y))
  }
  
  # Fill with NA start and end of season
  k <- m - (n %% m)
  ks <- s-1
  ke <- m - ((n+ks) %% m)
  ynt <- c(rep(NA,times=ks),as.vector(ynt),rep(NA,times=ke))
  ns <- length(ynt)/m
  ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
  colnames(ynt) <- paste("p",1:m,sep="")
  rownames(ynt) <- paste("s",1:ns,sep="")
  
  # If h>0 then produce forecasts of seasonality 
  g <- NULL
  if (type=="mean"){
    # Calculate the seasonality as the overall mean
    if (!is.null(w)){
      ynt <- colWins(ynt,w)
    }
    season <- colMeans(ynt, na.rm=TRUE)
    if (h>0){
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h)])
    } else {
      f.season <- NULL
    }
    i.season <- rep(season,ns)
    i.season <- i.season[(ks+1):(ns*m-ke)]+y*0
  }
  if (type=="median"){
    # Calculate the seasonality as the overall median
    season <- array(NA,c(1,m))
    for (si in 1:m){
      season[si] <- median(ynt[,si], na.rm=TRUE)
    }
    if (h>0){
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-ke+1):(m-ke+h)])
    } else {
      f.season <- NULL
    }
    i.season <- rep(season,ns)
    i.season <- i.season[(ks+1):(ns*m-ke)] + y*0
  }
  if (type=="pure.seasonal"){
    # Seasonality is modelled with a pure seasonal smoothing
    sout <- decomp.opt.sfit(ynt=ynt,costs="MSE",n=n,m=m)
    g <- sout$g
    if (h>0){
      season <- sout$season
      f.season <- rep(season,h+m*2)
      f.season <- as.vector(f.season[(m-k+1):(m-k+h)])
      f.season <- rep(season, h %/% m + 1)[1:h]
    } else {
      f.season <- NULL
    }
    season.sample <- matrix(t(ynt),ncol=1)        
    season.sample <- season.sample[!is.na(season.sample)]
    i.season <- as.vector(decomp.fun.sfit(season.sample,g,n,m)$ins) + y*0
  }
  
  # Convert f.season to ts object
  if (class(y)[1] == "ts" && h>0){
    s <- end(y)
    if (s[2]==m){
      s[1] <- s[1]+1
      s[2] <- 1
    } else {
      s[2] <- s[2]+1
    }
    f.season <- ts(f.season,start=s,frequency=m)  
  } 
  
  if (decomposition == "multiplicative"){
    resid <- y - (trend*i.season)
  } else {
    resid <- y - (trend+i.season)
  }
  
  # Produce plots
  if (outplot == TRUE){
    # Write down the default values of par
    def.par <- par(no.readonly = TRUE);
    par(mfrow=c(4,1),mar=c(0,2,0,0),oma=c(2,2,2,2))
    
    # Series
    plot(1:n,as.vector(y),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),ylim=yminmax.y,xaxs = "i")
    if (decomposition == "multiplicative"){
      lines(1:n,trend*i.season,col="red",lty=1)
    } else {
      lines(1:n,trend+i.season,col="red",lty=1)
    }
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.y,yminmax.y[2:1]),border=NA,col="gray93")
    }
    mtext("Data",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    legend("topleft",c("Data","Reconstructed"),lty=c(1,1),lwd=c(1,1),col=c("black","red"),cex=0.8,bty="n",horiz=TRUE)
    
    # Trend
    plot(1:n,as.vector(trend),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.y,xaxs="i")
    mtext("Trend",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.y,yminmax.y[2:1]),border=NA,col="gray93")
    }
    
    # Season
    plot(1:n,as.vector(i.season),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.s,xaxs="i")
    mtext("Season",side=2,cex=0.8,padj=-2.5)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.s,yminmax.s[2:1]),border=NA,col="gray93")
      lines((n+1):(n+h),f.season,col="blue")
    }
    if (decomposition == "multiplicative"){
      lines(1:(n+h),rep(1,n+h),lty=2,col="grey")
    } else {
      lines(1:(n+h),rep(0,n+h),lty=2,col="grey")
    }
    
    # Irregular
    yminmax.i = yminmax.y-mean(yminmax.y)
    plot(1:n,as.vector(resid),type="l",xlab="",xaxt='n',ylab="",yaxt='n',xlim=c(1,n+h),lty=1,ylim=yminmax.i,xaxs="i")
    mtext("Irregular",side=2,cex=0.8,padj=-2.5)
    text(1,yminmax.i[2],paste("RMSE:",round(sqrt(mean(resid^2)),2)),pos=4,cex=0.8)
    axis(2,cex.axis=0.8)
    if (h>0){
      polygon(c(n+1,n+1,n+h,n+h),c(yminmax.i,yminmax.i[2:1]),border=NA,col="gray93")
    }
    lines(1:n,rep(0,n),lty=2,col="grey")
    
    par(def.par) 
  }
  
  return(list(trend=trend,season=i.season,irregular=resid,
              f.season=f.season,g=g))
  
}


decomp.opt.sfit <- function(ynt,costs,n,m){
  # Optimise pure seasonal model and predict out-of-sample seasonality
  g0 <- c(0.001,colMeans(ynt,na.rm=TRUE))       # Initialise seasonal model
  season.sample <- matrix(t(ynt),ncol=1)        # Transform back to vector
  season.sample <- season.sample[!is.na(season.sample)]
  opt <- optim(par=g0, decomp.cost.sfit, method = "Nelder-Mead", season.sample=season.sample, 
               cost=costs, n=n, m=m, control = list(maxit = 2000))
  g <- opt$par
  season <- decomp.fun.sfit(season.sample,g,n,m)$outs
  return(list("season"=season,"g"=g))
}

decomp.fun.sfit <- function(season.sample,g,n,m){
  # Fit pure seasonal model
  s.init <- g[2:(m+1)]
  season.fit <- c(s.init,rep(NA,n))
  for (i in 1:n){
    season.fit[i+m] <- season.fit[i] + g[1]*(season.sample[i] - season.fit[i])
  }
  return(list("ins"=season.fit[1:n],"outs"=season.fit[(n+1):(n+m)]))  
}

decomp.cost.sfit <- function(g,season.sample,cost,n,m){
  # Cost function of pure seasonal model
  err <- season.sample-decomp.fun.sfit(season.sample,g,n,m)$ins
  err <- decomp.cost.err(err,cost)
  if (g[1]<0 | g[1]>1){
    err <- 9*10^99
  }
  return(err)   
}

decomp.cost.err <- function(err,cost){
  # Cost calculation
  if (cost == "MAE"){
    err <- mean(abs(err))
  }
  if (cost == "MdAE"){
    err <- median(abs(err))
  }
  if (cost == "MSE"){
    err <- mean((err)^2)
  }
  if (cost == "MdSE"){
    err <- median((err)^2)
  }
  return(err)
}

mseastest <- function(y,m=NULL,type=c("pearson","spearman","kendall"),cma=NULL,
                      sn=1, alpha=0.05,outplot=c(0,1,2))
{
  # Multiplicative seasonality test
  # 
  # Inputs:
  #   y             Time series vector (can be ts object)
  #   m             Seasonal period. If y is a ts object then the default is its frequency  
  #   type          Test is based on:
  #                   - "pearson" correlation [Default]
  #                   - "spearman" correlation
  #                   - "kendall" correlation
  #   cma           Level/trend of time series. If not given a central moving average is calculated
  #   sn            Seasonal periods of decreasing magnitude to consider for the test. Default = 1.
  #   alpha         Significance level. Default = 0.05
  #   outplot       Produce a plot of seasonal elements and level
  #                   0 - No plot [Default]
  #                   1 - Scatterplot
  #                   2 - Time series plot
  #
  # Outputs:
  #   out$is.multiplicative   If TRUE the test found evidence of multiplicative seasonality
  #   out$statistic           The test statistic
  #   out$pvalue              P-value of the test
  #
  # Example:
  #   mseastest(referrals)
  #
  # Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Defaults
  type <- type[1]
  outplot <- outplot[1]
  
  # Get m (seasonality)
  if (is.null(m)){
    if (class(y) == "ts"){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }
  
  # Test sn
  if (sn>m){
    sn <- m
  }
  
  # Calculate CMA if not provided
  if (is.null(cma)){
    cma <- cmav(y,ma=m)
  }
  
  # Remove trend from time series
  seas <- (y - cma)
  
  # Convert to seasonal matrix
  n <- length(y)
  k <- m - (n %% m)
  seas <- c(as.vector(seas),rep(NA,times=k))
  ns <- length(seas)/m
  seas <- matrix(seas,nrow=ns,ncol=m,byrow=TRUE)
  colnames(seas) <- paste("m",1:m,sep="")
  
  # Find size of seasonality and adjust for direction (of mean size)
  mag <- colMeans(seas, na.rm = TRUE)
  idx <- order(abs(mag),decreasing=TRUE)
  ssign <- (mag<0)*-2+1
  seas <- seas*matrix(rep(ssign,ns),nrow=ns,ncol=m,byrow=TRUE)
  idx <- idx[1:sn]
  
  # Measure correlation for sn seasonal periods
  cor.size <- array(NA,c(sn,1))
  cor.pvalue <- array(NA,c(sn,1))
  for (i in 1:sn){
    X <- seas[,idx[i]]
    X <- X[!is.na(X)]
    Y <- cma[seq(idx[i],n,m)]
    Y <- Y[!is.na(X)]
    if (length(X)<3){
      stop("Not enough seasons to test for multiplicative seasonality.")
    }
    test <- cor.test(X,Y,method=type,alternative="greater")
    cor.size[i] <- test$estimate
    cor.pvalue[i] <- test$p.value
  }
  
  # Aggregate cases
  p <- median(cor.pvalue)
  c <- median(cor.size)
  
  # Reset p-value for negative correlations
  if (c <= 0){
    p <- 1
  }
  
  # Do test
  if (p <= alpha){
    is.multiplicative <- TRUE
    H <- "Multiplicative"
  } else {
    is.multiplicative <- FALSE
    H <- "Additive"
  }
  
  # Plot if requested
  if (outplot == 1){
    plot.title = paste(H, " seasonality (pval: ",round(p,3),")",sep="")
    xmin <- as.vector(seas[,idx[1:sn]])
    xmin <- xmin[!is.na(xmin)]
    xmax <- max(xmin)
    xmin <- min(xmin)
    xminmax <- c(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin))
    ymin <- as.vector(mapply(Vectorize(seq),idx[1:sn],idx[1:sn]+m*ns,m))
    ymin <- cma[ymin[ymin<n]]
    ymax <- max(ymin)
    ymin <- min(ymin)
    yminmax <- c(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
    cmp <- rainbow(sn)
    X <- seas[,idx[1]]
    X <- X[!is.na(X)]
    Y <- cma[seq(idx[1],n,m)]
    Y <- Y[!is.na(X)]
    plot(X, Y, col=cmp[1], pch=20, xlab="Seasonal value", ylab="Level value",
         main=plot.title, xlim=xminmax, ylim=yminmax)
    if (sn>1){
      text(X[1],Y[1],"s1",col=cmp[1], pos=3, cex=0.7)
      for (i in 2:sn){
        X <- seas[,idx[i]]
        X <- X[!is.na(X)]
        Y <- cma[seq(idx[i],n,m)]
        Y <- Y[!is.na(X)]
        points(X,Y,col=cmp[i],pch=20)
        text(X[1],Y[1],paste("s",i,sep=""),col=cmp[i], pos=3, cex=0.7)
      }
    }
  }
  
  if (outplot == 2){
    cmp <- rainbow(sn)
    plot.title = paste(H, " seasonality (pval: ",round(p,3),")",sep="")
    plot(1:n,as.vector(y),type="l",main=plot.title,ylab="",xlab="Period")
    lines(as.vector(cma),lty=1,col="black",lwd=2)
    midx <- seq(idx[1],n,m)
    points(midx,y[midx],col=cmp[1],pch=16)
    if (sn>1){
      if (y[midx[1]]>cma[midx[1]]){
        text(midx[1],y[midx[1]],"s1",col=cmp[1], pos=3, cex=0.7)
      }
      for (i in 2:sn){
        midx <- seq(idx[i],n,m)  
        points(midx,(y)[midx],col=cmp[i],pch=16)
        if (y[midx[1]]>cma[midx[1]]){
          text(midx[1],y[midx[1]],paste("s",i,sep=""),col=cmp[i], pos=3, cex=0.7)
        } else {
          text(midx[1],y[midx[1]],paste("s",i,sep=""),col=cmp[i], pos=1, cex=0.7)
        }
      }
    }
  }
  
  return(list("is.multiplicative"=is.multiplicative,"statistic"=c,"pvalue"=p))
  
}

seasplot <- function(y,m=NULL,s=NULL,trend=NULL,colour=NULL,alpha=0.05,
                     outplot=c(1,0,2,3,4,5),decomposition=c("multiplicative","additive"),
                     cma=NULL,labels=NULL,...)
{
  # Seasonal plots and crude trend/season tests
  #
  # Inputs:
  #   y             Time series vector (can be ts object)
  #   m             Seasonal period. If y is a ts object then the default is its frequency  
  #   s             Starting period in the season. If y is a ts object then default is read
  #   trend         If TRUE then a trend is assumed and is removed using CMA
  #                 If FALSE then no trend is assumed
  #                 If NULL then trend is identified and removed if found
  #   colour        Single colour override for plots
  #   alpha         Significance level for statistical tests (kpss and friedman)
  #   outplot       Provide plot output:
  #                   0 - None
  #                   1 - Seasonal diagramme
  #                   2 - Seasonal boxplots
  #                   3 - Seasonal subseries
  #                   4 - Seasonal distribution
  #                   5 - Seasonal density
  #   decomposition Type of seasonal decomposition: "multiplicative" or "additive".
  #   cma           Input precalculated level/trend for the analysis. Overrides trend=NULL.
  #   labels        External labels for the seasonal periods. Use NULL for default. 
  #                 If length(labels) < m, then this input is ignored.
  #   ...           Additional arguments can be passed to plotting functions. For example use 
  #                 main="" to replace the title.
  #
  # Outputs:
  #   List with the following elements:
  #     season        Matrix of (detrended) seasonal elements
  #     season.exist  TRUE/FALSE results of friedman test
  #     season.pval   Friedman test p-value
  #     trend         CMA estimate or NULL if trend == FALSE
  #     trend.exist   TRUE/FALSE result of kpss test
  #     trend.pval    kpss (approximate) p-value
  #
  # Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Defaults
  decomposition <- match.arg(decomposition,c("multiplicative","additive"))
  outplot <- outplot[1]
  
  # Get m (seasonality)
  if (is.null(m)){
    if (any(class(y) == "ts")){
      m <- frequency(y)
    } else {
      stop("Seasonality not defined (y not ts object).")
    }
  }
  
  # Get starting period in seasonality if available
  if (is.null(s)){
    if (any(class(y) == "ts")){
      s <- start(y)
      s <- s[2]
      # Temporal aggregation can mess-up s, so override if needed
      if (is.na(s)){s<-1}
    } else {
      s <- 1
    }
  } 
  
  # Make sure that labels input is fine
  if (!is.null(labels)){
    if (length(labels) < m){
      labels  <- NULL
    } else {
      labels <- labels[1:m]
    }
  }
  
  if (is.null(labels)){
    labels <- paste(1:m)
  }
  
  n <- length(y)
  
  if ((decomposition == "multiplicative") && (min(y)<=0)){
    decomposition <- "additive"
  }
  
  # Override trend if cma is given
  if (!is.null(cma)){
    trend <- NULL
    if (n != length(cma)){
      stop("Length of series and cma do not match.")
    }
  }
  
  # Calculate CMA
  if ((is.null(trend) || trend == TRUE) && (is.null(cma))){
    cma <- cmav(y=y,ma=m,fill=TRUE,outplot=FALSE)
  }
  
  # Test for changes in the CMA (trend)
  if (is.null(trend)){
    trend.pval <- coxstuart(cma)$p.value 
    trend.exist <- trend.pval <= alpha/2
    trend <- trend.exist
  } else {
    trend.exist <- NULL
    trend.pval <- NULL
  }
  
  if (trend == TRUE){
    if (decomposition == "multiplicative"){
      ynt <- y/cma  
    } else {
      ynt <- y - cma
    }
    title.trend <- "(Detrended)"
  } else {
    ynt <- y
    title.trend <- ""
    cma <- NULL
  }
  
  ymin <- min(ynt)
  ymax <- max(ynt)
  ymid <- median(ynt)
  
  # Fill with NA start and end of season
  k <- m - (n %% m)
  ks <- s-1
  ke <- m - ((n+ks) %% m)
  ynt <- c(rep(NA,times=ks),as.vector(ynt),rep(NA,times=ke))
  ns <- length(ynt)/m
  ynt <- matrix(ynt,nrow=ns,ncol=m,byrow=TRUE)
  colnames(ynt) <- labels
  rownames(ynt) <- paste("s",1:ns,sep="")
  
  # Check seasonality with Friedman
  if (m>1 && (length(y)/m)>=2){
    season.pval <- friedman.test(ynt)$p.value
    season.exist <- season.pval <= alpha
    if (season.exist==TRUE){
      title.season <- "Seasonal"
    } else {
      title.season <- "Nonseasonal"
    }
  } else {
    season.pval <- NULL
    season.exist <- NULL
    title.season <- "Nonseasonal"
  }
  
  # Produce plots
  if (outplot != 0){
    yminmax <- c(ymin - 0.1*(ymax-ymin),ymax + 0.1*(ymax-ymin))
    yminmax <- c(-1,1)*max(abs(ymid-yminmax))+ymid
    if (is.null(season.pval)){
      plottitle <- paste(title.trend, "\n", title.season,sep="")
    } else {
      plottitle <- paste(title.trend, "\n", title.season,
                         " (p-val: ",round(season.pval,3),")",sep="")
    }
    # Allow user to override plot defaults
    args <- list(...)
    if (!("main" %in% names(args))){
      addtitle <- TRUE
    } else {
      addtitle <- FALSE
    }
    if (!("xlab" %in% names(args))){
      args$xlab <- "Period"
    }
    if (!("ylab" %in% names(args))){
      args$ylab <- ""
    }
    if (!("yaxs" %in% names(args))){
      args$yaxs <- "i"
    }
    if (!("ylim" %in% names(args))){
      args$ylim <- yminmax
    }
    # Remaining defaults
    args$x <- args$y <- NA
  }
  
  if (outplot == 1){
    # Conventional seasonal diagramme
    if (is.null(colour)){
      cmp <- colorRampPalette(brewer.pal(9,"YlGnBu")[4:8])(ns)
    } else {
      cmp <- rep(colour,times=ns)
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    if (addtitle){
      args$main <- paste0("Seasonal plot ",main=plottitle)
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    for (i in 1:ns){
      lines(ynt[i,],type="l",col=cmp[i])
    }
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
    legend("topleft",c("Oldest","Newest"),col=c(cmp[1],cmp[ns]),lty=1,bty="n",lwd=2,cex=0.7)
    axis(1,at=1:m,labels=labels)
  } 
  if (outplot == 2){
    # Seasonal boxplots
    if (is.null(colour)){
      cmp <- brewer.pal(3,"Set3")[1]
    } else {
      cmp <- colour
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (addtitle){
      args$main <- paste0("Seasonal boxplot ",main=plottitle)
    }
    args$x <- ynt
    args$col <- cmp
    # Produce plot
    do.call(boxplot,args)
    lines(c(0,m+1),c(ymid,ymid),col="black",lty=2)
  }
  if (outplot == 3){
    # Subseries plots
    if (is.null(colour)){
      cmp <- brewer.pal(3,"Set1")
    } else {
      cmp <- colour
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m*ns)
    }
    if (addtitle){
      args$main <- paste0("Seasonal subseries ",main=plottitle)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    lines(c(1,ns),median(ynt[,1],na.rm=TRUE)*c(1,1),col=cmp[1],lwd=2)
    for (i in 1:m){
      lines((1+(i-1)*ns):(ns+(i-1)*ns),ynt[,i],type="o",col=cmp[2],pch=20,cex=0.75)
      lines(c((1+(i-1)*ns),(ns+(i-1)*ns)),median(ynt[,i],na.rm=TRUE)*c(1,1),col=cmp[1],lwd=2)
      if (i < m){
        lines(c(1,1)*i*ns+0.5,yminmax,col="gray")
      }
    }
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    axis(1,at=seq(0.5+(ns/2),0.5+m*ns-ns/2,ns),labels=labels)
  }
  if (outplot == 4){
    # Seasonal distribution
    if (is.null(colour)){
      cmp <- "royalblue"
    }
    qntl <- matrix(NA,nrow=9,ncol=m)
    for (i in 1:m){
      qntl[,i] <- quantile(ynt[!is.na(ynt[,i]),i], c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1))
    }
    # Additional plot options
    if (!("xlim" %in% names(args))){
      args$xlim <- c(1,m)
    }
    if (addtitle){
      args$main <- paste0("Seasonal distribution ",main=plottitle)
    }
    if (!("xaxs" %in% names(args))){
      args$xaxs <- "i"
    }
    args$xaxt <- "n"
    # Produce plot
    do.call(plot,args)
    polygon(c(1:m,rev(1:m)),c(qntl[7,],rev(qntl[1,])),col=gray(0.8),border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[6,],rev(qntl[2,])),col="lightblue",border=NA)
    polygon(c(1:m,rev(1:m)),c(qntl[5,],rev(qntl[3,])),col="skyblue",border=NA)
    lines(1:m,qntl[4,],col=cmp,lwd=2)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    legend("topleft",c("Median","25%-75%","10%-90%","MinMax"),col=c(cmp,"skyblue","lightblue",gray(0.8)),lty=1,bty="n",lwd=2,cex=0.7)
    axis(1,at=1:m,labels=labels)
    box()
  }
  if (outplot == 5){
    dnst <- matrix(NA,nrow=m,ncol=512)
    for (i in 1:m){
      tmp <- density(ynt[!is.na(ynt[,i]),i], bw = "SJ", n = 512, from = yminmax[1], to = yminmax[2])
      dnst[i,] <- tmp$y
      if (i == 1){
        llc <- tmp$x
      }
    }
    cmp <- c(rgb(0,0,0,0), colorRampPalette((brewer.pal(9,"Blues")[3:9]))(100))
    # Additional plot options
    if (addtitle){
      args$main <- paste0("Seasonal density ",main=plottitle)
    }
    args$xaxt <- "n"
    args$col <- cmp
    args$x <- 1:m
    args$y <- llc
    args$z <- dnst
    # Produce plot
    do.call(image,args)
    lines(colMeans(ynt,na.rm=TRUE),type="o",lty=1,bg=brewer.pal(3,"Set1")[1],pch=21,cex=1.1,lwd=2)
    lines(c(0,m*ns+1),c(ymid,ymid),col="black",lty=2)
    box()
    axis(1,at=1:m,labels=labels)
  }
  
  out <- list(season=ynt,season.exist=season.exist,season.pval=season.pval,
              trend=cma,trend.exist=trend.exist,trend.pval=trend.pval)
  out <- structure(out,class="classseas")
  return(out)
  
}

summary.classseas <- function(object,...){
  print(object)
}

print.classseas <- function(x,...){
  cat("Results of statistical testing\n")
  if (!is.null(x$trend.exist)){
    cat(paste0("Evidence of trend: ",x$trend.exist, "  (pval: ",round(x$trend.pval,3),")\n"))
  } else {
    cat("Presence of trend not tested.\n")
  }
  cat(paste0("Evidence of seasonality: ",x$season.exist, "  (pval: ",round(x$season.pval,3),")\n"))
}
