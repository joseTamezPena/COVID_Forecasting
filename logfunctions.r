

logisticcdf <- function(days,ro,to,rp=NULL) 
{ 
  if (is.null(rp)) 
  {
     rp = 0.45*ro
  }
  cdf <- 1.0/(1.0+exp(ro*(days-to)));
  postto <- days > (to-21)
  if (any(postto))
  {
    cdfp <- 1.0/(1.0+exp(rp*(days[postto]-to)));
    cdf3 <- 1.0/(1.0+exp(1.0*ro*(days[postto]-(to-21))));
    w = sqrt(2.0*(1.0-cdf3))
    cdf[postto] <- w*cdf[postto] + (1.0-w)*cdfp
  }
  return(cdf);
}

logisticpdf <- function(days,ro,to,rp=NULL) 
{ 
  if (is.null(rp))
  {
    rp <- 0.45*ro
  }
  pdf <- -ro*exp(ro*(days-to))/((1.0+exp(ro*(days-to)))^2); 
  postto <- days > (to-21)
  if (any(postto))
  {
    cdf3 <- 1.0/(1.0+exp(1.0*ro*(days[postto]-(to-21))));
    pdf2 <- -rp*exp(rp*(days[postto]-to))/((1.0+exp(rp*(days[postto]-to)))^2); 
    w = sqrt(2.0*(1.0-cdf3))
    pdf[postto] <- w*pdf[postto] + (1.0-w)*pdf2
  }
  return(pdf);
}

logisitcfit <- function(data,ro,to,gratio=2,adjini=1,lowf=1/4,daysrange=c(1:nrow(data)))
{
  lastObs <- daysrange[length(daysrange)];
#  cat(lastObs,":")
  if (lowf > 0)
  {
    data <- data[c(1:lastObs),];
    data$fatalities <- supsmu(c(1:lastObs),data$fatalities[c(1:lastObs)])$y;
    data$newfatalities <- 0.25*data$newfatalities[c(1:lastObs)] + 0.75*c(data$newfatalities[1],(data$fatalities[2:lastObs]-data$fatalities[1:lastObs-1]));
    data$newfatalities <- runmed(data$newfatalities,5)
    data$newfatalities <- supsmu(c(1:lastObs),data$newfatalities[c(1:lastObs)])$y
  }
#  print(daysrange)
  filterCDF <- data$fatalities;
  filtercdf <- data$newfatalities;
  data <- data[daysrange,];
  lastObs <- nrow(data);
  cdfestimate <- NULL;
  pdfestimate <- NULL;
  firscdftestimation <- cdfestimate;
  firspdftestimation <- pdfestimate;
  adjust <- 1.0;
  fro <- ro;
  fto <- to;
  rpo <- 0.5*fro
  opLeft <- data$fatalities[lastObs];
  adjusto <- 0.0;
  loops <- 0;
  defecitRatio <- 1.0;
  accAdjust <- adjini;
  error <- TRUE;

  data$fatalities <- adjini*data$fatalities;
  data$newfatalities <- adjini*data$newfatalities;
#  cat("Start:",adjusto,":",adjust,":",accAdjust,":",gratio,":",abs(log(accAdjust)),":",log(gratio),"\n")
  
  while ((abs((adjusto - accAdjust)/accAdjust) >= 0.03) && (abs(log(accAdjust)) <= abs(log(gratio))))
  {
    loops <- loops + 1;
    pLeft <- 1.0-data$fatalities[lastObs];
    if (data$fatalities[lastObs] > 0.5) 
    {
      pLeft <- 1.0 - pLeft;
    }
#    cat(sprintf("%8.5f : pLeft= %5.3f, ro= %8.3f to=%8.3f",abs((adjusto - accAdjust)/accAdjust),pLeft,fro,fto),"\n")
    adjusto <- accAdjust;
    cdfestimate <- try(nls(fatalities ~ logisticcdf(days, ro, to),
                           data = data,
                           start=list(ro = fro,to = fto),
    ));
    if (!inherits(cdfestimate, "try-error"))
    {
        smo <- try(summary(cdfestimate))
        if (!inherits(smo, "try-error"))
        {
          error <- FALSE;
          fro <- smo$coefficients[1,1];
          fto <- smo$coefficients[2,1];
          if (fto < 0.10*to)
          {
            fto <- 0.10*to;
          }
          if (fto > 10*to)
          {
            fto <- 10*to;
          }
#         cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f rpo=%8.3f",pLeft,fro,fto,rpo),"\n")
          pdfestimate <- try(nls(newfatalities ~ logisticpdf(days, ro, to),
                                 data = data,
  #                               control=list(warnOnly = TRUE),
                                 start=list(ro = fro ,to = fto),
                                 ))
          if (!inherits(pdfestimate, "try-error"))
          {
            psmo <- try(summary(pdfestimate))
            if (!inherits(psmo, "try-error"))
            {
              fro <- 0.50*fro+0.50*psmo$coefficients[1,1];
              fto <- 0.90*fto+0.10*psmo$coefficients[2,1];
              nleft <- 1.0-logisticcdf(data$days[lastObs], fro, fto);
              if (data$fatalities[lastObs] > 0.5) 
              {
                nleft <- 1.0 - nleft;
                adjust <- (0.01 + nleft)/(0.01 + pLeft);
              }
              else
              {
                adjust <- (0.01 + pLeft)/(0.01 + nleft);
              }
              if (adjust > 1.2)
              {
                adjust <- 1.2;
              }
              if (adjust < 0.8)
              {
                adjust <- 0.8;
              }
              accAdjust <- accAdjust*adjust;
              data$fatalities <- adjust*data$fatalities;
              data$newfatalities <- adjust*data$newfatalities;
#              cat(sprintf("%8.5f : (%5.3f,%5.3f) Adjust= %5.3f, ro= %8.3f to=%8.3f",abs((adjusto - accAdjust)/accAdjust),pLeft,nleft,accAdjust,fro,fto),"\n")
            }
          }
      }
    }
    if (loops > 100)
    {
      adjusto <- accAdjust
    }
    if (loops == 1)
    {
      firscdftestimation <- cdfestimate;
      firspdftestimation <- pdfestimate;
    }
  }

  
  rpo <- 0.45*fro
#  cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f rpo=%8.3f",pLeft,fro,fto,rpo),"\n")
  pdfestimatef <- try(nls(newfatalities ~ logisticpdf(days, fro, fto,rpo),
                         data = data,
#                         start=list(rpo = rpo,fto=fto),
                        start=list(rpo = rpo),
  ))
  if (!inherits(pdfestimatef, "try-error"))
  {
    psmo <- try(summary(pdfestimatef))
    if (!inherits(psmo, "try-error"))
    {
      pdfestimate <- pdfestimatef
      rpo <- psmo$coefficients[1,1]
      if (rpo > 0.25*fro)
      {
        rpo <- 0.45*fro
      }
      if (rpo < fro)
      {
        rpo <- fro
      }
      #      fto <- psmo$coefficients[2,1]
    }
  }
  
#  cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f rpo=%8.3f",pLeft,fro,fto,rpo),"\n")

  if (opLeft > 0)
  {
    defecitRatio <- logisticcdf(data$days[lastObs], fro, fto,rpo)/opLeft;
  }
#  cat(sprintf("Defecit= %5.3f, ro= %8.3f to=%8.3f  rpo=%8.3f",defecitRatio,fro,fto,rpo),"\n")
  models <- list(CDF=cdfestimate,
                 pdf=pdfestimate,
                 ro=fro,
                 to=fto,
                 rpo=rpo,
                 firstcdf = firscdftestimation,
                 firstpdf = firspdftestimation,
                 defecitRatio=defecitRatio,
                 adjust = accAdjust,
                 filterCDF = filterCDF,
                 filterpdf = filtercdf
  )
  class(models) <- "logisticFit";
  if (error) class(models) <- append(class(models),"try-error")
  return (models);
}

bootstraplogisitcfit <- function(data,inifit,ratiorange=1.75,n=200,daysrange=c(1:nrow(data)))
{
  toestimations <- numeric(n);
  roestimations <- numeric(n);
  rpoestimations <- numeric(n);
  defratios <- numeric(n);
  optGain <- ratiorange^(runif(n, -1,1));

  for (bsamples in 1:n)
  {
    toestimations[bsamples] <- inifit$to;
    roestimations[bsamples] <- inifit$ro
    rpoestimations[bsamples] <- inifit$rpo;
    defratios[bsamples] <- inifit$defecitRatio;
#    bootresample <- sample(nrow(data),nrow(data),TRUE)
#    bootresample <- bootresample[order(bootresample)];
#    print(bootresample)
#    ndata <- data[bootresample,];
    ndata <- data;
    maxgain <- min(2.0/max(ndata$fatalities),optGain[bsamples]);
    ndata$fatalities <- maxgain*ndata$fatalities;
    ndata$newfatalities <- maxgain*ndata$newfatalities;
    iniadjs <- inifit$adjust;
#    iniadjs <- 1.0;
    nfit <- try(logisitcfit(ndata,inifit$ro,inifit$to,1.5*inifit$defecitRatio,iniadjs,daysrange=daysrange))
    if (!inherits(nfit, "try-error"))
    {
      if ((nfit$ro < 0) && (nfit$to>0) && (nfit$rpo < 0) )
      {
        toestimations[bsamples] <- nfit$to;
        roestimations[bsamples] <- nfit$ro;
        rpoestimations[bsamples] <- nfit$rpo;
        defratios[bsamples] <- nfit$defecitRatio*maxgain;
      }
    }
  }
  result <- list(to=toestimations,ro=roestimations,rpo=rpoestimations,defratios=defratios,gains=optGain)
  return(result)
}
