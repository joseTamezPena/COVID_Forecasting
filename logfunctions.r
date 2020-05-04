logisticcdf <- function(days,ro,to) 
{ 
  cdf <- 1.0/(1.0+exp(ro*(days-to))); 
  return(cdf);
}

logisticpdf <- function(days,ro,to) 
{ 
  pdf <- -ro*exp(ro*(days-to))/((1.0+exp(ro*(days-to)))^2); 
  return(pdf);
}

logisitcfit <- function(data,ro,to,gratio=2,adjini=1,lowf=1/4,daysrange=c(1:nrow(data)))
{
  lastObs <- daysrange[length(daysrange)];
  if (lowf > 0)
  {
    data <- data[c(1:lastObs),];
    data$fatalities <- lowess(data$fatalities[c(1:lastObs)],f=lowf)$y;
    data$newfatalities <- 0.25*data$newfatalities[c(1:lastObs)] + 0.75*c(data$newfatalities[1],(data$fatalities[2:lastObs]-data$fatalities[1:lastObs-1]));
    data$newfatalities <- runmed(data$newfatalities,5)
    data$newfatalities <- lowess(data$newfatalities,f=lowf)$y
  }
#  print(daysrange)
  data <- data[daysrange,];
  lastObs <- nrow(data);
  filterCDF <- data$fatalities;
  filtercdf <- data$newfatalities;
  cdfestimate <- NULL;
  pdfestimate <- NULL;
  firscdftestimation <- cdfestimate;
  firspdftestimation <- pdfestimate;
  adjust <- 1.0;
  fro <- ro;
  fto <- to; 
  opLeft <- data$fatalities[lastObs];
  adjusto <- 0.0;
  loops <- 0;
  defecitRatio <- 1.0;
  accAdjust <- adjini;
  error <- FALSE;

  data$fatalities <- adjini*data$fatalities;
  data$newfatalities <- adjini*data$newfatalities;
#  cat(adjusto,":",adjust,":",accAdjust,":",gratio,":",abs(log(accAdjust)),":",log(gratio),"\n")
  
  while ((abs(adjusto - adjust) > 0.01) && (abs(log(accAdjust)) <= abs(log(gratio))))
  {
    adjusto <- 1.0;
    loops <- loops + 1;
    pLeft <- 1.0-data$fatalities[lastObs];
    if (data$fatalities[lastObs] > 0.5) 
    {
      pLeft <- 1.0 - pLeft;
    }
#    cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f",pLeft,fro,fto),"\n")
    cdfestimate <- try(nls(fatalities ~ logisticcdf(days, ro, to),
                           data = data,
#                           control=list(warnOnly = TRUE),
                           start=list(ro = fro,to = fto),
    ));
    if (!inherits(cdfestimate, "try-error"))
    {
      smo <- try(summary(cdfestimate))
      if (!inherits(smo, "try-error"))
      {
        fro = smo$coefficients[1,1];
        fto = smo$coefficients[2,1];
        if (fto < 0.5*to)
        {
          fto <- to;
        }
        if (fto > 10*to)
        {
          fto <- to;
        }
#         cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f",pLeft,fro,fto),"\n")
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
            fro = psmo$coefficients[1,1];
            fto = psmo$coefficients[2,1];
            if (fto < 0.5*to)
            {
              fto <- to;
            }
            if (fto > 10*to)
            {
              fto <- to;
            }
            nleft <- 1.0-logisticcdf(data$days[lastObs], fro, fto);
            adjust <- (0.5 + pLeft)/(0.5 + nleft);
            if (data$fatalities[lastObs] > 0.5) 
            {
              nleft <- 1.0 - nleft;
              adjust <- (0.25 + nleft)/(0.25 + pLeft);
            }
            else
            {
              adjust <- (0.25 + pLeft)/(0.25 + nleft);
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
 #                 cat(sprintf("(%5.3f,%5.3f) Adjust= %5.3f, ro= %8.3f to=%8.3f",pLeft,nleft,accAdjust,fro,fto),"\n")
          }
        }
      }
      else
      {
        error <- TRUE;
      }
    }
    else
    {
      error <- TRUE;
    }
    if (loops > 1000)
    {
      adjusto <- adjust
    }
    if (loops == 1)
    {
      firscdftestimation <- cdfestimate;
      firspdftestimation <- pdfestimate;
    }
  }
  if (opLeft > 0)
  {
    defecitRatio <- logisticcdf(data$days[lastObs], fro, fto)/opLeft;
  }
#    cat(sprintf("Defecit= %5.3f, ro= %8.3f to=%8.3f",defecitRatio,fro,fto),"\n")
  
  models <- list(CDF=cdfestimate,
                 pdf=pdfestimate,
                 ro=fro,
                 to=fto,
                 firstpdf = firscdftestimation,
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

bootstraplogisitcfit <- function(data,inifit,ratiorange=1.25,n=500,daysrange=c(1:nrow(data)))
{
  toestimations <- numeric(n);
  roestimations <- numeric(n);
  optGain <- ratiorange^(runif(n, -1,1));

  for (bsamples in 1:n)
  {
    toestimations[bsamples] <- inifit$to;
    roestimations[bsamples] <- inifit$ro
    bootresample <- sample(nrow(data),nrow(data),TRUE)
    bootresample <- bootresample[order(bootresample)];
#    print(bootresample)
    ndata <- data[bootresample,];
    maxgain <- min(0.9/max(ndata$fatalities),optGain[bsamples]);
    ndata$fatalities <- maxgain*ndata$fatalities;
    ndata$newfatalities <- maxgain*ndata$newfatalities;
    nfit <- try(logisitcfit(ndata,inifit$ro,inifit$to,1.2*inifit$adjust,inifit$adjust,daysrange=daysrange))
    if (!inherits(nfit, "try-error"))
    {
      toestimations[bsamples] <- nfit$to;
      roestimations[bsamples] <- nfit$ro
    }
  }
  result <- list(to=toestimations,ro=roestimations,gains=optGain)
  return(result)
}
