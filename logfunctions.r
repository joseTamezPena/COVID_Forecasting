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
  opLeft <- data$fatalities[lastObs];
  adjusto <- 0.0;
  loops <- 0;
  defecitRatio <- 1.0;
  accAdjust <- adjini;
  error <- TRUE;

  data$fatalities <- adjini*data$fatalities;
  data$newfatalities <- adjini*data$newfatalities;
#  cat(adjusto,":",adjust,":",accAdjust,":",gratio,":",abs(log(accAdjust)),":",log(gratio),"\n")
  
  while ((abs((adjusto - accAdjust)/accAdjust) > 0.025) && (abs(log(accAdjust)) <= abs(log(gratio))))
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
#                           control=list(warnOnly = TRUE),
                           start=list(ro = fro,to = fto),
    ));
    if (!inherits(cdfestimate, "try-error"))
    {
      smo <- try(summary(cdfestimate))
      if (!inherits(smo, "try-error"))
      {
        error <- FALSE;
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
            fro = 0.50*fro+0.50*psmo$coefficients[1,1];
            fto = 0.90*fto+0.10*psmo$coefficients[2,1];
            if (fto < 0.5*to)
            {
              fto <- to;
            }
            if (fto > 10*to)
            {
              fto <- to;
            }
            nleft <- 1.0-logisticcdf(data$days[lastObs], fro, fto);
            if (data$fatalities[lastObs] > 0.5) 
            {
              nleft <- 1.0 - nleft;
              adjust <- (0.1 + nleft)/(0.1 + pLeft);
            }
            else
            {
              adjust <- (0.1 + pLeft)/(0.1 + nleft);
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
#            cat(sprintf("%8.5f : (%5.3f,%5.3f) Adjust= %5.3f, ro= %8.3f to=%8.3f",abs((adjusto - accAdjust)/accAdjust),pLeft,nleft,accAdjust,fro,fto),"\n")
          }
        }
      }
    }
    if (loops > 100)
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
#  cat(sprintf("Defecit= %5.3f, ro= %8.3f to=%8.3f",defecitRatio,fro,fto),"\n")
  
  models <- list(CDF=cdfestimate,
                 pdf=pdfestimate,
                 ro=fro,
                 to=fto,
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

bootstraplogisitcfit <- function(data,inifit,ratiorange=1.5,n=500,daysrange=c(1:nrow(data)))
{
  toestimations <- numeric(n);
  roestimations <- numeric(n);
  optGain <- ratiorange^(runif(n, -1,1));

  for (bsamples in 1:n)
  {
    toestimations[bsamples] <- inifit$to;
    roestimations[bsamples] <- inifit$ro
#    bootresample <- sample(nrow(data),nrow(data),TRUE)
#    bootresample <- bootresample[order(bootresample)];
#    print(bootresample)
#    ndata <- data[bootresample,];
    ndata <- data;
    maxgain <- min(2.0/max(ndata$fatalities),optGain[bsamples]);
    ndata$fatalities <- maxgain*ndata$fatalities;
    ndata$newfatalities <- maxgain*ndata$newfatalities;
    iniadjs <- 0.7+0.3*inifit$defecitRatio;
#    iniadjs <- 1.0;
    nfit <- try(logisitcfit(ndata,inifit$ro,inifit$to,1.2*inifit$defecitRatio,iniadjs,daysrange=daysrange))
    if (!inherits(nfit, "try-error"))
    {
      toestimations[bsamples] <- nfit$to;
      roestimations[bsamples] <- nfit$ro
    }
  }
  result <- list(to=toestimations,ro=roestimations,gains=optGain)
  return(result)
}
