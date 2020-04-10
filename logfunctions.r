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

logisitcfit <- function(data,ro,to)
{
  cdfestimate <- NULL;
  pdfestimate <- NULL;
  firscdftestimation <- cdfestimate;
  firspdftestimation <- pdfestimate;
  adjust <- 1.0;
  lastObs <- nrow(data);
  fro <- ro;
  fto <- to; 
  opLeft <- data$fatalities[lastObs];
  adjusto <- 0.0;
  loops <- 0;
  defecitRatio <- 1.0;
  accAdjust <- 1.0;
  error <- FALSE;
  
  while ((abs(adjusto - adjust) > 0.005) && (accAdjust < 1.25) && (accAdjust > 0.75))
  {
    adjusto <- 1.0;
    loops <- loops + 1;
    pLeft <- 1.0-data$fatalities[lastObs];
    #    cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f",pLeft,fro,fto),"\n")
    cdfestimate <- try(nls(fatalities ~ logisticcdf(days, ro, to),
                           data = data,
                           start=list(ro = fro,to = fto),
    ));
    if (!inherits(cdfestimate, "try-error"))
    {
      smo <- summary(cdfestimate)
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
      #     cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f",pLeft,fro,fto),"\n")
      pdfestimate <- try(nls(newfatalities ~ logisticpdf(days, ro, to),
                             data = data,
                             start=list(ro = fro ,to = fto),
                             control=list(warnOnly=TRUE)))
      if (!inherits(pdfestimate, "try-error"))
      {
        psmo <- summary(pdfestimate)
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
        adjust <- (0.25 + pLeft)/(0.25 + nleft);
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
        #      cat(sprintf("(%5.3f,%5.3f) Adjust= %5.3f, ro= %8.3f to=%8.3f",pLeft,nleft,accAdjust,fro,fto),"\n")
      }
    }
    else
    {
      error <- TRUE;
    }
    if (loops >50)
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
                 firstpdf = firscdftestimation,
                 firstpdf = firspdftestimation,
                 defecitRatio=defecitRatio,
                 adjust = accAdjust
  )
  class(models) <- "logisticFit";
  if (error) class(models) <- append(class(models),"try-error")
  return (models);
}