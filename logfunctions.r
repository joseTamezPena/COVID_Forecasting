

logisticcdf <- function(days,ro,to,alpha=NULL,lro=NULL,lto=NULL) 
{ 
  ef <- ro*(days-to);
  ef[ef > 20] <- 20;
  ef[ef < -20] <- -20;
  cdf <- 1.0/(1.0+exp(ef));
  if (!is.null(alpha))
  {
    ef <- lro*(days-lto);
    ef[ef > 20] <- 20;
    ef[ef < -20] <- -20;
    cdf2 <- 1.0/(1.0+exp(ef));
    cdf <- alpha*cdf + (1.0-alpha)*cdf2;
  }
  return(cdf);
}

logisticpdf <- function(days,ro,to,alpha=NULL,lro=NULL,lto=NULL) 
{ 
  ef <- ro*(days-to);
  ef[ef > 20] <- 20;
  ef[ef < -20] <- -20;
  pdf <- -ro*exp(ef)/((1.0+exp(ef))^2); 
  if (!is.null(alpha))
  {
    ef <- lro*(days-lto);
    ef[ef > 20] <- 20;
    ef[ef < -20] <- -20;
    pdf2 <- -lro*exp(ef)/((1.0+exp(ef))^2); 
    pdf <- alpha*pdf + (1.0-alpha)*pdf2;
  }
  return(pdf);
}

logisitcfit <- function(data,ro,to,gratio=2,adjini=1,daysrange=c(1:nrow(data)),lro=NULL,lto=NULL,alpha=NULL)
{
  lastObs <- daysrange[length(daysrange)];
  olastObs <- lastObs
#  print(data$fatalities[olastObs])
  #  cat(lastObs,":")
  data <- data[c(1:lastObs),];
  data$fatalities <- supsmu(c(1:lastObs),data$fatalities[c(1:lastObs)])$y;
  data$newfatalities <- 0.25*data$newfatalities[c(1:lastObs)] + 0.75*c(data$newfatalities[1],(data$fatalities[2:lastObs]-data$fatalities[1:lastObs-1]));
  data$newfatalities <- runmed(data$newfatalities,5)
  data$newfatalities <- supsmu(c(1:lastObs),data$newfatalities[c(1:lastObs)])$y
  #  print(daysrange)
  filterCDF <- data$fatalities;
  filtercdf <- data$newfatalities;
  allData <- data
  firsthalf <- data[-daysrange,]
#  print(allData$fatalities[olastObs])
  cdfestimate <- NULL;
  pdfestimate <- NULL;
  firscdftestimation <- cdfestimate;
  firspdftestimation <- pdfestimate;
  adjust <- 1.0;
  fro <- ro;
  fto <- to;
  adjusto <- 0.0;
  loops <- 0;
  defecitRatio <- 1.0;
  accAdjust <- adjini;
  error <- TRUE;
  firstdays <- cbind(days= c(-10000,-1000,-10),fatalities=c(0.0,0.0,0.0),newfatalities=c(0.0,0.0,0.0))
  lasdays <- cbind(days= c(1000,5000,10000),fatalities=c(1.0,1.0,1.0),newfatalities=c(0.0,0.0,0.0))
  allData$fatalities <- adjini*allData$fatalities;
  allData$newfatalities <- adjini*allData$newfatalities;
  fAlldata <- rbind(firstdays,allData,lasdays)
  nls.control(maxiter = 150);
  #  cat("Start:",adjusto,":",adjust,":",accAdjust,":",gratio," Ladd:",abs(log(accAdjust)),"LGr:",log(gratio),"ro:",fro,"to:",fto,"\n")
  oro <- fro;
  oto <- fto;
  {
#    if (length(daysrange) < 0.85*olastObs) gratio <- 1.5*gratio;
    if (is.null(alpha))
    {
      while ((abs((adjusto - accAdjust)/accAdjust) >= 0.01) && (abs(log(accAdjust)) <= abs(log(gratio))))
      {
        loops <- loops + 1;
        pLeft <- 1.0-allData$fatalities[lastObs];
        if (allData$fatalities[lastObs] > 0.5) 
        {
          pLeft <- 1.0 - pLeft;
        }
  #      cat(sprintf("%8.5f : pLeft= %5.3f, ro= %8.3f to=%8.3f",abs((adjusto - accAdjust)/accAdjust),pLeft,fro,fto),"\n")
        adjusto <- accAdjust;
        fAlldata <- rbind(firstdays,allData,lasdays)
        cdfestimate <- try(nls(fatalities ~ logisticcdf(days, ro, to),
                                    data = fAlldata,
                                    start=list(ro = fro,to = fto),
                                    algorithm="plinear"
                                    #                             control = list(maxiter = 150, warnOnly = TRUE)
          ));
        pdfestimate <- cdfestimate
        if (!inherits(cdfestimate, "try-error"))
        {
            smo <- try(summary(cdfestimate))
            if (!inherits(smo, "try-error"))
            {
              if (smo$coefficients[1,1] < 0)
              {
                error <- FALSE;
                fro <- 0.5*fro+0.5*smo$coefficients[1,1];
                fto <- 0.5*fto+0.5*smo$coefficients[2,1];
                if (fro > -0.0001)
                {
                  fro <- -0.0001;
                }
                if (fto < 0.10*to)
                {
                  fto <- 0.10*to;
                }
                if (fto > 10*to)
                {
                  fto <- 10*to;
                }
  #               cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f ",pLeft,fro,fto),"\n")
#                fAlldata <- rbind(firstdays,allData[c(as.integer(2*nrow(allData)/3):olastObs),],lasdays)
                # cdfestimate <- try(nls(fatalities ~ logisticcdf(days, ro, to),
                #                        data = fAlldata,
                #                        start=list(ro = fro ,to = fto),
                #                       algorithm="plinear",
                # ))
                # if (!inherits(cdfestimate, "try-error"))
                # {
                #   psmo <- try(summary(cdfestimate))
                #   if (!inherits(psmo, "try-error"))
                #   {
                #     if (psmo$coefficients[1,1] < 0)
                #     {
                #       fro <- 0.125*fro+0.875*psmo$coefficients[1,1];
                #       fto <- 0.5*fto+0.5*psmo$coefficients[2,1];
                #     }
                #   }
                # }
                nleft <- 1.0-logisticcdf(allData$days[olastObs], fro, fto);
                if (allData$fatalities[olastObs] > 0.5) 
                {
                  nleft <- 1.0 - nleft;
                  adjust <- (0.01 + nleft)/(0.01 + pLeft);
                }
                else
                {
                  adjust <- (0.01 + pLeft)/(0.01 + nleft);
                }
                if (adjust > 1.25)
                {
                  adjust <- 1.25;
                }
                if (adjust < 0.9)
                {
                  adjust <- 0.9;
                }
                if ((accAdjust >= 0.9) || (adjust > 1.0))
                {
                  accAdjust <- accAdjust*adjust;
                  allData$fatalities <- adjust*allData$fatalities;
                  allData$newfatalities <- adjust*allData$newfatalities;
#                  cat(sprintf("%8.5f : (%5.3f,%5.3f) Adjust= (%5.3f,%5.3f), ro= %8.3f to=%8.3f",abs((adjusto - accAdjust)/accAdjust),pLeft,nleft,accAdjust,abs(log(accAdjust)),fro,fto),"\n")
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
    }
    print(accAdjust)
    adjusto <- 0.0;
    oro <- fro;
    oto <- fto;
    if (is.null(alpha))
    {
      lro <- fro
      lto <- fto
    }
    padjust <- accAdjust
    opLeft <- allData$fatalities[olastObs];
    if (opLeft > 0)
    {
      defecitRatio <- accAdjust*logisticcdf(allData$days[olastObs], oro, oto,alpha,lro,lto)/opLeft;
    }
    if (length(daysrange) < 0.85*olastObs)
    {
#      accAdjust <- adjini;
      pLeft <- 1.0 - allData$fatalities[olastObs];
#      cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f ",pLeft,fro,fto),"\n")
      oro <- 1.5*fro;
      oto <- 0.5*fto;
      if (is.null(alpha))
      {
        ffro <- oro;
        ffto <- oto;
        fro <- 1.5*fro;
        fto <- min(2.0*fto,olastObs);
        alphao <- 0.15+0.85*allData$fatalities[daysrange[1]]
        alpha <- alphao;
      }
      else
      {
        alphao <- alpha
        ffro <- ro;
        ffto <- to;
        fro <- lro;
        fto <- lto;
      }
#      gratio = gratio/2.0;
#      print(gratio)
      totchanges <- 1.0;
      slopeDeltaRatio <- 100.0;
      while ((totchanges > 0.0001) || ((abs((adjusto - accAdjust)/accAdjust) >= 0.005) && (abs(log(accAdjust)) <= abs(log(gratio)))))
      {
        aalpha <- alpha
        afro <- ffro;
        afto <- ffto;
        aro <- fro;
        ato <- fto;
        oro <- ffro
        oto <- ffto
        lro <- fro
        lto <- fto
        totchanges <- 0;
#        alpha <- 0.8*alpha + 0.15*allData$fatalities[daysrange[1]]+0.05*alphao;
        loops <- loops + 1;
        pLeft <- 1.0 - allData$fatalities[olastObs];
        if (pLeft < 0.5) 
        {
          pLeft <- 1.0 - pLeft;
        }
#        cat(sprintf("%8.5f : pLeft= %5.3f, ro= %8.3f to=%8.3f, lro= %8.3f lto=%8.3f",abs((adjusto - accAdjust)/accAdjust),pLeft,ffro,ffto,fro,fto),"\n")
#        adjusto <- accAdjust;
#        print(firsthalf$days)
        firsthalf <- rbind(firstdays,allData[-daysrange,],lasdays)
        cdfestimate <- try(nls(newfatalities ~ logisticpdf(days, oro, oto,alpha,lro,lto),
                               data = firsthalf,
                               start=list(oro = ffro,oto = ffto),
                               algorithm="plinear",
#                               control = list(maxiter = 150, warnOnly = TRUE)
                              control = list(maxiter = 150)
        ));
#        pdfestimate <- cdfestimate
        if (!inherits(cdfestimate, "try-error"))
        {
          smo <- try(summary(cdfestimate))
          if (!inherits(smo, "try-error"))
          {
            if (smo$coefficients[1,1] < 0)
            {
              adjusto <- 1.1*accAdjust
              error <- FALSE;
              ffro <- 0.5*ffro+0.5*smo$coefficients[1,1];
              ffto <- 0.5*ffto+0.5*smo$coefficients[2,1];
#              print(c(smo$coefficients[1,1],smo$coefficients[2,1]))
              
              if (ffro > -0.0001)
              {
                ffro <- -0.0001;
              }
              if (ffto < 0.10*to)
              {
                ffto <- 0.10*to;
              }
              if (ffto > 10*to)
              {
                ffto <- 10*to;
              }
              #               cat(sprintf("pLeft= %5.3f, ro= %8.3f to=%8.3f ",pLeft,fro,fto),"\n")
              fdata <- rbind(firstdays,allData[daysrange,],lasdays)
              pdfestimate <- try(nls(newfatalities ~ logisticpdf(days, ffro, ffto,alpha,lro,lto),
                                     data = fdata,
                                     start=list(lro = fro ,lto = fto),
                                     algorithm="plinear",
#                                     control = list(maxiter = 150, warnOnly = TRUE)
                                      control = list(maxiter = 150)
              ))
              if (!inherits(pdfestimate, "try-error"))
              {
                psmo <- try(summary(pdfestimate))
                if (!inherits(psmo, "try-error"))
                {
                  if (psmo$coefficients[1,1] < 0)
                  {
                    fro <- 0.75*fro+0.25*psmo$coefficients[1,1];
                    fto <- 0.5*fto+0.5*psmo$coefficients[2,1];
                    fdata <- rbind(firstdays,allData,lasdays)
                    cdfestimate <- try(nls(fatalities ~ logisticcdf(days, ffro, ffto,alpha,fro,fto),
                                           data = fdata,
                                           start=list(alpha=alpha,fro=fro,fto=fto),
                                           algorithm="plinear",
                                           control = list(maxiter = 150, warnOnly = TRUE)
                    ))
                    if (!inherits(cdfestimate, "try-error"))
                    {
                      psmo <- try(summary(cdfestimate))
                      if (!inherits(psmo, "try-error"))
                      {
                        alpha = 0.750*alpha + 0.250*psmo$coefficients[1,1]
                        fro = 0.750*fro + 0.250*psmo$coefficients[2,1]
                        fto = 0.750*fto + 0.250*psmo$coefficients[3,1]
#                        print(c(alpha,psmo$coefficients[1,1]))
                      }
                    }
                  }
                }
              }
              nleft <- 1.0-logisticcdf(allData$days[olastObs], ffro, ffto,alpha,fro,fto);
              if (allData$fatalities[olastObs] > 0.5) 
              {
                nleft <- 1.0 - nleft;
                adjust <- (0.01 + nleft)/(0.01 + pLeft);
              }
              else
              {
                adjust <- (0.01 + pLeft)/(0.01 + nleft);
              }
              if (adjust > 1.1)
              {
                adjust <- 1.1;
              }
              if (adjust < 0.9)
              {
                adjust <- 0.9;
              }
              adjusto <- accAdjust;
              totchanges <- abs(ffro-afro) + abs(ffto-afto) + abs(alpha-aalpha) + abs(fro-aro) + abs(fto-ato);
              if (((accAdjust >= 0.95) || (adjust > 1.0)) && (abs(log(accAdjust)) <= abs(log(gratio))))
              {
                accAdjust <- accAdjust*adjust;
                allData$fatalities <- adjust*allData$fatalities;
                allData$newfatalities <- adjust*allData$newfatalities;
              }
              slopeDeltaRatioa <- slopeDeltaRatio
              slopeDeltaRatio <- abs(allData$newfatalities[olastObs]-logisticpdf(allData$days[olastObs],ffro, ffto,alpha,fro,fto))/(allData$newfatalities[olastObs]+0.0001)
              if (slopeDeltaRatio < 0.01)
              {
                totchanges <- 0;
                adjusto <- accAdjust;
              }
#              cat(loops,":",slopeDeltaRatio, sprintf(" %8.5f : (%5.3f,%5.3f) Adjust= (%5.3f,%5.3f), alpha= %8.3f ro = %8.3f to=%8.3f (lro= %8.3f lto=%8.3f) %8.6f",abs((adjusto - accAdjust)/accAdjust),pLeft,nleft,accAdjust,abs(log(accAdjust)),alpha,ffro,ffto,fro,fto,totchanges),"\n")
            }
          }
        }
        if (loops > 100)
        {
          adjusto <- accAdjust
          totchanges <- 0
        }
        if (loops == 1)
        {
          firscdftestimation <- cdfestimate;
          firspdftestimation <- pdfestimate;
        }
      }
      oro <- ffro
      oto <- ffto
      lro <- fro
      lto <- fto
      opLeft <- allData$fatalities[olastObs];

      if (opLeft > 0)
      {
        defecitRatio <- accAdjust*logisticcdf(allData$days[olastObs], oro, oto,alpha,lro,lto)/opLeft;
      }
    }
#    cat(olastObs,":",defecitRatio,":",logisticcdf(allData$days[olastObs], oro, oto,alpha,lro,lto),sprintf("%5.3f, pLeft= %5.3f, ro= %8.3f to=%8.3f alpha %8.3f ,lro= %8.3f lto=%8.3f ",opLeft,pLeft,fro,fto,alpha,lro,lto),"\n")
  }
  
#   cat(defecitRatio,sprintf(" oLeft= %5.3f, ro= %8.3f to=%8.3f ",opLeft,fro,fto),"\n")
  
#  cat(sprintf("%5.3f,Defecit= %5.3f, alpha = %8.3f ro= %8.3f to=%8.3f, lro=%8.3f,lto=%8.3f,",accAdjust,defecitRatio,alpha,oro,oto,lro,lto),"\n")

  models <- list(CDF=cdfestimate,
                 pdf=pdfestimate,
                 ro=oro,
                 to=oto,
                 lro=lro,
                 lto=lto,
                 alpha=alpha,
                 firstcdf = firscdftestimation,
                 firstpdf = firspdftestimation,
                 defecitRatio=defecitRatio,
                 adjust = accAdjust,
                 filterCDF = filterCDF,
                 filterpdf = filtercdf
  )
  class(models) <- "logisticFit";
#  print(error)
  if (error) class(models) <- append(class(models),"try-error")
  return (models);
}

bootstraplogisitcfit <- function(data,inifit,ratiorange=1.5,n=200,daysrange=c(1:nrow(data)))
{
  toestimations <- numeric(n);
  roestimations <- numeric(n);
  ltoestimations <- numeric(n);
  lroestimations <- numeric(n);
  alphaestimations <- numeric(n);
  defratios <- numeric(n);
#  print(c(inifit$ro,inifit$to))
  optGain <- ratiorange^(runif(n, -1,1));

  for (bsamples in 1:n)
  {
    toestimations[bsamples] <- inifit$to;
    roestimations[bsamples] <- inifit$ro
    ltoestimations[bsamples] <- inifit$to;
    lroestimations[bsamples] <- inifit$ro;
    if (!is.null(inifit$alpha))
    {
      ltoestimations[bsamples] <- inifit$lto;
      lroestimations[bsamples] <- inifit$lro;
    }
    defratios[bsamples] <- inifit$defecitRatio;
    bootresample <- sample(nrow(data),nrow(data)-5,TRUE)
    bootresample <- c(1:5,bootresample[order(bootresample)]);
#    print(bootresample)
    ndata <- data[bootresample,];
#    ndata <- data;
    maxgain <- min(2.0/max(ndata$fatalities),optGain[bsamples]);
    ndata$fatalities <- maxgain*ndata$fatalities;
    ndata$newfatalities <- maxgain*ndata$newfatalities;
    iniadjs <- inifit$adjust;
#    iniadjs <- 1.0;
    nfit <- try(logisitcfit(ndata,inifit$ro,inifit$to,1.25*iniadjs,iniadjs,daysrange=daysrange,lro=inifit$lro,lto=inifit$lto,alpha=inifit$alpha))
    if (!inherits(nfit, "try-error"))
    {
      if ((nfit$ro < 0) && (nfit$to > 0) && (nfit$to < 1000))
      {
        cat(".")
        toestimations[bsamples] <- nfit$to;
        roestimations[bsamples] <- nfit$ro;
        defratios[bsamples] <- nfit$defecitRatio*maxgain;
        if (!is.null(nfit$alpha))
        {
          ltoestimations[bsamples] <- nfit$lto;
          lroestimations[bsamples] <- nfit$lro;
          alphaestimations[bsamples] <- nfit$alpha;
        }
        else
        {
          ltoestimations[bsamples] <- nfit$to;
          lroestimations[bsamples] <- nfit$ro;
        }
      }
    }
  }
  result <- list(to=toestimations,ro=roestimations,lto=ltoestimations,lro=lroestimations,alpha=alphaestimations,defratios=defratios,gains=optGain)
  return(result)
}
