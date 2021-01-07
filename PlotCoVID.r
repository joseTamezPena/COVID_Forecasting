plotCovid <- function(data,feature,mainName,totalEsperado,startDate,currentdate)
{
    dataPerc <-  data/totalEsperado;

    lastobs <- length(dataPerc)
    datsetchange <- c(dataPerc[1],dataPerc[2:lastobs]-dataPerc[1:(lastobs-1)])
    datasetperc <- as.data.frame(cbind(days=c(1:lastobs),fatalities = dataPerc,newfatalities = datsetchange))

    colorcode <- rep("gray",lastobs);
    colorcode[1:(lastobs-14)] <- "Black";
    plot(datsetchange*totalEsperado,col=colorcode,
         main=mainName,
         xlab="Día desde primer evento",
         ylab="Eventos Observados",
         cex.axis=0.75,
         cex.lab=0.65
    );
    
    sdatsetchange <- runmed(datsetchange,5)

    
    mindaytoincludeID <- which.max(datsetchange);
        
    startDate <- as.Date(startDate)
    lastDataobs <- startDate + lastobs - 1;
    mindaytoinclude <- lastDataobs - 14;
    mindaytoinclude <- max(mindaytoincludeID,as.numeric(mindaytoinclude - startDate) + 1);
    print(mindaytoinclude)
    #    print(startDate)
    maxplotdata <- 3*lastobs;
    dateAxis <- c(startDate,startDate+c(1:maxplotdata-1))
#    print(dateAxis)
    
    mxchange <- 0.8*max(sdatsetchange);
    lastday <- as.integer(mean(c(1:lastobs)[(sdatsetchange > mxchange)]) + 0.5) - 3;
#    lastday <- lastobs-14;
#    peakDate <- 0.25*lastday + 0.75*which.max(sdatsetchange)
    peakDate <- as.integer(mean(c(1:lastobs)[(sdatsetchange > mxchange)]) + 0.5)
    print(c(lastobs,lastday,peakDate))
    if (lastday < mindaytoinclude)
    {
#      print(c(lastobs,lastday,lastobs-lastday))
      lastday <- mindaytoinclude
    }
    print(c(lastobs,lastday,peakDate,lastobs-lastday))
    print(lastday)
    
    stdate <- max(1,peakDate-2*lastobs/3)
    daysrange <- c(1:lastday)
    print(daysrange)
    cdffitglobalo <- logisitcfit(datasetperc[c(1:lastday),],-0.035,peakDate,10,daysrange=daysrange)
    daysrange <- c(stdate:min(lastday,peakDate+62))
    print(daysrange)
    cdffitglobal <- logisitcfit(datasetperc[c(1:lastday),],cdffitglobalo$ro,peakDate,2,daysrange=daysrange)
    lines(cdffitglobal$filterpdf*totalEsperado,lty=1,col="red",lwd=3)
    #    adjsini <- 0.5*(cdffitglobal$adjust + 1.0)
#    cdffitglobal <- logisitcfit(datasetperc[c(1:lastday),],cdffitglobal$ro,cdffitglobal$to,2.50*adjsini,adjini=adjsini,daysrange=daysrange)
    
    cat(cdffitglobal$adjust,"\n")
    adjsini <- cdffitglobal$adjust
    cat(adjsini,"\n")
    roo <- cdffitglobal$ro
    daysrange <- c(as.integer(3*lastday/4):lastday)
    if (cdffitglobal$to < (lastobs-24))
    {
      peakDate <- 0.5*(cdffitglobal$to + which.max(cdffitglobal$filterpdf))
      if ((datsetchange[lastday] > 0.001) && (datsetchange[lastday] > 0.25*datsetchange[min(c(0.5*lastday,as.integer(cdffitglobal$to)))]) )
      {
        cat("Tail")
        daysrange <- c(max(1,as.integer(cdffitglobal$to-60),as.integer(3*lastday/4)):lastday)
        peakDate <- 0.15*cdffitglobalo$to + 0.85*lastobs
        roo <-  0.75*cdffitglobalo$ro
        adjsini <- 0.95*cdffitglobalo$adjust
      }
      else
      {
        daysrange <- c(max(c(1,(cdffitglobal$to-60))):lastday)
      }
    }
    print(daysrange)
    #    lastday <- (lastobs-14)
    print(c(adjsini,roo,peakDate))
    cdffit <- cdffitglobalo;
#    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,1.2*cdffitglobal$adjust,adjini=cdffitglobal$adjust,daysrange=daysrange)
    cdffit <- logisitcfit(datasetperc[c(1:lastday),],roo,peakDate,1.10*adjsini,adjini=0.75*adjsini,daysrange=daysrange)
    #    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,3,daysrange=daysrange)
    cat(cdffit$adjust,"\n")
    print(c(cdffit$to,cdffit$ro,cdffit$lto,cdffit$lro,cdffit$alpha))
    
    lines(cdffitglobalo$filterpdf*totalEsperado,lty=1,col="pink",lwd=3)
    lines(logisticpdf(c(1:maxplotdata),cdffitglobal$ro,cdffitglobal$to)/cdffitglobal$defecitRatio*totalEsperado,lty=5,col="blue",lwd=1)

    lines(logisticpdf(c(1:maxplotdata),cdffit$ro,cdffit$to,cdffit$alpha,cdffit$lro,cdffit$lto)/cdffit$defecitRatio*totalEsperado,lty=3,col="green",lwd=3)    
    
    legend("topleft",
           legend = c("Observado","Suavizado","Estimación Global","Ultima Estimación"),
           col = c(1,"red","blue","green"),
           pch = c(1,NA,NA,NA),
           lty = c(NA,1,5,3),
           lwd = c(NA,3,1,3),
           cex=0.5)
    
        hostCDFIT <- logisticcdf(c(1:maxplotdata),cdffit$ro,cdffit$to,cdffit$alpha,cdffit$lro,cdffit$lto)/cdffit$defecitRatio
#    hostCDFIT <- logisticcdf(c(1:maxplotdata),cdffit$ro,cdffit$to)
    plot(hostCDFIT,
         ylim=c(0,1),
         type="l",
         lty=2,
         ylab=paste("Fracción esperada de",feature),
         xaxt="none",
         xlab="",
         main=mainName,
         sub=as.character(currentdate),
         cex.axis=0.75,
         cex.lab=0.65
    )
    lines(dataPerc,lty=1,col="gray",lwd=3)
    lines(cdffitglobalo$filterCDF,lty=1,col="red")
    pdfpredictions <- logisticpdf(c(1:maxplotdata),cdffit$ro,cdffit$to,cdffit$alpha,cdffit$lro,cdffit$lto)/cdffit$defecitRatio;
    print(c(cdffit$ro,cdffit$to))
    lines(30*pdfpredictions,lty=2,col="blue")
    atpeak <- as.integer(cdffit$to+0.5)
    atpeak2 <- as.integer(cdffit$lto+0.5)
    estimadoActual <- totalEsperado*sum(pdfpredictions[1:lastobs]);
    abline(v=atpeak,col="pink",lty=2)
    abline(v=atpeak2,col="red",lty=2)
    maxhostpital <- startDate + atpeak
    
    text(1,0.75,
         paste("Total de",feature,":",sprintf("Reportado: %6.0f, Estimado: %6.0f",max(data),estimadoActual)),
         pos=4,
         cex=0.5)
    text(1,0.70,
         paste("Total estimado", 
               sprintf("%6.0f",totalEsperado/cdffit$defecitRatio),
               "en",
               as.character(dateAxis[maxplotdata])
         ),
         pos=4,
         cex=0.5)
    text(1,0.65,
         paste("Pico dia de nuevas",feature,":",sprintf("%5.0f",max(pdfpredictions*totalEsperado))),
         pos=4,
         cex=0.5)

    text(1,0.60,
         paste("Pico semana ",feature,":",sprintf("%5.0f",sum(pdfpredictions[(atpeak-3):(atpeak+3)])*totalEsperado)),
         pos=4,
         cex=0.5)
    
    cat("Ro:",cdffit$ro," to:",cdffit$to,"\n")
    bootsamples <- 25
    bootest <- bootstraplogisitcfit(datasetperc[c(1:lastday),],cdffit,n=bootsamples,daysrange=daysrange)

    torange <- quantile(bootest$to,probs = c(0.025,0.5,0.975))

    meanto <- median(bootest$to)
    meanro <- median(bootest$ro)
    meanlto <- median(bootest$lto)
    meanlro <- median(bootest$lro)
    meanalpha <- median(bootest$alpha)
    meandr <- median(bootest$defratios)
    cat("medro:",meanro," medto:",meanto,"\n")
    for (n in 1:bootsamples)
    {
      if ((bootest$to[n] > 17) && (bootest$lto[n] > bootest$to[n]) && (bootest$ro[n] < -0.001))
      {
        lines(c(lastobs:maxplotdata),logisticcdf(c(lastobs:maxplotdata),bootest$ro[n],bootest$to[n],bootest$alpha[n],bootest$lro[n],bootest$lto[n])/bootest$defratios[n],lwd=1,lty=8,col="light gray")
      }
      
    }
    lines(c(lastobs:maxplotdata),logisticcdf(c(lastobs:maxplotdata),meanro,meanto,meanalpha,meanlro,meanlto)/meandr,lwd=2,lty=8,col="green")
    
    daystopeakrange <- startDate + as.integer(torange+0.5);
    text(1,0.55,
         paste("Día pico entre",as.character(daystopeakrange[1]),"y",as.character(daystopeakrange[3])),
         pos=4,
         cex=0.5)
    
    
    axis(1,at=c(1:length(dateAxis)),labels=dateAxis,las=2,cex.axis=0.4)
    
    legend("topright",
           legend = c("Total Estimado","Total Observado","Datos utilizados","Ingreso por día","linea de pico","Estimado"),
           col = c(1,"gray","red","blue","pink","green"),
           lty = c(2,1,1,2,2),
           lwd = c(1,3,1,1,1),
           cex=0.5)
    
    z <- c(0:10)/100;
    axis(4, at=30*z,labels=round(z,digits=2),
         col.axis="blue", las=2, cex.axis=0.6, tck=-.01)
    mtext("Fracción de casos nuevos", cex=0.30,side=4, line=3, cex.lab=0.2, col="black",las=3)
    
    
    medianpdfpredictions <- logisticpdf(c(1:(2*maxplotdata)),meanro,meanto,meanalpha,meanlro,meanlto)/meandr;
    
    result <- list(name=mainName,
                   daysRange=daystopeakrange,
                   peak=max(pdfpredictions*totalEsperado),
                   medianPred=medianpdfpredictions,
                   lastNew=cdffitglobalo$filterpdf[lastday])
    return (result)
    
    
}