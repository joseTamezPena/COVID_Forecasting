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
    
    mxchange <- 0.8*max(datsetchange);
    lastday <- as.integer(mean(c(1:lastobs)[(datsetchange > mxchange)]) + 0.5) - 3;
    peakDate <- 0.5*(lastday+ which.max(datsetchange))
    print(c(lastobs,lastday))
    if (lastday < mindaytoinclude)
    {
#      print(c(lastobs,lastday,lastobs-lastday))
      lastday <- mindaytoinclude
    }
    print(c(lastobs,lastday,lastobs-lastday))
  #    print(lastday)
    
    daysrange <- c(1:lastday)
    cdffitglobal <- logisitcfit(datasetperc[c(1:lastday),],-0.075,peakDate,1.5,daysrange=daysrange)
#    adjsini <- 0.5*(cdffitglobal$adjust + 1.0)
#    cdffitglobal <- logisitcfit(datasetperc[c(1:lastday),],cdffitglobal$ro,cdffitglobal$to,2.50*adjsini,adjini=adjsini,daysrange=daysrange)
    
    daysrange <- c((lastday - 17):lastday)
    cat(cdffitglobal$adjust,"\n")
    adjsini <- 0.5*(cdffitglobal$adjust + 1.0)
    cat(adjsini,"\n")
    if (peakDate < (lastday-17))
    {
      daysrange <- c((peakDate-17):lastday)
    }
#    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,1.2*cdffitglobal$adjust,adjini=cdffitglobal$adjust,daysrange=daysrange)
    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,2.50*adjsini,adjini=adjsini,daysrange=daysrange)
    #    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,3,daysrange=daysrange)
    cat(cdffit$adjust,"\n")
    
    lines(cdffitglobal$filterpdf*totalEsperado,lty=1,col="red",lwd=3)
    lines(logisticpdf(c(1:maxplotdata),cdffitglobal$ro,cdffitglobal$to)/cdffitglobal$defecitRatio*totalEsperado,lty=5,col="blue",lwd=1)    
    lines(logisticpdf(c(1:maxplotdata),cdffit$ro,cdffit$to)/cdffit$defecitRatio*totalEsperado,lty=3,col="green",lwd=3)    
    
    legend("topleft",
           legend = c("Observado","Suavizado","Global Estimation","Last 17 Days"),
           col = c(1,"red","blue","green"),
           pch = c(1,NA,NA,NA),
           lty = c(NA,1,5,3),
           lwd = c(NA,3,1,3),
           cex=0.5)
    
        hostCDFIT <- logisticcdf(c(1:maxplotdata),cdffit$ro,cdffit$to)/cdffit$defecitRatio
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
    lines(cdffitglobal$filterCDF,lty=1,col="red")
    pdfpredictions <- logisticpdf(c(1:maxplotdata),cdffit$ro,cdffit$to)/cdffit$defecitRatio;
    lines(30*pdfpredictions,lty=2,col="blue")
    atpeak <- as.integer(cdffit$to+0.5)
    estimadoActual <- totalEsperado*sum(pdfpredictions[1:lastobs]);
    abline(v=atpeak,col="pink",lty=2)
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
    
    bootest <- bootstraplogisitcfit(datasetperc,cdffit,n=500,daysrange=daysrange)
    torange <- quantile(bootest$to,probs = c(0.025,0.5,0.975))
    daystopeakrange <- startDate + as.integer(torange+0.5);
    text(1,0.55,
         paste("Día pico entre",as.character(daystopeakrange[1]),"y",as.character(daystopeakrange[3])),
         pos=4,
         cex=0.5)
    
    
    axis(1,at=c(1:length(dateAxis)),labels=dateAxis,las=2,cex.axis=0.4)
    
    legend("topright",
           legend = c("Total Estimado","Total Observado","Datos utilizados","Ingreso por día","linea de pico"),
           col = c(1,"gray","red","blue","pink"),
           lty = c(2,1,1,2,2),
           lwd = c(1,3,1,1,1),
           cex=0.5)
    
    z <- c(0:10)/100;
    axis(4, at=30*z,labels=round(z,digits=2),
         col.axis="blue", las=2, cex.axis=0.6, tck=-.01)
    mtext("Fracción de casos nuevos", cex=0.30,side=4, line=3, cex.lab=0.2, col="black",las=3)
    result <- list(name=mainName,daysRange=daystopeakrange,peak=max(pdfpredictions*totalEsperado))
    return (result)
    
    
}