plotCovid <- function(data,feature,mainName,totalEsperado,startDate,currentdate)
{
    dataPerc <-  data/totalEsperado;

    lastobs <- length(dataPerc)
    datsetchange <- c(dataPerc[1],dataPerc[2:lastobs]-dataPerc[1:(lastobs-1)])
    datasetperc <- as.data.frame(cbind(days=c(1:lastobs),fatalities = dataPerc,newfatalities = datsetchange))
    
    startDate <- as.Date(startDate)
    lastDataobs <- startDate + lastobs - 1;
    mindaytoinclude <- lastDataobs - 14;
    mindaytoinclude <- as.numeric(mindaytoinclude - startDate) + 1
    print(mindaytoinclude)
    #    print(startDate)
    maxplotdata <- 3*lastobs;
    dateAxis <- c(startDate,startDate+1,c((startDate+1):(startDate - 1 + maxplotdata)))
#    print(dateAxis)
    
    mxchange <- 0.8*max(datsetchange);
    lastday <- as.integer(mean(c(1:lastobs)[(datsetchange > mxchange)]) + 0.5) - 3;
    print(c(lastobs,lastday))
    if (lastday < mindaytoinclude)
    {
#      print(c(lastobs,lastday,lastobs-lastday))
      lastday <- mindaytoinclude
    }
    print(c(lastobs,lastday,lastobs-lastday))
  #    print(lastday)
    
    daysrange <- c(1:lastday)
    cdffitglobal <- logisitcfit(datasetperc[c(1:lastday),],-0.075,lastobs,2,daysrange=daysrange)
    daysrange <- c((lastday - 21):lastday)
    adjsini <- 0.5*(cdffitglobal$adjust + 1.0)
#    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,1.2*cdffitglobal$adjust,adjini=cdffitglobal$adjust,daysrange=daysrange)
    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,1.25*adjsini,adjini=adjsini,daysrange=daysrange)
    #    cdffit <- logisitcfit(datasetperc,cdffitglobal$ro,cdffitglobal$to,3,daysrange=daysrange)
    
    hostCDFIT <- logisticcdf(c(1:maxplotdata),cdffit$ro,cdffit$to)/cdffit$defecitRatio
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
    
    abline(v=atpeak,col="pink",lty=2)
    maxhostpital <- startDate + atpeak
    
    text(1,0.75,
         paste("Total de",feature,":",sprintf("%6.0f",max(data))),
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
    
    
}