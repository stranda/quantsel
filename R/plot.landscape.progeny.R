"landscape.plot.progeny" <-
function(Rland)
  {
    if (is.landscape(Rland))
      {
        parents <- Rland$individuals[Rland$individuals[,3]==0,1:9]
        progeny <- Rland$individuals[Rland$individuals[,3]==1,1:9]
        ymin <- min(Rland$individuals[,5])
        ymax <- max(Rland$individuals[,5])
        xmin <- min(Rland$individuals[,4])
        xmax <- max(Rland$individuals[,4])
        plot(parents[,4],parents[,5],ylim=c(ymin,ymax),xlim=c(xmin,xmax),pch=19,xlab="x",ylab="y")
        points(progeny[,4],progeny[,5],pch=21)
      }
  }

