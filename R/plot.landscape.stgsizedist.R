"landscape.plot.stgsizedist" <-
function(Rland)
  {
    if (is.landscape(Rland))
      {
        plot(table(Rland$individuals[,1]+1),
                   main=c("Frequency distribution of demographic stage sizes"),
                   xlab=c("Demographic stage"),
                   ylab=c("Number of individuals"),
             )
      }
  }

