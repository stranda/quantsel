"landscape.plot.popsizedist" <-
function(Rland)
  {
    if (is.landscape(Rland))
      {
        plot(table(landscape.populations(Rland)),
             main="Frequency distribution of population sizes",
             xlab=c("Population"),
             ylab=c("Number of individuals"),
             )
      }
  }

