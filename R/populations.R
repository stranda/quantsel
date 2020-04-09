"landscape.populations" <-
function(Rland)
  {
    if (is.landscape(Rland,verb=F))
      {
        (Rland$individuals[,1]%/%Rland$intparam$stages)+1
      }
  }

