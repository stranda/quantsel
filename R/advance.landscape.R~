"landscape.advance" <-
function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland=landscape.coerce(Rland)
        .Call("advance_landscape",Rland,PACKAGE = "kernelPop2")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

