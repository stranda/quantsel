"landscape.survive" <-
function(Rland, seed=-1)
  {
    if (is.landscape(Rland))
      {
        if (!(seed<0))
          {
            set.seed(seed)
          }
        Rland=landscape.coerce(Rland)
        .Call("survive_landscape",Rland,PACKAGE = "quantsel")
      }
    else
      {
        print("Rland not a landscape object...exiting")
      }
  }

