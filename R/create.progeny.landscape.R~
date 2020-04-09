"landscape.create.progeny" <-
function(Rland,gen=1,seed=-1)
  {
    if (seed<0)
      {
        seed <- 100000 * runif(1);
      }
    if (is.landscape(Rland))
      {
        if (gen>4) gen<-4;
        .Call("create_progeny",Rland,gen,seed,PACKAGE = "kernelPop2")
      }
  }

