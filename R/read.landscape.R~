"landscape.read" <-
function(fn = "filename")
  {
    if (file.exists(fn))
      {
        .Call("read_landscape",fn,PACKAGE = "kernelPop2")
      }
    else
      {
        print (paste("Filename: ",fn,"does not appear to exist"))
        NULL
      }
  }

