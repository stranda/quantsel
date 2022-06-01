"landscape.read" <-
function(fn = "filename")
  {
    if (file.exists(fn))
      {
        .Call("read_landscape",fn,PACKAGE = "quantsel")
      }
    else
      {
        print (paste("Filename: ",fn,"does not appear to exist"))
        NULL
      }
  }

