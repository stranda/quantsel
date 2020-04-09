"landscape.test.function" <-
function()
  {
    .Call("test",PACKAGE = "quantsel")
  }

landscape.l2w <-
    function(Rland,numind)
{
    .Call("l2w",Rland,numind,PACKAGE = "quantsel")
}
