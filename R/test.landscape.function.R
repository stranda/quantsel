#' tests whatever c code is compiled in
#' 
#' currently generates 10,000 random pulls from a pdf
#' 
#' 
#' @return vector
#' @keywords misc
#' @export landscape.test.function
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
