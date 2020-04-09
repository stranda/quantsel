"landscape.new.switchparam" <-
function(rland, re=0,rd=0,mp=1,sc=1)
{
  rland$switchparam <- list(re,rd,mp)
  names(rland$switchparam) <- c( "randepoch","randdemo","multp")
  rland
}

