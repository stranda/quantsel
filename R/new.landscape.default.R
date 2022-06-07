"landscape.new.default" <-
function()
{
  rland <- landscape.new.empty()
  rland <- landscape.new.intparam(rland)
  rland <- landscape.new.floatparam(rland)
  rland <- landscape.new.switchparam(rland)
  rland
}

