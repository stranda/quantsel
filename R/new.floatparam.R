#probably needs a good bit of error checking written in
#
landscape.new.floatparam <-
function(rland,s=0,seedscale=c(10,10),seedshape=c(10,10),seedmix=1,
         pollenscale=c(2,10),pollenshape=c(2,6),pollenmix=1,asp=1,
         mindens=1e-25)
{
  rland$floatparam <- list(selfing=s,seedmu=seedscale[1],seedshape=seedshape[1],
                           seedmu2=seedscale[2],seedshape2=seedshape[2],seedmix=seedmix[1],
                           pollenmu=pollenscale[1],pollenshape=pollenshape[1],
                           pollenmu2=pollenscale[2],pollenshape2=pollenshape[2],
                           pollenmix=pollenmix,aspect=asp,mindens=mindens)
  rland
}

