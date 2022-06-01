l=library(quantsel)
#library(holoSimCell)
library(ggplot2)
library(dplyr)


###the following four lines read in an environment containing a holosim suitability surface
### and then they rename the one object to 'hsl' and delete the object (they are stored in
### holosimcell as 'landscape' and that is a confusing name
### there will be a problem if an object called landscape already exisits
surf='naive.rda'
load(paste0(system.file(package="holoSimCell"),"/extdata/landscapes/",surf))
hsl=landscape
rm(landscape)



rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=hsl$details$ncells, s=2, totgen=21000,maxland=3e6)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,seedscale=c(40,1500),
                                  seedshape=c(1,500),seedmix=c(0.1),
                                  pollenscale=c(40,1250),pollenshape=c(1,100),
                                  pollenmix=0.2 , asp=0.5)


S <- matrix(c(
    0.5   ,   0,
    0.15,   0.85
), byrow=T, nrow = 2)
  R <- matrix(c(
      0,  4,
      0,   0
  ), byrow=T, nrow = 2)
M <- matrix(c(
      0, 0,
      0, 1
), byrow=T, nrow = 2)

rland <- landscape.new.local.demo(rland,S,R,M)

S <- matrix(0,ncol = (rland$intparam$habitats*rland$intparam$stages),
            nrow = (rland$intparam$habitats*rland$intparam$stages))

R <- S
M <- S

k=calc.k.vec(200, hsl,1)

e=rep(0,rland$intparam$habitat)

locs=raster2popcrds(hsl$sumrast)

rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                             carry=k,
                             extinct=e,
                             leftx=locs[,2],
                             rightx=locs[,3],
                             boty=locs[,5],
                             topy=locs[,4],
                             maxland=c(min(locs[,2]),min(locs[,5]),max(locs[,3]),max(locs[,4])))

nl = 10
for (i in 1:nl)
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)

#rland <- landscape.nophen(rland)
initpopsize <- 100
inits <- matrix(50,ncol=rland$intparam$habitats,nrow=2)
inits[1,1233] <- initpopsize
inits[2,1233] <- initpopsize

rland <- landscape.carry(landscape.new.individuals(rland,c(inits)))

print(dim(rland$individuals))

l=landscape.simulate(rland,1)


gen=50


for (i in 1:gen)
{
    print(dim(l$individuals))
    if (dim(l$individuals)[1]>0) l.old=l
    l=landscape.simulate(l,1)
    print(i)
    print(dim(l$individuals))
#    print(landscape.allelefreq(l) )
}
#dev.off()

