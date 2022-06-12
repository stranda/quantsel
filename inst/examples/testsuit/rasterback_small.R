library(quantsel)
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

habdiagonal = 144935



rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=hsl$details$ncells, s=2, totgen=21000,maxland=3e6)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,
                                  seedscale=c(0.0001,0.06)*habdiagonal,
                                  seedshape=c(1,habdiagonal*0.1),seedmix=c(0.001),
                                  pollenscale=c(0.001,0.08)*habdiagonal,
                                  pollenshape=c(1,habdiagonal*0.05),
                                  pollenmix=0.2 , asp=1)


S <- matrix(c(
    0.07   ,   0,
    0.03,   0.85
), byrow=T, nrow = 2)

R <- matrix(c(
      0,   15,
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

carry=600

k=rep(carry,rland$intparam$habitat) #need some sort of carrying capacity to initalize.  replaced in first gen

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
initpopsize <- 1000
inits <- matrix(0,ncol=rland$intparam$habitats,nrow=2)
refuge="ALL"

  if(refuge == "PA") {
    pops <- c(999, 1000, 998, 1050, 948)
  } else if(refuge == "TX") {
    pops <- c(525, 526, 524, 576, 474)
  } else if(refuge == "GA") {
    pops <- c(535, 536, 534, 586, 484)
  } else if(refuge == "ALL") {
    pops <- c(999, 1000, 998, 1050, 948, 525, 526, 524, 576, 474,535, 536, 534, 586, 484 )
  }

#pops=c(167,166,165,114,115,116,216,217,218) + 1
inits[1,pops] <- initpopsize
inits[2,pops] <- initpopsize

l <- landscape.carry(landscape.new.individuals(rland,c(inits)))

pdf("landscape_small.pdf",width=20,height=20)
landscape.plot.locations(l,label=T)

rasterlayers=c(1:10) #low numbers are early, 1:33 simulates 990 years
genperlayer=30

for (i in rasterlayers)
{
    ##change the carrying capacitites every change in rasterlayer
    l$demography$epochs[[1]]$Carry = calc.k.vec(carry, hsl,i)
    print(dim(l$individuals))
    l=landscape.simulate(l,genperlayer)
    saveRDS(file="landscape-latest.RDS", l)
    landscape.plot.locations(l,label=T)
}

dev.off()
