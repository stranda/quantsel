library(kernelPop2)
library(ggplot2)
library(dplyr)
source("helpers.R")
source("analysis.R")

### this script maks a 1024 population grid and 
### populates the entire thing, then three traits evolve


gapprop = 0



rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=2, s=2,np=0,totgen=20000,maxland=3e5)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,seedscale=c(40,290),
                                  seedshape=c(1,300),seedmix=c(0.12),
                                  pollenscale=c(40,100),pollenshape=c(1,10),
                                  pollenmix=0.2 , asp=0.5)


S <- matrix(c(
    0,     0,
    0.8, 0.0
), byrow=T, nrow = 2)
  R <- matrix(c(
      0,  12,
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


locs=cbind(lft=c(1,20001),bot=c(1,1),rgt=c(20000,40000),top=c(20000,20000))

lfts=1
k=round((0.06 * (sqrt((locs[,3]-locs[,1])*(locs[,4]-locs[,2])))))
e=rep(gapprop,rland$intparam$habitat)
rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                             carry=k,
                             extinct=e,
                             leftx=locs[,1],
                             rightx=locs[,3],
                             boty=locs[,2],
                             topy=locs[,4],
                             maxland=c(min(locs[1]),min(locs[2]),max(locs[3]),max(locs[4])))


for (i in 1:16)
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)


expmat <- matrix(c(
    0.5,0,0.5,0,
    0.5,0,0.5,0,
    0.5,0,0.5,0,
    0.5,0,0,0.5,
    0,1,0,0,
    0,1,0,0,
    0,1,0,0,
    0,1,0,0,
    0.5,0,0.5,0,
    0.5,0,0.5,0,
    0.5,0,0.5,0,
    0,0,1,0,
    0,0,0,1,
    0,0,0,1,
    0.5,0,0,0.5,
    0.5,0,0,0.5
                   ),byrow=T,ncol=4)
hsq <- c(1,1,1,1)
rland <- landscape.new.expression(rland,expmat=expmat*0.125,hsq=hsq) #0.125 -> 4 diploid loci, up to 8 alelle additive doses
rland <- landscape.new.gpmap(rland,
                             ## 4 cols 5 rows.  Cols correspond to
                             ##phenotype (-1 is none), curvature, range of effect on vital rate
                             ##phenotypes are in C indexing so, add 1 to compare to pehnotypes above
                             matrix(c(-1,   0,    .1, 1,  #short scale #no selection
                                       1, 0.01, .1, 1,#long scale  #no selection
                                      -1,   0,    .1, 1, #long shape
                                       -1,   1,    .1, 1,   #mixture   #phenotype 2
                                      -1,   0,    .1, 1 ),
                                    ncol=4,byrow=T),
                             
                             matrix(c( 2,   0.01,  0.2, -1, #survival
                                      0,     0.01,  0.2,  1, #reproduction #phenotype 1
                                       3,       1,  0.005,  -1), #density tolerance 
                                    ncol=4,byrow=T))


initpopsize <- k[1]/2
inits <- matrix(initpopsize,ncol=rland$intparam$habitats,nrow=2)

rland <- landscape.new.plasticity(rland)
rland <- landscape.new.phenohab(rland)

rland <- landscape.new.individuals(rland,c(inits))


if (FALSE)
    {
rland$individuals[landscape.populations(rland)==pop1,c(10,11,12,13)]=1L
rland$individuals[landscape.populations(rland)==pop2,c(10,11,12,13)]=2L

rland$individuals[landscape.populations(rland)==pop1,c(14,15,16,17)]=2L
rland$individuals[landscape.populations(rland)==pop2,c(14,15,16,17)]=1L
}
#rland$individuals[,5] <- 4500+floor(rland$individuals[,5]/10)
###############################

print(names(rland))

rland$plasticity[1,3]=1.1
rland$plasticity[2,1]=1.2

l=landscape.simulate(rland,1)

landscape.plot.phenotypes(l,1)

print(names(l))

print(l$plasticity)
print(l$phenohab)
