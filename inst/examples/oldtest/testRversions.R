library(kernelPop2)
library(ggplot2)
library(dplyr)
source("helpers.R")
source("newlandscape.R")
source("analysis.R")

### this script maks a 250 population grid and 
###


gapprop = 0



rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=10, s=3,np=0,totgen=20000,maxland=3e5)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,seedscale=c(40,290),
                                  seedshape=c(1,300),seedmix=c(0.12),
                                  pollenscale=c(40,100),pollenshape=c(1,10),
                                  pollenmix=0.2 , asp=0.5)


S <- matrix(c(
    0.1, 0, 0,
    0.4, 0, 0,
    0.1, 0.4, 0
), byrow=T, nrow = 3)

R <- matrix(c(
      0,  5, 9,
      0,   0, 0,
      0,   0, 0
  ), byrow=T, nrow = 3)

M <- matrix(c(
    0, 0, 0,
    0, 1, 1,
    0, 1, 1
), byrow=T, nrow = 3)

rland <- landscape.new.local.demo(rland,S,R,M)

S <- matrix(0,ncol = (rland$intparam$habitats*rland$intparam$stages),
            nrow = (rland$intparam$habitats*rland$intparam$stages))

R <- S
M <- S

rights <- floor(seq(0,10000,length=11))
tops <-   floor(seq(0,10000,length=2))

locs=NULL
for (i in 1:(length(tops)-1))
    locs <- rbind(locs,
                  data.frame(lft=c(rights[-1]-diff(rights[])+1),
                             bot=0,
                             rgt=rights[-1],
                             top=tops[i+1]))



k=rep(5000,10)
e=rep(0,10)
rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                             carry=k,
                             extinct=e,
                             leftx=locs[,1],
                             rightx=locs[,3],
                             boty=locs[,2],
                             topy=locs[,4],
                             maxland=c(min(locs[1]),min(locs[2]),max(locs[3]),max(locs[4])))


for (i in 1:2)
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)

expmat <- matrix(c(
    1,0,
    0,1
    ),byrow=T,ncol=2)
hsq <- c(1,1)
rland <- landscape.new.expression(rland,expmat=expmat*0.5,hsq=hsq)
rland <- landscape.new.gpmap(rland,
                             matrix(c(-1,0,1,0, #short scale
                                      -1,-0.5,1,0, #long scale
                                      -1,0,1,0, #long shape
                                       -1,-0.5,1,0,     #mixture
                                      -1,0,1,0),ncol=4,byrow=T),
                             matrix(c(-1,0,1,0,
                                      -1,0,1,0,
                                      -1,0.5,-0.4,0    #reproduction
                                      ),ncol=4,byrow=T))

initpopsize <- c(c(0,5000,5000),rep(c(0,5000,5000),(rland$intparam$habitats-1)))
rland <- landscape.new.individuals(rland,initpopsize)


l=rland


