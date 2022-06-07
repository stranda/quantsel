library(kernelPop2)
library(dplyr)
source("helpers.R")


rland <- NULL
  rland <- landscape.new.empty()
  rland <- landscape.new.intparam(rland, h=50, s=2,np=0,totgen=20000)
  rland <- landscape.new.switchparam(rland,mp=0)
  rland <- landscape.new.floatparam(rland,s=0,seedscale=c(300,1000),
                                    seedshape=c(1,500),seedmix=c(0.2),
                                    pollenscale=c(50,200),pollenshape=c(1,1),
                                    pollenmix=0.1 , asp=0.5)


  S <- matrix(c(
                  0.25, 0,
                0.15, 0.02
                ), byrow=T, nrow = 2)
  R <- matrix(c(
                0, 200,
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
  
rights <- floor(seq(500,100000,length=50))
locs <- as.matrix(data.frame(lft=c(0,501,rights[c(-1:-2)]-diff(rights[-1])+1),
                   bot=rep(0,50),
                   rgt=rights,
                   top=10000))
  
  rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                     carry=(0.50 * (sqrt((locs[,3]-locs[,1])*(locs[,4]-locs[,2])))),
                     extinct=rep(0.05,rland$intparam$habitat),
                     leftx=locs[,1],
                     rightx=locs[,3],
                     boty=locs[,2],
                     topy=locs[,4],
                     maxland=c(min(locs[1]),min(locs[2]),max(locs[3]),max(locs[4])))
  

for (i in 1:16)
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)


expmat <- matrix(c(1,0,0,
                   1,0,0,
                   1,0,0,
                   1,0,0,
                   0,1,0,
                   0,1,0,
                   0,1,0,
                   0,1,0,
                   0,0,1,
                   0,0,1,
                   0,0,1,
                   0,0,1
                   ),byrow=T,ncol=3)
  hsq <- c(1,1,1)
  rland <- landscape.new.expression(rland,expmat=expmat*0.125,hsq=hsq)
  rland <- landscape.new.gpmap(rland,
                             matrix(c(-1,0,1,
                                       0,0,1,
                                      -1,0,1,
                                       1,0,1,     #mixture
                                      -1,0,1),ncol=3,byrow=T),
                             matrix(c(-1,0,1,
                                      -1,0,1,
                                       -1,0.5,-1    #reproduction
                                      ),ncol=3,byrow=T))
initpopsize <- 5000
rland <- landscape.new.individuals(rland,c(initpopsize,initpopsize,rep(0,98)))

rland$individuals[,5] <- 4500+floor(rland$individuals[,5]/10)
###############################


l=rland
landscape.plot.locations(l)

locs <- landscape.generate.locations(npop=250,
                                     xrange=c(0,100000),yrange=c(0,10000),
                                     sizexkernel=c(1000,65),sizeykernel=c(1000,65)
                                     )

gen=100
sumlst=list()[1:gen]

for (i in 1:gen)
{
    print(dim(l$individuals))
 #   l=landscape.kill.locs(l,locs)
    print(dim(l$individuals))
    if (dim(l$individuals)[1]>0) l.old=l
    l=landscape.simulate(l,1)
    landscape.plot.locations(l)
    print(i)
    print(dim(l$individuals))
#    print(landscape.allelefreq(l) )
    print(colMeans(landscape.phenotypes.c(l)))
    sumlst[[i]] <- landscape.phenosummary(l)
    sumlst[[i]]$gen=i
}

sumdf <- do.call(rbind,sumlst)

p = sumdf %>%
    ggplot(aes(y=gen,x=pop,z=phen1_mean))+
    geom_tile(aes(fill = phen1_mean)) +
    geom_contour() +
    scale_fill_gradientn(colors = rev(cm.colors(100)))
