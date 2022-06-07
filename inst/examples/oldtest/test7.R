library(kernelPop2)
library(ggplot2)
library(dplyr)
source("helpers.R")

### this script maks a 250 population grid and 
###


gapprop = 0



rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=256, s=2,np=0,totgen=20000)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,seedscale=c(300,3000),
                                  seedshape=c(1,1000),seedmix=c(0.2),
                                  pollenscale=c(50,200),pollenshape=c(1,1),
                                  pollenmix=0.1 , asp=0.5)


S <- matrix(c(
    0, 0,
    0.8, 0.0
), byrow=T, nrow = 2)
  R <- matrix(c(
      0, 12,
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

rights <- floor(seq(00,40000,length=17))
tops <-   floor(seq(0,40000,length=17))
locs=NULL
for (i in 1:(length(tops)-1))
{
    locs <- rbind(locs,
                  data.frame(lft=c(rights[-1]-diff(rights[])+1),
                             bot=rep(tops[i+1]-(tops[2]-tops[1]),16),
                             rgt=rights[-1],
                             top=tops[i+1]))
}

lfts=which(locs$lft==1)
k=(0.90 * (sqrt((locs[,3]-locs[,1])*(locs[,4]-locs[,2]))))
e=rep(gapprop,rland$intparam$habitat)
k[128] <- k[128]*2
rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                             carry=k,
                             extinct=e,
                             leftx=locs[,1],
                             rightx=locs[,3],
                             boty=locs[,2],
                             topy=locs[,4],
                             maxland=c(min(locs[1]),min(locs[2]),max(locs[3]),max(locs[4])))


for (i in 1:12)
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)


expmat <- matrix(c(1,0,0,
                   1,0,0,
                   0,0,0,
                   0,0,0,
                   0,1,0,
                   0,1,0,
                   0,0,0,
                   0,0,0,
                   0,0,1,
                   0,0,1,
                   0,0,0,
                   0,0,0
                   ),byrow=T,ncol=3)
hsq <- c(1,1,1)
rland <- landscape.new.expression(rland,expmat=expmat*0.25,hsq=hsq)
rland <- landscape.new.gpmap(rland,
                             matrix(c(-1,0,1, #short scale
                                      2,-0.5,1, #long scale
                                      -1,0,1, #long shape
                                      1,-0.5,1,     #mixture
                                      -1,0,1),ncol=3,byrow=T),
                             matrix(c(-1,0,1,
                                      -1,0,1,
                                      0,-0.5,1    #reproduction
                                      ),ncol=3,byrow=T))
initpopsize <- 5000
inits <- matrix(0,ncol=rland$intparam$habitats,nrow=2)
inits[1:2,120] <- initpopsize  
rland <- landscape.new.individuals(rland,c(inits))

#rland$individuals[,5] <- 4500+floor(rland$individuals[,5]/10)
###############################


l=rland
landscape.plot.phenotypes(l,1)

locs <- landscape.generate.locations(npop=250,
                                     xrange=c(0,100000),yrange=c(0,10000),
                                     sizexkernel=c(1000,65),sizeykernel=c(1000,65)
                                     )


phens=c(1,2,3) #represented as 0 in c++
gen=200
sumlst=list()[1:gen]

#pdf(paste0("gaps_",gapprop,".pdf"), width=15,height=7.5)

for (i in 1:gen)
{
    print(dim(l$individuals))
 #   l=landscape.kill.locs(l,locs)
    print(dim(l$individuals))
    if (dim(l$individuals)[1]>0) l.old=l
    l=landscape.simulate(l,1)
    if ((i %% 1)==0)
    {
        par(mfrow=c(1,length(phens)))
        for (phen in phens)
            landscape.plot.phenotypes(l,phen)
        par(mfrow=c(1,1))
    }
    print(i)
    print(dim(l$individuals))
#    print(landscape.allelefreq(l) )
    print(colMeans(landscape.phenotypes.c(l)))
    sumlst[[i]] <- data.frame(landscape.phenosummary(l))
    sumlst[[i]]$gen=i
}
#dev.off()

sumdf <- do.call(rbind,sumlst)

save(file=paste0("gap_",gapprop,"_res.rda"),sumlst,sumdf)


p = sumdf %>%
    ggplot(aes(y=gen,x=pop,z=phen1_mean))+
    geom_tile(aes(fill = phen1_mean)) +
    geom_contour() +
    scale_fill_gradientn(colors = rev(cm.colors(100)))
