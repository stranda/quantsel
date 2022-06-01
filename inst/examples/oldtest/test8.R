library(kernelPop2)
library(ggplot2)
library(dplyr)
source("helpers.R")
source("analysis.R")

### this script maks a 250 population grid and 
###


gapprop = 0



rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=1024, s=2,np=0,totgen=20000,maxland=3e5)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,seedscale=c(40,290),
                                  seedshape=c(1,300),seedmix=c(0.12),
                                  pollenscale=c(40,100),pollenshape=c(1,10),
                                  pollenmix=0.2 , asp=0.5)


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

rights <- floor(seq(0,40000,length=33))
tops <-   floor(seq(0,40000,length=33))
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
k=(0.2 * (sqrt((locs[,3]-locs[,1])*(locs[,4]-locs[,2]))))
e=rep(gapprop,rland$intparam$habitat)
k[525:526] <- k[525:526]*2
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
    1,0,0,0,
    1,0,0,0,
    1,0,0,0,
    1,0,0,0,
    0,1,0,0,
    0,1,0,0,
    0,1,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,1,0,
    0,0,1,0,
    0,0,1,0,
    0,0,0,1,
    0,0,0,1,
    0,0,0,1,
    0,0,0,1
                   ),byrow=T,ncol=4)
hsq <- c(1,1,1,1)
rland <- landscape.new.expression(rland,expmat=expmat*0.125,hsq=hsq)
rland <- landscape.new.gpmap(rland,
                             matrix(c(-1,0,1,0, #short scale
                                       2,-0.5,1,0, #long scale
                                      -1,0,1,0, #long shape
                                       0,-0.5,1,0,     #mixture
                                      -1,0,1,0),ncol=4,byrow=T),
                             matrix(c(-1,0,1,0,
                                      -1,0,1,0,
                                      1,0.5,-0.4,0    #reproduction
                                      ),ncol=4,byrow=T))

initpopsize <- 10000
inits <- matrix(0,ncol=rland$intparam$habitats,nrow=2)
inits[1:2,c(528:529,560:561)] <- initpopsize  
rland <- landscape.new.individuals(rland,c(inits))

#rland$individuals[,5] <- 4500+floor(rland$individuals[,5]/10)
###############################

l=rland
landscape.plot.phenotypes(l,1)

locs <- landscape.generate.locations(npop=1024,
                                     xrange=c(0,40000),yrange=c(0,40000),
                                     sizexkernel=c(400,65),sizeykernel=c(400,65)
                                     )

phens=c(1,2,3,4) #represented as 0:3 in c++
gen=250
sumlst=list()[1:ceiling(1+gen/5)]

#pdf(paste0("gaps_",gapprop,".pdf"), width=15,height=7.5)
slc=1
for (i in 1:gen)
{
    print(dim(l$individuals))
 #   l=landscape.kill.locs(l,locs)
    print(dim(l$individuals))
    if (dim(l$individuals)[1]>0) l.old=l
    l=landscape.simulate(l,1)
    if ((i %% 1)==0)
    {
        par(mfrow=c(2,2))
        for (phen in phens)
            landscape.plot.phenotypes(l,phen)
        par(mfrow=c(1,1))
    }

    if ((i %% 5)==0)
    {
        print("summarizing")
        sumlst[[slc]]$pheno <- data.frame(landscape.phenosummary(l))
        print("summarizing genos")
        sumlst[[slc]]$geno <- landscape.gensummary(l)
        sumlst[[slc]]$neighbor <- landscape.neighborhood(l)
        sumlst[[slc]]$gen=i
        slc <- slc+1
    }
    
    print(i)
    print(dim(l$individuals))
#    print(landscape.allelefreq(l) )
    print(colMeans(landscape.phenotypes.c(l)))

}
#dev.off()


save(file=paste0("gap_",gapprop,"_res.rda"),sumlst)

