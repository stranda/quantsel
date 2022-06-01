library(quantsel)
library(ggplot2)
library(dplyr)

### this script maks a 1024 population grid and 
### populates the entire thing, then three traits evolve


gapprop = 0 #proportion of gaps put in landscape



rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=1024, s=2,np=4,totgen=20000,maxland=3e5)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,seedscale=c(40,1500),
                                  seedshape=c(1,500),seedmix=c(0.05),
                                  pollenscale=c(40,1250),pollenshape=c(1,100),
                                  pollenmix=0.2 , asp=0.5)


S <- matrix(c(
    0   ,   0,
    0.75, 0.0
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
k=rep(200,rland$intparam$habitat)
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


expmat <- matrix(c(  #16rows for 16 loci, 4 cols for 4 phenotypes
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
rland <- landscape.new.expression(rland,
                                  expmat=expmat*0.125, #0.125 -> 1 diploid locus per phen, 
                                  hsq=hsq) #up to 8 alelle additive doses, when summed across 4 loci.  
rland <- landscape.new.gpmap(rland,
                             ## 4 cols 5 rows.  Cols correspond to phenotype effects on fit components
                             ##for each phenotype (0 is no effect, 4 phenotypes in this example)
                             ##phenotypes are in C indexing so, add 1 to compare to pehnotypes above
                             matrix(c(0,   0,   0, 0,   #short scale #no selection
                                      0,   0,   0, 0,   #long scale 
                                      0,   0,   0, 0,   #long shape  #no selection
                                      0,   0,   0, 0,   #mixture   #phenotype 2   #no selection
                                      0,   0,   0, 0,   #not used  #no selection
                                      0,   0,   0, 0,   #survive  
                                      1,   0,   0, 0,   #reproduce 
                                      0,   0,   0, 0    #density tolerance  #no selection
                                      ),
                                    ncol=4,byrow=T)
                             )
                             
rland <- landscape.new.plasticity(rland)

rland <- landscape.new.phenohab(rland)

##reproduction gradient on pheno 1
ph <- matrix(rep(c(1,2,0.05,0),  #alpha=1, beta=2, range around 1 = 0.05 and use original sign (=0)
, rland$intparam$habitats),nrow=rland$intparam$habitats,ncol=4,byrow=T)

rland <- landscape.new.phenohab(rland,7,ph=ph)
if (FALSE)
    {
##survive gradient on pheno 2
ph <- matrix(c(
    1,2,0.1,0,  #alpha=1, beta=2, range around 1 = 0.05 and use original sign (=0)
    1,2,0.1,0
),nrow=2,ncol=4,byrow=T)

rland <- landscape.new.phenohab(rland,6,ph=ph)

##dispersal distance on phenotype 3
ph <- matrix(c(
    1,2,0.3,0,  #alpha=1, beta=2, range around 1 = 0.05 and use original sign (=0)
    1,2,0.3,0
),nrow=2,ncol=4,byrow=T)

rland <- landscape.new.phenohab(rland,2,ph=ph)
    }

initpopsize <- 50
inits <- matrix(initpopsize,ncol=rland$intparam$habitats,nrow=2)

rland <- landscape.new.individuals(rland,c(inits))

l=landscape.simulate(rland,1)

landscape.plot.phenotypes(l,1,F)

locs <- landscape.generate.locations(npop=1024,
                                     xrange=c(0,40000),yrange=c(0,40000),
                                     sizexkernel=c(400,65),sizeykernel=c(400,65)
                                     )

phens=c(1,2,3,4) #represented as 0 in c++
gen=100
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
        par(mfrow=c(2,2))
        for (phen in phens)
            landscape.plot.phenotypes(l,phen,F)
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
