

library(quantsel)
library(ggplot2)
library(dplyr)
library(tidyr)

### this script makes a 2 population grid,
### populates the entire thing,
### then four phenotypes,
###           each determined by a single locus
###           evolve through drift
### however, phenotype 1 has higher values in population 1 than 2 independent of genotype
###
gapprop = 0

rland <- NULL
rland <- landscape.new.empty()
rland <- landscape.new.intparam(rland, h=2, s=2,np=0,totgen=20000,maxland=3e5)
rland <- landscape.new.switchparam(rland,mp=0)
rland <- landscape.new.floatparam(rland,s=0,seedscale=c(500,2000),
                                  seedshape=c(1,300),seedmix=c(0.05),
                                  pollenscale=c(40,100),pollenshape=c(1,10),
                                  pollenmix=0.2 , asp=0.5)

S <- matrix(c(
    0,     0,
    0.8, 0.0
), byrow=T, nrow = 2)
  R <- matrix(c(
      0,  18,
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
k=rep(2500,rland$intparam$habitat)
e=rep(gapprop,rland$intparam$habitat)
rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                             carry=k,
                             extinct=e,
                             leftx=locs[,1],
                             rightx=locs[,3],
                             boty=locs[,2],
                             topy=locs[,4],
                             maxland=c(min(locs[1]),min(locs[2]),max(locs[3]),max(locs[4])))


for (i in 1:4)
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)


expmat <- matrix(c(  #4rows for 4 loci, 4 cols for 4 phenotypes
    1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
    ),byrow=T,ncol=4)
hsq <- c(1,1,1,1)
rland <- landscape.new.expression(rland,
                                  expmat=expmat*0.5,
                                  hsq=hsq) #0.5 -> 1 diploid locus per phen, up to 2 alelle additive doses
rland <- landscape.new.gpmap(rland,
                             ## 4 cols 5 rows.  Cols correspond to phenotype effects on fit components
                             ##for each phenotype (0 is no effect, 4 phenotypes in this example)
                             ##phenotypes are in C indexing so, add 1 to compare to pehnotypes above
                             matrix(c(0,   0,   0, 0,   #short scale #no selection
                                      0,   0,   0, 0,   #long scale  #no selection
                                      0,   0,   0, 0,   #long shape  #no selection
                                      0,   0,   0, 0,   #mixture   #phenotype 2   #no selection
                                      0,   0,   0, 0,   #not used  #no selection
                                      0,   0,   0, 0,   #survive  #no selection
                                      0,   0,   0, 0,   #reproduce  #no selection
                                      0,   0,   0, 0    #density tolerance  #no selection
                                      ),
                                    ncol=4,byrow=T)
                             )
                             
rland <- landscape.new.plasticity(rland,
                                  matrix(c(
                                      1.5, 1, 1, 1,
                                      0.5, 1, 1, 1
                                      ),nrow=2,ncol=4,byrow=T)) #two habitats, four phenotypes

rland <- landscape.new.phenohab(rland)

initpopsize <- 10
inits <- matrix(initpopsize,ncol=rland$intparam$habitats,nrow=2)

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

phens=c(1,2,3,4) #represented as 0 in c++
gen=100
sumlst=list()[1:gen]

l <- landscape.simulate(rland,1)

for (i in 1:gen)
{
    print(dim(l$individuals))
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
    print(i)
    print(dim(l$individuals))
#    print(landscape.allelefreq(l) )
    print(colMeans(landscape.phenotypes.c(l)))
    sumlst[[i]] <- list(phenosum=cbind(gen=i,data.frame(landscape.phenosummary(l))),
                        afreq=cbind(gen=i,landscape.allelefreq(l)))
                        
    sumlst[[i]]$gen=i
}

par(mfrow=c(1,1))
phenosumdf <- do.call(rbind,lapply(sumlst,function(x){x$phenosum}))
afreqdf <- do.call(rbind,lapply(sumlst,function(x){x$afreq}))
                      

phenosumdf <- pivot_longer(phenosumdf,cols=ends_with(c("mean","sd"))) %>% mutate(variable=name) %>% select(-name)

save(file=paste0("testPlasticity_res.rda"),sumlst,phenosumdf,afreqdf)

    phenosumdf %>% filter(!grepl("sd",variable)) %>%
    ggplot(aes(x=gen,y=value,color=variable))+
        geom_point()+geom_smooth()+
        facet_wrap(~pop)
