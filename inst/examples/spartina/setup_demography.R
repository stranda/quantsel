spartina.landscape <- function(nloc=1,
                               nphen=0,
                               K=1000,
                               seedmix=0.1) #just the landscape setup, demography and dispersal
    {
        rland <- NULL
        rland <- landscape.new.empty()
        rland <- landscape.new.intparam(rland, h=50, s=2,np=nphen,totgen=20000,maxland=3e5)
        rland <- landscape.new.switchparam(rland,mp=0)
        rland <- landscape.new.floatparam(rland,s=0,seedscale=c(40,100),
                                          seedshape=c(1,50),seedmix=seedmix,
                                          pollenscale=c(40,100),pollenshape=c(1,200),
                                          pollenmix=0.1 , asp=1)
        
        
        S <- matrix(c(
            0   ,   0,
            0.25, 0.02
        ), byrow=T, nrow = 2)
        R <- matrix(c(
            0,  7,
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
        
        rights <- rep(200,50)
        lefts <- rep(0,50)
        bots <-   floor(seq(0,980,length=50))
        tops <-   c(bots[2:length(bots)],1000)
        locs <- data.frame(lft=lefts,
                           bot=bots+1,
                           rgt=rights,
                           top=tops)

        
        k=rep(K,rland$intparam$habitat)
        e=rep(0,rland$intparam$habitat)
        rland <- landscape.new.epoch(rland,S=S,R=R,M=M,
                                     carry=k,
                                     extinct=e,
                                     leftx=locs[,1],
                                     rightx=locs[,3],
                                     boty=locs[,2],
                                     topy=locs[,4],
                                     maxland=c(min(locs[1]),min(locs[2]),max(locs[3]),max(locs[4])))
        for (i in 1:nloc)
            rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.00,transmission=0,numalleles=2)
        rland
    }
