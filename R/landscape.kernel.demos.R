


###tries out the different kernels.  Plots offspring and mothers
###also compares the dispersal distances to what we intended
landscape.pollen.kernel.demo <- function(s=0,pmn1=30,pshp1=1.01,
                                         pmn2=600,pshp2=200,
                                         pmix=0.75,plot=T)
{

  ##the landscape produced by this function is inteded to be an example of how different stages within populations can have different
  ##seed dispersal parameters (can also differ among populations too)
  ## this is set up as a demostration, using landscape.simulate will not work, but landscape.reproduce will produce the seeds and
  ##seed dispersal shadow
  landscape.stage.kernel.test<-function (s, sdkernel, polkernel, polscale=c(90,600), polshape=c(1,200), polmix=1,
                                         offspring=200)
    {
      rland <- NULL
      rland <- landscape.new.empty()
      rland <- landscape.new.intparam(rland, h = 1, s = 4, np = 0, 
                                      totgen = 20000)
      rland <- landscape.new.switchparam(rland, mp = 1)
      rland <- landscape.new.floatparam(rland, s = s, pollenscale = polscale,
                                        pollenshape = polshape, pollenmix=polmix,
                                        asp = 1)
      
      ##all non-seed stages survive forever.
      ##seeds (stage 0 die every year)
      S <- matrix(c(  0, 0, 0,   0,
                    0, 0, 0,   0,
                    0, 0, 0,   0,
                    0.01, 0, 0,   0
                    ), byrow = T, nrow = rland$intparam$stages)
      
      ##lower number stages produce less offspring
      ##
      R <- matrix(c(0, offspring, offspring, offspring,
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0
                    ), byrow = T, nrow = rland$intparam$stages)
      
      ##
      ## every parent has the same male gamete output.
      M <- matrix(c(0, 0, 0, 0,
                    0, 1, 1, 1,
                    0, 1, 1, 1,
                    0, 1, 1, 1
                    ), byrow = T, nrow = rland$intparam$stages)
      rland <- landscape.new.local.demo(rland, S, R, M)
      S <- matrix(rep(0, (rland$intparam$stages * rland$intparam$habitat)^2), nrow = rland$intparam$stages * 
                  rland$intparam$habitat)
      R <- matrix(rep(0, (rland$intparam$stages * rland$intparam$habitat)^2), nrow = rland$intparam$stages * 
                  rland$intparam$habitat)
      M <- matrix(rep(0, (rland$intparam$stages * rland$intparam$habitat)^2), nrow = rland$intparam$stages * 
                  rland$intparam$habitat)
      
      
      if ((dim(sdkernel)[1]==dim(S)[1])&&(dim(sdkernel)[2]==6))
        {
          rland <- landscape.new.epoch(rland, S = S, R = R, M = M, 
                                       carry = 20000,
                                       extinct = rep(0.05, rland$intparam$habitat),
                                       seed.kernels = sdkernel,
                                       pollen.kernels=polkernel,
                                       leftx = 0, rightx = 4000, boty = 0, topy = 4000,
                                       maxland = c(0,0,4000,4000))
        }
      else
        {
          stop("dimensions of the kernel object does not correspond to the landscape")
        }
      
      rland <- landscape.new.locus(rland,type=0, ploidy=1, mutationrate=0.005, numalleles=3, frequencies = c(0.2, 0.2, 0.6))
      rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.001,transmission=0,numalleles=5)
      ##print(unlist(rland$floatparam))
      
      rland <- landscape.new.individuals(rland, c(1500,0,0,0))
      rland
    }
#########end included function to create landscapes

#######begin main function
  ##this is the object that describe pollen kernels for different stages in a landscape. This expects 4 stages
  ##with the first stage non-reproducing.  I'm using an exponential, weibull and mixed distribution for the four stages
  ##column 1 of this matrix.  The other columns are scale and shape parameters (two sets because some distributions are mixtures)
  ##and a mixting parameter in the last col.
  
  sk <- matrix(c( 0,   0, 0,   0,   0,    0,    #stage 0 does not produce
                                                #seed so parameters not needed
                  1, 100, 0,   0,   0,    0,    #exponential with mean 100
                  2, 200, 2,   0,   0,    0,    #weibull with scale 200 and shape 2   
                  3,  50, 2, 300, 300, 0.25      #mixture with weibull scale 50, shape 2,
                                                #normal mean 300, sd 200, mixture 75% weibull
                 ), ncol=6,byrow=T)

  ###all mixed pollen kernel with exponential short, normal long
  pk <- matrix(c( 3,    pmn1, pshp1,  pmn2,   pshp2,    pmix,  
                  3,    pmn1, pshp1,   pmn2,   pshp2,    pmix,  
                  3,    pmn1, pshp1,   pmn2,   pshp2,    pmix,  
                  3,    pmn1, pshp1,   pmn2,   pshp2,    pmix   
                 ), ncol=6,byrow=T)


###actually create the landscape and simulate a few gens
  
  l <- landscape.stage.kernel.test(s=s,sdkernel=sk,polkernel=pk,offspring=109)
  print(l$floatparams$selfing)

  l <- landscape.simulate(l,3)

  # the rest of this function is just fancy plotting
if (is.landscape(l)&&plot==T)
  {
    mtr <- apply(l$ind[l$ind[,1]==0,],1,function(x,l){intersect(which(x[6]==l$ind[,4]),which(x[7]==l$ind[,5]))},l)
    xmi <- min(c(l$ind[l$ind[,1]==0,4]))
    xma <- max(c(l$ind[l$ind[,1]==0,4]))
    ymi <- min(c(l$ind[l$ind[,1]==0,5]))
    yma <- max(c(l$ind[l$ind[,1]==0,5]))
    pd=pollination.dist(l)
    par(mfrow=c(2,1))
    pd=pollination.dist(l)
    psq=seq(0,max(pd),length=50)
    lv=l$demography$epochs[[1]]$pollenkern[dim(l$demography$epochs[[1]]$pollenkern)[1],]
    hist(pd,freq=F,breaks=40,xlab="Pollination distance",main="Pollen cloud and predicted density",sub=paste("Realized selfing rate: ",round(realized.selfing(l),2)))
    points(psq,dmixed(psq,lv[2],lv[3],lv[4],lv[5],lv[6]),type="l")
    plot(l$ind[l$ind[,1]==0,5]~l$ind[l$ind[,1]==0,4],ylab="Y-coordinate",xlab="X-coordinate",
         pch=16,col=unclass(factor(unlist(mtr))),ylim=c(ymi,yma),xlim=c(300,3700))
    par(mfrow=c(1,1))
    l
  } else {NULL}
}

