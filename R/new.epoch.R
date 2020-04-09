
"landscape.new.epoch" <- function(rland,S=NULL,R=NULL,M=NULL,epochprob=1,startgen=0,extinct=NULL,carry=NULL,localprob=NULL,
                                  pollen.kernels=NULL,
                                  seed.kernels=NULL,
                                  leftx=NULL, rightx=NULL, boty=NULL, topy=NULL,
                                  maxland=c(0,0,10000,10000))
{
  popw=0.1 * sum(maxland[c(1,3)])/2
  poph=0.1 * sum(maxland[c(2,4)])/2
  if (is.null(rland$demography$epochs))
    {
      rland$demography$epochs <- list(NULL)
      epochnum <- 1
    }
  else
    {
      epochnum <- length(rland$demography$epochs) + 1
      rland$demography$epochs[[epochnum]] <- list(RndChooseProb=NULL,StartGen=NULL,Extinct=NULL,
                                                  Carry=NULL,Localprob=NULL,S=NULL,R=NULL,M=NULL)
    }
  
  rland$demography$epochs[[epochnum]]$RndChooseProb <- epochprob
  rland$demography$epochs[[epochnum]]$StartGen <- startgen
  if (is.null(extinct))
    {
      extinct <- rep(0, rland$intparam$habitats)
    }  
  if (length(extinct) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$Extinct <- extinct
    }
  else
    {
      stop("Wrong size for extinct vector", dim((extinct)))
    }
  
  if (is.null(carry)) 
    {
      carry <- rep(1000, rland$intparam$habitats)
    }
  if (length(carry) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$Carry <- carry
    }
  else
    {
      stop("Wrong size for carry vector")
    }

  numdem <- length(rland$demography$localdem)
  if (is.null(localprob)) 
    {
      localprob <- rep(1/numdem, numdem)
    }
  if (length(localprob == numdem))
    {
      rland$demography$epochs[[epochnum]]$Localprob <- localprob
    }
  else
    {
      stop("Wrong size for localprob vector")
    }
  
  if (is.null(leftx)|is.null(rightx)|is.null(topy)|is.null(boty))
    {
      leftx <- rep(seq(maxland[1],maxland[3]-popw,length=ceiling(rland$intparam$habitats/2)),2)
      if (length(leftx)>rland$intparam$habitats)
        {
          leftx <- leftx[-length(leftx)]
        }
      rightx <- leftx+popw
      boty <- rep(seq(maxland[2],maxland[4]-poph,length=ceiling(rland$intparam$habitats/2)),2)
      if (length(leftx)>rland$intparam$habitats)
        {
          boty <- boty[-length(boty)]
        }
      topy <- boty+poph
    }  
  if (length(leftx) == rland$intparam$habitats)
    {
      rland$demography$epochs[[epochnum]]$leftx  <- leftx
      rland$demography$epochs[[epochnum]]$rightx <- rightx
      rland$demography$epochs[[epochnum]]$boty   <- boty
      rland$demography$epochs[[epochnum]]$topy   <- topy
    }
  else
    {
      stop("Wrong size for coords vector", dim((extinct)))
    }

  matrixsize <- rland$intparam$habitats * rland$intparam$stages
  
  if (is.null(S))
    {
      S <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.null(R))
    {
      R <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.null(M))
    {
      M <- matrix(rep(0, matrixsize * matrixsize),ncol=matrixsize,nrow=matrixsize)
    }

  if (is.nsquare(S,matrixsize) && 
      is.nsquare(R,matrixsize) && 
      is.nsquare(M,matrixsize))
    {
      rland$demography$epochs[[epochnum]]$S <- S
      rland$demography$epochs[[epochnum]]$R <- R
      rland$demography$epochs[[epochnum]]$M <- M
      rland$intparam$numepochs <- length(rland$demography$epochs)
    }
  else
    {
      stop("S, R, and M matricies are not the correct size")
    }

  if(is.null(seed.kernels))
    {# a kernel matrix has the same number of rows as a side of the landscape matrix and
      #six columns.  The columns correspond to 1 kernel, 2-3 scale and shape parameters of the kernel
      #4-5 scale and shape parameters for the 2nd distribution in a mixture and 6 mixing parameter (proportion of the first dist)
      seed.kernels <- cbind(rep(3,matrixsize),rep(rland$floatparam$seedmu,matrixsize),
                            rep(rland$floatparam$seedshape,matrixsize),rep(rland$floatparam$seedmu2,matrixsize),
                            rep(rland$floatparam$seedshape2,matrixsize),rep(rland$floatparam$seedmix,matrixsize))
    }
  if ((dim(seed.kernels)[2]!=6)||(dim(seed.kernels)[1]!=rland$intparam$habitats*rland$intparam$stages))
    {
      stop("Kernel matrix of incorrect dimension or landscape.intparam not run yet")
    }

  if(is.null(pollen.kernels))
    {# a kernel matrix has the same number of rows as a side of the landscape matrix and
      #six columns.  The columns correspond to 1 kernel, 2-3 scale and shape parameters of the kernel
      #4-5 scale and shape parameters for the 2nd distribution in a mixture and 6 mixing parameter (proportion of the first dist)
      pollen.kernels <- cbind(rep(3,matrixsize),rep(rland$floatparam$pollenmu,matrixsize),rep(rland$floatparam$pollenshape,matrixsize),rep(rland$floatparam$pollenmu2,matrixsize),rep(rland$floatparam$pollenshape2,matrixsize),rep(rland$floatparam$pollenmix,matrixsize))
    }
  
  if ((dim(pollen.kernels)[2]!=6)||(dim(pollen.kernels)[1]!=rland$intparam$habitats*rland$intparam$stages))
    {
      stop("Kernel matrix of incorrect dimension or landscape.intparam not run yet")
    }
  
  rland$demography$epochs[[epochnum]]$seedkern <- seed.kernels
  rland$demography$epochs[[epochnum]]$pollenkern <- pollen.kernels
  rland
}

