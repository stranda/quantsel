#"landscape.Gst" <-
#function(Rland,verb=F)
#  {
#    pl <- landscape.ploidy(Rland)
#    strt <-  9 + cumsum(pl) - pl + 1  
#    stp <- strt + (pl - 1)
#    rl <- matrix(0, nrow = Rland$intparam$habitats, ncol = length(Rland$loci))
#    afrq <- allelefreq.landscape(Rland)
#    afrq$Frqsq <- afrq$Freq^2
#
#    homframe <- aggregate(afrq$Frqsq,by=list(pop=afrq$pop,loc=afrq$loc),sum)
#    homframe <- aggregate(homframe$x,by=list(loc=homframe$loc),mean,na.rm=T)
#    homframe$loc <- as.numeric(as.character(homframe$loc))
#    Jss <- homframe[order(homframe$loc),2]
#    Jts <- sapply(1:length(Rland$loci),function(x,ind,lv,pl)
#                  {
#                    sum((table(as.numeric(ind[,9+which(lv==x)]))/(pl[x]*dim(ind)[1]))^2)
#                  },
#                  ind = Rland$individuals,
#                  lv = landscape.locusvec(Rland),
#                  pl=landscape.ploidy(Rland)
#                  )
#    Jts2 <- aggregate(afrq$
#    Hss = 1 - Jss
#    Hts = 1 - Jts
#    Dst <- Hts - Hss
#    Gst <- Dst/Hts
#    Gst
#  }

"landscape.Fst" <-
function(rland,verb=F)
  {
    aft <- landscape.allelefreq(rland,T)
    Fst <- matrix(0,length(landscape.ploidy(rland)),length(aft[1,,1]))
    for (i in 1:length(Fst[,1]))
      {
        Fst[i,] <- ((colSums(sweep(aft[,,i],2,colMeans(aft[,,i]))^2)/(rland$intparam$habitats))/(colMeans(aft[,,i])*(1-colMeans(aft[,,i]))))
      }
    if (verb)
      {
        print(paste("Populations:",rland$intparam$habitats,"Loci:",rland$intparam$locusnum))
        print(paste("Mean per locus:"))
        print(rowMeans(Fst,na.rm=T))
        print(paste("Overall mean:",mean(rowMeans(Fst,na.rm=T))))
        print("===")
      }
    Fst
  }

landscape.Fst.pairwise <- function(rland)
  {
    retmat <- matrix(NA,nrow=rland$intparam$habitat,ncol=rland$intparam$habitat)
    for (i in 1:rland$intparam$habitat)
      for (j in 1:i-1)
        {
          if (length(landscape.populations(rland))>0)
            {
              if ((max(landscape.populations(rland)==i)>0)&(max(landscape.populations(rland)==j)>0))
                {
                  compland <- rland
                  compland$individuals <- compland$individuals[landscape.populations(compland) %in% c(i,j),]
                  retmat[i,j] <- mean(rowMeans(landscape.Fst(compland),na.rm=T))
                }
            }
        }
    retmat
  }

landscape.pairwise.dist <- function(rland)
  {
    xcoords <- (rland$demography$epochs[[1]]$leftx + rland$demography$epochs[[1]]$rightx)/2
    ycoords <- (rland$demography$epochs[[1]]$topy + rland$demography$epochs[[1]]$boty)/2
    as.matrix(dist(cbind(xcoords,ycoords),diag=T,upper=T))
  }

mantel <- function(dist1,dist2,perm=1000)
  {
    Z <- sum(dist1*dist2)
    rearrange <- rep(0,perm)
    for (i in 1:perm)
      {
        rearrange[i] <- sum(sample(dist1,length(dist1),replace=F)*dist2)
      }
    sum(rearrange>Z)/perm
  }
