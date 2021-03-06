#' Calculate demographic parameters
#' 
#' Calculate demographic parameters from a landscape: CURRENTLY BROKEN!
#' 
#' 
#' @param rland the Rmetasim landscape object
#' @return A list of length populations+1.  The first 1..populations elements
#' are lists comprised of lambda, the equilibrium stage-structure, the actual
#' stage structure, a chi^2 value for the test of difference between predicted
#' and actual, and an estimate of significance for that test.  The last element
#' of the main list is the same as the previous ones except it refers to the
#' entire landscape
#' @keywords misc
#' @export landscape.demography
"landscape.demography" <-
function(Rland)
{
  #this routine is not ready for use.
  rl <- vector("list",length=Rland$intparam$habitats+1)
  A.full <- Rland$demography$epoch[[Rland$intparam$currentepoch+1]]$R + Rland$demography$epoch[[Rland$intparam$currentepoch+1]]$S
  oblanddist <- table(Rland$individuals[,1])
  rlcnt <- 1
  for (i in 0:(Rland$intparam$habitats-1))
    {

      loclist <- NULL
      strt <- (i*Rland$intparam$stages)+1
      len <- Rland$intparam$stages-1
      obdist <- c(oblanddist[(strt):(strt+len)])
      A <- NULL
      A <- A.full[(strt):(strt+len),(strt):(strt+len)]
      es <- eigen(A)
      l <- es$values
      statedist <- es$vectors[,(which(Re(l)==max(Re(l))))]
      statedist <- (statedist/sum(statedist))
      exdist <- statedist*sum(obdist)
      chisq <- sum(((obdist-exdist)^2)/exdist)
      df <- length(obdist)-1
      prob <- 1-pchisq(chisq,df)
      loclist <- list(lambda=l[which(Re(l)==max(Re(l)))],statedist=statedist,obdist=obdist,exdist=exdist,chisq=chisq,df=df,prob=prob)
      rl[[rlcnt]] <- loclist
      rlcnt <- rlcnt+1
    }
  A <- NULL
  es <- NULL
  A <- A.full
  es <- eigen(A)
  l <- es$values
  print(A)
  statedist.1 <- es$vectors[,(which(Re(l)==max(Re(l))))]
  statedist.2 <- (statedist.1/sum(statedist.1))
  exdist <- statedist.2*sum(oblanddist)
print("exdist")
  
  print(exdist)
print("oblanddist")
  print(oblanddist)
  chisq <- sum(((oblanddist-exdist)^2)/exdist)
  df <- length(oblanddist)-1
  prob <- 1-pchisq(chisq,df)

  loclist <- list(lambda=l[which(Re(l)==max(Re(l)))],statedist=statedist.2,obdist=oblanddist,exdist=exdist,chisq=chisq,df=df,prob=prob)
  rl[[length(rl)]] <- loclist 
  rl
}

