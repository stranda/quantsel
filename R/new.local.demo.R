"landscape.new.local.demo" <-
function(rland,S,R,M)
{
  if (is.null(rland$demography$localdem))
    {
      rland$demography$localdem <- list(NULL)
      demonum <- 1
    }
  else
    {
      demonum <- length(rland$demography$localdem) + 1
      rland$demography$localdem[[demonum]] <- list(LocalS=NULL,LocalR=NULL,LocalM=NULL)
    }

  if (is.nsquare(S,rland$intparam$stages) && 
      is.nsquare(R,rland$intparam$stages) && 
      is.nsquare(M,rland$intparam$stages))
    {
      rland$demography$localdem[[demonum]]$LocalS <- S    
      rland$demography$localdem[[demonum]]$LocalR <- R    
      rland$demography$localdem[[demonum]]$LocalM <- M
      rland$intparam$numdemos <- length(rland$demography$localdem) 
    }
  else
    {
      stop("Matricies do not conform to stages set in intparam!")
    }
  rland
}

