"landscape.rm.lambda" <-
function(Rland)
{
  rl <- NULL
  for (i in 1:length(Rland$demography$localdem))
    {
       A <- Rland$demography$localdem[[i]]$LocalR + Rland$demography$localdem[[i]]$LocalS
       l <- eigen(A)$values
       rl <- c(rl,sort(Re(l[Im(l)==0]),decreasing=T)[1])
     }
  rl
}

