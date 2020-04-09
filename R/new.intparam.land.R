"landscape.new.intparam" <-
function(rland,h=1,s=1,cg=0,ce=0,totgen=1000,maxland=200000,np=0)
{
  rland$intparam <- list(h,s,0,0,cg,ce,totgen,0,maxland,np)
  names(rland$intparam) <- c("habitats","stages","locusnum","numepochs","currentgen","currentepoch","totalgens","numdemos","maxlandsize","nphen")
  rland
}

