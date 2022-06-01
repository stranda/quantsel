"new.intparam" <-
function(h=1,s=2,l=1,ne=1,cg=0,ce=0,totgen=1,nd=1,maxland=200000,np=0,rd=0)
{
  rl <- list(h,s,l,ne,cg,ce,totgen,nd,maxland,np,rd)
  names(rl) <- c("habitats","stages","locusnum","numepochs","currentgen","currentepoch","totalgens","numdemos","maxlandsize","nphen","randdemo")
  rl
}

