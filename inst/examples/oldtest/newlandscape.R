message("adding some libraries in newlandscape.R")
library(compiler)
#library(data.table)
library(inline)
library(Rcpp)
enableJIT(3)



##
## attempting fast vectorized R instead of C++
##
landscape.survive.Rversion <- function(rland)
{

    ##while there are a lot of for loops in this function, each only runs for a few iterations
    ##each click could do 100,000 operations so, there is a pretty good efficiency here.  Could be
    ##improved no doubt.
    
    S <- rland$demography$localdem[[1]]$LocalS
      
    stgs <- rland$individuals[,1] %% rland$intparam$stages  #convert democlass to local stage
    
    cm <- sm  <-  matrix(0,nrow=dim(rland$individuals)[1],ncol=dim(S)[2]) ##one col per from stage
    ##(will contain trans probs)
    for (i in 1:dim(S)[2])  #populate the sm matrix rows with columns of trans probs 
        sm[stgs==(i-1),] <- matrix(rep(S[,i],sum(stgs==(i-1))),nrow=sum(stgs==(i-1)),byrow=T)
    

    for (i in 1:dim(S)[2]) #make each row a cumulative sum
        if (i==1) cm[,i] <- sm[,i] else cm[,i] <- cm[,(i-1)]+sm[,i]

    
    rn <- runif(length(stgs))
    ns <- rep(-1,length(rn))
    for (i in 1:dim(S)[2]) #truck through cols and identify the column by slot in cum prob dist
        ns <- ifelse(((ns<0) &(rn <= cm[,i])),i-1,ns) #need to convert these local stages to global

    rland$individuals[,1] <- (rland$individuals[,1]-stgs) + ns #globalizing
    rland$individuals <- rland$individuals[ns>=0,] #get rid of dead inds.

    
    rland                                           
}


FindMale <- function(rland,ind) #returns rowname with male id in rland, based on female
                                #location and pollen kernels;
                               
{
    M <- rland$demography$localdem[[1]]$LocalM
    Mlong <- as.data.frame(cbind((which(M>0,arr.ind=T)-1),M[which(M>0,arr.ind=T)]))
    names(Mlong)[3] <- "maleCont"
    possible.males <- which(rland$individuals[,1]%in%(Mlong$col-1))

    if (ind>dim(rland$individuals)[1])
    {
        print(ind)
        print(dim(rland$individuals))
        stop("bad ind identifier, must refer to row in land individs")
        }
    allcoords <- rland$individuals[possible.males,4:5]
    allcoords[,1] <- allcoords[,1] - rland$individuals[ind,c(4)] #translate all the coords
    allcoords[,2] <- allcoords[,2] - rland$individuals[ind,c(5)] #so the 'ind' is at 0,0
    d2m <- data.frame(maleids=possible.males,
                      dist=sqrt(allcoords[,1]^2+allcoords[,2]^2)) #dist from origin=dist from ind
#    d2m <- d2m[d2m$dist!=0,]
    kern <- rland$demography$epochs[[1]]$pollenkern[(rland$individuals[ind,1]+1),]
    d <- rdispfunc(1,kern)
    d2m$maleids[which(abs(d2m$dist-d)==min(abs(d2m$dist-d)))]
}



##
## slower
fm2.dst.cpp <- '
          
          Rcpp::NumericVector rv = Rcpp::clone<Rcpp::NumericVector>(ac);
          std::transform(rv.begin(),rv.end(),rv.begin(),::sqrt);
          return rv;
          
'
#fm2.dst <- cxxfunction(signature(ac="numeric"),fm2.dst.cpp,plugin="Rcpp")



FindMale2 <- function(individuals,ind,Mlong,possible.males,pollenkern) #returns rowname with male id in rland, based on female
                                #location and pollen kernels;
{

    acrd <- individuals[possible.males,4:5]
    asq <- ( acrd-c(individuals[ind,c(4,5)]))^2

    dst <- sqrt(asq[,1]+asq[,2])
    kern <- pollenkern[(individuals[ind,1]+1),]
    d <- rdispfunc(1,kern)
    ret <- possible.males[which.min(abs(dst-d))]
    ret[1]
}


AllMales <- function(rland,reprows) #gonna need optimization along with FindMale
{
    sapply(reprows,function(x) FindMale(rland,x))
}


AllMales2 <- function(rland,reprows) #gonna need optimization along with FindMale
{
    M <- rland$demography$localdem[[1]]$LocalM
    Mlong <- as.data.frame(cbind((which(M>0,arr.ind=T)-1),M[which(M>0,arr.ind=T)]))
    names(Mlong)[3] <- "maleCont"
    possible.males <- which(rland$individuals[,1]%in%(Mlong$col))
    pollenkern=rland$demography$epochs[[1]]$pollenkern
    vapply(reprows,function(ind) FindMale2(rland$individuals,ind,Mlong,possible.males,pollenkern),integer(1)) 
}


InCircle.R <- function(circle,points)
{
    points <- as.matrix(points)
    cx=circle[1]
    cy=circle[2]
    r=circle[3]
    x=points[,1]
    y=points[,2]
    ifelse((x-cx)*(x-cx) + (y-cy)*(y-cy)<=r*r,TRUE,FALSE)
}

sourceCpp("FindMales3.cpp")

###uses data.table
###tries to figure out subsets of fathers data appropriate for each mother

AllMales3 <- function(rland,reprows,dens=0.01)
{
#    print("running allmales3")
    M <- rland$demography$localdem[[1]]$LocalM
    Mlong <- as.data.frame(cbind((which(M>0,arr.ind=T)-1),M[which(M>0,arr.ind=T)]))
    names(Mlong)[3] <- "maleCont"
    pollenkern=rland$demography$epochs[[1]]$pollenkern
    
    individuals <- (rland$individuals)

    #moms <- individuals[reprows,1:5] below I just use the righthand side of this
    
    dadclass <- apply(t(rmultinom(n=length(reprows),size=1,
                                  prob=Mlong$maleCont[Mlong$row==(unlist(individuals[reprows[1],1])%%rland$intparam$stages)])),1,which.max)
    
    v <- table(dadclass)
    rdists <- lapply(1:length(v),function(i){list(idx=1,
                                                  dists=rdispfunc(v[i],
                                                                  pollenkern[as.numeric(names(v)[i])+1,]))})

    mrdist=rep(-1,length(reprows))
    for (i in 1:length(rdists))  {  mrdist[dadclass==i]  <-  rdists[[i]]$dists  }

###
### use the c++ approach
###
#### use a series of rings to identify ok dads for each mom
#### this can be sped up in in C++ for sure

###   c(FindMales3(reprows,as.matrix(individuals),mrdist,1,0.1,0.1)[,2])
    c(FindMales4(reprows,
                 as.matrix(rland$individuals),
                 landscape.populations(rland),
                 rep(1.5*rland$demography$epochs[[1]]$pollenkern[2,4],dim(rland$individuals)[1]),
                 c(rland$intparam$stages),
                 cbind(rland$demography$epochs[[1]]$leftx,
                            rland$demography$epochs[[1]]$boty,
                            rland$demography$epochs[[1]]$rightx,
                            rland$demography$epochs[[1]]$topy),
                 rland$demography$epochs[[1]]$pollenkern
                 )[,2])

 
}

AllMales4 <- function(rland,Rlong,dens=0.01)
{

    side=4

    xrange=c(min(rland$demography$epoch[[1]]$leftx),max(rland$demography$epoch[[1]]$rightx))
    yrange=c(min(rland$demography$epoch[[1]]$boty),max(rland$demography$epoch[[1]]$topy))

    vert=diff(seq(xrange[1],xrange[2],length=side+1))[1]
    horiz=diff(seq(yrange[1],yrange[2],length=side+1))[1]

    locs <- NULL
    for (i in 1:side)
        for (j in 1:side)
        {
            locs <- rbind(locs,
                          c(xleft=xrange[1]+(i-1)*horiz,
                            boty=yrange[1]+(j-1)*vert,
                            xright=xrange[1]+(i-1)*horiz+horiz,
                            topy=yrange[1]+(j-1)*vert+vert)
                          )
        }
    
    
    M <- rland$demography$localdem[[1]]$LocalM
    Mlong <- as.data.frame(cbind((which(M>0,arr.ind=T)-1),M[which(M>0,arr.ind=T)]))
    names(Mlong)[3] <- "maleCont"
    pollenkern=rland$demography$epochs[[1]]$pollenkern
    
###
### use the c++ approach
###
#### use a series of rings to identify ok dads for each mom
#### this can be sped up in in C++ for sure
    orig=1:dim(rland$individuals)[1]
###   c(FindMales3(reprows,as.matrix(individuals),mrdist,1,0.1,0.1)[,2])
    do.call(rbind,lapply(1:dim(locs)[1],function(l)
    {
        loc=locs[l,]
        sv1 <- (rland$individuals[,4]>c(loc[1]))&(rland$individuals[,4]<c(loc[3]))&
            (rland$individuals[,5]>c(loc[2]))&(rland$individuals[,5]<c(loc[4]))
        inds <- rland$individuals[sv1,]
        reprows <- which((c(as.matrix(inds[,1]))%%rland$intparam$stages) %in% Rlong[,2])
        oids  <- orig[sv1]
#        print(loc)
#        print(dim(inds)[1])
        if (dim(inds)[1]>0)
        {
            dads=c(FindMales4(reprows,
                         as.matrix(inds),
                         rep(1.5*pollenkern[2,4],dim(inds)[1]),
                         pollenkern
                         )[,2])
            cbind(oids[reprows],oids[dads])
        }
    }))
    
}


FindMales4_R <- function(reprows){NULL}


rdispfunc <- function(n,krn)
    {
      if (krn[1]==1)
        {rweibull(n,scale=krn[2],shape=1)}
      if (krn[1]==2)
        {rweibull(n,scale=krn[2],shape=krn[3])}
      if (krn[1]==3)
        {
          if (krn[6]<1)
            {
              p <- 1*(runif(n)<krn[6])
              p*rweibull(n,scale=krn[2],shape=krn[3])+
                (1-p)*rnorm(n,mean=krn[4],sd=krn[5])
            } else {
              rweibull(n,scale=krn[2],shape=krn[3])
            }
        } else {stop("undefined kernel")}
    }


#
# fast reproduction using vectorized functions
# (trying to avoid kernelPop)
#
landscape.reproduce.Rversion <- function(rland)
{
    start=Sys.time()
    ##put the repro matrix in a nice form
    habs <- rland$intparam$habitats
    stgs <- rland$intparam$stages
    R <- rland$demography$localdem[[1]]$LocalR
    Rlong <- as.data.frame(cbind((which(R>0,arr.ind=T)-1),R[which(R>0,arr.ind=T)]))
    colnames(Rlong)[3] <- "rate"
    for (i in 1:dim(Rlong)[1])
    {
        Rlong <- unique(rbind(Rlong,cbind(row=Rlong$row[i],col=seq(Rlong$col[i],habs*stgs,by=stgs),rate=Rlong$rate[i])))
    }
    
    if (length(unique(Rlong[,2]))!=length(Rlong[,2])) {stop("should be only one juvenile class per species")}
  

################# find males:  look for all males, rank them based on criteria and return one at random
    ##this df holds indices to mothers rows, their offspring size and fathers rows.
    ## the AllMales call could be slow...

 #   print("about to about to assemble rdf, calls AllMales2")


 #   print(length(reprows))
    partners <- AllMales4(rland,Rlong)
    fathers <- partners[,2]
    mothers <- partners[,1]
#print(Sys.time()-start)
###calculate the number of offspring per mother

###use the reproduction phenotype here to control the_mean_ of the distribution
### rather than change individual values
    repind = rland$individuals[mothers,]
    rclass <- data.frame(table(repind[,1]))
    names(rclass)[1] <- "col"
    rclass <- merge(rclass,Rlong)
    
    noff <- rep(0,dim(repind)[1])
    
    for (i in 1:dim(Rlong)[1])
    {
        noff[repind[,1]==Rlong$col[i]] <- rpois(sum(repind[,1]==Rlong$col[i]),Rlong$rate[i])
    }

    
    rdf <- data.frame(mother=mothers, noff=noff, father=fathers)
    rdf$noff[rdf$father==0] <- 0
    rdf <- rdf[rdf$noff>0,]
    
#    print("assembled rdf")
    ##sex.  there is no physical linkage in this setup  probably should be separate function
    ## but trying for simplicity and speed, also figure out offspring coords here...

    lv=landscape.locusvec(rland)

    offmat = cbind(mother=inverse.rle(list(lengths=rdf$noff,values=rdf$mother)),
                   father=inverse.rle(list(lengths=rdf$noff,values=rdf$father)))
    offmat= cbind(offmat,
                  stg=rland$individuals[offmat[,1],1]+1,
                  dist=0,
                  direction=runif(dim(offmat)[1],0,2*pi))
    
    seedkern=rland$demography$epochs[[1]]$seedkern
    
    for (i in 1:dim(seedkern)[1])
        if (sum(offmat[,3]==(i))>0)
            offmat[offmat[,3]==(i),4] <- abs(rdispfunc(sum(offmat[,3]==(i)),seedkern[i,]))
    
    offmat = cbind(offmat,
                   newx=rland$individuals[offmat[,1],4]+cos(offmat[,5])*offmat[,4],
                   newy=rland$individuals[offmat[,1],5]+sin(offmat[,5])*offmat[,4],
                   momx=rland$individuals[offmat[,1],4],
                   momy=rland$individuals[offmat[,1],5],
                   dadx=rland$individuals[offmat[,2],4],
                   dady=rland$individuals[offmat[,2],5]                   
                   )
    newind = Inheritance(offmat,as.matrix(rland$individuals),lv,rle(lv)$lengths);
    newind[,1]=0
    newind[,4:9]=offmat[,6:11]

    
  topy <- rland$demography$epochs[[1]]$topy
  boty <- rland$demography$epochs[[1]]$boty
  leftx <- rland$demography$epochs[[1]]$leftx
  rightx <- rland$demography$epochs[[1]]$rightx

#    print(dim(newind))

    inx=(newind[,4]>=min(leftx))&(newind[,4]<=max(rightx))
    iny=(newind[,5]>=min(boty))&(newind[,5]<=max(topy))

    if (sum(inx&iny)>0)
    {
#        print(sum(!(inx&iny)))
        newind <- newind[inx&iny,] #this removes individuals

    } 

#    print(dim(newind))
    
    rland$individuals <- as.matrix(rbind(rland$individuals,newind))
    rownames(rland$individuals) <- 1:dim(rland$individuals)[1]
    rland
}


###
###
landscape.advance.Rversion <-function(Rland)
  {
      if (is.landscape(Rland))
      {
          Rland$intparam$currentgen <- Rland$intparam$currentgen+1
          Rland
      }
      else
      {
          stop("Rland not a landscape object in landscape.advance...exiting")
      }
  }
###
### applies the randomly choose individuals within a pop and kill them 
### to reach carrying capacity.  so, density dependence with hard K.
###
landscape.carry.Rversion <-function(Rland)
  {
    if (is.landscape(Rland))
      {
          Ks <- data.frame(pops=1:Rland$intparam$habitats,
                           k=Rland$demography$epoch[[1]]$Carry)
          
          pops <- landscape.populations(Rland)
           
          tdf <- as.data.frame(table(pops)) %>% mutate(pops=as.numeric(as.character(pops))) %>%
                               arrange(pops) %>% left_join(Ks,"pops")
          tdf <- tdf[tdf$k/tdf$Freq<1,]
          kill <- rep(FALSE,length(pops))
          for (i in tdf$pops)
          {
              kill[pops==i] <- runif(sum(pops==i))>(tdf$k[tdf$pops==i]/tdf$Freq[tdf$pops==i])
          }
          
          Rland$individuals <- Rland$individuals[!kill,]
          Rland
      }
    else
      {
        print("Rland not a landscape object in landscape.carry...exiting")
      }
  }




###
### takes every location and kills everybody if a random number comes up
###
"landscape.extinct.Rversion" <- function(Rland)
  {
    if (is.landscape(Rland))
      {
          e <- Rland$demography$epochs[[1]]$Extinct
          kp <- which(e>runif(length(e))) #vec of pops to kill results
          pops <- landscape.populations(Rland)
          Rland$individuals <- Rland$individuals[!(pops%in%kp),]
          Rland
      }
    else
      {
        print("Rland not a landscape object in landscape.extinct...exiting")
      }
  }

###
### take geographic coordinates and find population.
### 
###
landscape.geopops <- function(rland)
{
    leftx <- rland$demography$epochs[[1]]$leftx
    rightx <- rland$demography$epochs[[1]]$rightx
    boty <- rland$demography$epochs[[1]]$boty
    topy <- rland$demography$epochs[[1]]$topy

    x <- rland$individuals[,4]
    y <- rland$individuals[,5]
    
    np <- rep(-1,dim(rland$individuals)[1])
    
    for (i in 1:length(leftx))
    {
        np[which( ((x>=leftx[i])&(x<=rightx[i])) &
               ((y>=boty[i])&(y<=topy[i])) )] <- i
    }   
    np
}

###adds new population numers to stages
landscape.popassign <- function(rland)
{
    npops <- landscape.populations(rland)
    gpops <- landscape.geopops(rland)
    inds <- rland$individuals
    rootstages <- unique((0:(rland$intparam$stages*rland$intparam$habitats)) - ((0:(rland$intparam$stages*rland$intparam$habitats)) %% rland$intparam$stages))
    
    if (sum(npops!=gpops)>0)
    {
        di = which(npops!=gpops)
        inds[di,1] <- rootstages[gpops] + (inds[di,1] %% rland$intparam$stages)
        rland$individuals <- inds
    }
    rland
    
}


###
### 
### 
###
landscape.simulate.Rversion <- function(rland,clicks=1)
{
    for ( i in 1:clicks )
        rland <- landscape.extinct.Rversion(rland) %>%
            landscape.survive.Rversion() %>%
            landscape.carry.Rversion() %>%
            landscape.reproduce.Rversion() %>%
            landscape.advance.Rversion()
    rland
}


