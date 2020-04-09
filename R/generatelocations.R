#
#
#routines to generate the locations of habitats in the landscape.  Should be able to specify
#distribution of populations and distribution of sizes approximately and boundaries of the
#landscape exactly
#

landscape.generate.locations <- function(npop=10,
                               xrange=c(0,15000),yrange=xrange,
                               sizexkernel=c(300,80),sizeykernel=sizexkernel,
                               boundaries=NULL
                               )
  {
    if(!is.null(boundaries))
      warning("boundaries not used at this point in time")

    lft <- runif(npop,min=xrange[1],max=xrange[2])
    bot <- runif(npop,min=yrange[1],max=yrange[2])
    rgt <- lft+abs(rnorm(npop,mean=sizexkernel[1],sd=sizexkernel[2]))
    top <- bot+abs(rnorm(npop,mean=sizeykernel[1],sd=sizeykernel[2]))
    overlap=TRUE
    while (overlap)
      {
        overlap=FALSE
        regen=rep(FALSE,npop)
        for (i in 1:npop)
          for(j in i:npop)
            if (i!=j)
              {
                if(rectangle.overlap(bot[i],lft[i],top[i],rgt[i],bot[j],lft[j],top[j],rgt[j]))
                  {
                    overlap <- T
                    regen[i] <- T
                  }
               }
        if (overlap)
          {
            lft[regen] <- runif(sum(regen),min=xrange[1],max=xrange[2])
            bot[regen] <- runif(sum(regen),min=yrange[1],max=yrange[2])
            rgt[regen] <- lft[regen]+abs(rnorm(sum(regen),mean=sizexkernel[1],sd=sizexkernel[2]))
            top[regen] <- bot[regen]+abs(rnorm(sum(regen),mean=sizeykernel[1],sd=sizeykernel[2]))
          }
      }
    cbind(lft,bot,rgt,top)
  }

line.intersect <- function(line1,line2)
  {
    ret=F
    if (((max(line1[,1])<=min(line2[,1]))|
         (min(line1[,1])>=max(line2[,1])))&
        ((max(line1[,2])<=min(line2[,2]))|
         (min(line1[,2])>=max(line2[,2]))))
      ret=FALSE
    else
      {
        #line 1 is horizontal
        line1=line1[order(line1[,1]),]
        line2=line2[order(line2[,2]),]
        if (((line1[1,1]<=line2[1,1])&(line1[2,1]>=line2[1,1]))&
            (((line1[1,2]>=line2[1,2])&(line1[2,1]<=line2[2,2]))))
             ret=TRUE
        #line 1 is vertical
        line1=line1[order(line1[,2]),]
        line2=line2[order(line2[,1]),]
        if (((line2[1,1]<=line1[1,1])&(line2[2,1]>=line1[1,1]))&
            (((line2[1,2]>=line1[1,2])&(line2[2,1]<=line1[2,2]))))
             ret=TRUE
      }
    ret
  }

rectangle.overlap <- function(b1,l1,t1,r1,b2,l2,t2,r2)
  {
    !(
        ((l1)<(l2) && (r1)<(l2))
        || ((l1)>(r2) && (r1)>(r2))
        || ((b1)<(b2) && (t1)<(b2))
        || ((b1)>(t2) &&
            (t1)>(t2))   )
  }
