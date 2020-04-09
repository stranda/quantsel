landscape.kill.locs <- function(l,locs)
{
    for (i in 1:dim(locs)[1])
    {
        kvec=!((locs[i,1]<l$individuals[,4])&(locs[i,3]>l$individuals[,4])&
                    (locs[i,2]<l$individuals[,5])&(locs[i,4]>l$individuals[,5]))
        l$individuals <- l$individuals[kvec,]
    }
    l
}


landscape.gridlocs <- function(l,nx=1,ny=1)
{
    spatial.slices <- function(l,nx=1,ny=1)
    {
        lft=min(l$demography$epochs[[1]]$leftx)
        rt=max(l$demography$epochs[[1]]$rightx)
        bot=min(l$demography$epochs[[1]]$boty)
        top=max(l$demography$epochs[[1]]$topy)
        
        xs = seq(lft,rt,length.out=nx+1)
        ys = seq(bot,top,length.out=ny+1)
        list(xs=xs,ys=ys)
    }

 
    sl = spatial.slices(l,nx=nx,ny=ny)
    locs=matrix(0,ncol=4,nrow=nx*ny)
    cnt=1
    for (i in 1:(length(sl$xs)-1))
        for (j in 1:(length(sl$ys)-1))
        {
            locs[cnt,] <- cbind(lft=sl$xs[i],bot=sl$ys[j],
                                rgt=sl$xs[[i+1]],top=sl$ys[j+1])
            cnt <- cnt+1
        }
    
    locs
}

landscape.density.reg <- function(l,nx=1,ny=1,indPunit=6.4e-5)
{
    rlocs=landscape.gridlocs(l,nx,ny)
    for (i in 1:dim(rlocs)[1])
    {
        inds=which(
            l$individuals[,4]>rlocs[i,1] & l$individuals[,5]>rlocs[i,2] &
            l$individuals[,4]<=rlocs[i,3] & l$individuals[,5]<=rlocs[i,4]
        )
        area=abs((rlocs[i,1]-rlocs[i,3])*(rlocs[i,2]-rlocs[i,4]))

        if (indPunit*area < length(inds))
        {
            print(paste("length inds:",length(inds)))

            
            numtokill=length(inds)-(indPunit*area)
                   
            sinds = sample(inds,floor(numtokill),replace=F)
            print(paste("numtokill",numtokill))
            print(paste("rloc:",i," num to kill", length(sinds)))

            l$individuals=l$individuals[-sinds,]
        }
    }
    l
}



landscape.plot.phenotypes <- function (rland, phen, annotate=T, palletfunc=terrain.colors) 
{
    if (is.landscape(rland, FALSE))
    {
        plot(1, 1, type = "n",
             xlim = c(min(rland$demography$epochs[[1]]$leftx), 
                      max(rland$demography$epochs[[1]]$rightx)),
             ylim = c(min(rland$demography$epochs[[1]]$boty), 
                      max(rland$demography$epochs[[1]]$topy)),
             xlab = "X coordinate", 
             ylab = "Y coordinate",
             main = paste("landscape state at gen", 
                          rland$intparam$currentgen,"for phenotype",phen))
        for (i in 1:rland$intparam$habitats) {
            rect(rland$demography$epochs[[1]]$leftx[i], rland$demography$epochs[[1]]$boty[i], 
                 rland$demography$epochs[[1]]$rightx[i], rland$demography$epochs[[1]]$topy[i], 
                 lwd = 2, border = "white",col="lightgrey")
        }
        if (length(landscape.populations(rland)) > 1) {
            ph <- landscape.phenotypes.c(rland)[,phen]
            sites <- landscape.populations(rland)
            afreq <- landscape.allelefreq(rland)
            icol <- floor(ph*10)
            points(rland$individuals[, c(4, 5)], type = "p", 
                   pch = 15 + (rland$individuals[, 1] - (rland$intparam$stages * 
                                                        (landscape.populations(rland) - 1))),
                   col = palletfunc(11)[icol + 1],
                   cex = 0.6)
        }
        for (i in 1:rland$intparam$habitats)
            if (annotate)
            {
                text(x=floor(rland$demography$epochs[[1]]$rightx[i]*0.99),
                     y=floor(rland$demography$epochs[[1]]$topy[i]*0.97),
                     labels=paste("Mean phenotype: ",round(mean(ph[sites==i]),3)),
                     adj=1
                     )
                
            }           


    }
    else {
        print("no landscape to plot")
    }
}


landscape.geodist <- function(rland,cmp=c("offspring","mother","father")[1:2])
{
    if (length(unique(cmp))==2)
        {
            if (cmp[1]=="offspring") crd1=rland$individuals[,4:5]
            if (cmp[1]=="mother") crd1=rland$individuals[,6:7]
            if (cmp[1]=="father") crd1=rland$individuals[,8:9]
            
            if (cmp[2]=="offspring") crd2=rland$individuals[,4:5]
            if (cmp[2]=="mother") crd2=rland$individuals[,6:7]
            if (cmp[2]=="father") crd2=rland$individuals[,8:9]

            d=sqrt((crd2[,1]-crd1[,1])^2 + (crd2[,2]-crd1[,2])^2)
            d
            } else { stop("need two individuals to calc distance with")}
    
}

#map between loci and phenotypes for plotting, etc
landscape.phenoloc <- function(l)
{
    m=l$expression$expmat!=0
    colnames(m) <- paste0("ph",1:dim(m)[2])
    rownames(m) <- paste0("l",1:dim(m)[1])
    m
}

