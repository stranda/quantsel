##
## here is a function that returns the coordinates of each population in a landscape
## we now depend on dplyr
##
landscape.popcoords <- function(rland)
{
    as.data.frame(data.frame(pop=landscape.populations(rland),
               x=rland$individuals[,4],
               y=rland$individuals[,5]) %>%
        group_by(pop) %>% summarize(x=mean(x),y=mean(y)))
    
}

##makes a data frame of phenotypic means and variances for each population phenotype comparison
landscape.phenosummary  <- function(l)
{
    if ((is.landscape(l))&&(l$intparam$nphen>0))
    {
        phens <- as.data.frame(cbind(pop=landscape.populations(l),landscape.phenotypes.c(l)))
        names(phens)[-1] <- paste0("phen",1:l$intparam$nphen)
        phensum <- left_join(group_by(phens,pop) %>% summarise_all(.funs=c("mean","sd")),
                             group_by(phens,pop) %>% summarise(n=n()))
        data.frame(full_join(data.frame(pop=1:l$intparam$habitats),phensum)%>%arrange(pop))
    }  else {NULL}
}


###this function takes a list of items that describe calcuation of summary statistics (essntially stats functions)
landscape.gensummary <- function(l,analyses=list(He=list(format=1, #1 is landscape 2 is genetics::genotype
                                                         name="ExpHet",
                                                         func=function(l){
                                                             as.data.frame(cbind(pop=unique(landscape.populations(l)),landscape.exp.het(l)))
                                                         }
                                                         ),
                                                 LD=list(format=1, #1 is landscape 2 is genetics::genotype
                                                         name="LD",
                                                         func=function(l){
                                                             as.data.frame(landscape.LD(l))
                                                         }
                                                         ))
                                 )
{
    if (is.landscape(l))
    {
        lst=lapply(analyses, function(x)
        {
            if (x$format==1) res = x$func(l)
            else if (x$format==2) res = x$func(g)
            else stop("specify a correct input format for the statistic")
            res$stat=x$name
            res
        })
        
    } else {NULL}
}

###the genetics package genotypes
### returns a list of genotypes
landscape2genotype <- function(l)
{
    lapply(which(landscape.ploidy(l)==2),function(loc)
    {
        tloc <- matrix(landscape.locus(loc,l)[,-1:-9],ncol=2)
        tloc2=cbind(ifelse((tloc[,1]>tloc[,2]),tloc[,2],tloc[,1]),
                    ifelse((tloc[,1]>tloc[,2]),tloc[,1],tloc[,2]))
        tloc2f=factor(paste(tloc2[,1],tloc2[,2],sep="/"))
        attr(tloc2f,"allele.names")=as.character(unique(c(tloc)))
        attr(tloc2f,"allele.map")=do.call(rbind,strsplit(levels(tloc2f),"/"))
        attr(tloc2f,"genotypeOrder")=as.character(c("1/1","1/2","2/1","2/2"))
        class(tloc2f) <- "genotype"
        class(tloc2f) <- append(class(tloc2f),"factor")
        tloc2f
    })
}


landscape.LD <- function(l)
{

    ret <- NULL
    ##rsq is standardized with the overall allele freqs
    afreqs <- matrix(NA,nrow=sum(landscape.ploidy(l)==2),ncol=2)
    locs <- list()
    pops=landscape.populations(l)
    for (i in which(landscape.ploidy(l)==2))  ###be careful because could have haploid low-numbered loci (probably should just toss the haploids
    {
        locs[[i]] <- landscape.locus(i,l)[,-1:-9]
        tbl= table(locs[[i]][,1])
        afreqs[i,] <- tbl/sum(tbl)
    }

    ret <- do.call(rbind,lapply(which(landscape.ploidy(l)==2),function(i)
    {
        rdf=NULL
        for (j in 1:i)
            if (j<i)
                if (((afreqs[i,1]>0)&(afreqs[i,1]<1))&((afreqs[j,1]>0)&(afreqs[j,1]<1))) #polymorphism
                {
                    df <- data.frame(pop=pops,
                                     ava=locs[[i]][,1],
                                     eva=locs[[i]][,1]==locs[[j]][,1],
                                     gva=locs[[i]][,1]>locs[[j]][,1],
                                     avb=locs[[i]][,2],
                                     evb=locs[[i]][,2]==locs[[j]][,2],
                                     gvb=locs[[i]][,2]>locs[[j]][,2]
                                     )
                    mns <- df%>%mutate(p11a=(ava==1)&(eva&(!gva)),
                                       p22a=(ava==2)&(eva&(!gva)),
                                       p12a=((!eva)&(!gva)),
                                       p21a=((!eva)&(gva)),
                                       p11b=(avb==1)&(evb&(!gvb)),
                                       p22b=(avb==2)&(evb&(!gvb)),
                                       p12b=((!evb)&(!gvb)),
                                       p21b=((!evb)&(gvb))) %>%
                        group_by(pop)%>%
                        summarise(p11a=mean(p11a),p12a=mean(p12a),p21a=mean(p21a),p22a=mean(p22a),
                                  p11b=mean(p11b),p12b=mean(p12b),p21b=mean(p21b),p22b=mean(p22b),
                                  n=n())
                    mns = mutate(mns,Dv=( (p11a*p22a)-(p12a*p21a)+(p11b*p22b)-(p12b*p21b) )/2)
                    rdf <- rbind(rdf,mutate(mns,
                                            D=Dv,
                                            rsq=(Dv/sqrt(afreqs[i,1]*afreqs[j,1]*afreqs[i,2]*afreqs[j,2]))^2,
                                            loc1=i,
                                            loc2=j)
                                 )
                }
        rdf
    }))
    ret[apply(ret[,2:5],1,function(v)max(v)<1),]
}




###
###
###
landscape.dist2origin <- function(l,origin=NULL)
{
    if (is.null(origin)) {origin <- sample(unique(landscape.populations(l)),1)}
    crd=data.frame(landscape.popcoords(l))
    rownames(crd) <- crd$pop
    d <- as.matrix(dist(crd[,c("x","y")],diag=T,upper=T))
    df=as.data.frame(cbind(pop=unique(landscape.populations(l)),
                  d=d[rownames(d)==as.character(origin),]))
}

landscape.neighborhood <- function(rland)
{
    mo  <- landscape.geodist(rland,cmp=c("offspring","father"))
    fo  <- landscape.geodist(rland,cmp=c("offspring","father"))
    mp2 <- ((mo+fo)/2)^2
    df  <- data.frame(pop=landscape.populations(rland),mp2=mp2) %>%
        group_by(pop)%>%summarise(sigma2 = mean(mp2,na.rm=T))
    area=((rland$demography$epochs[[1]]$rightx-rland$demography$epochs[[1]]$leftx) *
          (rland$demography$epochs[[1]]$topy-rland$demography$epochs[[1]]$boty))[as.numeric(as.character(df$pop))]
    df <- as.data.frame(df)
    df$D =  rland$demography$epochs[[1]]$Carry[as.numeric(as.character(df$pop))]/area
    
    df$Nn = c(4*pi*df[,"sigma2"]*df$D)
    df[,c("pop","Nn")]
}
