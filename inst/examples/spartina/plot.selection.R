###
### quick little function to take the parameters for selection and
###    display the selection gradient

plot.selection <- function(a=1,b=1,r=0,d=0)  #these are the four parameters in the phenohab object
{
    if ((a<1) | (b<1)) stop("a and b need to be >=1")
    if (r>2) stop("r must range between 0 and 2")
    dom <- seq(0,1,0.05)
    db <- dbeta(dom,a,b)
    ndb <- db/max(db,na.rm=T)
    pdb <- ndb*r + (1-(r/2))
    if (d==1)
    {
        pdb = -1*(pdb-1)+1
    }
    plot(x=dom,y=pdb,type="l",xlab="Phenotype value",ylab="relative fitness",main=paste0("Surface for a,b,r,s: ",a,", ",b,", ",r,", ",d),ylim=c(0,2))
    
}

