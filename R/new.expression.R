#probably needs a good bit of error checking written in
                                        #

###aindex is a vector of the aindices for each locus that confer an additive effect
###hsq is a vector of heritabilities (now always 1)
landscape.new.expression <-
function(rland,expmat,addstates=NULL,hsq=NULL)  
{
    if (is.null(hsq)) hsq <- rep(1,dim(expmat)[2])  #change if we change the heritabilities for each trait
                                        #if diff than 0, used to add noise
    
    if (is.null(addstates)) addstates <- rep(1,
                                           length=dim(expmat)[1])
    
    if (dim(expmat)[1]!=length(rland$loci)) {stop("expression mat should have nloc rows and nphen columns")}
    rland$expression <- list(expmat=expmat,addstates=addstates,hsq=hsq)
    ##set the nphen in intparam
    rland$intparam$nphen <- dim(expmat)[2]
    
  rland
}


#### This object controls the mapping between phenotypes and 
#### Fitness components It basically informs whether to implement
#### selection on a fitness component
#### Columns are: 1->nphen phenotypes
####              
#### The gpdemo rows correspond to: 1) ShortScale
####                                2) Long Scale
####                                3) Long Shape
####                                4) Mix
####                                5) Unused
####                                6) survival
####                                7) reproduction
####                                8) density tolerance
####
#### Each cell in this matrix corresponds to the proportional contribution of a phenotype to the
#### fitness component The rows will always be normalized (if they are not all nonzero)
####
landscape.new.gpmap <- function(rland,gpdemo=NULL)
{
    if ((rland$intparam$nphen<1)||(is.null(rland$expression)))
    {
        rland$gpmap=matrix(0,ncol=1,nrow=8)
    } else {
        if ((is.null(gpdemo))||(!is.matrix(gpdemo)))
        {
            gpdemo=matrix(0,nrow=8,ncol=rland$intparam$nphen)
        }
        rs <- rowSums(gpdemo)
        rland$gpmap <- t(sapply(1:length(rs),function(i){if(rs[i]>0) {gpdemo[i,]/rs[i]} else {gpdemo[i,]}}))
    }
    rland
}
