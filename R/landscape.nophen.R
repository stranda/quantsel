#' creates a dummy set of quantgen params to allow marker-only simulation (nphen has to be 0)
#' @param rland a quantsel landscape object.  Assumes no phneotypes
#' @export
landscape.nophen <- function(rland)
{
    if (rland$intparam$nphen!=0) stop("This landscape has nphen >0; do you want to make a dummy set?")
    nl=length(rland$loci)

    expmat <- matrix(rep(0,nl),  #16rows for 16 loci, 4 cols for 4 phenotypes
                     byrow=T,ncol=1)
    hsq <- c(1)
    rland <- landscape.new.expression(rland,
                                      expmat=expmat*0.125, #0.125 -> 1 diploid locus per phen, 
                                      hsq=hsq) #up to 8 alelle additive doses, when summed across 4 loci.  
    rland <- landscape.new.gpmap(rland,
                                 ## 1 cols 8 rows.  Cols correspond to phenotype effects on fit components
                                 ##for each phenotype (0 is no effect, 1 phenotypes in this example)
                                 ##phenotypes are in C indexing so, add 1 to compare to pehnotypes above
                                 matrix(c(0,     #short scale #no selection
                                          0,     #long scale 
                                          0,     #long shape  #no selection
                                          0,     #mixture   #phenotype 2   #no selection
                                          0,     #not used  #no selection
                                          0,     #survive  
                                          0,     #reproduce 
                                          0      #density tolerance  #no selection
                                          ),
                                        ncol=1,byrow=T)
                                 )
    
    rland <- landscape.new.plasticity(rland)
    
    rland <- landscape.new.phenohab(rland)
    rland
}
