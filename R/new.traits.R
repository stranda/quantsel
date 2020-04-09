###
### make a "plasticity landscape" that encodes the effect of environment on phenotypes
###
### create a landscape sub-object that holds matrices for each phenotype
### matrices have  rows=number of habitats
###                cols=Single col that has factor to multiply x phenotype
###
###  
landscape.new.plasticity <- function(rland,pm=NULL)
{
    print("running landscape.new.plasticity");
    if (length(rland$plasticity)<1)
    {
        rland$plasticity <- matrix(rep(1,rland$intparam$habitats*rland$intparam$nphen),ncol=rland$intparam$nphen)
    }
    
    if (ncol(rland$plasticity)==rland$intparam$nphen)
    {
        if (is.null(pm))
        {
            rland$plasticity <- matrix(rep(1,rland$intparam$habitats*rland$intparam$nphen),ncol=rland$intparam$nphen)
        } else { #is null
            rland$plasticity <- pm  #need error checking 
        }
    } else { #plasticiy list length
            stop ("mismatch between plasticity list and nphen")        
    }
rland
}

###
### make a phenotype*habitat map onto fitness.  Effectively this means that
### fitness surfaces for each phenotype can be different for each habitat
### this allows certain phenotypes to be better matched to some habitats and
### not to others
###
### create a landscape sub-object that holds matrices for each fitness component
### matrices have  rows = number of habitats
###                cols = 1) param k from logistic
###                       2) param range to vary across
###                       3) direction of slope
###
###  
landscape.new.phenohab <- function(rland,fitcomp=1,ph=NULL)
{
    if (length(rland$phenohab)!=8)
    {
        rland$phenohab <- vector("list",8)
        for (i in 1:8)
            rland$phenohab[[i]] <- cbind(rep(1,rland$intparam$habitats), # a (must be>=1)
                                         rep(1,rland$intparam$habitats), # b(must be>=1)
                                         rep(0,rland$intparam$habitats), # range around 1 for multiplier
                                         rep(0,rland$intparam$habitats))   # (0 = beta) vs (1 = beta inverted)

    } else {
        
        if (fitcomp<=8)
        {
            if (is.null(ph))
            {
            rland$phenohab[[fitcomp]] <- cbind(rep(1,rland$intparam$habitats), # a (must be>=1)
                                         rep(1,rland$intparam$habitats), # b(must be>=1)
                                         rep(0,rland$intparam$habitats), # range around 1 for multiplier
                                         rep(0,rland$intparam$habitats)) # (0 = beta) vs (1 = beta inverted)
            }
            else #is null
            {
                rland$phenohab[[fitcomp]] <- ph
            }
        } else #plasticiy list length
        { 
            stop ("wrong length list of fitness components")        
        }
    }
    rland
}
     

