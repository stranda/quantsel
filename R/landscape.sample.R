#returns a landscape with a sample of populations and or individuals within populations
#np is the number of populations, if NULL then all populations, pvec gives the populations to sample
#ns is the number of individuals to sample
landscape.sample.old <- function(rland,np=NULL,ns=NULL,pvec=NULL)
  {
    
    if (!is.null(np))
      {
        pops <- sample(unique(landscape.populations(rland)),
                       ifelse(length(unique(landscape.populations(rland)))>np,np,
                              length(unique(landscape.populations(rland)))),
                       ,replace=F)
        rland$individuals <- rland$individuals[landscape.populations(rland) %in% pops,]
      }
    else if (!is.null(pvec))
      {
        rland$individuals <- rland$individuals[landscape.populations(rland) %in% pvec,]
      }
    if (!is.null(ns))
      {
        pops <- landscape.populations(rland)
        ptbl <- table(pops)
        if ((is.null(pvec) & is.null(np)))
          {
            names <- as.numeric(names(ptbl))
           } else {
             names <- as.numeric(names(which(ptbl>=np)))
           }
        rland$individuals <- rland$individuals[
                                               as.numeric(unlist(sapply(unique(names),
                                                                        function(x,pops,ns)
                                                                        {
                                                                          ss <- ifelse(length(which(pops==x))>ns,
                                                                                       ns,length(which(pops==x)))
                                                                          sample(which(pops==x),ss,replace=F)
                                                                        },pops=pops,ns=ns)))
                                               ,]
      }
    rland
  }





#' simulates sampling for genetics on the landscape
#' 
#' 
#' Randomly pulls a max of \code{ns} individuals from a max of \code{np}
#' populations and returns a landscape object that could be used for further
#' simulation, but is usually used for analyses and summary statistics
#' calculatiuons
#' 
#' 
#' @param rland landscape object
#' @param np number populations
#' @param ns number samples per population
#' @param pvec a vector of populations to sample.  Should be numbers from 1 to
#' number of habitats
#' @return landscape object
#' @keywords misc
#' @examples
#' 
#' 	l <- landscape.new.example()
#' 	l <- landscape.simulate(l,1)
#' 	l.samp <- landscape.sample(l,np=3,ns=24)
#' 	landscape.amova.pairwise(l.samp)
#' 
#' @export landscape.sample
landscape.sample <- function(rland,np=NULL,ns=NULL,pvec=NULL)
  {
    
    if (!is.null(np))
      {
        pops <- sample(unique(landscape.populations(rland)),
                       ifelse(length(unique(landscape.populations(rland)))>np,np,
                              length(unique(landscape.populations(rland)))),
                       ,replace=F)
        rland$individuals <- rland$individuals[landscape.populations(rland) %in% pops,]
      }
    else if (!is.null(pvec))
      {
        rland$individuals <- rland$individuals[landscape.populations(rland) %in% pvec,]
      }
    if (!is.null(ns))
      {
        pops <- landscape.populations(rland)
        ptbl <- table(pops)
        if ((is.null(pvec) & is.null(np)))
          {
            names <- as.numeric(names(ptbl))
           } else {
             names <- as.numeric(names(which(ptbl>=np)))
           }
        rland$individuals <- rland$individuals[
                                               as.numeric(unlist(sapply(unique(names),
                                                                        function(x,pops,ns)
                                                                        {
                                                                          ss <- ifelse(length(which(pops==x))>ns,
                                                                                       ns,length(which(pops==x)))
                                                                          sample(which(pops==x),ss,replace=F)
                                                                        },pops=pops,ns=ns)))
                                               ,]
      }
    rland
  }
