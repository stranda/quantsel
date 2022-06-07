"geneseq" <-
function(size)
  {
    if(size <= 0)
      {
        stop("Allelesize must be positive")
      }
    retval <- floor(runif(size)*4)
    retval[retval==0] <- 'A'
    retval[retval==1] <- 'T'
    retval[retval==2] <- 'G'
    retval[retval==3] <- 'C'
    paste(retval, sep="", collapse="")
  }

