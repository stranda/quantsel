makealleles <- function(type,numalleles,allelesize,frequencies,states)
{
  retval <- 0

  if(is.null(frequencies))
    {
      frequencies <- rep(1.0/numalleles, numalleles)
    }

  if(length(frequencies) != numalleles)
    {
      stop("Frequency list is not the right size")
    }
  
  if(type == 0 || type == 1)
    {
      retval <- vector("list", numalleles)
      for (x in 1:numalleles)
        {
          retval[[x]]$aindex <- as.integer(x)
          retval[[x]]$birth <- as.integer(0)
          retval[[x]]$prop <- frequencies[x]
          if (is.null(states))
            {
              retval[[x]]$state <- as.integer(x)
            } else
          {
            retval[[x]]$state <- as.integer(states[x])
          }
        }
    }
  else if(type == 2)
    {
      retval <- vector("list", numalleles)
      for (x in 1:numalleles)
        {
          retval[[x]]$aindex <- as.integer(x)
          retval[[x]]$birth <- as.integer(0)
          retval[[x]]$prop <- frequencies[x]
          if (is.null(states))
            {
              retval[[x]]$state <- geneseq(allelesize)
            } else {
              retval[[x]]$state <- states[x]
            }
        }
    }
  retval          
}
