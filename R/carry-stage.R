#
#
# function to implement carrying capacity on certain stages
# the stages are specified, and all death is implemented in these
# stages.  All that is necessary is to specify the stages in the
# local demography.
#
#
landscape.carry.stage <- function(land, stages=NULL)
  {
    if (!is.null(stages)) #just return 
      {
        if (length(stages)>land$intparam$stages) stop("more stages than present in the landscape")
        for (i in 1:land$intparam$habitats)
          {
            indrows <-
              which((landscape.populations(land)==i)&((land$individuals[,1]%%land$intparam$stages)%in%stages))
#            print(i)
#            print(length(indrows))
#            print(table(land$ind[indrows,1]))
#            print(table(land$ind[indrows,1]))
#            print(table((land$ind[indrows,1]%%5)))
            if (length(indrows)>land$demography$epochs[[1]]$Carry[i])
              {
                ss=length(indrows)-land$demography$epochs[[1]]$Carry[i]
                land$individuals <- land$individuals[-sample(indrows,ss,replace=F),]
              }
          }
      }
    land
  }
