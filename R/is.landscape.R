#' Test whether an object is a (fairly) legitimate landscape
#' 
#' Test whether a genuine landscape
#' 
#' 
#' @param Rland the Rmetasim landscape object
#' @param verb print why not a landscape
#' @param exact more strict
#' @keywords misc
#' @examples
#' 
#'   exampleland <- landscape.new.example()
#'   is.landscape(exampleland)
#'   rm(exampleland)
#' 
#' @export is.landscape
 is.landscape <- function (Rland = NULL, verb = TRUE, exact = FALSE) 
{
    ok <- TRUE
    if (!is.list(Rland)) {
        if (verb) {
            message("Landscape not a list")
        }
        ok <- FALSE
    }
    if (is.null(Rland$intparam)) {
        if (verb) {
            message("intparam not found.")
        }
        ok <- FALSE
    }
    if (is.null(Rland$plasticity) && (Rland$intparam$nphen>0)) {
        if (verb) {
            message("plasticity not found.")
        }
        ok <- FALSE
    }
    if (is.null(Rland$phenohab) && (Rland$intparam$nphen>0)) {
        if (verb) {
            message("phenohab not found.")
        }
        ok <- FALSE
    }
    if ((length(Rland$phenohab)!=8) && (Rland$intparam$nphen>0)) {
        if (verb) {
            message("length of phenohab no equal to number of fitness components")
        }
        ok <- FALSE
    }
    if (Rland$intparam$nphen>0)
        if (sum(sapply(1:8, function(i) {Rland$intparam$habitats!=nrow(Rland$phenohab[[i]])}))>0) {
            if (verb) {
                message("phenohab matrices dont match number of habitats")
            }
            ok <- FALSE
    }
    if (is.null(Rland$switchparam)) {
        if (verb) {
            message("switchparam not found.")
        }
        ok <- FALSE
    }
    if (is.null(Rland$floatparam)) {
        if (verb) {
            message("floatparam not found.")
        }
        ok <- FALSE
    }
    if (is.null(Rland$demography)) {
        if (verb) {
            message("demography not found.")
        }
        ok <- FALSE
    } else {
        for (i in 1:length(Rland$demo$localdem)) {
            if ((min(dim(Rland$demo$localdem[[i]]$LocalS) == 
                c(Rland$intparam$stages, Rland$intparam$stages)) == 
                0) || (min(dim(Rland$demo$localdem[[i]]$LocalR) == 
                c(Rland$intparam$stages, Rland$intparam$stages)) == 
                0) || (min(dim(Rland$demo$localdem[[i]]$LocalM) == 
                c(Rland$intparam$stages, Rland$intparam$stages)) == 
                0)) {
                if (verb) {
                  message("One or more of the local demography matrices is not of the correct dimensions")
                }
                ok <- FALSE
            }
            if (max(apply(Rland$demography$localdem[[i]]$LocalS, 
                2, sum)) > 1) {
                if (verb) {
                  message(paste("Local survival matrix", i, "has a column that sums to a number greater than one"))
                }
                ok <- FALSE
            }
        }
        for (i in 1:length(Rland$demo$epochs)) {
            if ((length(Rland$demo$epochs[[i]]$Extinct) != Rland$intparam$habitats) || 
                (length(Rland$demo$epochs[[i]]$Carry) != Rland$intparam$habitats) || 
                (length(Rland$demo$epochs[[i]]$Localprob) != 
                  Rland$intparam$numdemos) || (min(dim(Rland$demo$epochs[[i]]$S) == 
                c(Rland$intparam$stages * Rland$intparam$habitats, 
                  Rland$intparam$stages * Rland$intparam$habitats)) == 
                0) || (min(dim(Rland$demo$epochs[[i]]$R) == c(Rland$intparam$stages * 
                Rland$intparam$habitats, Rland$intparam$stages * 
                Rland$intparam$habitats)) == 0) || (min(dim(Rland$demo$epochs[[i]]$M) == 
                c(Rland$intparam$stages * Rland$intparam$habitats, 
                  Rland$intparam$stages * Rland$intparam$habitats)) == 
                0)) {
                if (verb) {
                  message(paste("One or more of the epoch paramters is of incorrect dimension in epoch", 
                    i))
                }
                ok <- FALSE
            }
            if (Rland$switchparam$randdemo == 1 || length(Rland$demography$localdem) == 
                1) {
                for (j in 1:length(Rland$demography$localdem)) {
                  strt <- seq(1, dim(Rland$demo$epochs[[i]]$S)[1], 
                    dim(Rland$demography$localdem[[j]]$LocalS)[1])
                  stp <- strt + (dim(Rland$demography$localdem[[j]]$LocalS)[1] - 
                    1)
                  slice <- cbind(strt, stp)
                  tmpS <- Rland$demography$epochs[[i]]$S
                  if (Rland$intparam$habitats != 1) 
                    for (l in 1:Rland$intparam$habitats) {
                      
                      tmpS[slice[l,1]:slice[l,2],slice[l,1]:slice[l,2]] <- Rland$demography$localdem[[j]]$LocalS
                    }
                  if (max(apply(tmpS, 2, sum)) > 1) {
                    if (verb) {
                      message(paste("Columns in the landscape S matrix associated with localdem", 
                        j, "total more than 1"))
                    }
                    ok <- FALSE
                  }
                }
            }
        }
    }
    if (is.null(Rland$loci)) {
        if (verb) {
            message("loci not found.")
        }
        ok <- FALSE
    } else {
        if (length(Rland$loci) != Rland$intparam$locusnum) {
            if (verb) {
                message("conflict between size of loci object and size specified in $intparam")
            }
            ok <- FALSE
        }
    }


    if (is.null(Rland$expression) && (Rland$intparam$nphen>0)) {
        if (verb) {
            message("gpmap not found.")
        }
        ok <- FALSE
    }
    if (Rland$intparam$nphen>0)
        if (is.null(Rland$gpmap)) {
            if (verb) {
                message("gpmap not found.")
            } 
            ok <- FALSE
        } else {
            if (dim(Rland$gpmap)[2]!=Rland$intparam$nphen) {
                if (verb) message("gpdemo does not have nphen columns")
                ok <- FALSE
            }
    }

    if (is.null(Rland$individuals)) {
        if (verb) {
            message("no individuals section found.")
        }
        ok <- FALSE
    } else {
        noind <- FALSE
        if (!(dim(Rland$individuals)[1] > 0)) {
            if (verb) {
                message("No individuals in this landscape: spontaneous generation is notallowed")
            }
            ok <- FALSE
            noind <- TRUE
        }
        if (!noind) {
            if (max(Rland$individuals[, 1] > ((Rland$intparam$habitats * 
                Rland$intparam$stages) - 1))) {
                if (verb) {
                  message("There are individuals in stages greater than demography allows")
                }
                ok <- FALSE
            }
            if ((!noind) & (dim(Rland$individuals)[2] != (9 + 
                sum(landscape.ploidy(Rland))))) {
                if (verb) {
                  message("The number of loci do not correspond to the number of columns in $individuals\nThis indicates a problem with locus number and/or ploidy")
                }
                ok <- FALSE
            }
            if ((!noind) & (max(Rland$individuals[, 3]) > Rland$intparam$currentgen)) {
                if (verb) {
                  message("There are individuals born in the future")
                }
                ok <- FALSE
            }
        }
    }
    return(ok)
}
