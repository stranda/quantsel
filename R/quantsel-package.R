

#' Code components to simulate a landscape
#' 
#' These functions can be used to construct custom simulations of landscapes.
#' Each conducts only a single generations worth of change
#' 
#' \code{landscape.advance()} merely advances the generation counter and
#' selects the new generations demographic conditions if such conditions can
#' vary.  The other functions implement carrying capacity, local extinction,
#' reproduction, and survival/growth, respectively.  The function
#' \code{landscape.simulate()} bundles the functionality of these components
#' into a single function (and executes it slightly faster all within linked
#' C++ code).
#' 
#' @aliases landscape.advance landscape.carry landscape.extinct
#' landscape.reproduce landscape.survive
#' @param rland the Rmetasim landscape object
#' @seealso \code{landscape.simulate}
#' @keywords misc
NULL



