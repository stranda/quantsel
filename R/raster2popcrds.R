#' Extract resolution and extents to convert a raster to quantsel popcoords
#' @param r raster input
#' @param zeroOrigin make sure that no coordinates were negative
#' @export
raster2popcrds <- function(r, zeroOrigin=TRUE)
{
    yl = raster::extent(r)[3:4]
    xl = raster::extent(r)[1:2]
    lefts = seq(xl[1],(xl[2]-raster::res(r)[1]),length=dim(r)[2])
    bots = seq(yl[1],(yl[2]-raster::res(r)[2]),length=dim(r)[1])
    if (zeroOrigin)
    {
        if (min(lefts)<0) lefts = abs(min(lefts))+lefts
        if (min(bots)<0) bots = abs(min(bots))+bots
    }
    rights = lefts+raster::res(r)[1]
    tops = bots+raster::res(r)[2]

    locs = data.frame(cell=1:raster::ncell(r),lft=NA,rgt=NA,top=NA,bot=NA)
    cnt=1
    for (i in 1:dim(r)[1])
        for (j in 1:dim(r)[2])
        {
            locs$lft[cnt]=lefts[j]
            locs$rgt[cnt]=rights[j]
            locs$top[cnt]=tops[i]
            locs$bot[cnt]=bots[i]
            cnt=cnt+1
        }
    locs
}
