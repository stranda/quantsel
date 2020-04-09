
linecirc <- function(x0,y0,x1,y1,xc,yc,r)
{
    found=0
###x first
    if ((xc<x0)&((xc+r)>x0)) found=1;  #circle bounds from left of cell
    if ((xc>x1)&((xc-r)<x1)) found=1;  #circle bounds from right of cell

    if ((yc<y0)&((yc+r)>y0)) found=1;  #circle bounds from bottom of cell
    if ((yc>y1)&((yc-r)<y1)) found=1;  #circle bounds from top of cell

}
