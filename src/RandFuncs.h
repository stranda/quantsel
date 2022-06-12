/**

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

Instantiates an object that implements a GNU Scientific Library RNG

set_seed() takes a single long and sets the seed.

*/

#ifndef RANDFUNCS_H
#define RANDFUNCS_H

/*includes
*/

#include <metasim.h>


using namespace std;

void getRNDstate();
void putRNDstate();

///This function picks a number from the lookup table
int PickMultinomial(std::vector< double > pv);

/// returns a number from 0 to n from a multinomial distribution.  *p points to an array of class proportions.  
/// n is the size of the class.
int multinomial(double *p, int ncat);

///returns an integer on the interval between 0 and maxval at random from a uniform dist.
int unirange(int maxval=1);
///returns an integer on the interval between 0 and maxval at random from a uniform dist.
int uniminmax(int minval, int maxval);
///returns a value between 0 and 1 from a uniform distribution
double uniform();

///returns a value from a Poisson Distribution with mean equal to mu.
int poisson(double mu);

////the next three functions return random pulls from different pdfs.  
//// they are useful for getting a random distance.
///
//returns a value from a weibull distribution
double rndweibull(double sc, double sh);  

//returns a value from a geometric distribution
double rndgeom(double sc, double sh);  

//returns a value from a mixed distribution
double rndmixed(double mu1, double mu2, double sd1, double sd2, double mix);  

///Takes a xy coordinate and returns a new xy coordinate derived from 
///choosing a direction uniformly and a distance based upon a negative exponential 
///distribution 
void rnegexp_xy(double ix, double iy, double mu, double aspect, double &newx, double &newy);

///weibull
void rweibull_xy(double ix, double iy, double sc, double sh, double aspect, double &newx, double &newy);

///geometric distribution ala Clark 1999 (also in Austerlitz 2004)
void rgeom_xy(double ix, double iy, double sc, double sh, double aspect, double &newx, double &newy);

///
///generates data
///
void rmixed_xy(double ix, double iy, double mu1, double mu2, double sd2, double mix, double aspect, double &newx, double &newy);
void rassym_mixed_xy(double ix, double iy, double nmu1, double nmu2, double mu1, double mu2, 
		     double nsh1, double nsh2, double sh1, double sh2, 
		     double nmix,double mix, double aspect, double &newx, double &newy);


double normal(double mu=0, double sd=1);
double weibull(double d, double sc, double sh);
double mixed(double dist, double sc1, double sc2, double sh1, double sh2, double mix);


/**/
#endif //randfuncs.h


/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
