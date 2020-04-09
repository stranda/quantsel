/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim
*/

/**
   
*/


/*includes
*/

#include <RandLib.h>
#include <R.h>
#include <Rmath.h>

/*
 */

RandLib::RandLib () 
{
  init();  
}


void RandLib::init()
{
    GetRNGstate();
}


RandLib::~RandLib () 
{
    PutRNGstate();
}


void RandLib::FreeDiscreteLookup()
{
}

int RandLib::CheckDiscreteLookup()
{
  return 1;
}

void RandLib::SetDiscreteLookup(double *p, int ncat) //returns a randomly chosen
                              //outcome from a multinomial distribution given in array
                              //p, of length n
     {
       double tot=0.0;
       int i;

       lp.resize(ncat);

       for (i=0;i<(ncat-1);i++)
	 {
	   tot = tot + p[i];
	   lp[i]=p[i];
	 }
       ///do different things if the passed vector totals to less than, equal to, or greater than 1
       if ((tot>=0) && (tot < 1)) //the total of the passed vector on [0,1)
	 {
	   lp[ncat-1]=1-tot;
	 }
       else if (tot==1)
	 {
	   lp[ncat-1]=0;
	 }
       else if (tot > 1)
	 {
	   if (tot > 1.5) ///the total of the vector is really big.  something is wrong.
	     {
	       cerr <<"In Randlib.cc, the total of a vector passed to multinomial is much greater than 1:  "<<tot<<endl;
	       assert(tot<1);
	     }
	   else //scale all the probs to sum to one
	     {
	       for (i=0;i<(ncat);i++)
		 {
		   lp[i]= p[i]/tot;
		 }
	       lp[ncat-1]=0;
	     }
	 }
     }

int RandLib::PickMultinomial()
{
  int i;
  int *resvec = new int[lp.size()];
  double *pvec = new double[lp.size()];
  for (i=0;i<int(lp.size());i++)
    {
      pvec[i]=lp[i];
    }
  rmultinom(1,pvec,int(lp.size()),resvec);
  i=0;
  while (resvec[i]<1)
    {
      i++;
    }
  delete [] pvec;
  delete [] resvec;
  return i;
}

int RandLib::multinomial(double *p, int ncat)
{
  int rv;

  SetDiscreteLookup(p,ncat);
  rv = PickMultinomial();
  return rv;
}


///this function should return numbers from the range inclusive  of the endpoints
int RandLib::unirange(int maxval)
{
  return uniminmax(0,maxval);
}

int RandLib::uniminmax(int minval, int maxval)
{
  int rv;
  rv=int(fround(runif(minval,maxval),0));
  return rv;
}

double RandLib::uniform()
{
  return runif(0.0,1.0);
}

int RandLib::poisson(double mu)
{
  int rv;
  //  rv=int(round(rpois(mu)));
  rv=int(fround(rpois(mu),0));
  return rv;
}


double RandLib::rndweibull(double sc, double sh)
{
  return rweibull(sh,sc);
}


double RandLib::rndgeom(double sc, double sh)
{
  return pow(2, 1/(1 - sh))*pow((pow(sc, 1 - sh)*uniform()*M_PI)/(-2 + sh), 1/(1 - sh)) - sc;
}

double RandLib::rndmixed(double mu1, double mu2, double sd1, double sd2, double mix)
{
    double mt = uniform();
    if (mt<mix)
      {
	return  rweibull(mu1,sd1);
      }
    else
      {
	return  rnorm(mu2,sd2);
      }

}

void RandLib::rnegexp_xy(double ix, double iy, double mu, double aspect, double &newx, double &newy)
{
  double dir = uniform()* 2 * M_PI;
  double dist = rexp(mu);
  double cflag = 1.0;
  double sflag = 1.0;
  if (dir>M_PI) {cflag=-1.0;}
  if (((dir>M_PI/2)&(dir<=M_PI))|((dir>1.5*M_PI)&(dir<2*M_PI))) {sflag=-1.0;}
  newy = (sin(dir)*aspect*dist)*sflag + iy;
  newx = pow(1-pow((sin(dir)*aspect),2),0.5)*dist*cflag + ix;
}


void RandLib::rweibull_xy(double ix, double iy, double sc, double sh, double aspect, double &newx, double &newy)
{
  double dir = uniform()* 2 * M_PI;
  double dist = rweibull(sh,sc);
  double cflag = 1.0, sflag = 1.0;
  if (dir>M_PI) {cflag=-1.0;}
  if (((dir>M_PI/2)&(dir<=M_PI))|((dir>1.5*M_PI)&(dir<2*M_PI))) {sflag=-1.0;}
  newy = (sin(dir)*aspect*dist)*sflag + iy;
  newx = pow(1-pow((sin(dir)*aspect),2),0.5)*dist*cflag + ix;

}

void RandLib::rgeom_xy(double ix, double iy, double sc, double sh, double aspect, double &newx, double &newy)
{
  double a = sc;
  double b = sh;
  double dir = uniform()* 2 * M_PI;
  double m = uniform();
  double dist = pow(2, 1/(1 - b))*pow((pow(a, 1 - b)*m*M_PI)/(-2 + b), 1/(1 - b)) - a;
  double dflag = 1.0;
  if (dir>M_PI) {dflag=-1.0;}
  newy = (sin(dir)*aspect*dist) + iy;
  newx = pow(1-pow((sin(dir)*aspect),2),0.5)*dist*dflag + ix;

}

void RandLib::rmixed_xy(double ix, double iy, double mu1, double mu2, double sd2, double mix,
			double aspect, double &newx, double &newy)
{
  double dir = uniform()* 2 * M_PI;
  double dist ;
  double mt = uniform();
  double dflag = 1.0;
  if (dir>M_PI) {dflag=-1.0;}

  /*
 if (mt<mix)
    {
      dist = rexp(mu1);
      newy = (sin(dir)**dist) + iy;
      newx = pow(1-pow((sin(dir)),2),0.5)*dist*dflag + ix;
    }
  else
    {
      dist = rnorm(mu2,sd2);
      newy = (sin(dir)*aspect*dist) + iy;
      newx = pow(1-pow((sin(dir)*aspect),2),0.5)*dist*dflag + ix;
    }
  */
  if (mt>=mix) //mix is the prop of ldd
    {
      dist = rexp(mu1);
      newy = (sin(dir)*dist) + iy;
      newx = (cos(dir)*dist) + ix;
    }
  else
    {
      dist = rnorm(mu2,sd2);
      newy = (sin(dir)*dist) + iy;
      newx = (cos(dir)*dist) + ix;
    }
}
void RandLib::rassym_mixed_xy(double ix, double iy, 
			     double nmu1, double nmu2, 
			     double mu1, double mu2, 
			     double nsh1, double nsh2, 
			     double sh1, double sh2, 
			     double nmix,double mix, 
			     double aspect, 
			     double &newx, double &newy)
{
  double dir = uniform() * 2 * M_PI;
  double dist ;
  double mt = uniform();
  double cflag = 1.0, sflag=1.0;
  if (dir>M_PI) {cflag=-1.0;}
  if (((dir>M_PI/2)&(dir<=M_PI))|((dir>1.5*M_PI)&(dir<2*M_PI))) {sflag=-1.0;}
  //  cerr << "mix "<<mix<<" nmix "<<nmix<<" mt "<< mt <<endl;
  //  cerr << "sh1 "<<sh1<<" mu1 "<<mu1<<endl;

  if ((dir>0) & (dir<=3.1415)) ///moving towards the right
    {
      if (mt<mix) //weibull
	{
	  //	  cerr << "heading right weibull"<<endl;
	  //cerr << "w";
	  dist = rweibull(sh1,mu1);
	}
      else //normal
	{
	  //	  cerr << "heading right norm"<<endl;
	  //cerr << "n"; 
	  dist = -1;
	  while (dist<0)
	    dist =  rnorm(mu2,sh2) ;
	}
    }
  else ///moving towards the left
    {
      if (mt<nmix) //exponetial
	{
	  //	  cerr << "heading left weibull"<<endl;
	  //cerr << "w";
	  dist = rweibull(sh1,mu1);
	}
      else //normal
	{
	  //	  cerr << "heading left norm"<<endl;
	  //cerr << "n";
	  dist = -1;
	  while (dist<0)
	    dist =  rnorm(mu2,sh2) ;
	}
    }
  newy = (sin(dir)*aspect*dist)*sflag + iy;
  newx = pow(1-pow((sin(dir)*aspect),2),0.5)*dist*cflag + ix;

  //  cerr <<"dist "<<dist<<" dir "<<dir<<" ix "<<ix<<" newx "<<newx<<" iy "<<iy<< " newy "<< newy << endl;
}

double RandLib::normal(double mu, double sd)
{
  return rnorm(mu,sd);
}

double RandLib::weibull(double d, double sc, double sh)
{
  return dweibull(d,sh,sc,0);
}

double RandLib::geom(double d, double sc, double sh)
{
  double a = sc;
  double b = sh;
  return ((-2 + b)*(-1 + b))/(2.*a*3.1415*pow(1 + d/a, b));
}

double RandLib::negexp(double dist, double mu)
{
  return dexp(dist,mu,0);
}

double RandLib::mixed(double dist, double sc1, double sc2, double sh1, double sh2, double mix)
{
  return (dweibull(dist,sh1,sc1,0)*mix + (1-mix)*dnorm(dist,sc2,sh2,0));
}

void RandLib::SetSeed(long int sd)
{
  //  gsl_rng_set(sd);
}


/**
Declaration  of a global RandLib 'RandLibObj';
 */

RandLib RandLibObj;


/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/



