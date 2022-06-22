/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the site object type.


In this file the transition matrix class is defined
*/


/*
includes
*/


#include <TransMat.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>


TransMat::TransMat (size_t s )
{
  size_t lt;
  size = s;
  tm.resize(s);
  for (lt=0; lt<s; lt++)
    {
      tm[lt].resize(s);
    }
}

TransMat::~TransMat ()
{
#ifdef RDEBUG
   cout << "destructing TransMat" << "\n";
#endif
   size_t i;
   for (i=0;i<tm.size();i++)
     {
       tm[i].resize(0);
     }
   tm.resize(0);
#ifdef RDEBUG
   cout << "finished destructing TransMat" << "\n";
#endif
}



void TransMat::SetSize(size_t sz) 
{
  size_t i;

  size = sz; 

#ifdef RDEBUG
  cerr << "Resizing a first dimension of a transmatrix of size"<<sz <<endl;
#endif

  tm.resize(sz);

#ifdef RDEBUG
  cerr << "Finished resizing the first dimension" <<endl;
#endif
  for (i=0;i<sz;i++)
    {
      tm[i].resize(sz);
    }
}

///Sets an entire TransMat of size s from a 2d array pointed to by a.  The colummns should represent from and rows, to
void TransMat::SetMat(TransMat a)
{
  size_t i,j;
  for (i=0;i<a.Size();i++)
    {
      for (j=0;j<a.Size();j++)
	{
	  SetElement(j,i,a.GetElement(j,i));
	}
    }
  /// don't have an error function here probably need to check for reasonable values.
}


int TransMat::RandomState(double adj, int frm)
{
  
  vector < double > p;
  p.resize(Size()+1);
  double s=0.0;
  size_t i=0;
  for (i=0; i<(p.size()-1); i++)
    {
      //    p.push_back(tm[i][frm]);
      p[i] = GetElement(frm,i) * adj;
      //      s = s+tm[i][frm];
      s = s+p[i];;
    }
  if (s<1.0) {p[p.size()-1]=1.0 - s;} else {p[p.size()-1]=0.0;}
  int rs = PickMultinomial(p);
  //  cerr << "p.back() "<<p.back()<<", p.size() "<<p.size()<<", i: "<<i<<", rs "<<rs<<endl;
  if (rs == int(p.size()-1)) {rs=-1;}
  return rs;
}

/*
This function assumes that the local, habitat-specific demography is the only demography that
matters.  It completely ignores "off-diagonal" elements. In other words, the output state is forced
to be within the same habitat as 'frm', or -1 which signifies "dead". 
This should only be called for survive-type applications

h is the number of habitats, need that to figure out in which habitat 'frm' is located
 */

int TransMat::RandomStateLocal(const double &adj, const int &frm, const size_t &h)
{
  //frm is the column of the transition matrix
  //adj is a factor that might be unique for an individual
  //  cerr << "in randomstatelocal "<<endl;
  vector < double > p(Size()+1,0);
  int stg = Size()/h;  //number of stages per habitat
  size_t strt=frm - (frm%stg);
  double s=0.0;
  size_t i=0;
    cerr << "strt "<<strt<<", frm "<<frm<<", h "<<h<<", stg"<<stg<<endl;
  for (i=strt; i<(strt+stg); i++)
    {
      p[i]=tm[i][frm];
      s = s+tm[i][frm];
    }
  if (s<1.0) {p[Size()-1]=(1.0 - s);} else {}
  int rs = PickMultinomial(p);
  if (rs == int(p.size()-1)) {rs=-1;} else {}
  return rs;
}
/******************  OLD VERSION
int TransMat::RandomStateLocal(const double &adj, const int &frm, const size_t &h)
{
  //frm is the column of the transition matrix
  //adj is a factor that might be unique for an individual
  //  cerr << "in randomstatelocal "<<endl;
  vector < double > p;
  int stg = Size()/h;  //number of stages per habitat
  size_t strt=frm - (frm%stg);
  double s=0.0;
  size_t i=0;

    cerr << "strt "<<strt<<", frm "<<frm<<", h "<<h<<", stg"<<stg<<endl;
  
  for (i=strt; i<(strt+stg); i++)
    {
      //      cerr << "i: "<<i<<endl;
      p.push_back(tm[i][frm]);
      //p[i] = GetElement(frm,i) * adj;
      s = s+tm[i][frm];
    }
  //  cerr << "s: "<<s<<endl;
  //  for (i=0;i<p.size();i++) {cerr << p[i]<<", " ;} cerr <<endl;
  if (s<1.0) {p.push_back(1.0 - s);} else {p.push_back(0.0);}
  //  for (i=0;i<p.size();i++) {cerr << p[i]<<", " ;} cerr <<endl;
  int rs = PickMultinomial(p);
  //  cerr << "p.back() "<<p.back()<<", p.size() "<<p.size()<<", i: "<<i<<", rs "<<rs<<endl;
  if (rs == int(p.size()-1)) {rs=-1;} else {rs = rs + strt;}
  //  cerr<<"about to return from randomstatelocal"<<endl;
  return rs;
}
****/

int TransMat::AnyFrom(size_t fs)
{
  size_t i=0;
  double tot = 0.0;
  SetFromState(fs);
  while ((tot==0)&&(i<size))
    {
      SetToState(i);
      tot = tot + Value();
      i++;
    }
  return (tot>0);
}

void TransMat::Diag()
{
  size_t i,j;
  for (i=0;i<tm.size();i++)
    for (j=0;j<tm.size();j++)
      if (i==j) 
	tm[i][j]=1.0;
      else
	tm[i][j]=0.0;
}

double TransMat::Lambda()
{
  int i,j,k, n, lwork, info;
  double *work, *wR, *wI, *left, *right, *xvals, tmp, maxval;
  char jobVL[1], jobVR[1];
  
  n = int(tm.size());
  xvals = new double[n*n];
  
  k=0;
  for (i=0; i<n; i++)
    for (j=0; j<n;j++)
      {
	xvals[k]=tm[j][i];
	k++;
      }
  
  //    vectors = !ov;
  jobVL[0] = jobVR[0] = 'N';
  left = right = (double *) 0;
  wR = new double[n];
  wI = new double[n];
  /* ask for optimal size of work array */
  lwork = -1;
#ifdef HAVE_LAPACK
  F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		  left, &n, right, &n, &tmp, &lwork, &info);
#else
  F77_CALL(rgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		  left, &n, right, &n, &tmp, &lwork, &info);
#endif
  if (info != 0)
    error("error code %d from Lapack routine dgeev", info);
  lwork = (int) tmp;
  work = new double[lwork];
#ifdef HAVE_LAPACK
  F77_CALL(dgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		  left, &n, right, &n, work, &lwork, &info);
#else
  F77_CALL(rgeev)(jobVL, jobVR, &n, xvals, &n, wR, wI,
		  left, &n, right, &n, work, &lwork, &info);
#endif
  if (info != 0)
    error("error code %d from Lapack routine dgeev", info);
  
  maxval=-1000000.0;
  for (i = 0; i < n; i++)
    {
      if (wI[i]==0.0)
	{
	  if (wR[i]>maxval) maxval=wR[i];
	}
      
    }
  delete work;
  delete wR;
  delete wI;
  delete xvals;

  return maxval;
}


//matrix addition operator
TransMat TransMat::operator+(TransMat TM)
{
  TransMat ret(TM.Size());
  size_t f,t;

  if (Size()==TM.Size())
    {
      for (f=0; f<ret.Size(); f++)
	for (t=0;t<ret.Size();t++)
	  {
	    ret.SetElement(f,t,(GetElement(f,t)+TM.GetElement(f,t)));
	  }
    }
  else
    {
      error("Matrices of different order in addition");
    }
  return ret;
}

//matrix multiplication operator
TransMat TransMat::operator*(TransMat TM)
{
  TransMat ret(TM.Size());
  size_t f,t,i;
  double accum;

  if (Size()==TM.Size())
    {
      for (f=0; f<ret.Size(); f++)
	for (t=0;t<ret.Size();t++)
	  {
	    accum=0;
	    for (i=0;i<ret.Size();i++)
	      {
		accum=accum+(GetElement(i,t)*TM.GetElement(f,i));
	      }
	    ret.SetElement(f,t,accum);
	  }
    }
  else
    {
      error("Matrices of different order in multiplication");
    }
  return ret;
}

ostream &operator<<(ostream &stream, TransMat & TM)
{
  size_t i,j;
  stream.precision(3);
  stream << TM.Size() << endl;
  for (i=0;i<TM.Size();i++)
    {
      TM.SetToState(i);
      for (j=0;j<TM.Size();j++)
	{
	  TM.SetFromState(j);
	  stream << TM.Value() << " ";
	}
      stream << endl;
    }
  return stream;
}

istream &operator>>(istream &stream, TransMat & TM)
{
  size_t i,j, n=199;
  stream >> n;
  TM.SetSize(n);
  for (j=0;j<n;j++)
    {
      for (i=0;i<n;i++)
	{
	  stream >> TM.tm[j][i] ;
	}
    }
  return stream;
}


ostream &operator<<(ostream & stream, DemoVec &d)
{
  int i,sz;
  sz = d.v.size();
  stream << sz << endl ;
  for (i=0;i<sz;i++)
    {
      stream << sz << " " ;
    }
  stream << endl ;
  return stream  ; 
}

istream &operator>>(istream & stream, DemoVec &d)
{
  
  int i,sz;
  int tmp;

  sz = d.size()   ;
  stream >> sz    ;
  d.resize(sz)    ;
  for (i=0;i<sz;i++)
    {
      stream >> tmp;
      d.Set(tmp,i) ;
    }

  return stream ; 
}



/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
