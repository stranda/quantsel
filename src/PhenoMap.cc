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


#include <PhenoMap.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RS.h>

PhenoMap::PhenoMap (size_t rows, size_t cols )
{
  size_t lt;
  pm.resize(rows);
  for (lt=0; lt<rows; lt++)
    {
      pm[lt].resize(cols);
    }
}

PhenoMap::~PhenoMap ()
{
#ifdef RDEBUG
   cout << "destructing PhenoMap" << "\n";
#endif
   size_t i;
   for (i=0;i<pm.size();i++)
     {
       pm[i].resize(0);
     }
   pm.resize(0);
#ifdef RDEBUG
   cout << "finished destructing PhenoMap" << "\n";
#endif
}



void PhenoMap::SetSize(size_t rows, size_t cols) 
{
  size_t i;


  //  cerr << "Resizing a first dimension of a PHENOMAP of size"<<rows <<endl;

  pm.resize(rows);

  //  cerr << "Finished resizing the first dimension" <<endl;

  //  cerr << "Resizing a second dimension of a PHENOMAP of size"<<cols <<endl;
  for (i=0;i<rows;i++)
    {
      pm[i].resize(cols);
    }
}



/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
