/**

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

The implementation of the demo class object

*/


/*includes
*/

#include <Democlass_space.h>

DemoClass_space::DemoClass_space ()
{
  maxind = 0;
  indcnt = 0;
}
DemoClass_space::~DemoClass_space ()
{
}

///add an individual to the data structure.  Returns the index to the individual.
///this method also adds the alleles to the allele table for each locus.  This maintains the allele
///frequency tables.
int DemoClass_space::AddIndividual (PackedIndividual_space & PI)
{
  int i;
  i = -1;
  if (UNUSED.empty())
    {
      I[maxind]=PI;
      i=maxind;
      maxind++;
    } //end     "if UNUSED is empty"
  else
    {
      i = UNUSED.back();
      UNUSED.pop_back();
      I[i]=PI;
    }
  return i;
}

void DemoClass_space::ClearClass (int t,AlleleLookTbl &Atbls)
{
  map<int,PackedIndividual_space,less <int> >::iterator iiter;
  iiter=I.begin();
  while (iiter!=I.end())
    {
      (*iiter).second.Death(t,Atbls);
      iiter++;
    }
  I.clear();
  UNUSED.clear();
  maxind=0;
}

int DemoClass_space::GetRandomIndex ()
{
  int f=0;
  int indx;

  while (f==0) // keep trying until an individual is found
    {
      indx = RandLibObj.unirange(maxind);
      if (I.find(indx)!=I.end())
	{
	  f=1;
	}
      else
	{
	  f=0;
	}
    }
  return indx;
}  

void DemoClass_space::RemoveRandomInd (int t,AlleleLookTbl &Atbls)
{
  while (!RemoveInd(RandLibObj.unirange(maxind),t,Atbls)) // keep trying until an individual is erased
    {
    }
}  

void DemoClass_space::CompressClass (double frac)
{
  vector <PackedIndividual_space> tvec;
  size_t i,sz;
  if (I.size()>0&&(I.size()<maxind*frac))
    {
      tvec.reserve(I.size());
      ResetIndividuals();
      i=0;
      for (nextind=I.begin();nextind!=I.end();nextind++)
	{
	  tvec.push_back((*nextind).second);
	}
      I.clear();
      UNUSED.clear();
      maxind=0;
      sz=tvec.size();
      for (i=0;i<sz;i++)
	{
	  I[i]=tvec[i];
	  maxind++;
	}
      ResetIndividuals();
    }
}

double DemoClass_space::GenLength (int t)
{
  double genoff, totoff;
  //  int indx;

  int lr, no;

  if (I.size()>0)
    {
      genoff=0.0;
      totoff=0.0;
      ResetIndividuals();
      do
	{
	  //	  indx = GetCurrentIndex();
	  lr = GetCurrentLastRep();
	  no =  GetCurrentNumOff();
	  genoff =+ ((t - lr) * no);
	  totoff =+ no;
	}
      while (!NextIndividual());
      if (totoff==0) 
	{
	  return 0;
	}
      else
	{
	  return double(genoff/totoff);
	}
    }
else
  {
    return 0.0;
  }
}

vector < PackedIndividual_space > DemoClass_space::ReturnAsVector ()
{
  size_t i;
  vector < PackedIndividual_space > dvec;
  map<int, PackedIndividual_space, less<int> >::iterator tmpind = nextind;

  ResetIndividuals();
  for (i=0;i<I.size();i++)
    {
      dvec.push_back((*nextind).second) ;
      NextIndividual();
    }
  nextind = tmpind;
  return dvec;
}


ostream &operator<<(ostream & stream, DemoClass_space & DC)
{
  size_t i;
  PackedIndividual_space tmpI;
  DC.ResetIndividuals();
  for (i=0;i<DC.I.size();i++)
    {
      stream << (*DC.nextind).second ;
      DC.NextIndividual() ;
    }
  return stream;
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
