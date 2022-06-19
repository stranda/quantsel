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

//copy constructor
DemoClass_space::DemoClass_space (const DemoClass_space& dc)
{
  //  cerr << "in DC_space copy constructor" <<endl;

  maxind = dc.maxind;
  indcnt = dc.indcnt;
  cl = dc.cl;
  I.insert(begin(dc.I),end(dc.I));
  for (size_t i=0; i<dc.UNUSED.size();i++) UNUSED[i]=dc.UNUSED[i];
  nextind = begin(I);
  //  cerr << "finished DC_space copy constructor, length: " <<I.size()<<endl;
}

DemoClass_space::~DemoClass_space ()
{
  I.clear();
  UNUSED.clear();
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
      indx = unirange(maxind);
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
  while (!RemoveInd(unirange(maxind),t,Atbls)) // keep trying until an individual is erased
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
