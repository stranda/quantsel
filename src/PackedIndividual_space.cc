/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

*/


/*includes
*/

#include <PackedIndividual_space.h>
#include <RandLib.h>

/**
Class PackedIndividual_space methods
 */

PackedIndividual_space::PackedIndividual_space (int c, int g, int nl)
{
  cl=c;
  sex=0;
  gen=g;
  assert(nl <= MAXLOCI);
  ///  SetLoci();
  Change(g-1);
  SetLastRep(0);
  SetNumOff(0);
}

PackedIndividual_space::~PackedIndividual_space ()
{
}

/// from a uniform distribution
int PackedIndividual_space::RandomizeClass(int numclass)
{
  return RandLibObj.unirange(numclass);
}

void  PackedIndividual_space::SetClass(int numclass)
{
  cl = numclass;
}

int  PackedIndividual_space::GetClass()
{
  int c;
  c = int(cl);
  return c;
}

void PackedIndividual_space::SetSex(int newsex)
{
  sex = newsex;
}

void PackedIndividual_space::SetGen(int newgen)
{
  gen = newgen;
}

int PackedIndividual_space::GetGen()
{
  return gen;
}


void PackedIndividual_space::SetLoci(AlleleLookTbl &Atbls)
{
  int i;
  for (i=0;i<MAXLOCI;i++)
    {
      if (i<int(Atbls.size()))
	{
	  PL[i]=Atbls[i]->getPloidy();
	}
      else
	{
	  PL[i]=0;
	}
    }
  nloc=int(Atbls.size());
}

void PackedIndividual_space::resetLoci(AlleleLookTbl &Atbls)
{
  int i,j;
  SetLoci(Atbls);
  for (i=0;i<MAXLOCI;i++)
    {
      for (j=0;j<MAXPLOIDY;j++)
	{
	  G[ ((i * MAXPLOIDY) + j) ]   = -1;
	}
    }
}


PackedIndividual_space PackedIndividual_space::MakeGamete(AlleleLookTbl &Atbls)
{
  int i;
  int lsize;
  //  int a=0,b=0,c=0;
  PackedIndividual_space pi;
  pi.resetLoci(Atbls);
  lsize = nloc;
  for (i=0;i<lsize;i++)
    {
      if (Atbls[i]->getTrans()==0)  //biparental inheritance
	{
	  //	  a++;
	  assert(Atbls[i]->getPloidy()==2);
	  pi.G[((i * MAXPLOIDY))] = G [((i * MAXPLOIDY) + RandLibObj.unirange(1))];
	}
      else if (Atbls[i]->getTrans()==1) //maternal inheritance
	{
	  //	  b++;
	  pi.G[((i * MAXPLOIDY))] = G [ i * MAXPLOIDY ];
	}
      else if (Atbls[i]->getTrans()==2) //Paternal inheritance. male parent
	{
	  //	  c++;
	  pi.G[((i * MAXPLOIDY))] = G [ i * MAXPLOIDY ]; 
	}
      else
	{
	  cerr << "Fell through all of the inheritance types in MakeGamete " << endl;
	  assert(1==0);
	}
      if (!(pi.G[((i * MAXPLOIDY))]>-1))
	{
	  //	  cerr << "a "<<a<< " b "<<b<<" c "<<c<<endl;
	    assert(pi.G[((i * MAXPLOIDY))]>-1);
	}
    }
  return pi;
}

int PackedIndividual_space::GetRandAlleleIndex(int l)
{
  int index;
  double ru;

  ru=RandLibObj.uniform();
  if (ru==1) {ru=0.999999999999;} //if ru=1 then index will equal ploidy below, instead of ploidy-1
  index = (int)floor(ru * (PL[l]));

  assert(G [((l * MAXPLOIDY) + index)]>=0);

  return G [((l * MAXPLOIDY) + index)];
}


int PackedIndividual_space::IsGenotypeSet()
{
  int j,s,i;
  s=1;
  for (i=0;i<nloc;i++)
    {
      for (j=0;j<PL[i];j++)
	{
	  if(G[ ((i * MAXPLOIDY) + j) ]<0)
	    {
	      s=0;
	    }
	}
    }

  return s;

}

void PackedIndividual_space::SetRandGenotype(AlleleLookTbl &Atbls)
{
  int j,i;


  for (i=0;i<nloc;i++)
    {
      for (j=0;j<PL[i];j++)
	{
	  G[ ((i * MAXPLOIDY) + j) ]   =  Atbls[i]->getRandAlleleIndex();
	}
    }
}





PackedIndividual_space  PackedIndividual_space::repro_sex(PackedIndividual_space & SO1, PackedIndividual_space & SO2, int t, AlleleLookTbl &Atbls)
{
  int i;
  int k, l;
  PackedIndividual_space pi(SO1), ti0, ti1;
  pi.resetLoci(Atbls);

  //  assert(SO2.IsGenotypeSet());

  ti0 = SO1.MakeGamete(Atbls);
  ti1 = SO2.MakeGamete(Atbls);

  for (i=0;i<nloc;i++)
    {
      k = ti0.GetAllele(i,0);
      assert(k>=0); 
      pi.G[((i * MAXPLOIDY) + 0)] = k;
      l = ti1.GetAllele(i,0);
      assert(l>=0); 
      pi.G[((i * MAXPLOIDY) + 1)] = l;

      //swap alleles so that diploid heterozygotes are sorted //if not then can tract haplotypes from parents
      //     if (PL[i]==2)
      //	{ 
      //	  if ((pi.G[((i * MAXPLOIDY) + 1)] >= 0) && (pi.G[((i * MAXPLOIDY) + 0)] > pi.G[((i * MAXPLOIDY) + 1)] ))
      //	    {
      //	      pi.swap_allele(i);
      //	    }
      //	}
    }
  return pi;
}

PackedIndividual_space  PackedIndividual_space::repro_asex(PackedIndividual_space & SO, int t)
{
  //  int i,j;

  //  int tmpi;

  //  int k;

  PackedIndividual_space pi(SO);
  //  pi.resetLoci(Atbls);

  return pi;
}

void PackedIndividual_space::Growth(AlleleLookTbl &Atbls)
{
  int j,i;
  for (i=0;i<nloc;i++)
    {
      for (j=0;j<PL[i];j++)
	{
	  Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + j) ]);
	}
    }
}

void PackedIndividual_space::Birth(int t, AlleleLookTbl &Atbls)
{
  int i;

  for (i=0;i<nloc;i++)
    {
      if (PL[i]==1)
	{
	  if (t>=0)
	    {
	      G[ ((i * MAXPLOIDY) + 0) ] = Atbls[i]->mutator(G[ ((i * MAXPLOIDY) + 0) ],t);
	    }
	  else
	    {
	      Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + 0) ]);
	    }
	}
      if (PL[i]==2)
	{
	  if (t>=0)
	    {
	      G[ ((i * MAXPLOIDY) + 0) ] = Atbls[i]->mutator(G[ ((i * MAXPLOIDY) + 0) ],t);
	      G[ ((i * MAXPLOIDY) + 1) ] = Atbls[i]->mutator(G[ ((i * MAXPLOIDY) + 1) ],t);
	    }
	  else
	    {
	      Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + 0) ]);
	      Atbls[i]->AddAlleleFreq(G[ ((i * MAXPLOIDY) + 1) ]);
	    }
	}
    }
}

void PackedIndividual_space::Death(int t, AlleleLookTbl &Atbls)
{
  size_t i, sz;
  sz = Atbls.size();
  for (i=0;i<sz;i++)
    {
      if (Atbls[i]->getPloidy()==1)
	{
	  Atbls[i]->KillAlleleCopy(G[ ((i * MAXPLOIDY) + 0) ],t);
	}
      if (Atbls[i]->getPloidy()==2)
	{
	  Atbls[i]->KillAlleleCopy(G[ ((i * MAXPLOIDY) + 0) ],t);
	  Atbls[i]->KillAlleleCopy(G[ ((i * MAXPLOIDY) + 1) ],t);
	}
    }
}


ostream & operator<<(ostream & stream, PackedIndividual_space &ind)
{
  int j,i,is;

  stream << ind.GetClass() << " " << ind.GetSex() << " " << ind.GetGen() << "  "<< ind.GetX() << "  " << ind.GetY()<<" " << ind.GetMX() << "  " << ind.GetMY()<<" " << ind.GetFX() << "  " << ind.GetFY()<<" ";

  for (j=0;j<ind.nloc;j++)
    {
      for (i=0;i<ind.PL[j]; i++)
	{
	  is=ind.GetAllele(j,i);
	  stream << is << " " ;
	}
      stream << "   " ;
    }
  stream << endl;
  return stream;
}

istream & operator>>(istream & stream, PackedIndividual_space &ind)
{
  int i, is, j;

  stream >> ind.cl >> ind.sex >> ind.gen >> ind.x >> ind.y >> ind.mx >> ind.my >> ind.fx >> ind.fy;

  for (j=0;j<ind.nloc;j++)
    {
      for (i=0;i<ind.PL[j]; i++)
	{
	  stream >> is ;
	  ind.SetAllele(j,i,is);
	}
    }

  return stream;
}

int PackedIndividual_space::operator==(PackedIndividual_space ti)
{
  short f=0;
  short i,j;
  if ((cl==ti.GetClass())&&(gen==ti.GetGen())&&(x==ti.GetX())&&(y==ti.GetY()))
  {
    f=1;
    for (i=0;i<MAXLOCI;i++)
      {
	for (j=0;j<MAXPLOIDY;j++)
	  {
	    if (G[(i*MAXPLOIDY) + j]!=ti.GetAllele(i,j))
	      {
		f=0;
		break;
	      }
	  }
      }
  }
  return (f==1);
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
