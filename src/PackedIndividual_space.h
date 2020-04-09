/*

$Modified: astrand$

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

This is the declaration of the site object type.

This is the object that defines an individual.

there is a base object: Individual



*/

#ifndef  PACKED_INDIVIDUAL_SPACE_H
#define PACKED_INDIVIDUAL_SPACE_H

/*
includes
*/

#include <metasim.h>
#include <TransMat.h>
#include <FastAllele.h>
//#include <FastSeqAllele.h>
#include <AlleleTbl.h>

using namespace std;

/**

The indivdual class is essentially abstract.  It mainly contains a
demographic class variable to pass to its descendents

 */
class PackedIndividual_space {

protected:
  /// xy coordinates of an individual in a population
  double x,y;  
  /// xy coordinates of mother of this individual
  double mx,my; 
  /// xy coordinates of father of this individual
  double fx,fy;
  /// time click in which the individual is born
  int gen;
  /// the time click the ind was last modified by survival.
  int changed;
  /// the time the individual last reproduced
  int lastrep;
  /// the number of offspring produced last gen
  int noff;
  ///demographic age or stage of the individual
  short cl;
  ///sex
  short sex;
  ///The number of loci
  int nloc;
  ///array of the ploidy for each locus
  short PL[MAXLOCI];
  ///linearized matrix of diploid genotypes of the ind.
  short G[MAXLOCI * MAXPLOIDY];

public:
  PackedIndividual_space(int c=0, int g=0, int nl=0);
  ~PackedIndividual_space();
  ///class stands for demographic class (across the entire landscape of habitats)
  ///randomly sets the class of an individual (for certain types of initialization)
  int RandomizeClass(int numclass=1);

  void SetClass(int newclass=0);
  int GetClass();

  void SetSex(int newsex=0);
  inline int GetSex() {  return sex;}

  void SetGen(int newgen=0);
  int GetGen();

  void SetLoci(AlleleLookTbl &Atbls);

  void resetLoci(AlleleLookTbl &Atbls);

  inline  double GetX()
    {
      return x;
    }

  inline void SetX(double nf)
    {
      x = nf;
    }

  inline  double GetY()
    {
      return y;
    }

  inline void SetY(double nf)
    {
      y = nf;
    }

  inline  double GetMX()
    {
      return mx;
    }

  inline void SetMX(double nf)
    {
      mx = nf;
    }

  inline  double GetMY()
    {
      return my;
    }

  inline void SetMY(double nf)
    {
      my = nf;
    }

  inline  double GetFX()
    {
      return fx;
    }

  inline void SetFX(double nf)
    {
      fx = nf;
    }

  inline  double GetFY()
    {
      return fy;
    }

  inline void SetFY(double nf)
    {
      fy = nf;
    }


  inline  int GetLastRep()
    {
      return lastrep;
    }

  inline void SetLastRep(int lr)
    {
      lastrep = lr;
    }
  inline  int GetNumOff()
    {
      return noff;
    }

  inline void SetNumOff(int nf)
    {
      noff = nf;
    }

  inline  int GetAllele(int l, int a)
    {
      return G[ ((l * MAXPLOIDY) + a) ];
    }

  inline void SetAllele(int l, int a, int al)
    {
      assert (al<MAXALLELES);
      G[ ((l * MAXPLOIDY) + a) ] = al;
    }

  inline int SumAlleles(int loc)
  {
    return (GetAllele(loc,0)+GetAllele(loc,1));
  }

  inline void WriteState(int l, int a, AlleleLookTbl &Atbls, ostream & stream = cout)
    {
      return Atbls[l]->WriteAlleleState(G [ ((l * MAXPLOIDY) + a) ], stream );
    }
  /*
  inline int GetIntState(int l, int a, AlleleLookTbl &Atbls)
  {
    return (Atbls[l]->getAllele(G [ ((l * MAXPLOIDY) + a) ])).GetState();
  }
  */

  inline void swap_allele(int l)
    {
      int tmpi;
      tmpi = G[((l * MAXPLOIDY) + 0)];
      G[((l * MAXPLOIDY) + 0)] = G[((l * MAXPLOIDY) + 1)];
      G[((l * MAXPLOIDY) + 1)] = tmpi;
    }
  ///sets the survival change flag to true
  inline void Change(int t)
    {
      changed = t;
    }
  inline int GetChanged()
    {
      return changed;
    }

  PackedIndividual_space MakeGamete(AlleleLookTbl &Atbls);


  int GetRandAlleleIndex(int l);
  void SetRandGenotype(AlleleLookTbl &Atbls);
  int IsGenotypeSet();


  /**
     Updates the allele tables.  If t>0 then also runs the mutation algorithms on each allele.
   */
  void Birth(int t,AlleleLookTbl &Atbls);
  /**
     Used for transferring individuals among classes
   */
  void Growth(AlleleLookTbl &Atbls);
  /**
     Updates the allele tables as an individual is removed from the population
   */
  void Death(int t, AlleleLookTbl &Atbls);

  PackedIndividual_space repro_sex(PackedIndividual_space & SO1, PackedIndividual_space & SO2, int t, AlleleLookTbl &Atbls);
  PackedIndividual_space repro_asex(PackedIndividual_space & SO, int t=0);


  friend ostream & operator<<(ostream & stream, PackedIndividual_space &ind);
  friend istream & operator>>(istream & stream, PackedIndividual_space &ind);

  size_t Sizeof();

  int operator==(PackedIndividual_space ti);

}; // end PackedIndividual_space


#endif /*PACKED_INDIVIDUAL*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
