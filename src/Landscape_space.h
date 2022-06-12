/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim
*/

#ifndef LANDSCAPE__SPACE_H
#define LANDSCAPE_SPACE_H



/*includes
*/

#include <metasim.h>
#include <utilities.h>
#include <PackedIndividual_space.h>
#include <Democlass_space.h>
#include <TransMat.h>
#include <RandFuncs.h>
#include <PhenoMap.h>
#include <string>
#include <vector>
#include <list>
#include <iostream> 
#include <fstream> 


using namespace std;




///Local Matrix information Class
/**

   This class implements a set of matrices that model demography within populations.
   There are three matrices all of type TransMat.

   1) a survival matrix
   2) a reproduction matrix
   3) a male gamete weight matrix (usually if a category can contribute gametes, it is 
      given a 1, but not necessarily) 

 */
class LocalMat {

  ///size of each of the matrices
  size_t sz;
  ///list of survival matrices local to a particular habitat.  
  TransMat Slocal;
  ///list of repro matrices local to a particular habitat
  TransMat Rlocal;
  ///list of male function matrices local to a particular habitat
  TransMat Mlocal;

public:

  LocalMat () {}
  ~LocalMat() {}
  
inline  double GetSlocalVal(size_t from, size_t to)
    {
      assert(from<Slocal.Size());
      assert(to<Slocal.Size());
      Slocal.SetFromState(from);
      Slocal.SetToState(to);
      return Slocal.Value();
    }

inline  double GetRlocalVal(size_t from, size_t to)
    {
      assert(from<Rlocal.Size());
      assert(to<Rlocal.Size());
      Rlocal.SetFromState(from);
      Rlocal.SetToState(to);
      return Rlocal.Value();
    }

inline  double GetMlocalVal(size_t from, size_t to)
    {
      assert(from<Mlocal.Size());
      assert(to<Mlocal.Size());
      Mlocal.SetFromState(from);
      Mlocal.SetToState(to);
      return Mlocal.Value();
    }



inline  void SetSlocalVal(size_t from, size_t to, double val)
    {
      assert(from<Slocal.Size());
      assert(to<Slocal.Size());
      assert(val>=0);
      Slocal.SetElement(from,to,val);
    }

inline  void SetRlocalVal(size_t from, size_t to, double val)
    {
      assert(from<Rlocal.Size());
      assert(to<Rlocal.Size());
      assert(val>=0);
      Rlocal.SetElement(from,to,val);
    }

inline  void SetMlocalVal(size_t from, size_t to, double val)
    {
      assert(from<Mlocal.Size());
      assert(to<Mlocal.Size());
      assert(val>=0);
      Mlocal.SetElement(from,to,val);
    }


  void  SetSize(size_t sz);

  friend ostream &operator<<(ostream & stream, LocalMat &lm);

  friend istream &operator>>(istream & stream, LocalMat &lm);



}; //end localmat



///Landscape_space Class
  /**
This is the declaration of the Landscape_space class.  The landscape class
implements the simulation, basically.

o It is built upon the concept of a 2d array of suitable sites.

o The landscape moves through time one year (generation) at a time.  

o every x generations a new epoch can take over, where the migration
  matrices, demography, everything can change.

o the migration matrix provides a value, m_ij which is the
  probablility that an individual will migrate from population i to
  population j for each pair of populations

*/


class Landscape_space {
protected:
  /// Individuals data structure.
  /**

     This data structure (I) holds all of the individuals in the
     landscape.

     Right now, it's implemented as a vector of lists.  Each of the
     lists refers to the individuals within a particular habitat.
     Each generation individuals will be added and deleted from the
     lists in each habitat, but the number of habitats will not
     change.

     These habitats are determined solely by the state of the
     individual.  So, a landscape with 3 habitats (p) and 3
     demographic stages (s) is a 9x9 matrix.  All individuals with
     states in columns 1-3 are in habitat 1, 4-6 in habitat 2, and 7-9
     in habitat 3.  Individuals don't know how to convert state into
     habitat, so this has to be a Landscape-level function.

  */
  /*
    
    Variable declarations

   */

  ///Title of landscape
  string title;

  /// Individuals
  vector< DemoClass_space > I;

  ///the following vector is a list of males that can function as fathers
  ///during the execution of Landscape_space::Reproduce()
  ///it saves recalculating the vector over and over for 
  ///each mother
  //  vector < PackedIndividual_space > valid_males;  //no longer a member recalc in Reproduce()

  ///similarly, occupied_subpops indicates which subpops have
  ///males in them and is reset each Reproduce() call
  //  vector < int > occupied_subpops;


  ///Lookup table of alleles at all loci
  AlleleLookTbl Atbls;

  /// number of potentially suitable habitats in the landscape.
  int nhab;

  /// number of demographic stages
  int s   ;

  /// number of genetic loci
  int nloc;
  ///number of phenotypes to express
  int nphen;

  ///number of columns in the landscape (assuming grid populations)
  int cols;

  ///number of rows in the landscape (assuming grid populations)
  int rows;
  
  ///aspect ratio (number times y mean relative to x mean)
  double asp;

  ///minimum density to accept as possible for pollination events
  double mindens;

  ///mean of seed dispersal curve
  double seed_nmu;
  double seed_nmu2;

  double seed_mu;
  double seed_mu2;

  ///shape of seed dispersal curve
  double seed_nshape;
  double seed_nshape2;

  double seed_shape;
  double seed_shape2;

  ///seed mixing parameter (proportion of distribution given by seed_mu and seed_shape.  
  ///1-seed_mix is the proportion of distribution with seed_mu2 and seed_shape2 as parameters)
  double seed_mix;
  double seed_nmix;

  ///mean of pollen dispersal curve
  double pollen_mu, pollen_mu2;
  ///shape of pollen dispersal curve
  double pollen_shape, pollen_shape2;
  ///mixing rate when appropriate
  double pollen_mix;
  /// selfing rate
  double self;

  /// the number of different types of migration and other demographic conditions (number of epochs)
  int nep ;

  ///number of different within-habitat demographies to choose from.
  ///This variable determines the length of the vector LM defined in
  ///this class.  It also defines the length of each element in the
  ///outer vector of the variable "demoProbVec"

  int ndemo;

  /// a flag to tell whether the demographies are chosen at random from demoProbVec
  int rdemo;

  /// a counter of the current epoch
  int e   ;

  /// when true, choose epochs at random at every advance of the clock.  Ignores epoch vector.

  int randepoch;

  /// total number of generations
  int ngen;

  /// the current generation number 
  int t   ;

  ///the percentage of deviation allowed when imposing carrying capacity
  double habdelta;

  ///The maximum number of individuals allowed in the landscape
  long maxlandsz;

  ///is multiple paternity allowed? 0=single father, 1=every child has a randomly selected father 
  int multiple_paternity ;

  /**
    the next several lines define vectors of "vital matrices" each of
    these vectors is nep long.  Therefore each"epoch" can have a
    drastically different migration model, number of sites allowed to
    be occupied, demography, extinction, etc..
   */

  /// a vector of epoch begin dates
  std::vector<int> epochs;
  /// a vector of probabilities of seeing the conditions of a particular epoch in any given year
  std::vector<double> epochprobs;

  ///Survival matrix
  std::vector<TransMat> S;
  ///Reproduction matrix 
  std::vector<TransMat> R;
  ///Probability of an individual contributing a gamete to another individual in the same or different class
  std::vector<TransMat> M;

  ///A vector of Local Matrix Types.  These represent the demography within populations and can be used to
  ///set the diagonal elements of S, R, and M.
  ///The length of this vector is the same as the number of different within-pop demographies possible

  std::vector<LocalMat> LM;

  ///This structure is a vector of length nep.  Each element is is a
  ///vector of probabilities of observing a particular local
  ///demography of type LocalMat. The length of each of these
  ///sub-vectors is the number of different possible local demographies

  std::vector<std::vector <double> > demoProbVec ;

  ///"non-demographic" annual extinction probabilities for each site
  std::vector<std::vector<double> > evec;
  /// carrying capacities for each site
  std::vector<std::vector<int> > kvec;

  /// a vector of population left x coordinates
  std::vector<std::vector<double> > popleftx;
  /// a vector of population right x coordinates
  std::vector<std::vector<double> > poprightx;
  /// a vector of population left y coordinates
  std::vector<std::vector<double> > poptopy;
  /// a vector of population right y coordinates
  std::vector<std::vector<double> > popboty;

  ///The next two vectors hold the center coordinates for the subpop grid
  ///these are calculated based on the extent of the full landscape 
  ///(max poprightx - min popleftx)/SubPopSize (similar for y)
  ///so only SubPopSize needs to be defined in a landscape.
  /// a vector of subpopulations.  The inner vector has length 2 and corresponds to the x and y coords
  ///std::vector<std::vector<std::vector<double> > > subpops;

  ///the dimension of subpops
  //double SubPopSize;



  ///matrix to describe pollen kernels and parameters (varys through time)
  std::vector<std::vector<std::vector<double> > > PK; 
  ///matrix to describe pollen kernels and parameters (varys through time)
  std::vector<std::vector<std::vector<double> > > SK;

  ///expression matrix for phenotypes.  Right now, it is constant through time
  std::vector<std::vector<double> > expmat;
  ///heritability vector for phenotypes.  Right now, it is constant through time
  std::vector<double>  hsq;
  ///ainidices vector for phenotypes.  Right now, it is constant through time
  std::vector<int>  addstates;

  ///map of which phenotypes go with which demographic traits
  std::vector<std::vector <double> > gpdemo;


  ///per-habitat plasticity 
  PhenoMap plasticity;

  ///per-habitat selection 
  std::vector< PhenoMap > phenohab;
  
  /**
     Table of allele frequencies rows are alleles and cols are habitats

   */  

public:

  ///Constructor
  Landscape_space (int h=1, int stg=2, int loc=1, int ep=1, int nd=1, int gn=2);
  ///Destructor
  virtual ~Landscape_space();

  /*
    set the initial parameters
  */

  inline void setseed_mu(double m=0, double m2=0, double nm=0, double nm2=0) {seed_mu=m;seed_mu2=m2;seed_nmu=nm;seed_nmu2=nm2;}
  inline void setseed_shape(double s=0, double s2=0,double ns=0, double ns2=0) {seed_shape=s;seed_shape2=s2;seed_nshape=ns;seed_nshape2=ns2;}
  inline void setseed_mix(double mx=1, double nmx=1) {seed_mix=mx;seed_nmix=nmx;}

  inline void setpollen_mu(double m=0, double m2=0) {pollen_mu=m; pollen_mu2=m2;}
  inline void setpollen_shape(double s=0, double s2=0) {pollen_shape=s;pollen_shape2=s2;}
  inline void setpollen_mix(double mx=1) {pollen_mix=mx;}

  inline void setnphen(int p=1){nphen=p;}
  inline void setrows(int r=0){rows=r;}
  inline void setcols(int c=0){cols=c;}

  inline void setaspect(double as=0){asp=as;}
  inline void setmindens(double md=0){mindens=md;}
  inline void setself(double slf=0) {self=slf;}
  inline void setmultp(int mp=1) {multiple_paternity=mp;}
  inline void setranddemo(int rd=1) {rdemo=rd;}
  inline void setgens(int gn=2) {ngen=gn;}
  inline void setCgen(int cg) {t=cg;}
  inline void setCepoch(int ce) {e=ce;}

  inline void setMaxLandSize(int mls=300000) {maxlandsz=mls;}
  inline void assignRandEpoch(int re=1) {randepoch=re;}

  inline void setRandEpoch() {randepoch=1;}
  inline void unsetRandEpoch() {randepoch=0;}
  inline void reserveclasses() {I.resize(s*nhab);}

         void setepochs(int ep=1);
         void setndemo(int nd=1);
         void sethabs(int h=1);
         void setexpression();
         void setgpmap();
         void setupplasticity();
         void setupphenohab(); 
         void setstages(int stg=2);
  inline void setloci() {nloc=Atbls.size();}
         void setepochprob(int ce, double prob);
         void setepochstart(int ce, int strt);

  //  inline void setSubPop(double sps=10) {SubPopSize = sps;}

  /*
    get values of the parameters
  */


  inline int gethabs() {return nhab;}
  inline int getstages() {return s;}

/// find the population based on the demographic stage
///
  inline int habfromstage(int stg) {return stg/getstages();}
  
  inline double getseed_mu() {return seed_mu;}
  inline double getseed_mu2() {return seed_mu2;}
  inline double getseed_nmu() {return seed_nmu;}
  inline double getseed_nmu2() {return seed_nmu2;}
  inline double getseed_shape() {return seed_shape;}
  inline double getseed_shape2() {return seed_shape2;}
  inline double getseed_nshape() {return seed_nshape;}
  inline double getseed_nshape2() {return seed_nshape2;}
  inline double getseed_mix() {return seed_mix;}
  inline double getseed_nmix() {return seed_nmix;}
  inline double getpollen_mu() {return pollen_mu;}
  inline double getpollen_shape() {return pollen_shape;}

  inline double getpollen_mu2() {return pollen_mu2;}
  inline double getpollen_shape2() {return pollen_shape2;}
  inline double getpollen_mix() {return pollen_mix;}

  inline double getaspect(){return asp;}
  inline double getmindens(){return mindens;}
  inline double getself() {return self;}
  inline int getmultp() {return multiple_paternity;}
  inline int getloci() {setloci(); return nloc;}
  inline int getepochs() {return nep;}
  inline int getCgen() {return t;}
  inline int getCepoch() {return e;}
  inline double getepochprob(int ce) {return epochprobs[ce];}
  inline int getepochstart(int ce) {return epochs[ce];}
  inline int getgens() {return ngen;}
  inline int getrandepoch() {return randepoch;}
  inline int getranddemo() {return rdemo;}
  inline int getndemo() {return ndemo;}
  inline int getMaxLandSize() {return maxlandsz;}

  inline int getnphen(){return nphen;}
  inline int getrows(){return rows;}
  inline int getcols(){return cols;}

  //  inline double getSubPop() {return SubPopSize;}

  void init(int h=1, int stg=2, int loc=0, int ep=1, int nd=1, int gn=1);

  void setupexpression();
  /**
    set the appropriate matrices
   */
  void setS(TransMat a, int ep=-1);
  void setR(TransMat a, int ep=-1);
  void setM(TransMat a, int ep=-1);

  inline void setSmatElement(int ep, int t, int f, double val)  {S[e].SetElement(f,t,val);}
  inline void setRmatElement(int ep, size_t t, size_t f, double val)  {R[e].SetElement(f,t,val);}
  inline void setMmatElement(int ep, size_t t, size_t f, double val)  {M[e].SetElement(f,t,val);}

  inline double getSmatElement(int ep, size_t t, size_t f)  {return S[e].GetElement(f,t);}
  inline double getRmatElement(int ep, size_t t, size_t f)  {return R[e].GetElement(f,t);}
  inline double getMmatElement(int ep, size_t t, size_t f)  {return M[e].GetElement(f,t);}


  inline void setLSmatElement(int d, size_t t, size_t f, double val) {LM[d].SetSlocalVal(f,t,val);}
  inline void setLRmatElement(int d, size_t t, size_t f, double val) {LM[d].SetRlocalVal(f,t,val);}
  inline void setLMmatElement(int d, size_t t, size_t f, double val) {LM[d].SetMlocalVal(f,t,val);}

  inline double getLSmatElement(int d, size_t t, size_t f)  {return LM[d].GetSlocalVal(f,t);}
  inline double getLRmatElement(int d, size_t t, size_t f)  {return LM[d].GetRlocalVal(f,t);}
  inline double getLMmatElement(int d, size_t t, size_t f)  {return LM[d].GetMlocalVal(f,t);}

  inline void setPollenkern(int e, size_t r, size_t c, double val) {PK[e][r][c]=val;}
  inline void setSeedkern(int e, size_t r, size_t c, double val) {SK[e][r][c]=val;}
  inline double getPollenkern(int e, int r, int c) {return PK[e][r][c];}
  inline double getSeedkern(int e, int r, int c) {return SK[e][r][c];}

  inline void setexpmatel(int l,int p, double val) {expmat[l][p]=val;}
  inline double getexpmatel(int l, int p){return expmat[l][p];}

  inline void setheritability(int p, double val){hsq[p]=val;}
  inline double getheritability(int p){return hsq[p];}

  inline void setaddstate(int l, int index){addstates[l]=index;}
  inline int getaddstate(int l){return addstates[l];}

  inline double getgpdemo(int resp, int phen) {return gpdemo[resp][phen];}

  inline void setgpdemo(int index, int coef, double val) {gpdemo[index][coef]=val;}

  double betamax(double a,double b);
  double getAdjDemo(int response, PackedIndividual_space Ind);
  double getAdjDisp(int response, PackedIndividual_space Ind); 
  double getAdjDemoDens(PackedIndividual_space Ind);
  inline double getplasticity(int hab, int phen){return plasticity.GetElement(hab,phen);}
  inline void   setplasticity(int hab, int phen, double val){plasticity.SetElement(hab,phen,val);}

  inline double getphenohab(int hab, int response, int param){return phenohab[response].GetElement(hab,param);}
  inline void   setphenohab(int hab, int response, int param, double v){phenohab[response].SetElement(hab,param,v);}

  
inline void ConstructDemoMatrix()
  {
    if (ndemo)
      {
	if (rdemo)
	  {
	    RandomlyConstructDemoMatrix();
	  }
	else
	  {
	    SequentiallyConstructDemoMatrix();
	  }
      }
  }
  /**

     This function selects sub-matrices from the vector of local
  matrices.  The sub-matrices are chosen in order, and are reused
  until all habitats have a matrix assigned to them.  The method takes
  these matrices and inserts them on the diagonal of the S, R, and M
  matrices.  This function DOES NOT ALTER THE OFF-DIAGONAL
  SUBMATRICES.  These are defined in the current epochs S,R, and M
  matrices.

     The current epoch determines which  S,R,M, and demoProbVec to choose.

   */
void SequentiallyConstructDemoMatrix();

  /**

     This function randomly selects sub-matrices from the vector of
     local matrices.  The sub-matrices are chosen from a multinomial
     distribution given by "demoProbVec".  The method takes these
     matrices and inserts them on the diagonal of the S, R, and M
     matrices.  This function DOES NOT ALTER THE OFF-DIAGONAL
     SUBMATRICES.  These are defined in the current epochs S,R, and M matrices.

     The current epoch determines which  S,R,M, and demoProbVec to choose.

   */
void RandomlyConstructDemoMatrix();


  /*
    set the population characteristic vectors
   */
  ///Set up the extinction vector(s).  'ev' is a vector of type double that correpond
  ///to extinction rates for each of the 0..nhab habitats. 'ep' is the epoch
void setextinct(int ep, double *ev);
  ///Get an extinction vector(s).  for each of the 0..nhab habitats. 'ep' is the epoch
void getextinct(int ep, double *ev);
  ///Set up the carrying capacity vector(s) 'cv' is a vector of type
  //int that correspond to extinction rates for each of the 0..nhab
  ///habitats. 'ep' is the epoch
void setk(int ep, int *cv);
  ///Get carry capacity vector(s).  0..nhab habitats. 'ep' is the epoch
void getk(int ep, int *cv);
  ///set the coordinates of the rectangular shaped populations
void setpoploc(int ep, double *lx, double *rx,double *topy, double *boty);
void getpoploc(int ep, double *lx, double *rx,double *topy, double *boty);

  std::vector<int> getSurroundingPops(int stg, double radius);
  
  ///takes an x and y coordinate and returns the population that contains that coordinate
  ///if found outside of any population, returns -1
  int  getpopulation(double x, double y);

  ///takes an x and y coordinate and returns the subpopulation that contains that coordinate
  ///if found outside of any population, returns -1
  ///int  getsubpopulation(double x, double y);

  ///Set up the local demography vector for a particular epoch
void setldemovector(int ep, double *dv);
  ///Get the local demography vector for a particular epoch
void getldemovector(int ep, double *dv);

void zeroextinct();
void zerok();

  ///This resets the pointer to members of the individual class 'cl' to the start of the list
  inline void resetStage(int cl)
  {
    I[cl].ResetIndividuals(); 
  }
  ///This function advances the pointer to inds in the demographic class 'cl'.
  ///it returns 1 if there is another ind to grab, otherwise it returns 0

  inline int advanceStagePtr(int cl)
  {
    return I[cl].NextIndividual();
  }
  ///This function tells the size of the demographic class pointed to by 'cl'
  inline int StageSize(int cl=0)
  {
    return I[cl].size();
  }

  inline PackedIndividual_space getNextInd(int cl)
  {
    return I[cl].GetCurrentIndividual();
  }

  inline void SetUpInd(PackedIndividual_space &ind)
  {
    ind.resetLoci(Atbls);
  }
  
  inline int addIndividual(PackedIndividual_space ind, int t)
  {
    ind.Birth(t,Atbls);
    return (I[ind.GetClass()].AddIndividual(ind));
  }
  ///returns a vector of phenotypes (additive model filtered through exp matrix)
  std::vector< double > IndividualPhenotype(PackedIndividual_space ind);
  /**
     returns all of the individual phenotypes
  */
  std::vector< double > Phenotypes();


  ///choose an epoch using some selection criteron
  void ChooseEpoch();
  ///randomly chooses an epoch, random epoch selection is in action
  void RandomlyChooseEpoch();
  ///Initialize the Populations
  void popsizeset(std::vector<int> & ps);

  /*
    functions that allow outside processes to interact with the allele lookup tbl
   */

  inline void Atbl_push_back(AlleleTbl * atp)
  {
    Atbls.push_back(atp);
  }

  inline int LocusGetClassType(int l)
  {
    return Atbls[l]->getClassType();
  }

  inline int LocusGetPloidy(int l)
  {
    return Atbls[l]->getPloidy();
  }

  inline int LocusGetTrans(int l)
  {
    return Atbls[l]->getTrans();
  }

  inline double LocusGetMutRate(int l)
  {
    return Atbls[l]->getMutationRate();
  }

  inline vector<int> LocusGetAindices(int l)
  {
    return Atbls[l]->getAindices();
  }

  inline void LocusGetAlleleRef(int l, int andx, Allele* ptr)
  {
    if (Atbls[l]->getClassType()==SEQALLELETBL)
      {
	Atbls[l]->getAlleleRef(andx,(dynamic_cast<SeqAllele *>(ptr)));
      }
    else if (Atbls[l]->getClassType()==INFALLELETBL)
      {
	Atbls[l]->getAlleleRef(andx,ptr);
      }
    else if (Atbls[l]->getClassType()==STEPALLELETBL)
      {
	Atbls[l]->getAlleleRef(andx,ptr);
      }
    else 
      {
	cerr <<"dont know what type of locus this is"<<endl;
	assert(1==0);
      }
  }

  /**

     This function takes a deomgraphic stage on a s*p x s*p matrix and returns the habitat that that stage belongs to


   */
int Habitat(int stage);


  /** 
      This function reports the population size of the habitat given
      by the integer argument.  If the argument is -1, then the total
      size of all habitats are reported
 */

int PopSize(int p=-1);




  /**

     Survive

     Survive might be better named because this is the method that
     implements survival, growth, and migration.
     
   */



void Survive();


/*
this function is supposed to take a vector of individuals and return a unique
set that has their subpop ids
 */
///  std::vector<int> uniquesubpops(vector < PackedIndividual_space > inds);

/*
This function takes a subpopulation id and a vector of individuals and returns the individuals that
are located in the subpopulation
 */
/// PackedIndividual_space  indsfromsubpop(int sp, vector <PackedIndividual_space> inds);

/** 

    This method sets the discrete lookup table in the global
    RandLibObj.  The table is set with probabilities that individual
    demographic classes will contribute a mate to a particular pi
    passed in.  

    ***Important: this function allocates a lookup table.  it must be
    ***freed at some point after the function is invoked by issuing the
    ***command: RandLibObj.FreeDiscreteLookup();

*/
vector < PackedIndividual_space > CalculateMaleGameteClassVector(PackedIndividual_space pi, vector< PackedIndividual_space > valid_males);


/** 

This function approximates the solution implemented in CalculateMaleGameteClassVector by 
contructing an ever expanding 'band' around each female at a distance determined by a random
pull from a pollen dispersal kernel.


    ***Important: this function allocates a lookup table.  it must be
    ***freed at some point after the function is invoked by issuing the
    ***command: RandLibObj.FreeDiscreteLookup();

*/
vector < PackedIndividual_space > CalculateMaleGameteClassVectorApproxDist(PackedIndividual_space pi, vector< PackedIndividual_space > valid_males);


  /**
     Reproduce:

     o traipses through all of the individuals 
     o Decides on the number of offspring produced by every individual.
     o Matches individuals to their mates (if necessary)
     
   */

void Reproduce();

  /**
     same as reproduce, but uses a different approach to movement of male gametes. This is _not_ the algorithm defined in the kernelPop paper in Molecular Ecology Res.
  **/
  //void Reproduce_approx();


  /** 

      extirpate will go through extinction vector pick random numbers
      from a binomial distribution and set the census size at those
      sites to zero

  */
void Extirpate(); 



  /**
     
     Adjust the size of a particular state down to a particular size
     
   */
void CarryState(size_t maxsz, int i);

  /**

     set habitat sizes down to carrying capacity

   */
void HabCarry(int k = -1);


  /**

     set habitat sizes down to carrying capacity by only removing individuals from the
     smallest category in each habitat

   */
void HabCarry_stg0(int k = -1);

  /**

     set overall landscape size down to "maxlandsz" carrying capacity

   */
void LandCarry();


  /**

     shuffles the individuals vector randomly.  Must be used between
     generations if the order in which individuals are chose is
     important.  Neglecting it could simulate inbreeding in the case of sexual reproduction


   */
void Randomize();
  

void Advance();

  /**
  This function takes the x,y coordinate, the demographic state, and returns new x and y coordinates for the propagule
  based on the internal seed kernel matrix
  */
  void new_propagule_xy(double ix, double iy, int cls, double aspect, double &x, double &y);

  double pollenKernelDensity(double dist, int i);

  //this function returns a distance based on the pollen kernel in place
  //this is used in the pollen approx functions and assumes that there 
  //is only one pollen kernel density (no difference among stages)
  double RandpollenKernelDensity();



  /**

     This function goes through and checks the number of copies of each allele present in the entire landscape.  If an allele is missing from the table, an error should be generated.  If there are no individuals with an allele, eliminate that allele from the table

   */
void GCAlleles();  

  ostream &WriteLoci(ostream &stream);

  friend ostream &operator<<(ostream & stream, Landscape_space &l);

  friend istream &operator>>(istream & stream, Landscape_space &l);


}; // end Landscape_space


#endif /*LANDSCAPE*/



class Landscape_space_statistics: public Landscape_space {

public:
Landscape_space_statistics (int h=1, int stg=2, int loc=1, int ep=1, int gn=2);
~Landscape_space_statistics ();

  /** Statistics 

      Reports various staistics about the state of the landscape.
      Sends the report to standard out.  Would be nice to somehow
      determine the nature of the output produced by Statistics with a
      run-time choice...

*/
void Statistics(ostream &streamout = cout);

  //Demographic Stats
  /**

     Return generation length by averaging the ages of reproductive individuals weighted by the number of offspring they produce

  */
double GenLength();



  //Methods to output  data to stream

  /**

     Sends an arlequin 2.0 project file(s) to the stream provided.  Takes numind  individuals from each habitat.

   */
void ArlequinDiploidOut(int numind = 200, ostream &streamout = cout);
void ArlequinHaploidOut(int numind = 200, ostream &streamout = cout);


  /**

     Sends a genepop 3.1 input files to the stream provided.  Takes
     numind individuals from each habitat. Includes both haploid and
     diploid data

   */
void GenepopOut(int numind = 200, ostream &streamout = cout);

  /**

     Sends a MicroRat input file to the stream provided.  Takes
     numind individuals from each habitat. Includes both haploid and
     diploid data

   */

void MicroRatOut(int numind = 200, ostream &streamout = cout);

  /**

     Sends a biosys input file to the stream provided.  Takes numind  individuals from each habitat.

   */
void BiosysDiploidOut(int numind = 200, ostream &streamout = cout);

  /**
     creates input files for Peter Beerli's "migrate" program
  */
void MigrateDiploidOut(int numind = 200, ostream &streamout = cout);

  /**
     creates input files for GDA by Paul Lewis
  */
void GdaOut(int numind = 200, ostream &streamout = cout);

  /**

     Sends a R input file to the stream provided.  Takes numind  individuals from each habitat.

   */
void ROut(int numind = 200, ostream &streamout = cout);

  /**
     returns the landscape as a vector of int to be used in some r popgen analyses.  Called by rmetasim
  */
  vector <int> Rmat(int numind=0);

};



/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; minor-mode: font-lock ***
;;; End: ***
*/
