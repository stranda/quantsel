/*

$Modified: astrand $

Copyright (C) 1999 Allan E. Strand

This file is part of Metasim

*/

#ifndef PHENOMAP_H
#define PHENOMAP_H
/*
includes
*/

#include <metasim.h>
#include <RandLib.h>
#include <iostream>

using namespace std;
/// Generic matrix to hold information for each locality

class PhenoMap {
private:

  std::vector< std::vector<double> > pm;

public:
  PhenoMap ( size_t rows=1, size_t cols=1 ) ;

  ~PhenoMap () ;

  /// Sets a matrix cell value (returns a zero if successful)
  inline int SetElement(size_t row, size_t col, double val) 
    { 
      return (pm[row][col] = val) > 0; 
    }
  
  
  ///returns value of element at r,c  
  inline double GetElement(size_t row, size_t col) 
    {
      return pm[row][col]; 
    } 

  ///Returns the extent of the matrix (assumes a square matrix)
  inline size_t nrows () 
    { 
      return pm.size(); 
    }
  ///Returns the extent of the matrix (assumes a square matrix)
  inline size_t ncols () 
    { 
      return pm[0].size(); 
    }
    ///Sets the size of the matrix: Resizes the matrix and sets "size"
  void SetSize(size_t rows, size_t cols);

};//end PhenoMap



#endif /* PHENOMAP */

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
