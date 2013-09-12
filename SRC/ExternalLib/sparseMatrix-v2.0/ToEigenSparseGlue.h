#ifndef _TOEIGENSPARSEGLUE_H_
#define _TOEIGENSPARSEGLUE_H_

#include "sparseMatrix.h"
#include <vector>
#include <eigen3/Eigen/Sparse>

namespace LSW_MATH{
  
  /**
   * @class ToEigenSparseGlue export a sparse matrix in barbic's code to a
   * sparsematrix in Eigen3 library.
   *  
   */
  template<class T>
  class ToEigenSparseGlue{

	typedef std::vector<Eigen::Triplet<T> > VecT;
	
  public:
	static Eigen::SparseMatrix<T> &convert
 	 (const BARBIC_LIB::SparseMatrix & in,Eigen::SparseMatrix<T> &out){

	  const int rows = in.GetNumRows();
	  const int cols = in.GetNumColumns();
	  const int n_zeros = in.GetNumEntries();
	  VecT out_triplet;
	  convert(in,out_triplet);
	  out.resize(rows,cols);
	  out.reserve(n_zeros);
	  out.setFromTriplets( out_triplet.begin(), out_triplet.end() );
	  return out;
	}

	static VecT &convert
 	 (const BARBIC_LIB::SparseMatrix & in,VecT &out_triplet){

	  const int rows = in.GetNumRows();
	  const int n_zeros = in.GetNumEntries();
	  out_triplet.clear();
	  if (out_triplet.capacity() <= (size_t)n_zeros){
		out_triplet.reserve(n_zeros*4/3);
	  }

	  for (int r = 0; r < rows; ++r){
	  	for (int j = 0; j < in.GetRowLength(r); ++j){
	  	  const int c = in.GetColumnIndex(r,j);
	  	  const T value = in.GetEntry(r,j);
		  out_triplet.push_back(Eigen::Triplet<T>(r,c,value));
	  	}
	  }
	  return out_triplet;
	}

	/**
	 * Only convert the lower part of the input matrix to the output matxix. This
	 * is especially useful when convert the stiffness and mass matrix to the
	 * Eigen::SparseMatrix and using aparck to compute the general eigen value
	 * problem, where only the lower part of the matrix is need. 
	 * 
	 */
	static Eigen::SparseMatrix<T> &convertLower
	(const BARBIC_LIB::SparseMatrix & in,Eigen::SparseMatrix<T> &out){

	  const int rows = in.GetNumRows();
	  const int cols = in.GetNumColumns();
	  const int n_zeros = in.GetNumEntries();
	  VecT out_triplet;
	  convertLower(in, out_triplet);
	  out.resize(rows,cols);
	  out.reserve(n_zeros);
	  out.setFromTriplets( out_triplet.begin(), out_triplet.end() );
	  return out;
	}

	static VecT &convertLower
	(const BARBIC_LIB::SparseMatrix & in,VecT &out_triplet){
	  
	  const int rows = in.GetNumRows();
	  const int n_zeros = in.GetNumEntries();
	  out_triplet.clear();
	  if (out_triplet.capacity() <= (size_t)n_zeros){
		out_triplet.reserve(n_zeros);
	  }
	  for (int r = 0; r < rows; ++r){
	  	for (int j = 0; j < in.GetRowLength(r); ++j){
	  	  const int c = in.GetColumnIndex(r,j);
	  	  const T value = in.GetEntry(r,j);
		  if(r >= c){
			out_triplet.push_back(Eigen::Triplet<T>(r,c,value));
		  }
	  	}
	  }
	  return out_triplet;
	}
  };
  
}//end of namespace

#endif /*_TOEIGENSPARSEGLUE_H_*/
