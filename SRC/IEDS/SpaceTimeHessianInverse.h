#ifndef _SPACETIMEHESSIANINVERSE_H_
#define _SPACETIMEHESSIANINVERSE_H_

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <SparseMatrixTools.h>
#include <Timer.h>
#include <assertext.h>
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;

namespace IEDS{
  
  /**
   * @class SpaceTimeHessianInverse provide some functions to calculate the
   * inverse of the hessian matrix H, which is produced by class
   * BaseSpaceTimeHessian.
   * 
   */
  class SpaceTimeHessianInverse{
	
  public:
	static void getModeBlocks(const SparseMatrix<double> &H, const int r, vector<MatrixXd> &Ms){

	  Ms.resize(r);
	  for (int i = 0; i < r; ++i){
		getModeBlock(H,r,i,Ms[i]);
	  }
	}
	static void getModeBlock(const SparseMatrix<double> &H, const int r, const int modeId, MatrixXd &M){

	  const int T = H.rows()/r;
	  std::set<int> reserve_rows_set;
	  for (int frame = 0; frame < T; ++frame){
		reserve_rows_set.insert(frame*r + modeId);
	  }
	  SparseMatrix<double> P;
	  genReshapeMatrix(H.rows(),reserve_rows_set,P,false);
	  const SparseMatrix<double> sparseM = P*(H*P.transpose());
	  M = sparseM;
	}
	static void inverse(const SparseMatrix<double> &H, const int r, SparseMatrix<double> &invH){

	  vector<MatrixXd> M;
	  Timer timer;
	  timer.start();
	  getModeBlocks(H,r,M);
	  timer.stop("getModeBlocks: ");

	  timer.start();
	  for (int i = 0; i < r; ++i){
		const MatrixXd invM = M[i].inverse();
		// cout << "cond(M): " << invM.norm()*M[i].norm() << endl;
		M[i] = invM;
	  }
	  timer.stop("inverse: ");

	  timer.start();
	  assembleFromModeBlocks(M,invH);
	  timer.stop("assemble: ");

	  assert_eq(invH.rows(),H.rows());
	  assert_eq(invH.cols(),H.cols());
	}
	static void assembleFromModeBlocks(const vector<MatrixXd> &Ms, SparseMatrix<double> &H){

	  if(Ms.size() <= 0){
		H.resize(0,0);
	  }else{
		const int T = Ms[0].rows();
		const int r = (int)Ms.size();
		const int nz = T*r*r;
		assert_gt(T,0);
		assert_gt(r,0);

		vector<Eigen::Triplet<double> > triplet;
		triplet.reserve(nz);
		for (int modeId = 0; modeId < r; ++modeId){
		  for (int i = 0; i < T; ++i){
			for (int j = 0; j < T; ++j){
			  const int row = i*r+modeId;
			  const int col = j*r+modeId;
			  const double value = Ms[modeId](i,j);
			  triplet.push_back( Eigen::Triplet<double>(row,col,value) );
			}
		  }
		}
		H.resize(T*r,T*r);
		H.setZero();
		H.reserve (nz);
		H.setFromTriplets( triplet.begin(), triplet.end() );
	  }
	}
  };
  
}//end of namespace

#endif /*_SPACETIMEHESSIANINVERSE_H_*/
