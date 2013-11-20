#ifndef _BASISUPDATING_H_
#define _BASISUPDATING_H_

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <boost/shared_ptr.hpp>
#include <ConMatrixTools.h>
#include <RS2Euler.h>
#include <SparseMatrixTools.h>
#include <ElasticForceTetFullStVK.h>
#include <DefGradOperator.h>
#include <RSCoordComp.h>
#include <MassMatrix.h>
#include <MatrixTools.h>
#include <MapMA2RS.h>
#include <assertext.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;
using namespace LSW_WARPING;
using namespace EIGEN3EXT;

namespace LSW_ANI_EDITOR{
  
  /// @class BasisUpdating given keyframe uk, update the basis B and W.
  class BasisUpdating{
	
  public:
	BasisUpdating(pTetMesh_const tet,const set<int>&fixed_nd,
				  const double tolB=1e-3,const double tolW=1e-3):
	  tetmesh(tet),fixed_nodes(fixed_nd),tolForUpdateB(tolB),tolForUpdateW(tolW){
	  precompute();
	}

	BasisUpdating(pTetMesh_const tet,const double tolB=1e-3,const double tolW=1e-3):
	  tetmesh(tet),tolForUpdateB(tolB),tolForUpdateW(tolW){
	  precompute();
	}

	// update basis
	void update(const MatrixXd&Uk,MatrixXd&B,MatrixXd&W,MatrixXd&subK,bool useOptSK=true);
	void update(const MatrixXd&Uk,MatrixXd &W,MatrixXd &subK,bool useOptSubK=true);

	// compute reduced keyframes
	void computeZk(const MatrixXd&Uk,const MatrixXd &B,const MatrixXd &W,MatrixXd &Zk);
	void computeZk(const MatrixXd&Uk,const MatrixXd &W,MatrixXd &Zk);

	static int updateB(const MatrixXd&Uk,MatrixXd &B,const double tolB);
	static int updateW(const MatrixXd&Pk,const DiagonalMatrix<double,-1>&M,
					   MatrixXd &W,const double tolB);

  protected:
	void precompute();
	void computeLinearP(const MatrixXd &U,MatrixXd &LinearP)const;
	template<class VECTOR_TMP>
	inline static double RMS_norm(const VECTOR_TMP &u){
	  assert_gt(u.size(),0);
	  return u.norm()/u.size();
	}
	template<class VECTOR_TMP,class MATRIX_TMP>
	inline static double RMS_norm2_M(const VECTOR_TMP&u,const MATRIX_TMP&M){
	  assert_gt(u.size(),0);
	  assert_eq(M.rows(),M.cols());
	  assert_eq(M.rows(),u.size());
	  return u.dot(M*u)/u.size();
	}

  private:
	pTetMesh_const tetmesh;
	set<int> fixed_nodes;
	double tolForUpdateB;
	double tolForUpdateW;
	SparseMatrix<double> K;
	DiagonalMatrix<double,-1> M;

	SparseMatrix<double> G, PG, rm_fixed_nodes, A;
	SparseQR<SparseMatrix<double>,COLAMDOrdering<int> > solver;
  };
  
  typedef boost::shared_ptr<BasisUpdating> pBasisUpdating;
  
}//end of namespace

#endif /*_BASISUPDATING_H_*/
