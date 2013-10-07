#ifndef _ROTATIONMATRIX_H_
#define _ROTATIONMATRIX_H_

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <boost/shared_ptr.hpp>
#include <RSWarperExt.h>
using namespace IEDS;

namespace LSW_WARPING{
  
  /**
   * @class RotationMatrix compute the rotational matrix of each node of a
   * tetrahedron mesh.
   * 
   */
  class RotationMatrix{
	
  public:
	RotationMatrix(pRSWarperExt warper, pTetMesh_const tet_mesh):
	  _warper(warper),_tet_mesh(tet_mesh){}

	void compute(const VectorXd &u, SparseMatrix<double> &R)const;

  protected:

  private:
	pRSWarperExt _warper;
	pTetMesh_const _tet_mesh;
  };
  
  typedef boost::shared_ptr<RotationMatrix> pRotationMatrix;
  
}//end of namespace

#endif /*_ROTATIONMATRIX_H_*/

