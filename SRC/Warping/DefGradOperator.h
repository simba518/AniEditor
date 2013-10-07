#ifndef _DEFGRADOPERATOR_H_
#define _DEFGRADOPERATOR_H_

#include <volumetricMesh.h>
#include <eigen3/Eigen/Sparse>
#include "ComputeGeAutoDiff.h"
#include <TetMesh.h>
using namespace Eigen;
using namespace UTILITY;

namespace LSW_WARPING{
  
  /**
   * @class DefGradOperator compute the discrete deformation gradient operator G
   * of a tetrahedron mesh, which is usually used to compute the deformation
   * gradient 
   * m = G*x or m = G*u + I. 
   * 
   * where m a vector with length of 9e, e.g. constructed by assembling all the
   * row-major 3x3 deformation gradient matrix of each tetrahedron. And x=r + u
   * is the deformed shape of the object, where r the rest shape and u is the
   * displacement.
   * 
   * @note it's the operator, not the deformation itself.
   * @note m is a 9x1 vector which represents a row-major 3x3 matrix, that is
   *                                                        |m0, m1, m2|
   * m = (m0,m1,m2,m3,m4,m5,m6,m7,m8), then the matrix is   |m3, m4, m5|
   *                                                        |m6, m7, m8|
   * @see the study note: Discrete-Gradient-Operator.tex
   */
  class DefGradOperator{

	typedef Eigen::Triplet<double> T;
	
  public:
	static bool compute(pTetMesh_const tet_mesh,SparseMatrix<double> &G);

	// assemble the gradient operator Ge of one tetrahedron to the tripletlist
	// of the full operator G.
	static inline void assembleGe2G(const int *elem_v,const int e, const double Ge[9][12], vector<T> &G_T){

	  assert(elem_v != NULL);
	  const int r0 = e*9;
	  for (int j = 0; j < 4; ++j){
    
		const int v = elem_v[j];
		for (int i = 0; i < 9; ++i){
		  if (Ge[i][j*3+0] != 0)
			G_T.push_back(  T(r0+i, v*3+0 ,Ge[i][j*3+0])  );

		  if (Ge[i][j*3+1] != 0)
			G_T.push_back(  T(r0+i, v*3+1 ,Ge[i][j*3+1])  );

		  if (Ge[i][j*3+2] != 0)
			G_T.push_back(  T(r0+i, v*3+2 ,Ge[i][j*3+2])  );
		}
	  }
	}

  protected:
	static inline void vertexOfTetEle(pTetMesh_const tet_mesh,int e,double V[4][3]){
	  const Vector4i &tet = tet_mesh->tets()[e];
	  const VVec3d &nodes = tet_mesh->nodes();
	  for (int i = 0; i < 4; ++i){
		const Vector3d &v = nodes[tet[i]];
		V[i][0] = v[0];
		V[i][1] = v[1];
		V[i][2] = v[2];
	  }
	}
	
  };
  
}//end of namespace

#endif /*_DEFGRADOPERATOR_H_*/
