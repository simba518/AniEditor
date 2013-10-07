#ifndef _NODEROTVEC_H_
#define _NODEROTVEC_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include <RSCoordComp.h>
#include <TetMesh.h>
using namespace UTILITY;

namespace LSW_WARPING{
  
  /**
   * @class NodeRotVec compute the rotational vector of each ndoe.
   * 
   */
  class NodeRotVec{
	
  public:
	static void compute(pTetMesh_const tetmesh, const SparseMatrix<double> &G, const VectorXd &p, VVec3d &w){

	  VectorXd y;
	  RSCoordComp::constructWithWarp(G,p,y);
	  RSCoordComp::RS2NodeRotVec(tetmesh,y,w);
	}
  };
  
}//end of namespace

#endif /*_NODEROTVEC_H_*/
