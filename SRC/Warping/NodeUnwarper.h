#ifndef _NODEUNWARPER_H_
#define _NODEUNWARPER_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include <RSCoordComp.h>
#include <NodeWarper.h>

namespace LSW_WARPING{
  
  /**
   * @class NodeUnwarper given a well deformed shape x, compute the rotation
   * vector and local displacement based on polar decomposition. These values
   * can be used as initial values to recover the final rotation vector w and
   * local displacement p as in: GPU-Friendly Shape Interpolation Based on
   * Trajectory Warping [CASA2005].
   * 
   */
  class NodeUnwarper{
	
  public:
	static void compute(pVolumetricMesh_const tetmesh, const SparseMatrix<double> &G, 
						const VectorXd &u, VectorV3 &w, VectorXd &p){
	  // compute w
	  VectorXd y;
	  RSCoordComp::constructWithoutWarp(G,u,y);
	  RSCoordComp::RS2NodeRotVec(tetmesh,y,w);

	  // compute p using w
	  p.resize(w.size()*3);
	  for (size_t i = 0; i < w.size(); ++i){
		const Matrix3d _R = NodeWarper::warpMat(w[i]);
		p.segment(i*3,3) = (_R.inverse())*u.segment(i*3,3);
	  }
	}
	
  };
  
  typedef boost::shared_ptr<NodeUnwarper> pNodeUnwarper;
  
}//end of namespace

#endif /*_NODEUNWARPER_H_*/
