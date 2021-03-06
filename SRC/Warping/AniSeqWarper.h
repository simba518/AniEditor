#ifndef _ANISEQWARPER_H_
#define _ANISEQWARPER_H_

#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <TetMesh.h>

using namespace std;
using namespace Eigen;
using namespace UTILITY;

namespace LSW_WARPING{
  
  /**
   * @class AniSeqWarper: Base class providing only interfaces for, given an
   * input animation sequence and some linear edits, construct the well deformed
   * shapes using modal warping approaches.
   * @see RSWarperExt.h and 
   */
  class AniSeqWarper{
	
  public:
	virtual void setTetMesh(pTetMesh_const mesh){
	  this->tetmesh = mesh;
	}

	// set the ref disp in full space with respect rest shape.
	virtual void setRefU(const vector<VectorXd> &ref_u) = 0;

	// set the fixed nodes to fixed some nodes all the frames.
	virtual void setFixedNodes(const vector<int> &c_nodes){}

	// set constrained nodes and con-disp for one frame.
	virtual void setConNodes(const vector<set<int> > &c_nodes,const VectorXd &bcen_uc){}

	virtual bool precompute() = 0;

	/** 
	 * @param p additional linear displacement with respect to the corresponding
	 * reference frame.
	 * @param frame_id the frame index, which is 0-indexed.
	 * @param u the warped displacement which respect to the rest shape, not the
	 * reference frame.
	 */
	virtual void warp(const VectorXd&p,int f_id,VectorXd&u)=0;

	pTetMesh_const getTetMesh()const{
	  return tetmesh;
	}

  protected:
	pTetMesh_const tetmesh;

  };
  
  typedef boost::shared_ptr<AniSeqWarper> pAniSeqWarper;
  
}//end of namespace

#endif
