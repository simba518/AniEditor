#ifndef _RSWARPEREXT_H_
#define _RSWARPEREXT_H_

#include <RS2Euler.h>
#include <AniSeqWarper.h>

namespace IEDS{
  
  /**
   * @class RSWarperExt a warping approach using the extention of the
   * RS-coordinates
   * 
   */
  class RSWarperExt: public LSW_WARPING::AniSeqWarper{
	
  public: 
	void setTetMesh(pTetMesh_const tet_mesh){
	  AniSeqWarper::setTetMesh(tet_mesh);
	  rs2euler.setTetMesh(tetmesh);
	}

	// set the reference displacement in full space with respect to rest shape.
	void setRefU(const vector<VectorXd> &ref_u);

	// set the fixed nodes to fixed some nodes all the frames.
	void setFixedNodes(const vector<int> &c_nodes);

	// set constrained nodes and con-displacement for one frame.
	void setConNodes(const vector<set<int> > &c_nodes,const VectorXd &bcen_uc);

	// initialize rs2euler and ref_y.
	bool precompute();

	// compute T^{-1}(u): warp the linear displacement p of frame with frame_id.
	void warp(const VectorXd &p,int frame_id, VectorXd &u);

	// compute T(u): compute the RS coordinates of u.
	void computeRSCoord(const VectorXd &u, VectorXd &RS_u,const bool linear_disp=false)const;

	// get G: deformation gradient operator.
	const SparseMatrix<double> &get_G()const{
	  return rs2euler.get_G();
	}

	// retrun the RS coordinates of the input sequcence
	const vector<VectorXd> &getInputRSCoord()const{
	  return ref_y;
	}
	const VectorXd &getInputRSCoord(const int frame_id)const;

	const LSW_WARPING::RS2Euler &getRS2Euler()const{
	  return rs2euler;
	}

  private:
	LSW_WARPING::RS2Euler rs2euler;
	vector<VectorXd> ref_y;
  };
  
  typedef boost::shared_ptr<RSWarperExt> pRSWarperExt;
  
}//end of namespace

#endif /*_RSWARPEREXT_H_*/
