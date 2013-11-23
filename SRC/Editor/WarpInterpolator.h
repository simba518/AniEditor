#ifndef _WARPINTERPOLATOR_H_
#define _WARPINTERPOLATOR_H_

#include <set>
#include <boost/shared_ptr.hpp>
#include <BaseInterpolator.h>
#include <RSWarperExt.h>
#include <RedRSWarperExt.h>
#include <ModalModeDisplayer.h>
using namespace std;
using namespace IEDS;
using namespace LSW_WARPING;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class WarpInterpolator base class for interpolation algorithm using warping.
   * 
   */
  class WarpInterpolator:public BaseInterpolator{
	
  public:
	WarpInterpolator(){
	  warper = pRSWarperExt(new RSWarperExt());
	  use_warp = false;
	}

	virtual bool init (const string init_filename);

	void useWarp(const bool use_warp){
	  this->use_warp = use_warp;
	}

	bool isUseWarp()const{
	  return use_warp;
	}

	virtual void setConGroups(const int f,const vector<set<int> >&g,const Matrix<double,3,-1>&uc){
	  vector<vector<set<int> > >group_rlst;
	  vector<Eigen::Vector3d> uc_rlst;
	  splitAllConstraints(g,uc,group_rlst,uc_rlst);
	  for (int i = 0; i < uc_rlst.size(); ++i){
		addConGroups(f,group_rlst[i],uc_rlst[i]);
	  }
	}

	virtual void setUc(const int frame_id,const Matrix<double,3,-1> &uc);

	virtual void removeAllPosCon(){
	  con_frame_id.clear();
	  con_nodes.clear();
	  uc.clear();
	}

	void removePartialCon(const int frame_id){
	  ///@todo
	  ERROR_LOG("undefined function.");
	}

	void setReducedEdits(const int frame_id, const VectorXd &z_i){
	  assert_in(frame_id,0,getT());
	  assert_eq(getT(),(int)delta_z.size());
	  delta_z[frame_id] = adjustToSubdim(z_i);
	}

	void setReducedEdits(const vector<VectorXd> &z){
	  assert_eq(getT(),(int)z.size());
	  this->delta_z = adjustToSubdim(z);
	}

	const vector<VectorXd>& getInterpSeq()const{
	  return getRefSeq();
	}
	const vector<VectorXd>& getDeltaInterpSeq()const{
	  return delta_z;
	}
	const vector<VectorXd>& getRefSeq()const{
	  return u_ref;
	}

	const VectorXd& getInterpU(const int frame_id);
	const VectorXd& getInputU(const int frame_id){
	  assert_in (frame_id, 0, getT()-1);
	  return u_ref[frame_id];
	}
	const VectorXd& getUnwarpedU(const int frame_id){
	  unwarped_full_u = u_ref[frame_id] + W * delta_z[frame_id]; ///@todo U^T
	  return unwarped_full_u;
	}

	// get dimensions
	int nonlinearSubDim()const{
	  return B.cols();
	}
	int linearSubDim()const{
	  return W.cols();
	}
	int reducedDim()const{
	  return linearSubDim();
	}
	int getT()const{
	  return (int)u_ref.size();
	}
	const VectorXd& getReducedEdits(const int frame_id)const{
	  assert_in(frame_id,0,getT());
	  assert_eq(getT(),(int)delta_z.size());
	  assert_eq(delta_z[frame_id].size(),reducedDim());
	  return delta_z[frame_id];
	}

	// more 
	pAniSeqWarper getWarper()const{
	  return warper;
	}
	const VectorXd &getWarpU(const int frame_id,
							 const vector<set<int> > &c_nodes,
							 const VectorXd &bcen_uc);

  protected:
	bool initWarper(JsonFilePaser &inf);
	bool validConstraints(const int frame_id)const{
	  const bool valid = editable(frame_id);
	  WARN_LOG_COND("the frame id of the constraints is invalid:"<<frame_id,valid);
	  return valid;
	}
	void addConGroups(const int f,const vector<set<int> >&g,const Matrix<double,3,-1>&uc);
	void removeConOfFrame(const int frame_id);
	bool loadUref(JsonFilePaser &json_f,vector<VectorXd> &u_ref)const;
	VectorXd adjustToSubdim(const VectorXd &z)const;
	vector<VectorXd> adjustToSubdim(const vector<VectorXd> &vz)const;
	
  protected:
	pRSWarperExt warper;
	pRedRSWarperExt nodeWarper;
	bool use_warp;

	vector<VectorXd> u_ref; // reference sequence in fullspace.
	vector<VectorXd> delta_z; // increamental displacements in linear subspace.

	VectorXd Lambda; // eigen values
	MatrixXd W; // linear basis of MA
	MatrixXd B; // nonlinear bais
  	VectorXd full_u; // displacements in full space.
  	VectorXd unwarped_full_u; // unwarped displacements in full space.

	vector<int> con_frame_id;   // frames that has position constraints.
	vector<vector<int> > con_nodes;   // con nodes of each con frames.
	vector<VectorXd> uc; // target displacements with respect to the rest.

	ModalModeDisplayer modalDisplayer;
  };
  
}//end of namespace

#endif /*_WARPINTERPOLATOR_H_*/
