#ifndef _MTLOPTINTERPOLATOR_H_
#define _MTLOPTINTERPOLATOR_H_

#include <set>
#include <boost/shared_ptr.hpp>
#include <BaseInterpolator.h>
#include <RSWarperExt.h>
#include <RedRSWarperExt.h>
#include <ConNodesOfFrame.h>
#include <ModalModeDisplayer.h>
#include <MtlOptEnergy.h>
#include <MtlOptIpopt.h>
using namespace std;
using namespace IEDS;
using namespace LSW_WARPING;
using namespace LSW_SIM;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class MtlOptInterpolator Interpolator using reduced RS method and elastic
   * material optimization.
   * 
   */
  class MtlOptInterpolator: public BaseInterpolator{
	
  public:
	MtlOptInterpolator();

	bool init (const string init_filename);

	void useWarp(const bool use_warp){
	  this->use_warp = use_warp;
	}

	bool isUseWarp()const{
	  return use_warp;
	}

	bool editable(const int frame_id)const{
	  return !(ctrlF->isKeyframe(frame_id));
	}

	void setAllConGroups(const set<pConNodesOfFrame> &newCons);

	void setConGroups(const int f,const vector<set<int> >&g,const VectorXd&uc){
	  addConGroups(f,g,uc);
	}

	void setUc(const int frame_id,const VectorXd &uc);

	void removeAllPosCon();

	bool interpolate ();

	void setReducedEdits(const int frame_id, const VectorXd &z_i){
	  assert_in(frame_id,0,getT());
	  assert_eq(getT(),(int)delta_z.size());
	  assert_eq(z_i.size(),reducedDim());
	  delta_z[frame_id] = z_i;
	}

	void setReducedEdits(const vector<VectorXd> &z){
	  assert_eq(getT(),(int)z.size());
	  this->delta_z = z;
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
	string getName()const{
	  return string("Reduced RS and Mtl Opt.");
	}

  protected:
	bool initWarper(JsonFilePaser &inf);
	bool validConstraints(const int frame_id)const{
	  const bool valid = editable(frame_id);
	  WARN_LOG_COND("the frame id of the constraints is invalid:"<<frame_id,valid);
	  return valid;
	}
	void addConGroups(const int f,const vector<set<int> >&g,const VectorXd&uc);
	void sortConstraints();
	void removeConOfFrame(const int frame_id);
	bool loadUref(JsonFilePaser &json_f,vector<VectorXd> &u_ref)const;
	
  private:
	// pCtrlForceEnergyNoWarp ctrlF;
	pCtrlForceEnergy ctrlF;
	pMtlOptEnergy mtlOpt;
	pNoConIpoptSolver ctrlFSolver;
	pNoConIpoptSolver mtlOptSolver;
	bool _optMtl;

	vector<VectorXd> u_ref; // reference sequence in fullspace.
	vector<VectorXd> delta_z; // increamental displacements in linear subspace.

	MatrixXd W; // linear basis of MA	
	MatrixXd B; // nonlinear bais
  	VectorXd full_u; // displacements in full space.
  	VectorXd unwarped_full_u; // unwarped displacements in full space.
	pRSWarperExt warper;
	pRedRSWarperExt nodeWarper;
	bool use_warp;

	// partial control
	vector<int> con_frame_id;   // frames that has position constraints.
	vector<vector<int> > con_nodes;   // con nodes of each con frames.
	vector<VectorXd> uc; // target displacements with respect to the rest.

	ModalModeDisplayer modalDisplayer;
  };
  
  typedef boost::shared_ptr<MtlOptInterpolator> pMtlOptInterpolator;
  
}//end of namespace

#endif /* _MTLOPTINTERPOLATOR_H_ */
