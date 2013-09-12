#ifndef _REDRSINTERPOLATOR_H_
#define _REDRSINTERPOLATOR_H_

#include <set>
#include <boost/shared_ptr.hpp>
#include <BaseInterpolator.h>
#include <FastIntuitiveAniEditor.h>
#include <RSWarperExt.h>
#include <RedRSWarperExt.h>
#include <ConNodesOfFrame.h>
#include <ModalModeDisplayer.h>
using namespace std;
using namespace IEDS;
using namespace LSW_WARPING;
using namespace LSW_SIM;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class RedRSInterpolator intuitive animation interpolator based on reduced
   * RS method.
   * 
   */
  class RedRSInterpolator: public BaseInterpolator{
	
  public:
	RedRSInterpolator();

	bool init (const string init_filename);

	// warping
	void useWarp(const bool use_warp){
	  this->use_warp = use_warp;
	}

	bool isUseWarp()const{
	  return use_warp;
	}

	bool editable(const int frame_id)const;

	void setAllConGroups(const set<pConNodesOfFrame> &newCons);

	void setConGroups(const int f,const vector<set<int> >&g,const VectorXd&uc);

	void setUc(const int frame_id, const VectorXd &uc);

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

	// return the reference one only
	const vector<VectorXd>& getInterpSeq()const{
	  return getRefSeq();
	}
	const vector<VectorXd>& getDeltaInterpSeq()const;
	const vector<VectorXd>& getRefSeq()const{
	  return u_ref;
	}

	const VectorXd& getInterpU(const int frame_id);
	const VectorXd& getInputU(const int frame_id);
	const VectorXd& getUnwarpedU(const int frame_id);

	// get dimensions
	int nonlinearSubDim()const{
	  return 0;
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
	  return string("Reduced RS");
	}

	void print(const string data_name)const;

  protected:
	bool initWarper(const string init_filename);
	bool validConstraints(const int frame_id)const;
	void addConGroups(const int f,const vector<set<int> >&g,const VectorXd&uc);
	void sortConstraints();
	void removeConOfFrame(const int frame_id);
	
  private:
	FastIntuitiveAniEditor anieditor;
	vector<VectorXd> u_ref; // reference sequence in fullspace.
	vector<VectorXd> delta_z; // increamental displacements in linear subspace.

	MatrixXd W; // linear basis of MA	
	MatrixXd B; ///@todo nonlinear bais
  	VectorXd full_u; // displacements in full space.
  	VectorXd unwarped_full_u; // unwarped displacements in full space.
	pRSWarperExt warper;
	pRedRSWarperExt nodeWarper;
	bool use_warp;

	// partial control
	vector<int> con_frame_id;   // frames that has position constraints
	vector<vector<int> > con_nodes;   // con nodes of each con frames
	vector<VectorXd> uc; // target displacements with respect to the rest.
	pVolumetricMesh_const tet_mesh;
	ModalModeDisplayer modalDisplayer;
  };
  
  typedef boost::shared_ptr<RedRSInterpolator> pRedRSInterpolator;
  
}//end of namespace

#endif /* _REDRSINTERPOLATOR_H_*/
