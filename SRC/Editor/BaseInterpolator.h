#ifndef _BASEINTERPOLATOR_H_
#define _BASEINTERPOLATOR_H_

#include <string>
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class BaseInterpolator an interface of the class for interpolating
   * animation sequence.
   * 
   * 1. Select one frame for edit, and set constrained nodes' groups and uc.
   * 2. Update uc for one frame. If full control is used, compute full control zk.
   * 3. Interpolate to satisfy the all constraints of all constrained frames.
   * 
   */
  class BaseInterpolator{
	
  public:
	/// read data from initfile, such as dampings and timestep.
	virtual bool init (const string init_filename) = 0;

	/// toggle use warping or not, nothing to do in default.
	virtual void useWarp(const bool use_warp){
	  
	}

	virtual bool isUseWarp()const{
	  return false;
	}

	/// update reference dense sequence. (it is "update", not set!)
	virtual void updateRefSeq (const vector<VectorXd> &q_ref){}

	/// update reference dense sequence using current interpolated result.
	virtual void updateRefSeq(){}

	/// set keyframes q_c with the frame indexes in keyframe_ids.
	virtual void setKeyframe (const vector<VectorXd> &q_c, 
							  const vector<int>& keyframe_ids){}

	virtual bool isKeyframe(const int frame_id)const{
	  return false;
	}

	virtual const VectorXd &getKeyframe(const int frame_id){
	  static VectorXd tempt_u;
	  return tempt_u;
	}

	/// Return true if the frame with frame_id can be edit. usually, only the
	/// frame in [0,T-1], except the boundary frames, can be edited.
	virtual bool editable(const int frame_id)const = 0;

	/// partial control (position constraints of selected frames)
	virtual void setConGroups(const int frame_id,const vector<set<int> >&group, 
							  const VectorXd &uc) = 0;
	virtual void setUc(const int frame_id, const VectorXd &uc) = 0;
	virtual void removeAllPosCon() = 0;

	/// set initial value of the control forces
	virtual void setInitCtrlForces(const VectorXd &w){}

	/// perform the interpolation process.
	virtual bool interpolate () = 0;

	virtual void printSolveInf()const{}

	/// return the interpolated displacements of frame frame_id in full space.
	virtual const VectorXd& getInterpU(const int frame_id) = 0;

	/// return the input displacements of frame frame_id in full space.
	virtual const VectorXd& getInputU(const int frame_id) = 0;

	/// return the unwarped displacements, used for manipulation. By default, it
	/// is equal to the interpolated displacements (when no warping used).
	virtual const VectorXd& getUnwarpedU(const int frame_id){
	  return getInterpU(frame_id);
	}

	virtual const VectorXd& getUforConstraint(const int frame_id){
	  return getInterpU(frame_id);
	}

	// return the additional edits of frame_id in subspace.
	virtual const VectorXd& getReducedEdits(const int frame_id)const = 0;

	// set the additional edits z_i to frame_id in subspace
	virtual void setReducedEdits(const int frame_id, const VectorXd &z_i){}

	// set all additional edits z in subspace
	virtual void setReducedEdits(const vector<VectorXd> &z){}

	/// return the total frame number.
	virtual int getT()const = 0;

	virtual int reducedDim()const = 0;

	// return the name of the interpolation method.
	virtual string getName()const = 0;
	
	// given the name of the private data, then save it to the file.
	virtual bool save(const string data_name, const string file_name)const{
	  cout << "unspported data saving operation." << endl;
	  return false;
	}

	// given the name of the private data, print it to the console.
	virtual void print(const string data_name)const{
	  cout << "unspported data print operation." << endl;
	}
	
  };
  
  typedef boost::shared_ptr<BaseInterpolator> pBaseInterpolator;
  typedef boost::shared_ptr<const BaseInterpolator> pBaseInterpolator_const;
  
}//end of namespace

#endif /*_BASEINTERPOLATOR_H_*/
