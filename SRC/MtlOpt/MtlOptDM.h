#ifndef _MTLOPTDM_H_
#define _MTLOPTDM_H_

#include <vector>
#include <set>
#include <utility>
#include <eigen3/Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <assertext.h>
#include <Log.h>
using namespace std;
using namespace Eigen;

namespace MTLOPT{
  
  /**
   * @class MtlOptDM data model for the algorithm of material optimization.
   * 
   */
  class MtlOptDM{
	
  public:
	MtlOptDM();

	void setTotalFrames(const int T);
	void setDamping(const double ak,const double am);
	void initVariables(const VectorXd &lambda0, const int rs);

	void fixHeadFrames(const int numFrames);
	void fixTailFrames(const int numFrames);
	void setMechanicConFrames(const int numFrames);

	void setConGroups(const int frame_id,const vector<int>&g,const VectorXd&pc);
	void setUc(const int frame_id, const VectorXd &pc);
	void removeAllPosCon();
	void removePartialCon(const int frame_id);
	void sortConGroups();

	void setKeyfames(const vector<int> &keyf, const MatrixXd &zk);

	int numConNodes()const;
	int reducedDim()const;
	int subFrames()const;

	void updateZ(const double *subZ);
	int expandS(const int rs);
	void decomposeK();

	void print()const;
	
  public:
	// input information
	int T;
	double alphaK0;
	double alphaM0;
	VectorXd Lambda0;

	// constraints
	int mechanic_con_frames;
	int T_begin; // begin of the editable frames,i.e [0,T_begin-1] is fixed.
	int T_end; // end of the editable frames, i.e [T_end+1,T-1] is fixed.
	vector<int> conFrames;
	vector<vector<int> > conNodes;
	vector<VectorXd> uc;
	vector<int> keyframes;
	MatrixXd Zk;
	
	// include all frames with fixed frames to be zero, e.g Z.cols()==T.
	MatrixXd K,Z,S; 
	VectorXd Lambda, Damping, alphaM, alphaK;
  };
  
  typedef boost::shared_ptr<MtlOptDM> pMtlOptDM;
  typedef boost::shared_ptr<const MtlOptDM> pMtlOptDM_const;
  
}//end of namespace

#endif /*_MTLOPTDM_H_*/
