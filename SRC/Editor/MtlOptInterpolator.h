#ifndef _MTLOPTINTERPOLATOR_H_
#define _MTLOPTINTERPOLATOR_H_

#include <JsonFilePaser.h>
#include <BaseInterpolator.h>
#include <MtlOptSolver.h>
#include <MtlOptSimulator.h>
using namespace UTILITY;
using namespace MTLOPT;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class MtlOptInterpolator Interpolator using reduced RS method and elastic
   * material optimization.
   * 
   */
  class MtlOptInterpolator: public BaseInterpolator{
	
  public:
	// set
	MtlOptInterpolator();

	bool init (const string init_filename);

	void setConGroups(const int frame_id,const vector<set<int> >&group, 
							  const Eigen::Matrix<double,3,-1> &pc);

	void setUc(const int frame_id, const Eigen::Matrix<double,3,-1> &pc);

	void removeAllPosCon(){
	  dataModel->removeAllPosCon();
	  mtlopt->getEcZ()->updateConstraints();
	  mtlopt->getEcS()->updateConstraints();
	}

	void removePartialCon(const int frame_id){
	  dataModel->removePartialCon(frame_id);
	  mtlopt->getEcZ()->updateConstraints();
	  mtlopt->getEcS()->updateConstraints();
	}

	void setKeyframes(const MatrixXd &U);

	void setKeyframes(const MatrixXd &U,const vector<int>& keyframe_ids);


	// solve
	bool interpolate ();

	// get
	const VectorXd& getInterpU(const int frame_id);

	const VectorXd& getInputU(const int frame_id){
	  assert_in (frame_id, 0, u_ref.size()-1);
	  return u_ref[frame_id];
	}

	const VectorXd& getReducedEdits(const int frame_id)const{

	  static VectorXd zi;/// @bug
	  assert_in(frame_id,0,getT()-1);
	  if (frame_id < dataModel->Z.cols()){
		zi = dataModel->Z.col(frame_id);
	  }else{
		zi = simulator->getZ( frame_id - dataModel->Z.cols() );
	  }
	  return zi;
	}

	bool isUseWarp()const{
	  return true;
	}

	bool editable(const int frame_id)const{
	  return (frame_id >= dataModel->T_begin && frame_id <= dataModel->T_end );
	}

	int getT()const{
	  return dataModel->T+simulator->totalFrames();
	}

	int reducedDim()const{
	  return dataModel->reducedDim();
	}

	string getName()const{
	  if (mtlopt->optimizeMtl()){
		return string("Mtl Opt.");
	  }
	  return string("No Mtl Opt.");
	}

	void printSolveInf()const;

  protected:
	bool initWarper(JsonFilePaser &inf);
	bool loadUref(JsonFilePaser &json_f,vector<VectorXd> &u_ref)const;
	void sortCubPoints(vector<int> &cubPoints, vector<double> &cubWeights);
	void getAllRedDisp(MatrixXd &Z)const;
	void getAllRedCtrlF(MatrixXd &Ef)const;

  private:
	pMtlOptDM dataModel;
	pReducedRS reducedRS;
	ReducedRSUnwarp unwarp;
	pMtlOptSolver mtlopt;
	pMtlOptSimulator simulator;

	vector<VectorXd> u_ref;
	VectorXd full_u_i;
	string save_S_to;
	MatrixXd scaled_W;
 };
  
  typedef boost::shared_ptr<MtlOptInterpolator> pMtlOptInterpolator;
  
}//end of namespace

#endif /* _MTLOPTINTERPOLATOR_H_ */
