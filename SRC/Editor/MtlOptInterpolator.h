#ifndef _MTLOPTINTERPOLATOR_H_
#define _MTLOPTINTERPOLATOR_H_

#include <WarpInterpolator.h>
#include <MtlOptEnergy.h>
#include <MtlOptIpopt.h>

namespace LSW_ANI_EDITOR{
  
  /**
   * @class MtlOptInterpolator Interpolator using reduced RS method and elastic
   * material optimization.
   * 
   */
  class MtlOptInterpolator: public WarpInterpolator{
	
  public:
	MtlOptInterpolator();
	bool init (const string init_filename);
	bool editable(const int frame_id)const{
	  return !(_ctrlF->isKeyframe(frame_id));
	}
	void setKeyframe(const vector<VectorXd> &keyZ, 
					 const vector<int>& keyframes){
	  _ctrlF->setKeyframes(adjustToSubdim(keyZ),keyframes);
	}
	void removeAllPosCon(){
	  WarpInterpolator::removeAllPosCon();
	  _ctrlF->clearPartialCon();
	}
	bool interpolate();
	string getName()const{
	  return string("Reduced RS and Mtl Opt.");
	}

  private:
	pCtrlForceEnergy _ctrlF;
	pMtlOptEnergy _mtlOpt;
	pNoConIpoptSolver _ctrlFSolver;
	pNoConIpoptSolver _mtlOptSolver;
	bool _optMtl;
  };
  
  typedef boost::shared_ptr<MtlOptInterpolator> pMtlOptInterpolator;
  
}//end of namespace

#endif /* _MTLOPTINTERPOLATOR_H_ */
