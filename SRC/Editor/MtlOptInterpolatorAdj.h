#ifndef _MTLOPTINTERPOLATORADJ_H_
#define _MTLOPTINTERPOLATORADJ_H_

#include <WarpInterpolator.h>
#include <MtlOptEnergyAdj.h>
#include <MtlOptIpopt.h>

namespace LSW_ANI_EDITOR{
  
  /**
   * @class MtlOptInterpolatorAdj Interpolator using reduced RS method and elastic
   * material optimization.
   * 
   */
  class MtlOptInterpolatorAdj: public WarpInterpolator{
	
  public:
	MtlOptInterpolatorAdj();
	bool init (const string init_filename);
	bool editable(const int frame_id)const{
	  return !(ctrlF->isKeyframe(frame_id));
	}
	void setKeyframe(const vector<VectorXd> &keyZ, 
					 const vector<int>& keyframes){
	  ctrlF->setKeyframes(adjustToSubdim(keyZ),keyframes);
	}
	void removeAllPosCon(){
	  WarpInterpolator::removeAllPosCon();
	  ctrlF->clearPartialCon();
	}
	bool interpolate();
	string getName()const{
	  return string("Reduced RS and Mtl Opt.");
	}

  private:
	pCtrlForceEnergy ctrlF;
	pMtlOptEnergy mtlOpt;
	pNoConIpoptSolver ctrlFSolver;
	pNoConIpoptSolver mtlOptSolver;
	bool _optMtl;
  };
  
  typedef boost::shared_ptr<MtlOptInterpolatorAdj> pMtlOptInterpolatorAdj;
  
}//end of namespace

#endif /* _MTLOPTINTERPOLATORADJ_H_ */
