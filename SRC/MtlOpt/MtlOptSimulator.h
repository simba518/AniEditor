#ifndef _MTLOPTSIMULATOR_H_
#define _MTLOPTSIMULATOR_H_

#include "MtlOptDM.h"
#include <HarmonicOscillator.h>

namespace MTLOPT{

  /**
   * @class MtlOptSimulator simulate a sequence using the optimized material.
   * 
   */
  class MtlOptSimulator{
	
  public:
	MtlOptSimulator(pMtlOptDM dm):dm(dm){}
	void setTotalFrames(const int T){
	  assert_ge(T,0);
	  Z.resize(T);
	  assert(dm);
	  for (int i = 0; i < Z.size(); ++i){
		Z[i].resize(dm->reducedDim());
		Z[i].setZero();
	  }
	}

	virtual void simulate() = 0;

	const VectorXd &getZ(const int frameId){
	  assert_in(frameId,0,Z.size()-1);
	  assert(dm);
	  assert_eq(Z[frameId].size(),dm->reducedDim());
	  return Z[frameId];
	}
	const vector<VectorXd> &getZ()const{
	  return Z;
	}
	int totalFrames()const{
	  return Z.size();
	}

	void checkDimension()const{
	  assert(dm);
	  assert_ge(dm->Z.cols(),3);
	  assert_eq(dm->Z.cols(),dm->T);
	}

  protected:
	pMtlOptDM dm;	
	vector<VectorXd> Z;
  };
  typedef boost::shared_ptr<MtlOptSimulator> pMtlOptSimulator;

  /**
   * @class MtlOptSimulator simulate a sequence using the optimized material
   * with implicit integration.
   * 
   */
  class MtlOptSimulatorImplicit:public MtlOptSimulator{
	
  public:
	MtlOptSimulatorImplicit(pMtlOptDM dm):MtlOptSimulator(dm){}
	void simulate(){

	  checkDimension();
	  if (0 == Z.size()){
		WARN_LOG("the simulated sequence length is 0");
		return;
	  }

	  const VectorXd &z0 = dm->Z.col(dm->T-2);
	  const VectorXd &z1 = dm->Z.col(dm->T-1);
	  if (Z.size() > 0)
		forward(z0,z1,Z[0]);

	  if (Z.size() > 1)
		forward(z1,Z[0],Z[1]);

	  for (int i = 1; i < Z.size()-1; ++i)
		forward(Z[i-1],Z[i],Z[i+1]);
	}

  protected:
	void forward(const VectorXd &z0,const VectorXd &z1,VectorXd &z2)const{

	  assert(dm);
	  const VectorXd &lambda = dm->Lambda;
	  const VectorXd &d = dm->Damping;
	  const VectorXd s = (2.0f*VectorXd::Ones(d.size())-lambda+d).asDiagonal()*z1 - z0;

	  z2.resize(lambda.size());
	  for (int i = 0; i < z2.size(); ++i){
		const double d1 = d[i]+1.0f;
		assert_gt(d1,0.0f);
		z2[i] = s[i]/d1;
	  }
	}
  };
  typedef boost::shared_ptr<MtlOptSimulatorImplicit> pMtlOptSimulatorImplicit;

  /**
   * @class MtlOptSimulator simulate a sequence using the optimized material
   * with analytic simulator.
   * 
   */
  class MtlOptSimulatorAnalytic:public MtlOptSimulator{
	
  public:
	MtlOptSimulatorAnalytic(pMtlOptDM dm):MtlOptSimulator(dm){}
	void simulate(){
	  
	  checkDimension();
	  if (0 == Z.size()){
		WARN_LOG("the simulated sequence length is 0");
		return;
	  }

	  assert_gt(Z.size(),3);
	  const VectorXd &lambda = dm->Lambda;
	  const VectorXd &d = dm->Damping;
	  const VectorXd &z0 = dm->Z.col(dm->T-2);
	  const VectorXd &z1 = dm->Z.col(dm->T-1);
	  const VectorXd v0 = z1-z0;
	  UTILITY::HarmonicOscillatorSet<double> oscillator(lambda,d,z0,v0);
	  const MatrixXd allZ = oscillator.generateSequence<MatrixXd>(0,1.0f,Z.size()+2);
	  assert_eq(allZ.cols(),Z.size()+2);
	  assert_eq(allZ.rows(),Z[0].size());
	  for (int i = 2; i < allZ.cols(); ++i){
		Z[i-2] = allZ.col(i);
	  }
	}
  };
  typedef boost::shared_ptr<MtlOptSimulatorAnalytic> pMtlOptSimulatorAnalytic;
  
}//end of namespace

#endif /*_MTLOPTSIMULATOR_H_*/
