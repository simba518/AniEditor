#ifndef _MTLOPTSOLVER_H_
#define _MTLOPTSOLVER_H_

#include "CtrlForcesEnergy.h"
#include "CtrlForcesEnergyAD.h"
#include "ReducedRS.h"
#include "PosConEnergy.h"
#include "KeyframeConEnergy.h"
#include "BasisSelRegEnergy.h"
#include <alglib/optimization.h>

namespace MTLOPT{
  
  /**
   * @class MtlOptSolver the solver for the material optimization algorithm.
   * 
   */
  class MtlOptSolver{
	
  public:
	MtlOptSolver(pMtlOptDM dm,pReducedRS_const rs,const int mIt_=50,
				 const double tol_=1e-12,const double maxFreq = 4.0f);

	void setMaxIt(const int maxIt_, const int maxItz, const int maxIts){
	  assert_gt(maxIt_,0);
	  assert_gt(maxIt_z,0);
	  assert_gt(maxIt_s,0);
	  maxIt = maxIt_;
	  maxIt_z = maxItz;
	  maxIt_s = maxIts;
	}
	void setTol(const double tol_, 
				const double tolzf, 
				const double tolzg, 
				const double tolsf,
				const double tolsg){
	  assert_gt(tol_,0.0f);
	  assert_gt(tolzf,0.0f);
	  assert_gt(tolzg,0.0f);
	  assert_gt(tolsf,0.0f);
	  assert_gt(tolsg,0.0f);
	  tol = tol_;
	  tol_z_f = tolzf;
	  tol_z_g = tolzg;
	  tol_s_f = tolsf;
	  tol_s_g = tolsg;
	}
	void setMaxFrequency(const double maxFreq){
	  assert_gt(maxFreq,0.0f);
	  maxFrequency = maxFreq;
	}
	void setR(const int r_min, const int r_max){
	  assert_ge(r_max, r_min);
	  this->r_min = r_min;
	  this->r_max = r_max;
	}
	void setMtlOpt(const bool mtl){
	  mtlOpt = mtl;
	}
	bool optimizeMtl()const{
	  return mtlOpt;
	}

	virtual void updateConstraints(){

	  assert(EcZ);
	  assert(EcS);
	  EcZ->updateConstraints();
	  EkZ->updateConstraints();
	  EcS->updateConstraints();
	}

	bool solve();

	pMtlOptDM getDataModel(){
	  return dataModel;
	}
	pReducedRS_const getReducedRS(){
	  return reducedRS;
	}

	pCtrlForcesEnergyZ getEwZ(){
	  return EwZ;
	}
	pCtrlForcesEnergyLambdaDamping getEwLambdaDamping(){
	  return EwLambdaDamping;
	}
	pPosConEnergyZ getEcZ(){
	  return EcZ;
	}
	pPosConEnergyS getEcS(){
	  return EcS;
	}
	pBasisSelRegEnergy getEsS(){
	  return EsS;
	}
	pKeyframeConEnergyZ getEkZ(){
	  return EkZ;
	}
	pKeyframeConEnergyS getEkS(){
	  return EkS;
	}
	pCtrlForcesEnergyK_AD getEwK(){
	  return EwK;
	}
	pMechanicalEnergy getEmZ(){
	  return EmZ;
	}

	const alglib::minlbfgsstate& getLbfgsState()const{
	  return state;
	}

	inline void setFunEw( double f){
	  fun_Ew = f;
	}
	inline void setFunEc( double f){
	  fun_Ec = f;
	}
	inline void setFunEk( double f){
	  fun_Ek = f;
	}
	inline void setFunEs( double f){
	  fun_Es = f;
	}

	inline double funEw()const{
	  return fun_Ew;
	}
	inline double funEc()const{
	  return fun_Ec;
	}
	inline double funEk()const{
	  return fun_Ek;
	}
	inline double funEs()const{
	  return fun_Es;
	}

	inline double funValue()const{

	  const double f = fun_Ew + fun_Ec + fun_Es + fun_Ek;
	  assert_ge(f,0.0f);
	  return f;
	}

  protected:
	void computeInitialS();
	void hierSolve();
	bool solve(const int maxIt, const double tol);
	virtual void optimizeZ();
	virtual void optimizeK();
	virtual void expandS(const int rs){
	  dataModel->expandS(rs);
	}

	void optimizeS_lm();
	void optimizeS_lbfgs();
	void optimizeS_cg();
	void optimizeS_lbfgs_num();

	void optimizeLambda();
	void optimizeDamping();

	void optimizeLambdaDamping();
	void optimizeLambdaDampingDiag();
	void optimizeLambdaDamping_ipopt();
	void optimizeLambdaDamping_alglib();
	void optimizeLambdaDamping_mosec();

	void optimizeK_sym();
	void optimizeAtA_lbfgs();
	void optimizeAtA_lm();
	void optimizeAtA_ipopt();

	bool alglibErrorReport(const int code)const;

  protected:
	bool mtlOpt;
	int maxIt, maxIt_s, maxIt_z;
	double tol, tol_s_f, tol_s_g, tol_z_f, tol_z_g;
	double maxFrequency;
	int r_min;
	int r_max;

	pMtlOptDM dataModel;
	pReducedRS_const reducedRS;

	pCtrlForcesEnergyZ EwZ;
	pCtrlForcesEnergyLambdaDamping EwLambdaDamping;
	pCtrlForcesEnergyLambda EwLambda;
	pCtrlForcesEnergyDamping EwDamping;
	pPosConEnergyZ EcZ;
	pPosConEnergyS EcS;
	pKeyframeConEnergyZ EkZ;
	pKeyframeConEnergyS EkS;
	pBasisSelRegEnergy EsS;
	pCtrlForcesEnergyK_AD EwK;
	pMechanicalEnergy EmZ;

	double fun_Ew, fun_Ec, fun_Es, fun_Ek;

	alglib::minlbfgsstate state;
	alglib::minlbfgsreport rep;
  };
  typedef boost::shared_ptr<MtlOptSolver> pMtlOptSolver;

  /**
   * @class MtlOptSolverIpOpt solve the material opt using ipopt.
   * 
   */
  class MtlOptSolverIpOpt:public MtlOptSolver{

  public:
	MtlOptSolverIpOpt(pMtlOptDM dm,pReducedRS_const rs,const int mIt_=50,
					  const double tol_=1e-12,const double maxFreq = 4.0f):
	  MtlOptSolver(dm, rs, mIt_, tol_, maxFreq){}
	
  protected:
	virtual void optimizeZ();
  };
  
}//end of namespace

#endif /*_MTLOPTSOLVER_H_*/

