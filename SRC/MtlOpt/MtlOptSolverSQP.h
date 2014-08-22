#ifndef _MTLOPTSOLVERSQP_H_
#define _MTLOPTSOLVERSQP_H_

#include "MtlOptSolver.h"
#include "SQPext.h"
#include <Timer.h>

namespace MTLOPT{

  class MtlOptSqpFunction: public SQPfun{

  public: 
	MtlOptSqpFunction(MtlOptSolver *mtl):mtlopt(mtl){
	  hessIsUpdated = false;
	}

	void updateHessian(){

	  FUNC_TIMER();
	  
	  pCtrlForcesEnergyZ EwZ = mtlopt->getEwZ();
	  const SparseMatrix<double> &HLower = EwZ->getHessian();
	  Hfull = HLower.selfadjointView<Lower>();

	  if (mtlopt->getEmZ()){
		pMechanicalEnergy EmZ = mtlopt->getEmZ();
		const SparseMatrix<double> &HLowerM = EmZ->getHessian();
		assert_eq(HLowerM.rows(), HLower.rows());
		assert_eq(HLowerM.cols(), HLower.cols());
		Hfull += HLowerM.selfadjointView<Lower>();
	  }

	  if (mtlopt->getDataModel()->keyframes.size()>0){
		Hfull += mtlopt->getEkZ()->getHess();
	  }

	  Hfull.makeCompressed();
	  H_JtJ = Hfull;
	  hessIsUpdated = true;
	}

	size_t dim(void) const {
	  const int r = mtlopt->getDataModel()->reducedDim();
	  const int Ts = mtlopt->getDataModel()->subFrames();
	  return r*Ts;
	}

	int prepare_for_new_x(const double *x){

	  FUNC_TIMER();

	  // mtlopt->getDataModel()->updateZ(x);
	  const MatrixXd &Z = mtlopt->getDataModel()->Z;
	  assert(&Z(0,mtlopt->getDataModel()->T_begin) == x);
	  mtlopt->getEcZ()->updateYc(x);
	  return 0;
	}

	int val(const double *x, double &v) const{

	  FUNC_TIMER();

	  pCtrlForcesEnergyZ EwZ = mtlopt->getEwZ();
	  pMechanicalEnergy EmZ = mtlopt->getEmZ();
	  pPosConEnergyZ EcZ = mtlopt->getEcZ();
	  pKeyframeConEnergyZ EkZ = mtlopt->getEkZ();

	  double fun1 = EwZ->fun();
	  if (EmZ){
		fun1 += EmZ->fun();
	  }
	  mtlopt->setFunEw(fun1);

	  double fun = 0.0f;
	  if (mtlopt->getDataModel()->conFrames.size() > 0){
		EcZ->funAddGrad(&fun,NULL);
	  }
	  mtlopt->setFunEc(fun);
	  mtlopt->setFunEk(EkZ->fun(x));

	  v = mtlopt->funValue();

	  return 0;
	}

	int gra(const double *x, double *grad_f) {

	  TRACE_FUN();
	  FUNC_TIMER();

	  pCtrlForcesEnergyZ EwZ = mtlopt->getEwZ();
	  const SparseMatrix<double> &HLower = EwZ->getHessian();
	  const VectorXd g = HLower.selfadjointView<Lower>()*(Map<VectorXd>(const_cast<double*>(x),dim()));
	  memcpy(grad_f,&g[0],sizeof(double)*g.size());
	  if (mtlopt->getDataModel()->conFrames.size() > 0){
		mtlopt->getEcZ()->funAddGrad(NULL, grad_f);
	  }
	  mtlopt->getEkZ()->gradAdd(x,grad_f);

	  return 0;
	}

	int hes(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx,double alpha) {

	  TRACE_FUN();
	  FUNC_TIMER();
	  const int r =  varHess(x,nnz,format,h,ptr,idx,alpha);
	  return r;
	}

  protected:
	int constHess(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx,double alpha);
	int varHess(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx,double alpha);

  private:
	MtlOptSolver *mtlopt;
	SparseMatrix<double> Hfull;
	SparseMatrix<double> JtJfull;
	SparseMatrix<double> H_JtJ;
	bool hessIsUpdated;
  };
  typedef boost::shared_ptr<MtlOptSqpFunction> pMtlOptSqpFunction;

  /**
   * @class MtlOptSolverSQP solve the material opt using SQP from zju library.
   * 
   */
  class MtlOptSolverSQP:public MtlOptSolver{

  public:
	MtlOptSolverSQP(pMtlOptDM dm,pReducedRS_const rs,const int mIt_=50,
					const double tol_=1e-12,const double maxFreq = 4);
	
	void updateConstraints(){
	  MtlOptSolver::updateConstraints();
	  sqp_solver.clear_hes();
	}

  protected:
	void optimizeZ();
	void expandS(const int rs){
	  MtlOptSolver::expandS(rs);
	  sqp_solver.clear_hes();
	}

  private:
	pMtlOptSqpFunction sqp_fun;
	SQPext sqp_solver;
	boost::property_tree::ptree solver_setting;
  };
  
}//end of namespace

#endif /*_MTLOPTSOLVERSQP_H_*/
