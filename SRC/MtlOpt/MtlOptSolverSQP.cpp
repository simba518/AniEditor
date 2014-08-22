#include "MtlOptSolverSQP.h"
using namespace MTLOPT;

MtlOptSolverSQP::MtlOptSolverSQP(pMtlOptDM dm,pReducedRS_const rs,const int mIt_,
								 const double tol_,const double maxFreq ):
  MtlOptSolver(dm, rs, mIt_, tol_, maxFreq){

  sqp_fun = pMtlOptSqpFunction(new MtlOptSqpFunction(this));
}

void MtlOptSolverSQP::optimizeZ(){

  FUNC_TIMER();

  solver_setting.put<int>("iter.value",maxIt_z);
  solver_setting.put<string>("linear_solver/type.value","direct");
  solver_setting.put<string>("linear_solver/name.value","cholmod");
  solver_setting.put<string>("customerlog.value","inner-it-fun(Z) ");

  solver_setting.put<double>("epsg.value", tol_z_g);
  solver_setting.put<double>("epsf.value", tol_z_f);
  solver_setting.put<double>("epsx.value", 0.0f);

  // prepare
  EwZ->reset(true);
  EcZ->reset();
  if (EmZ) EmZ->reset(dataModel->mechanic_con_frames);

  sqp_fun->updateHessian();

  // solve
  double *x = &(dataModel->Z(0,dataModel->T_begin));
  sqp_solver.set_f(sqp_fun.get());
  sqp_solver.solve(x,solver_setting);

  // get results
  EcZ->updateYc(x);
  const double fun1 = EwZ->fun();
  setFunEw(fun1);
  double fun = 0.0f;
  EcZ->funAddGrad(&fun,NULL);
  setFunEc(fun);
  setFunEk(EkZ->fun(x));
  DEBUG_LOG("norm(Z): "<<dataModel->Z.norm()<<endl);
}

int MtlOptSqpFunction::constHess(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx,double alpha){

  format = 1;

  if(h == 0 && ptr == 0 && idx == 0) { // query nnz
	assert(hessIsUpdated);
	assert_gt(Hfull.nonZeros(),0);
	nnz = Hfull.nonZeros();
	return 0;
  }

  if(h == 0 && ptr != 0 && idx != 0) { // query pattern

	assert_eq(nnz, Hfull.nonZeros());
	assert_gt(nnz,0);
	assert_eq(sizeof(int32_t),sizeof(int));
	memcpy(ptr, Hfull.outerIndexPtr(), (Hfull.cols()+1)*sizeof(int));
	memcpy(idx, Hfull.innerIndexPtr(), nnz*sizeof(int));
	return 0;
  }

  if(h != 0 && ptr != 0 && idx != 0) { // accumulate

	assert_gt(nnz,0);
	assert_eq(Hfull.nonZeros(), nnz);
	if (hessIsUpdated){
	  memcpy(h, Hfull.valuePtr(), nnz*sizeof(double));
	  hessIsUpdated = false;
	  return 0;
	}
	return -1;
  }

  return 0;
}

int MtlOptSqpFunction::varHess(const double *x, size_t &nnz, size_t &format, double *h, int32_t *ptr, int32_t *idx,double alpha){

  format = 1;

  if(h == 0 && ptr == 0 && idx == 0) { // query nnz	
	assert(hessIsUpdated);
	assert_gt(Hfull.nonZeros(),0);	
	if (mtlopt->getDataModel()->conFrames.size() > 0){
	  pPosConEnergyZ EcZ = mtlopt->getEcZ();
	  EcZ->initJtJfullStruct(JtJfull);
	  H_JtJ = Hfull + JtJfull;
	}else{
	  assert_eq(H_JtJ.nonZeros(), Hfull.nonZeros());
	}

	nnz = H_JtJ.nonZeros();
	return 0;
  }

  if(h == 0 && ptr != 0 && idx != 0) { // query pattern

	assert_eq(nnz, H_JtJ.nonZeros());
	assert_gt(nnz,0);
	assert_eq(sizeof(int32_t),sizeof(int));
	memcpy(ptr, H_JtJ.outerIndexPtr(), (H_JtJ.cols()+1)*sizeof(int));
	memcpy(idx, H_JtJ.innerIndexPtr(), nnz*sizeof(int));
	return 0;
  }

  if(h != 0 && ptr != 0 && idx != 0) { // accumulate
	assert_gt(nnz,0);
	if (mtlopt->getDataModel()->conFrames.size() > 0){
	  static MatrixXd J;
	  pPosConEnergyZ EcZ = mtlopt->getEcZ();
	  EcZ->jacbian(J);
	  EcZ->assembleJtJfull(J,JtJfull);
	  H_JtJ = JtJfull + Hfull;
	  hessIsUpdated = false;
	}
	assert_eq(H_JtJ.nonZeros(), nnz);
	memcpy(h, H_JtJ.valuePtr(), nnz*sizeof(double));
	return 0;
  }

  return 0;
}
