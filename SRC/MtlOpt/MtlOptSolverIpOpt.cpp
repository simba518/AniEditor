#include "MtlOptSolver.h"
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"
#include <MatrixIO.h>
#include <Timer.h>
using namespace Ipopt;
using namespace MTLOPT;

class OptSpaceTimeZ_NLP: public TNLP{

public:
  OptSpaceTimeZ_NLP(MtlOptSolver *mtl):mtlopt(mtl){}

  bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
					Index& nnz_h_lag, IndexStyleEnum& index_style){
	
	const int r = mtlopt->getDataModel()->reducedDim();
	const int Ts = mtlopt->getDataModel()->subFrames();
	n = r*Ts;
	m = 0;
	nnz_jac_g = 0;
	nnz_h_lag = (Ts-1)*3*r;/// @bug constraints not used
	index_style = C_STYLE;
	return true;
  }

  bool get_bounds_info(Index n, Number* x_l, Number* x_u,
					   Index m, Number* g_l, Number* g_u){
	checkDimension(n);
	assert_eq(m,0);
	for (int i = 0; i < n; ++i){
	  x_l[i] = -2e19;
	  x_u[i] = 2e19;
	}
	return true;
  }

  bool get_starting_point(Index n, bool init_x, Number* x,
						  bool init_z, Number* z_L, Number* z_U,
						  Index m, bool init_lambda,
						  Number* lambda){

	checkDimension(n);
	assert_eq(m,0);
	if(init_x){

	  const int r = mtlopt->getDataModel()->reducedDim();
	  const int Ts = mtlopt->getDataModel()->subFrames();
	  const MatrixXd &Z = mtlopt->getDataModel()->Z;

	  assert_ge(Z.size(),r*Ts);
	  memcpy( x, &Z(0,mtlopt->getDataModel()->T_begin),sizeof(double)*n );
	}
	return true;
  }

  bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){

	checkDimension(n);

	pCtrlForcesEnergyZ EwZ = mtlopt->getEwZ();
	pPosConEnergyZ EcZ = mtlopt->getEcZ();

	if(new_x){
	  mtlopt->getDataModel()->updateZ(x);
	  EcZ->updateYc(x);
	}

	const double fun1 = EwZ->fun();
	mtlopt->setFunEw(fun1);

	double fun = 0.0f;
	EcZ->funAddGrad(&fun,NULL);
	mtlopt->setFunEc(fun);

	obj_value = mtlopt->funValue();
	return true;
  }

  bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){

	checkDimension(n);

	pCtrlForcesEnergyZ EwZ = mtlopt->getEwZ();
	pPosConEnergyZ EcZ = mtlopt->getEcZ();

	if(new_x){
	  mtlopt->getDataModel()->updateZ(x);
	  EcZ->updateYc(x);
	}
	
	const VectorXd &subZV = Map<VectorXd>(const_cast<double*>(x),n);
	const SparseMatrix<double> &HLower = EwZ->getHessian();
	const VectorXd g = HLower.selfadjointView<Lower>()*subZV;
	assert_eq(n,g.size());
	memcpy(grad_f,&g[0],sizeof(double)*g.size());

	EcZ->funAddGrad(NULL, grad_f);
	return true;
  }

  bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
	checkDimension(n);
	assert_eq(m,0);
	return true;
  }

  bool eval_jac_g(Index n, const Number* x, bool new_x,
				  Index m, Index nele_jac, Index* iRow, Index *jCol,
				  Number* values){
	checkDimension(n);
	assert_eq(m,0);
	return true;
  }

  bool eval_h(Index n, const Number* x, bool new_x,
			  Number obj_factor, Index m, const Number* lambda,
			  bool new_lambda, Index nele_hess, Index* iRow,
			  Index* jCol, Number* values){
	
	checkDimension(n);
	assert_eq(m,0);
	
	/// @bug hess of constraints not used.
	const SparseMatrix<double> &HLower = mtlopt->getEwZ()->getHessian();
	assert_eq(nele_hess, HLower.nonZeros());
		
	if(NULL == values){

	  int p = 0;
	  for (int k=0; k<HLower.outerSize(); ++k){
		for (SparseMatrix<double>::InnerIterator it(HLower,k); it; ++it){
		  iRow[p] = it.row();
		  jCol[p] = it.col();
		  p++;
		}
	  }
	  assert_eq(p, nele_hess);
	}else{
	  memcpy(values,HLower.valuePtr(),nele_hess*sizeof(double));
	}

	return true;
  }

  void finalize_solution(SolverReturn status,
						 Index n, const Number* x, const Number* z_L, const Number* z_U,
						 Index m, const Number* g, const Number* lambda,
						 Number obj_value,
						 const IpoptData* ip_data,
						 IpoptCalculatedQuantities* ip_cq){
	
	checkDimension(n);
	assert_eq(m,0);
	assert_le(obj_value, mtlopt->funValue());
	DEBUG_LOG("inner-it-min(fun(Z)): "<< obj_value);
  }

protected:
  void checkDimension(const int n)const{

	assert_eq(sizeof(Number),sizeof(double));
	assert_eq(sizeof(Index),sizeof(int));
	assert(mtlopt);
	assert(mtlopt->getDataModel());
	const int r = mtlopt->getDataModel()->reducedDim();
	const int Ts = mtlopt->getDataModel()->subFrames();	
	assert_eq_ext(r*Ts,n," r = "<<r<<", Ts = "<<Ts);
	assert_gt(n,0);
  }

private:
  OptSpaceTimeZ_NLP(const OptSpaceTimeZ_NLP&){}

  OptSpaceTimeZ_NLP& operator=(const OptSpaceTimeZ_NLP&){}

private:
  MtlOptSolver *mtlopt;

};

void MtlOptSolverIpOpt::optimizeZ(){

  EwZ->reset(true);
  EcZ->reset();

  SmartPtr<TNLP> mynlp = new OptSpaceTimeZ_NLP(this);
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 1e-3);
  app->Options()->SetIntegerValue("max_iter",30);
  app->Options()->SetStringValue("mu_strategy", "adaptive");

  ApplicationReturnStatus status;
  status = app->Initialize();
  ERROR_LOG_COND("Error during initialization!",status == Solve_Succeeded);

  status = app->OptimizeTNLP(mynlp);
  ERROR_LOG_COND("The problem FAILED!",status == Solve_Succeeded);
}

class OptSpaceTimeAtA_NLP: public TNLP{

public:
  OptSpaceTimeAtA_NLP(MtlOptSolver *mtl):mtlopt(mtl){}

  bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
					Index& nnz_h_lag, IndexStyleEnum& index_style){
	
	const int r = mtlopt->getDataModel()->reducedDim();
	n = r*r;
	m = 0;
	nnz_jac_g = 0;
	nnz_h_lag = (n+1)*n/2;
	index_style = C_STYLE;
	return true;
  }

  bool get_bounds_info(Index n, Number* x_l, Number* x_u,
					   Index m, Number* g_l, Number* g_u){
	checkDimension(n);
	assert_eq(m,0);
	for (int i = 0; i < n; ++i){
	  x_l[i] = -2e19;
	  x_u[i] = 2e19;
	}
	return true;
  }

  bool get_starting_point(Index n, bool init_x, Number* x,
						  bool init_z, Number* z_L, Number* z_U,
						  Index m, bool init_lambda,
						  Number* lambda){

	checkDimension(n);
	assert_eq(m,0);
	if(init_x){
	  const MatrixXd &K = mtlopt->getDataModel()->K;
	  memcpy( x, &K(0,0), sizeof(double)*n );
	  const int r = K.rows();
	  for (int i = 0; i < r; ++i){
		for (int j = 0; j < r; ++j){
		  if ( i == j){
			assert_ge_ext(x[i*r+j],0.0f,"i="<<i);
			x[j*r+i] = sqrt(x[j*r+i]);
		  }else{
			assert_eq_ext(x[i*r+j],0.0f,"i="<<i<<",j="<<j);
		  }
		}
	  }
	}
	return true;
  }

  bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){

	checkDimension(n);
	pCtrlForcesEnergyK_AD EwK = mtlopt->getEwK();
	const double fun = EwK->fun(x, n);
	mtlopt->setFunEw(fun);
	obj_value = mtlopt->funValue();
	return true;
  }

  bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){

	checkDimension(n);
	pCtrlForcesEnergyK_AD EwK = mtlopt->getEwK();
	EwK->grad(x,n,grad_f);
	return true;
  }

  bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
	checkDimension(n);
	assert_eq(m,0);
	return true;
  }

  bool eval_jac_g(Index n, const Number* x, bool new_x,
				  Index m, Index nele_jac, Index* iRow, Index *jCol,
				  Number* values){
	checkDimension(n);
	assert_eq(m,0);
	return true;
  }

  bool eval_h(Index n, const Number* x, bool new_x,
			  Number obj_factor, Index m, const Number* lambda,
			  bool new_lambda, Index nele_hess, Index* iRow,
			  Index* jCol, Number* values){
	
	checkDimension(n);
	assert_eq(m,0);
	assert_eq(nele_hess, (n+1)*n/2);
		
	if(NULL == values){
	  int p = 0;
	  for (int i = 0; i < n; ++i){
		for (int j = 0; j <= i; ++j){
		  iRow[p] = j;
		  jCol[p] = i;
		  p ++;
		}
	  }
	  assert_eq(p, nele_hess);
	}else{
	  static MatrixXd hess;
	  pCtrlForcesEnergyK_AD EwK = mtlopt->getEwK();
	  hess.resize(n,n);
	  EwK->hessian(x,n,&hess(0,0));
	  int p = 0;
	  for (int i = 0; i < n; ++i){
		for (int j = 0; j <= i; ++j){
		  values[p] = hess(i,j);
		  p ++;
		}
	  }
	}

	return true;
  }

  void finalize_solution(SolverReturn status,
						 Index n, const Number* x, const Number* z_L, const Number* z_U,
						 Index m, const Number* g, const Number* lambda,
						 Number obj_value,
						 const IpoptData* ip_data,
						 IpoptCalculatedQuantities* ip_cq){
	
	checkDimension(n);
	assert_eq(m,0);
	assert_le(obj_value, mtlopt->funValue());
	DEBUG_LOG("inner-it-min(fun(Z)): "<< obj_value);

	// compute K = At*A, decompose K = Ut \Lambda U
	pMtlOptDM dataModel = mtlopt->getDataModel();
	const int r = dataModel->reducedDim();
	const MatrixXd &A = Map<MatrixXd>(const_cast<double*>(x), r, r);
	dataModel->K = A.transpose()*A;
	dataModel->decomposeK();
  }

protected:
  void checkDimension(const int n)const{

	assert_eq(sizeof(Number),sizeof(double));
	assert_eq(sizeof(Index),sizeof(int));
	assert(mtlopt);
	assert(mtlopt->getDataModel());
	const int r = mtlopt->getDataModel()->reducedDim();
	assert_eq_ext(r*r,n," r = "<<r);
	assert_gt(n,0);
  }

private:
  OptSpaceTimeAtA_NLP(const OptSpaceTimeAtA_NLP&){}

  OptSpaceTimeAtA_NLP& operator=(const OptSpaceTimeAtA_NLP&){}

private:
  MtlOptSolver *mtlopt;

};

void MtlOptSolver::optimizeAtA_ipopt(){

  FUNC_TIMER();

  EwK->reset();

  SmartPtr<TNLP> mynlp = new OptSpaceTimeAtA_NLP(this);
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 1e-3);
  app->Options()->SetIntegerValue("max_iter",30);
  app->Options()->SetStringValue("mu_strategy", "adaptive");

  ApplicationReturnStatus status;
  status = app->Initialize();
  ERROR_LOG_COND("Error during initialization!",status == Solve_Succeeded);

  status = app->OptimizeTNLP(mynlp);
  ERROR_LOG_COND("The problem FAILED!",status == Solve_Succeeded);

}

class OptSpaceTimeLambdaDamping_NLP:public TNLP{

public:
  OptSpaceTimeLambdaDamping_NLP(MtlOptSolver *mtl):mtlopt(mtl){}

  bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
					Index& nnz_h_lag, IndexStyleEnum& index_style){
	
	const int r = mtlopt->getDataModel()->reducedDim();
	n = r*2;
	m = 0;
	nnz_jac_g = 0;
	nnz_h_lag = r*3;
	index_style = C_STYLE;
	pCtrlForcesEnergyLambdaDamping EwLambdaDamping = mtlopt->getEwLambdaDamping();
	EwLambdaDamping->hessianLower(iRow,iCol,value);
	assert_eq(value.size(),iRow.size());
	assert_eq(value.size(),iCol.size());
	assert_eq(value.size(),nnz_h_lag);
	return true;
  }

  bool get_bounds_info(Index n, Number* x_l, Number* x_u,
					   Index m, Number* g_l, Number* g_u){
	checkDimension(n);
	assert_eq(m,0);
	for (int i = 0; i < n/2; ++i){
	  x_l[i] = 0;
	  x_u[i] = 4.0f;
	  // x_u[i] = 2e19;
	}
	for (int i = n/2; i < n; ++i){
	  x_l[i] = 0;
	  x_u[i] = 2e19;
	}
	return true;
  }

  bool get_starting_point(Index n, bool init_x, Number* x,
						  bool init_z, Number* z_L, Number* z_U,
						  Index m, bool init_lambda,
						  Number* lambda){

	checkDimension(n);
	assert_eq(m,0);
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	memcpy(x,&(mtlopt->getDataModel()->Lambda[0]),n/2*sizeof(double));
	memcpy(&x[n/2],&(mtlopt->getDataModel()->Damping[0]),n/2*sizeof(double));

	return true;
  }

  bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){

	checkDimension(n);
	pCtrlForcesEnergyLambdaDamping EwLambdaDamping = mtlopt->getEwLambdaDamping();
	const double fun = EwLambdaDamping->fun(x);

	mtlopt->setFunEw(fun);
	obj_value = mtlopt->funValue();
	return true;
  }

  bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){

	checkDimension(n);
	pCtrlForcesEnergyLambdaDamping EwLambdaDamping = mtlopt->getEwLambdaDamping();
	EwLambdaDamping->grad(x,grad_f);
	return true;
  }

  bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
	checkDimension(n);
	assert_eq(m,0);
	return true;
  }

  bool eval_jac_g(Index n, const Number* x, bool new_x,
				  Index m, Index nele_jac, Index* iRow, Index *jCol,
				  Number* values){
	checkDimension(n);
	assert_eq(m,0);
	return true;
  }

  bool eval_h(Index n, const Number* x, bool new_x,
			  Number obj_factor, Index m, const Number* lambda,
			  bool new_lambda, Index nele_hess, Index* iRow,
			  Index* jCol, Number* values){
	
	checkDimension(n);
	assert_eq(m,0);
		
	if(NULL == values){
	  assert_eq(this->iRow.size(),nele_hess);
	  assert_eq(this->iCol.size(),nele_hess);
	  memcpy(iRow, &(this->iRow[0]), nele_hess*sizeof(int));
	  memcpy(jCol, &(this->iCol[0]), nele_hess*sizeof(int));
	}else{
	  assert_eq(nele_hess, value.size());
	  memcpy(values, &(this->value[0]), nele_hess*sizeof(double));
	}

	return true;
  }

  void finalize_solution(SolverReturn status,
						 Index n, const Number* x, const Number* z_L, const Number* z_U,
						 Index m, const Number* g, const Number* lambda,
						 Number obj_value,
						 const IpoptData* ip_data,
						 IpoptCalculatedQuantities* ip_cq){
	
	checkDimension(n);
	assert_eq(m,0);
	assert_le(obj_value, mtlopt->funValue());
	DEBUG_LOG("inner-it-min(fun(lambda,damping)): "<< obj_value);

	pMtlOptDM dataModel = mtlopt->getDataModel();
	const int r = dataModel->reducedDim();
	memcpy( &(dataModel->Lambda[0]), x, r*sizeof(double));
	memcpy( &(dataModel->Damping[0]), &x[r], r*sizeof(double));
  }

protected:
  void checkDimension(const int n)const{

	assert_eq(sizeof(Number),sizeof(double));
	assert_eq(sizeof(Index),sizeof(int));
	assert(mtlopt);
	assert(mtlopt->getDataModel());
	const int r = mtlopt->getDataModel()->reducedDim();
	assert_eq_ext(r*2,n," r = "<<r);
	assert_gt(n,0);
  }

private:
  OptSpaceTimeLambdaDamping_NLP(const OptSpaceTimeLambdaDamping_NLP&){}

  OptSpaceTimeLambdaDamping_NLP& operator=(const OptSpaceTimeLambdaDamping_NLP&){}

private:
  MtlOptSolver *mtlopt;
  vector<int> iRow;
  vector<int> iCol;
  vector<double> value;
};

void MtlOptSolver::optimizeLambdaDamping_ipopt(){
  
  FUNC_TIMER();

  EwLambdaDamping->reset();

  SmartPtr<TNLP> mynlp = new OptSpaceTimeLambdaDamping_NLP(this);
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", tol_z_g);
  app->Options()->SetIntegerValue("max_iter",maxIt_z);
  app->Options()->SetStringValue("mu_strategy", "adaptive");

  ApplicationReturnStatus status;
  status = app->Initialize();
  ERROR_LOG_COND("Error during initialization!",status == Solve_Succeeded);

  status = app->OptimizeTNLP(mynlp);
  ERROR_LOG_COND("The problem FAILED!",status == Solve_Succeeded);
}
