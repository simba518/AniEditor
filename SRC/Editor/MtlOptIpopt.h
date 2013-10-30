#ifndef _MTLOPTIPOPT_H_
#define _MTLOPTIPOPT_H_

#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <MtlOptEnergy.h>
using namespace LSW_ANI_EDITOR;
using namespace Ipopt;

class FunGradNLP: public TNLP{
	
public:
  FunGradNLP(pBaseFunGrad energy):_energy(energy){}
  ~FunGradNLP(){}

  bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
					Index& nnz_h_lag, IndexStyleEnum& index_style){
	assert(_energy);
	n = _energy->dim();
	m = 0;
	nnz_jac_g = 0;
	nnz_h_lag = 0;
	index_style = C_STYLE;
	return true;
  }

  bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u){
	assert(_energy);
	assert_eq(_energy->dim(),n);
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
	assert(_energy);
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);
	assert_eq(_energy->dim(),n);
	_energy->init(x,n);
	return true;
  }

  bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value){
	assert(_energy);
	assert_eq(_energy->dim(),n);
	obj_value = _energy->fun(x);
	return true;
  }

  bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f){
	assert(_energy);
	assert_eq(_energy->dim(),n);
	_energy->grad(x,grad_f);
	return true;
  }

  bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g){
	return true;
  }

  bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values){
	return true;
  }

  void finalize_solution(SolverReturn status,
						 Index n,const Number* x,const Number* z_L,const Number* z_U,
						 Index m, const Number* g, const Number* lambda,
						 Number obj_value,
						 const IpoptData* ip_data,
						 IpoptCalculatedQuantities* ip_cq){
	assert(_energy);
	assert_eq(_energy->dim(),n);
	_energy->setRlst(x,obj_value);
  }

private:
  FunGradNLP(const FunGradNLP&);
  FunGradNLP& operator=(const FunGradNLP&);
	
private:
  pBaseFunGrad _energy;

};

typedef boost::shared_ptr<FunGradNLP> pFunGradNLP;

class NoConIpoptSolver{
  
public:
  NoConIpoptSolver(pBaseFunGrad energy){
	_mynlp = new FunGradNLP(energy);
	_app = IpoptApplicationFactory();
	_app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	_app->Options()->SetIntegerValue("print_level",0);
  }
  void setTol(const double tol){
	assert_gt(tol,0.0f);
	_app->Options()->SetNumericValue("tol",tol);
  }
  void setMaxIt(const int maxIt){
	assert_gt(maxIt,0);
	_app->Options()->SetIntegerValue("max_iter",maxIt);
  }
  void setPrintLevel(const int level){
	_app->Options()->SetIntegerValue("print_level",level);
  }
  bool initialize(){
	_status = _app->Initialize();
	return (_status == Solve_Succeeded);
  }

  bool solve(){
	_status = _app->OptimizeTNLP(_mynlp);
	return (_status == Solve_Succeeded);
  }

  const ApplicationReturnStatus &getStatus()const{
	return _status;
  }
  
private:
  SmartPtr<TNLP> _mynlp;
  SmartPtr<IpoptApplication> _app;
  ApplicationReturnStatus _status;
};

typedef boost::shared_ptr<NoConIpoptSolver> pNoConIpoptSolver;

#endif /* _MTLOPTIPOPT_H_ */
