#include <alglib/stdafx.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <alglib/optimization.h>
#include <Timer.h>
#include <IOTools.h>
#include <MatrixIO.h>
#include "MtlOptSolver.h"
#include "AuxTools.h"
using namespace alglib;
using namespace MTLOPT;

void FunGradOptZ(const real_1d_array &subZ, double &func, real_1d_array &grad,void*ptr) {

  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pMtlOptDM dataModel = mtl->getDataModel();
  pCtrlForcesEnergyZ EwZ = mtl->getEwZ();
  pPosConEnergyZ EcZ = mtl->getEcZ();

  const int r = dataModel->reducedDim();
  const int Ts = dataModel->subFrames();

  assert_eq(subZ.length(),r*Ts);
  assert_gt(subZ.length(),0);
  
  { // Ew(z)

	dataModel->updateZ(&subZ[0]);
	const double fun = EwZ->fun();
	mtl->setFunEw(fun);
  
	const VectorXd &subZV = Map<VectorXd>(const_cast<double*>(&subZ[0]),subZ.length());
	const SparseMatrix<double> &HLower = EwZ->getHessian();
	const VectorXd g = HLower.selfadjointView<Lower>()*subZV;

	assert_eq(grad.length(),g.size());
	memcpy(&grad[0],&g[0],sizeof(double)*g.size());
  }

  { // Ec(z)
	EcZ->updateYc(&subZ[0]);
	double fun = 0.0f;
	EcZ->funAddGrad(&fun, &grad[0]);
	mtl->setFunEc(fun);
  }

  { // Ek(z)
	assert(false);
  }

  func = mtl->funValue();
  DEBUG_LOG("inner-it-fun(Z): "<<func);
}

void FunOptS(const real_1d_array &S, double &func, void *ptr){

  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pPosConEnergyS EcS = mtl->getEcS();
  pBasisSelRegEnergy EsS = mtl->getEsS();

  { // Es(S)
  	const int rw = mtl->getDataModel()->S.rows();
  	assert_gt(rw,0);
  	const int rs = S.length()/rw;
  	assert_eq(rw*rs, S.length());
  	const double fun = EsS->fun(&S[0], rw, rs);
  	mtl->setFunEs(fun);
  }

  { // Ec(S)
  	EcS->updateYc(&S[0]);
  	double fun = 0.0f;
  	EcS->funAddGrad(&fun, NULL);
	mtl->setFunEc(fun);
  }

  { // Ek(S)
	assert(false);
  }

  func = mtl->funValue();
  DEBUG_LOG("inner-it-fun(S): "<<func);
}

void FunGradOptS(const real_1d_array &S, double &func, real_1d_array &grad,void*ptr) {

  FUNC_TIMER();

  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pPosConEnergyS EcS = mtl->getEcS();
  pKeyframeConEnergyS EkS = mtl->getEkS();
  pBasisSelRegEnergy EsS = mtl->getEsS();

  { // Es(S)
  	const int rw = mtl->getDataModel()->S.rows();
  	assert_gt(rw,0);
  	const int rs = S.length()/rw;
  	assert_eq(rw*rs, S.length());

  	const double fun = EsS->fun(&S[0], rw, rs);
  	EsS->grad(&S[0],&grad[0], rw, rs);

	mtl->setFunEs(fun);
  }

  { // Ec(S)
	double fun = 0.0f;
	if (mtl->getDataModel()->conFrames.size() > 0){
	  EcS->updateYc(&S[0]);
	  EcS->funAddGrad(&fun, &grad[0]);
	}
	mtl->setFunEc(fun);
  }

  { // Ek(S)
	const double fun = EkS->fun(&S[0]);
	EkS->gradAdd(&S[0],&grad[0]);
	mtl->setFunEk(fun);
  }

  func = mtl->funValue();

  {// report the iterations
	static int iterations = 1;
	if(0 == mtl->getLbfgsState().c_ptr()->repiterationscount && iterations != 0){
	  iterations = 0;
	  // DEBUG_LOG("inner-it-fun-grad(S): "<<func);
	}else if (mtl->getLbfgsState().c_ptr()->repiterationscount > iterations){
	  DEBUG_LOG("inner-it-fun-grad(S): "<<func);
	  iterations++;
	  DEBUG_LOG("iterations counter: "<<iterations);
	}
  }

}

void FunGradHessOptS(const real_1d_array &S, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr){
  
  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pPosConEnergyS EcS = mtl->getEcS();
  pBasisSelRegEnergy EsS = mtl->getEsS();

  VectorXd Hs;
  { // Es(S)
  	const int rw = mtl->getDataModel()->S.rows();
  	assert_gt(rw,0);
  	const int rs = S.length()/rw;
  	assert_eq(rw*rs, S.length());

  	const double fun = EsS->fun(&S[0], rw, rs);
  	EsS->grad(&S[0],&grad[0], rw, rs);
	mtl->setFunEs(fun);
	EsS->hessian(Hs,rw,rs);
  }

  { // Ec(S)
  	EcS->updateYc(&S[0]);
  	double fun = 0.0f;
  	EcS->funAddGrad(&fun, &grad[0]);
	mtl->setFunEc(fun);

	MatrixXd J;
	EcS->jacbian(J);
	MatrixXd JtJ = J.transpose()*J;
	assert_eq(JtJ.rows(), Hs.size());
	JtJ += Hs.asDiagonal();
	assert_eq(hess.rows(), JtJ.rows());
	assert_eq(hess.cols(), JtJ.cols());
	memcpy( &hess(0,0), &JtJ(0,0), JtJ.size()*sizeof(double) );
  }

  { // Ek(S)
	assert(false);
  }

  func = mtl->funValue();

  DEBUG_LOG("inner-it-fun-grad-hess(S): "<<func);
}

void FunOptK(const real_1d_array &X, double &func, void *ptr){
  
  FUNC_TIMER();

  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pCtrlForcesEnergyK_AD EwK = mtl->getEwK();

  { // Ew(K(x))
	const double r = mtl->getDataModel()->reducedDim();
	assert_eq(X.length(), r*r);
	assert_gt(r,0);

  	const double fun = EwK->fun(&X[0],X.length());
	mtl->setFunEw(fun);
  }
  func = mtl->funValue();

  DEBUG_LOG("inner-it-fun(Lambda): "<<func);
}

void FunGradOptK(const real_1d_array &X, double &func, real_1d_array &grad,void*ptr) {

  FUNC_TIMER();

  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pCtrlForcesEnergyK_AD EwK = mtl->getEwK();

  { // Ew(K(x))
	const double r = mtl->getDataModel()->reducedDim();
	assert_eq(X.length(), r*r);
	assert_eq(grad.length(), r*r);
	assert_gt(r,0);

  	const double fun = EwK->fun(&X[0],X.length());
  	EwK->grad(&X[0],X.length(),&grad[0]);
	mtl->setFunEw(fun);
  }
  func = mtl->funValue();

  DEBUG_LOG("inner-it-fun-grad(Lambda): "<<func);
  DEBUG_LOG("inner-it-grad(Lambda): "<<UTILITY::norm2(grad,grad.length()));
}

void FunGradHessOptK(const real_1d_array &X, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr){
 
  FUNC_TIMER();

  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pCtrlForcesEnergyK_AD EwK = mtl->getEwK();

  { // Ew(K(x))
	const double r = mtl->getDataModel()->reducedDim();
	assert_eq(X.length(), r*r);
	assert_eq(grad.length(), r*r);
	assert_gt(r,0);
	assert_eq(hess.cols(), r*r);
	assert_eq(hess.cols(), r*r);

  	const double fun = EwK->fun(&X[0], X.length());
	mtl->setFunEw(fun);

  	EwK->grad(&X[0], X.length(), &grad[0]);
	EwK->hessian(&X[0], X.length(), &hess(0,0));
  }
  func = mtl->funValue();

  DEBUG_LOG("inner-it-fun-grad(Lambda): "<<func);
  DEBUG_LOG("inner-it-grad(Lambda): "<<UTILITY::norm2(grad,grad.length()));
}

void FunGradOptLambdaDamping(const real_1d_array &X, double &func, real_1d_array &grad, void *ptr){

  assert(ptr);
  MtlOptSolver *mtl = (MtlOptSolver *)ptr;
  pCtrlForcesEnergyLambdaDamping EwLambdaDamping = mtl->getEwLambdaDamping();

  { // Ew(Lambda,Damping)
	const double r = mtl->getDataModel()->reducedDim();
	assert_eq(X.length(), r*2);
	assert_eq(grad.length(), r*2);
	assert_gt(r,0);
  	const double fun = EwLambdaDamping->fun(&X[0]);
  	EwLambdaDamping->grad(&X[0],&grad[0]);
	mtl->setFunEw(fun);
  }
  func = mtl->funValue();

  DEBUG_LOG("inner-it-fun-grad(LambdaDamping): "<<func);
  DEBUG_LOG("inner-it-grad(LambdaDamping): "<<UTILITY::norm2(grad,grad.length()));
  
}

MtlOptSolver::MtlOptSolver(pMtlOptDM dm,pReducedRS_const rs,const int mIt_,
						   const double tol_,const double maxFreq):
  maxIt(mIt_),tol(tol_),maxFrequency(maxFreq),dataModel(dm),reducedRS(rs){

  assert_gt(maxIt,0);
  assert_gt(tol,0.0f);
  assert_gt(maxFrequency,0.0f);

  EwZ = pCtrlForcesEnergyZ(new CtrlForcesEnergyZ(dm));
  EwLambda = pCtrlForcesEnergyLambda(new CtrlForcesEnergyLambda(dm));
  EwDamping = pCtrlForcesEnergyDamping(new CtrlForcesEnergyDamping(dm));
  EwLambdaDamping = pCtrlForcesEnergyLambdaDamping(new CtrlForcesEnergyLambdaDamping(dm));

  EcZ = pPosConEnergyZ(new PosConEnergyZ(dm,reducedRS));
  EcS = pPosConEnergyS(new PosConEnergyS(dm,reducedRS));
  EkZ = pKeyframeConEnergyZ(new KeyframeConEnergyZ(dm));
  EkS = pKeyframeConEnergyS(new KeyframeConEnergyS(dm));
  EsS = pBasisSelRegEnergy(new BasisNormalizeColEnergy());
  // EsS = pBasisSelRegEnergy(new BasisNormalizeEnergy());
  // EsS = pBasisSelRegEnergy(new BasisSelRegEnergy());
  EwK = pCtrlForcesEnergyK_AD(new CtrlForcesEnergyK_AD(dm));

  mtlOpt = true;
  maxIt_s = 30;
  maxIt_z = 30;
  tol_s_f = 1e-6;
  tol_s_g = 1e-6;
  tol_z_f = 1e-6;
  tol_z_g = 1e-6;

  fun_Ew = 0.0f;
  fun_Ec = 0.0f;
  fun_Ek = 0.0f;
  fun_Es = 0.0f;
  r_min = 2;
  r_max = 2;
}

void MtlOptSolver::computeInitialS(){

  const double p = EsS->getPenalty();
  EsS = pBasisSelRegEnergy(new BasisSelRegEnergy());
  EsS->setPenalty(p);
  optimizeZ();
  solve(1,tol);
  JacobiSVD<MatrixXd> svd(dataModel->S, ComputeThinU | ComputeThinV);
  DEBUG_LOG("Its singular values are:" << endl << svd.singularValues() << endl);
  assert_eq(svd.matrixU().rows(),dataModel->S.rows());
  assert_eq(svd.matrixU().cols(),dataModel->S.cols());
  dataModel->S = svd.matrixU();
  cout<< "initial S: " << dataModel->S << endl;
  /// @bug reset Z, Lambda and d?
  EsS = pBasisSelRegEnergy(new BasisNormalizeColEnergy());
  EsS->setPenalty(p);
}

void MtlOptSolver::hierSolve(){
  
}

bool MtlOptSolver::solve(){

  DEBUG_LOG("Z0.norm(): " << dataModel->Z.norm() << endl);
  if (!mtlOpt){
	optimizeZ();
	return true;
  }

  bool succ = false;
  assert_eq(r_min, dataModel->S.cols());

  // computeInitialS();

  int step = 1;
  for (int rs = r_min; rs <= r_max; rs += step){

	INFO_LOG("update rs = "<< rs);

	const double fk = funValue();	
	const MatrixXd Zk = dataModel->Z;
	const MatrixXd Sk = dataModel->S;
	expandS(rs);
	const MatrixXd &S = dataModel->S;
	fun_Es = EsS->fun( &S(0,0), S.rows(), S.cols() );

	optimizeZ();
	succ = solve(maxIt,tol);
	
	const double fk1 = funValue();
	if (rs > r_min && fk1 > fk){

	  assert_gt(dataModel->Z.rows(), Zk.rows());
	  assert_gt(dataModel->S.cols(), Sk.cols());
	  dataModel->Z.setZero();
	  dataModel->Z.topRows(Zk.rows()) = Zk;
	  dataModel->S.setZero();
	  dataModel->S.leftCols(Sk.cols()) = Sk;
	  step *= 2;
	  cout << "step = " << step << endl;
	}
	if (step != 1 && rs+step>r_max && rs != r_max){
	  step = r_max-rs;
	}
	if (rs == r_max && r_min != r_max){
	  fun_Es = EsS->fun( &S(0,0), S.rows(), S.cols() );
	  succ = solve(1,tol);
	}

	if (!succ) break;
  }
  return succ;
}

bool MtlOptSolver::solve(const int maxIt, const double tol){

  DEBUG_LOG("outer-it-fun(Z): "<<funValue());
  for (int it = 0; it < maxIt; ++it){

	const double oldFunZ = funValue();
	optimizeLambdaDampingDiag();
	DEBUG_LOG("outer-it-fun(Lambda): "<<funValue());

	optimizeS_lbfgs();
	DEBUG_LOG("outer-it-fun(S): "<<funValue());

	const double oldFunLambda = funValue();
	optimizeZ();
	DEBUG_LOG("outer-it-fun(Z): "<<funValue());

	if(oldFunZ - funValue() <= tol)  break;

  }

  DEBUG_LOG("optimal function value: "<<funValue());

  return true;
}

void MtlOptSolver::optimizeZ(){

  // update Ew(Z), Ec(Z) and Ek(Z).
  EwZ->reset(true);
  EcZ->reset();

  // set optimization parameters.
  const double epsg = 1e-3;
  const double epsf = 0.0f;
  const double epsx = 0.0f;
  const ae_int_t maxits = 3000;

  // set initial value
  static real_1d_array subZ;
  const int r = dataModel->reducedDim();
  const int Ts = dataModel->subFrames();
  const MatrixXd &Z = dataModel->Z;
  assert_ge(Z.size(),r*Ts);

  subZ.setlength(r*Ts);
  memcpy( &subZ[0], &Z(0,dataModel->T_begin), sizeof(double)*subZ.length() );

  // create solver
  minlbfgsstate state;
  minlbfgsreport rep;
  minlbfgscreate(1, subZ, state);
  minlbfgssetcond(state, epsg, epsf, epsx, maxits);

  // solve and get results
  alglib::minlbfgsoptimize(state, FunGradOptZ, NULL, (void*)this);
  minlbfgsresults(state, subZ, rep);
  dataModel->updateZ(&subZ[0]);

}

void MtlOptSolver::optimizeK(){
  
  // optimizeAtA_lm();
  // optimizeAtA_ipopt();
  optimizeK_sym();
}

void MtlOptSolver::optimizeK_sym(){
  
  FUNC_TIMER();

  EwK->reset();

  const int r = dataModel->reducedDim();
  const int n = r*(r+1)/2;
  VectorXd x(n), g(n);
  MatrixXd H(n,n);
  x.setZero();

  /// @bug for printing, should be removed
  fun_Ew=EwK->fun(&x[0], x.size());
  DEBUG_LOG( "inner-it-fun(K): "<<funValue() );
  
  EwK->grad(&x[0], x.size(), &g[0]);
  EwK->hessian(&x[0], x.size(), &H(0,0));
  g *= -1.0f;
  x = H.llt().solve(g);

  dataModel->K.resize(r,r);
  for (int i = 0; i < r; ++i){
    for (int j = 0; j < r; ++j){
	  const int p = CtrlForcesEnergyK_AD::symIndex(i,j);
	  assert_in(p,0,x.size()-1);
	  dataModel->K(i,j) = x[p];
	}
  }
  assert_eq(dataModel->K,dataModel->K.transpose());

  dataModel->decomposeK();
  for (int i = 0; i < dataModel->Lambda.size(); ++i){
	assert_eq(dataModel->K(i,i), dataModel->Lambda[i]);
    if (dataModel->Lambda[i] < 0){
	  dataModel->Lambda[i] = 0.0f;
	  dataModel->K(i,i) = 0.0f;
	}
  }

  fun_Ew=EwK->fun(&x[0], x.size());
  DEBUG_LOG( "inner-it-fun(K): "<<funValue() );
  
}

void MtlOptSolver::optimizeS_lm(){

  // update Ec(S)
  EcS->reset();
  EkS->reset();

  // set optimization parameters.
  const double epsg = 0.001f;
  const double epsf = 0.0f;
  const double epsx = 0.0f;
  const ae_int_t maxits = 100;

  // set initial value
  static real_1d_array S;
  S.setlength(dataModel->S.size());
  assert_gt(dataModel->S.size(),0);
  memcpy( &S[0], &(dataModel->S(0,0)), sizeof(double)*S.length() );

  // create solver
  minlmstate state;
  minlmreport rep;
  minlmcreatefgh(S, state);
  minlmsetcond(state, epsg, epsf, epsx, maxits);

  // solve and get results
  minlmoptimize(state, FunOptS, FunGradOptS, FunGradHessOptS, NULL, ((void*)this));
  minlmresults(state, S, rep);

  memcpy( &(dataModel->S(0,0)), &S[0], sizeof(double)*S.length() );

  alglibErrorReport(int(rep.terminationtype));
}

void MtlOptSolver::optimizeS_lbfgs(){

  FUNC_TIMER();

  // update Ec(S)
  EcS->reset();
  EkS->reset();

  // set optimization parameters.
  const double epsg = tol_s_g;
  const double epsf = tol_s_f;
  const double epsx = 0.0f;
  const ae_int_t maxits = maxIt_s;

  // set initial value
  static real_1d_array S;
  S.setlength(dataModel->S.size());
  assert_gt(dataModel->S.size(),0);
  memcpy( &S[0], &(dataModel->S(0,0)), sizeof(double)*S.length() );

  // create solver
  state = alglib::minlbfgsstate();
  minlbfgscreate(3, S, state);
  minlbfgssetcond(state, epsg, epsf, epsx, maxits);

  {
    // // scale 
	// static int it = 0;
	// it ++;
	// if (it > 5){
	// 	real_1d_array scale;
	// 	scale.setlength(S.length());
	
	// 	// scale using z, not help
	// 	// for (int j = 0; j < dataModel->S.cols(); ++j){
	// 	//   const double nj = dataModel->Z.row(j).norm();
	// 	//   assert_gt(nj,0.0f);
	// 	//   for (int i = 0; i < dataModel->S.rows(); ++i){
	// 	// 	scale[j*dataModel->S.rows()+i] = 1.0f/nj;
	// 	//   }
	// 	// }

	// 	// // scale using lambda0, not help
	// 	// for (int i = 0; i < dataModel->S.rows(); ++i){
	// 	//   for (int j = 0; j < dataModel->S.cols(); ++j){
	// 	// 	scale[j*dataModel->S.rows()+i] = 1.0f/EsS->getInitialLambda()[i];
	// 	//   }
	// 	// }

	// 	// minlbfgssetscale(state, scale);
	// 	// minlbfgssetprecscale(state);
	// }

	// solve and get results
  }

  getLbfgsState().c_ptr()->repiterationscount=0;
  alglib::minlbfgsoptimize(state, FunGradOptS, NULL, (void*)(this));
  DEBUG_LOG("iterations of (S): "<<getLbfgsState().c_ptr()->repiterationscount);
  // assert_gt(getLbfgsState().c_ptr()->repiterationscount,0);
  DEBUG_LOG("inner-it-fun-grad(S): "<<funValue());

  minlbfgsresults(state, S, rep);
  memcpy( &(dataModel->S(0,0)), &S[0], sizeof(double)*S.length() );

  alglibErrorReport(int(rep.terminationtype));
  // cout<< "S iterations count: " << rep.iterationscount << endl;
  // DEBUG_LOG("S^t: " << dataModel->S.transpose() << endl);
  // DEBUG_LOG("norm(S): " << dataModel->S.norm() << endl);
}

void MtlOptSolver::optimizeS_cg(){
  
  // update Ec(S)
  EcS->reset();
  EkS->reset();

  // set optimization parameters.
  const double epsg = tol_s_g;
  const double epsf = tol_s_f;
  const double epsx = 0.0f;
  const ae_int_t maxits = maxIt_s;

  // set initial value
  static real_1d_array S;
  S.setlength(dataModel->S.size());
  assert_gt(dataModel->S.size(),0);
  memcpy( &S[0], &(dataModel->S(0,0)), sizeof(double)*S.length() );

  // create solver
  mincgstate state;
  mincgreport rep;
  mincgcreate(S, state);
  mincgsetcond(state, epsg, epsf, epsx, maxits);

  {
    // // scale 
	// static int it = 0;
	// it ++;
	// if (it > 5){
	//   real_1d_array scale;
	//   scale.setlength(S.length());
	
	//   // scale using z, not help
	//   for (int j = 0; j < dataModel->S.cols(); ++j){
	// 	const double nj = dataModel->Z.row(j).norm();
	// 	assert_gt(nj,0.0f);
	// 	for (int i = 0; i < dataModel->S.rows(); ++i){
	// 	  scale[j*dataModel->S.rows()+i] = 1.0f/nj;
	// 	}
	//   }

	// 	// scale using lambda0, not help
	// 	// for (int i = 0; i < dataModel->S.rows(); ++i){
	// 	//   for (int j = 0; j < dataModel->S.cols(); ++j){
	// 	// 	scale[j*dataModel->S.rows()+i] = 1.0f/EsS->getInitialLambda()[i];
	// 	//   }
	// 	// }

	// 	// mincgsetscale(state, scale);
	// 	// mincgsetprecscale(state);
	// }
  }

  // solve and get results
  alglib::mincgoptimize(state, FunGradOptS, NULL, (void*)(this));
  mincgresults(state, S, rep);
  memcpy( &(dataModel->S(0,0)), &S[0], sizeof(double)*S.length() );

  alglibErrorReport(int(rep.terminationtype));
}

void MtlOptSolver::optimizeS_lbfgs_num(){

  FUNC_TIMER();
  
  // update Ec(S)
  EcS->reset();
  EkS->reset();

  // set optimization parameters.
  const double epsg = tol_s_g;
  const double epsf = tol_s_f;
  const double epsx = 0.0f;
  const ae_int_t maxits = maxIt_s;
  double diffstep = 1.0e-6;

  // set initial value
  static real_1d_array S;
  S.setlength(dataModel->S.size());
  assert_gt(dataModel->S.size(),0);
  memcpy( &S[0], &(dataModel->S(0,0)), sizeof(double)*S.length() );

  // create solver
  minlbfgsstate state;
  minlbfgsreport rep;
  minlbfgscreatef(3, S, diffstep, state);
  minlbfgssetcond(state, epsg, epsf, epsx, maxits);

  alglib::minlbfgsoptimize(state, FunOptS, NULL, (void*)(this));
  minlbfgsresults(state, S, rep);
  memcpy( &(dataModel->S(0,0)), &S[0], sizeof(double)*S.length() );

  alglibErrorReport(int(rep.terminationtype));
}

void MtlOptSolver::optimizeAtA_lbfgs(){
  
  FUNC_TIMER();

  // update Ew(K)
  EwK->reset();

  // set optimization parameters.
  const double epsg = tol_s_g;
  const double epsf = tol_s_f;
  const double epsx = 0.0f;
  const ae_int_t maxits = maxIt_s;
  double diffstep = 1.0e-6;

  // set initial value
  static real_1d_array K_sqrt;
  K_sqrt.setlength(dataModel->K.size());
  assert_gt(K_sqrt.length(),0);
  memcpy( &K_sqrt[0], &(dataModel->K(0,0)), sizeof(double)*K_sqrt.length() );

  const int r = dataModel->reducedDim();
  for (int i = 0; i < r; ++i){
    for (int j = 0; j < r; ++j){
	  if ( i == j){
		assert_ge_ext(K_sqrt[i*r+j],0.0f,"i="<<i);
		K_sqrt[j*r+i] = sqrt(K_sqrt[j*r+i]);
	  }else{
		assert_eq_ext(K_sqrt[i*r+j],0.0f,"i="<<i<<",j="<<j);
	  }
	}
  }

  // create solver
  minlbfgsstate state;
  minlbfgsreport rep;
  minlbfgscreate(3, K_sqrt, state);
  minlbfgssetcond(state, epsg, epsf, epsx, maxits);

  alglib::minlbfgsoptimize(state, FunGradOptK, NULL, (void*)(this));
  minlbfgsresults(state, K_sqrt, rep);

  // compute K = At*A, decompose K = Ut \Lambda U
  const MatrixXd &A = Map<MatrixXd>(&K_sqrt[0], r, r);
  dataModel->K = A.transpose()*A;
  dataModel->decomposeK();

  alglibErrorReport(int(rep.terminationtype));
}

void MtlOptSolver::optimizeAtA_lm(){

  FUNC_TIMER();

  // update Ew(K)
  EwK->reset();

  // set optimization parameters.
  const double epsg = tol_s_g;
  const double epsf = tol_s_f;
  const double epsx = 0.0f;
  const ae_int_t maxits = maxIt_s;
  double diffstep = 1.0e-6;

  // set initial value
  static real_1d_array K_sqrt;
  K_sqrt.setlength(dataModel->K.size());
  assert_gt(K_sqrt.length(),0);
  memcpy( &K_sqrt[0], &(dataModel->K(0,0)), sizeof(double)*K_sqrt.length() );

  const int r = dataModel->reducedDim();
  for (int i = 0; i < r; ++i){
    for (int j = 0; j < r; ++j){
	  if ( i == j){
		assert_ge_ext(K_sqrt[i*r+j],0.0f,"i="<<i);
		K_sqrt[j*r+i] = sqrt(K_sqrt[j*r+i]);
	  }else{
		assert_eq_ext(K_sqrt[i*r+j],0.0f,"i="<<i<<",j="<<j);
	  }
	}
  }

  // create solver, and solve
  minlmstate state;
  minlmreport rep;
  minlmcreatefgh(K_sqrt, state);
  minlmsetcond(state, epsg, epsf, epsx, maxits);
  minlmoptimize(state, FunOptK, FunGradOptK, FunGradHessOptK, NULL, ((void*)this));
  minlmresults(state, K_sqrt, rep);

  // compute K = At*A, decompose K = Ut \Lambda U
  const MatrixXd &A = Map<MatrixXd>(&K_sqrt[0], r, r);
  dataModel->K = A.transpose()*A;
  dataModel->decomposeK();

  alglibErrorReport(int(rep.terminationtype));
}

void MtlOptSolver::optimizeLambda(){
  
  assert(EwLambda);
  EwLambda->reset();
  static VectorXd g0;
  static VectorXd H;

   /// @bug for printing, should be removed
  fun_Ew=EwLambda->fun(dataModel->Lambda);
  DEBUG_LOG( "inner-it-fun(Lambda): "<<funValue() );

  EwLambda->grad(g0);
  EwLambda->hessian(H);
  assert_eq(H.size(),g0.size());
  assert_eq(dataModel->Lambda.size(),g0.size());

  for (int i = 0; i < g0.size(); ++i){

	if (0 != H[i]){

	  const double la = -g0[i]/H[i];
	  if (la < 0.0f){
		dataModel->Lambda[i] = 0.0f;
	  }else if(la > maxFrequency){
		dataModel->Lambda[i] = maxFrequency;
	  }else{
		dataModel->Lambda[i] = la;
	  }
	}
  }


  dataModel->Damping = dataModel->alphaM+dataModel->alphaK.asDiagonal()*dataModel->Lambda;

  DEBUG_LOG( "maxFrequency: "<< maxFrequency);
  DEBUG_LOG( "Lambda: "<< dataModel->Lambda.transpose() );
  fun_Ew=EwLambda->fun(dataModel->Lambda);
  DEBUG_LOG( "inner-it-fun(Lambda): "<<funValue() );
}

void MtlOptSolver::optimizeDamping(){
  
  assert(EwDamping);
  EwDamping->reset();
  static VectorXd g0;
  static MatrixXd H;

  /// @bug for printing, should be removed
  fun_Ew=EwDamping->fun(dataModel->alphaM, dataModel->alphaK);
  DEBUG_LOG( "inner-it-fun(Damping): "<<funValue() );

  const int r = dataModel->alphaK.size();
  EwDamping->grad(g0);
  assert_eq(r*2,g0.size());
  EwDamping->hessian(H);

  VectorXd am_ak = H.fullPivLu().solve(-g0);
  for (int i = 0; i < am_ak.size(); ++i){
	if (am_ak[i] < 0.0f) am_ak[i] = 0.0f;
  }
  dataModel->alphaM = am_ak.head(r);
  dataModel->alphaK = am_ak.tail(r);

  dataModel->Damping = dataModel->alphaM+dataModel->alphaK.asDiagonal()*dataModel->Lambda;

  DEBUG_LOG( "Damping_M: "<< dataModel->alphaM.transpose() );
  DEBUG_LOG( "Damping_K: "<< dataModel->alphaK.transpose() );
  fun_Ew=EwDamping->fun(dataModel->alphaM, dataModel->alphaK);
  DEBUG_LOG( "inner-it-fun(Damping): "<<funValue() );
}

void MtlOptSolver::optimizeLambdaDamping(){

  FUNC_TIMER();
  EwLambdaDamping->reset();  

  /// @bug for printing, should be removed
  // fun_Ew=EwLambdaDamping->fun(dataModel->Lambda, dataModel->Damping);
  // DEBUG_LOG( "inner-it-fun(Lambda): "<<funValue() );
  
  SparseMatrix<double> HLower;
  EwLambdaDamping->hessianLower(HLower);
  VectorXd g0;
  EwLambdaDamping->grad(g0);
  
  /// @bug solve sparse directly
  const SparseMatrix<double> hs = HLower.selfadjointView<Lower>();
  const MatrixXd H = hs;
  VectorXd lambda_damping = H.fullPivLu().solve(-g0);
  
  for (int i = 0; i < lambda_damping.size(); ++i){
	if (lambda_damping[i] < 0.0f) lambda_damping[i] = 0.0f;
  }

  const int r = dataModel->reducedDim();
  dataModel->Lambda = lambda_damping.head(r);
  dataModel->Damping = lambda_damping.tail(r);
  
  fun_Ew=EwLambdaDamping->fun(dataModel->Lambda, dataModel->Damping);
  DEBUG_LOG( "inner-it-fun(Lambda): "<<funValue() );

}

void MtlOptSolver::optimizeLambdaDampingDiag(){
  
  FUNC_TIMER();
  EwLambdaDamping->reset();

  /// @bug for printing, should be removed
  // fun_Ew=EwLambdaDamping->fun(dataModel->Lambda, dataModel->Damping);
  // DEBUG_LOG( "inner-it-fun(Lambda): "<<funValue() );

  VectorXd H_lambda, H_damping;
  EwLambdaDamping->hessianLambda(H_lambda);
  EwLambdaDamping->hessianDamping(H_damping);
  VectorXd g0;
  EwLambdaDamping->grad(g0);

  const int r = dataModel->reducedDim();
  for (int i = 0; i < r; ++i){

	if (0 != H_lambda[i]){
	  const double la = -g0[i]/H_lambda[i];
	  if (la < 0.0f){
		dataModel->Lambda[i] = 0.0f;
	  }else if(la > maxFrequency){
		dataModel->Lambda[i] = maxFrequency;
	  }else{
		dataModel->Lambda[i] = la;
	  }
	}

	if (0 != H_damping[i]){
	  const double d = -g0[i+r]/H_damping[i];
	  if (d < 0.0f){
		dataModel->Damping[i] = 0.0f;
	  }else{
		dataModel->Damping[i] = d;
	  }
	}
  }

  fun_Ew=EwLambdaDamping->fun(dataModel->Lambda, dataModel->Damping);
  DEBUG_LOG( "inner-it-fun(Lambda): "<<funValue() );
}

void MtlOptSolver::optimizeLambdaDamping_alglib(){

  FUNC_TIMER();
  EwLambdaDamping->reset();

  double epsg = 1e-13;
  double epsf = 0.0f;
  double epsx = 0.0f;
  ae_int_t maxits = maxIt_z;

  const int r = dataModel->reducedDim();
  real_1d_array x;
  real_1d_array bndl;
  real_1d_array bndu;
  x.setlength(r*2);
  bndl.setlength(r*2);
  bndu.setlength(r*2);

  for (int i = 0; i < r; ++i){

    x[i]  = dataModel->Lambda[i];
    x[i+r]= dataModel->Damping[i];

	bndl[i] = 0;
	bndl[i+r] = 0;

	bndu[i] = maxFrequency;
	bndu[i+r] = 2e19;
  }

  minbleicstate state;
  minbleicreport rep;
  minbleiccreate(x, state);
  minbleicsetbc(state, bndl, bndu);
  minbleicsetcond(state, epsg, epsf, epsx, maxits);
  alglib::minbleicoptimize(state, FunGradOptLambdaDamping, NULL, ((void*)this));
  minbleicresults(state, x, rep);

  for (int i = 0; i < r; ++i){
    dataModel->Lambda[i] = x[i];
    dataModel->Damping[i] = x[i+r];
  }

  alglibErrorReport(int(rep.terminationtype));

}

bool MtlOptSolver::alglibErrorReport(const int code)const{
  
  switch(code){
  case -7:cout << "gradient verification failed.\n"; break;
  case -3:cout << "inconsistent constraints.\n"; break;
  case 1:cout << "relative function improvement is no more than EpsF.\n"; break;
  case 2:cout << "relative step is no more than EpsX.\n"; break;
  case 4:cout << "gradient norm is no more than EpsG.\n"; break;
  case 5:cout << "MaxIts steps was taken.\n"; break;
  case 7:cout << "stopping conditions are too stringent,further improvement is impossible.\n"; break;
  default: cout << "the stop code is unkown: " << code << endl;
  }
  return code > 0;
}
