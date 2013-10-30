#include "MtlOptEnergy.h"
using namespace LSW_ANI_EDITOR;

void CtrlForceEnergy::precompute(){
  
  const int r = reducedDim();
  _G.resize(r*2,2);
  _S.resize(r*2,1);
  for (int i = 0; i < r; ++i){

	const double lambda = _Lambda[i];
	const double d = _diagD[i];
	const double gamma = 1.0f/(_h*(d+_h*lambda)+1.0f);

	const int i2 = i*2;
	const int i2_1 = i2+1;
	_G(i2,0) = gamma;
	_G(i2,1) = -_h*lambda*gamma;
	_G(i2_1,0) = _h*gamma;
	_G(i2_1,1) = 1.0f-_h*_h*lambda*gamma;

	_S(i2,0) = _h*gamma;
	_S(i2_1,0) = _h*_h*gamma;
  }
}

double CtrlForceEnergy::fun(const double *x){
  
  const VectorXd &w = Map<VectorXd>(const_cast<double*>(x),dim());
  double objValue = w.norm();
  objValue = objValue*objValue;

  forward(w,_V,_Z);
  static VectorXd ui;
  for (size_t i = 0; i < _conFrames.size(); ++i){
	const int f = _conFrames[i];
    const VectorXd &z = _Z.col(f);
	_warper->warp(z,f,_conNodes[i],ui);
	assert_eq(_uc[i].size(),ui.size());
	const double n = (_uc[i]-ui).norm();
	objValue += _penaltyCon*n*n;
  }
  return objValue*0.5f;
}

void CtrlForceEnergy::grad(const double *x,double *g){
  
  const VectorXd &w = Map<VectorXd>(const_cast<double*>(x),dim());
  forward(w,_V,_Z);
  compute_pEpz(_pEpz);
  adjoint(_r);

  static MatrixXd G;
  G.resize(reducedDim(),getT()-1);
  memcpy(&G(0,0),x,dim()*sizeof(double));
  for (int i = 0; i < reducedDim(); ++i){
	const int i2 = i*2;
	const int i2_1 = i2+1;
    G.row(i) += _S(i2,0)*_r.row(i2);
    G.row(i) += _S(i2_1,0)*_r.row(i2_1);
  }
  memcpy(g,&G(0,0),dim()*sizeof(double));
}

void CtrlForceEnergy::adjoint(MatrixXd &R)const{
  
  const int r = reducedDim();
  const int T = getT();
  R.resize(r*2,T);
  R.col(T-1).setZero();
  if (_confIndex[T-1] >= 0)
	R.block(r,T-1,r,1) = _pEpz.col(_confIndex[T-1]);

  for (int f = T-2; f >= 0; --f){
	const int f_1 = f+1;
	for (int i = 0; i < r; ++i){
	  const int i2 = i*2;
	  const int i2_1 = i2 + 1;
	  const double g0 = _G(i2,0);
	  const double g1 = _G(i2,1);
	  const double g2 = _G(i2_1,0);
	  const double g3 = _G(i2_1,1);
	  R(i2,f) = g0*R(i2,f_1) + g2*R(i2_1,f_1);
	  R(i2_1,f) = g1*R(i2,f_1) + g3*R(i2_1,f_1);
	}
	if (_confIndex[f] >= 0)
	  R.block(r,f,r,1) += _pEpz.col(_confIndex[f]);
  }
}

void CtrlForceEnergy::compute_pEpz(MatrixXd &pEpz)const{
  
  if (_conFrames.size() > 0){
	const int r = reducedDim();
	pEpz.resize(r,_conFrames.size());
	for (size_t i = 0; i < _conFrames.size(); ++i){
	  
	  const int f = _conFrames[i];
	  const VectorXd &z = _Z.col(f);
	  /// @todo
	  // _warper->warp(z,f,_conNodes[i],ui);
	  // assert_eq(_uc[i].size(),ui.size());
	  // const double n = (_uc[i]-ui).norm();
	  // objValue += _penaltyCon*n*n;
	}
  }
}

void CtrlForceEnergy::forward(const VectorXd &w,MatrixXd &V,MatrixXd &Z){

  const int T = getT();
  const int r = reducedDim();
  V.resize(r,T);
  Z.resize(r,T);
  V.col(0) = _v0;
  Z.col(0) = _z0;

  for (int i = 0; i < r; ++i){

	const int i2 = i*2;
	const int i2_1 = i2+1;
	const double g0 = _G(i2,0);
	const double g1 = _G(i2,1);
	const double g2 = _G(i2_1,0);
	const double g3 = _G(i2_1,1);
	const double s0 = _S(i2,0);
	const double s1 = _S(i2_1,0);
	for (int f = 1; f < T; ++f){
	  const int f1 = f - 1;
	  const double wi = w[f1*r+i];
	  V(i,f) = g0*V(i,f1) + g1*Z(i,f1) + s0*wi;
	  Z(i,f) = g2*V(i,f1) + g3*Z(i,f1) + s1*wi;
	}
  }
}
