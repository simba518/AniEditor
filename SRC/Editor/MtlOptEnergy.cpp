#include <set>
#include <MtlOptEnergy.h>
#include <SparseMatrixTools.h>
using namespace std;
using namespace LSW_ANI_EDITOR;

void CtrlForceEnergy::precompute(){

  const int r = reducedDim();
  if (_subZ.rows() != r || _subZ.cols() != getT()-_keyframes.size()){
	_subZ.resize(r,getT()-_keyframes.size());
	_subZ.setZero();
  }

  const double _1h = 1.0f/_h;
  const double _1h2 = _1h*_1h;
  _Ma = _diagD*_1h;
  _Mb = _Lambda-_Ma;
  for (int i = 0; i < r; ++i){
    _Ma[i] += _1h2;
	_Mb[i] -= 2.0f*_1h2;
  }
  
  _H.resize(r);
  const double c = _1h2;
  for (int i = 0; i < r; ++i){

    const double a = _Ma[i];
	const double b = _Mb[i];
	const double c2 = c*c;
	const double b2 = b*b;
	const double a2 = a*a;
	const double bc = b*c;
	const double ac = a*c;
	const double ab = a*b;
	
	_H[i](0,0) = c2;
	_H[i](0,1) = bc;
	_H[i](0,2) = ac;

	_H[i](1,0) = b2+a2;
	_H[i](1,1) = c2+b2;
	_H[i](1,2) = bc+ab;
	
	_H[i](2,0) = ab;
	_H[i](2,1) = a2;
	_H[i](2,2) = c2+b2+a2;
  }
}

void CtrlForceEnergy::keyframeChanged(){

  set<int> keyf;
  for (int i = 0; i < _keyframes.size(); ++i){
    keyf.insert(_keyframes[i]);
  }
  if (_keyframes.size() > 0){
	EIGEN3EXT::genReshapeMatrix(getT(),keyf,_insertKeyZ);
  }else{
	_insertKeyZ.resize(0,0);
  }
  if (_Z.rows() == reducedDim() && _Z.cols() == getT()){
	if (_keyframes.size() > 0){
	  _subZ = _Z*_insertKeyZ.transpose();
	}else{
	  _subZ = _Z;
	}
  }else{
	_subZ.resize(reducedDim(),getT()-_keyframes.size());
	_subZ.setZero();
  }
}

void CtrlForceEnergy::addKeyframes(MatrixXd &Z)const{

  if (_keyframes.size() > 0){
	assert_eq(_insertKeyZ.rows(),Z.cols());
	assert_eq(_insertKeyZ.cols(),getT());
	Z = Z*_insertKeyZ;
	for (int i = 0; i < _keyframes.size(); ++i){
	  const int f = _keyframes[i];
	  Z.col(f) = _keyZ[i];
	}
  }
}

double CtrlForceEnergy::fun(const double *x){

  const int T = getT();
  const int r = reducedDim();
  _Z = Map<MatrixXd>(const_cast<double*>(x),r,T-_keyframes.size());
  addKeyframes(_Z);
  const MatrixXd &Z0 = _Z.block(0,0,r,T-2);
  const MatrixXd &Z1 = _Z.block(0,1,r,T-2);
  const MatrixXd &Z2 = _Z.block(0,2,r,T-2);
  const double objValue = (_Ma.asDiagonal()*Z2+_Mb.asDiagonal()*Z1+Z0*(1.0f/(_h*_h))).norm();
  return objValue*objValue*0.5f;
}

void CtrlForceEnergy::grad(const double *x,double *g){
  
  const int T = getT();
  const int r = reducedDim();
  _Z = Map<MatrixXd>(const_cast<double*>(x),r,T-_keyframes.size());
  addKeyframes(_Z);

  MatrixXd G(r,T);
  for (int i = 0; i < r; ++i){

	const Matrix3d &H = _H[i];
	const VectorXd &z = _Z.row(i).transpose();
	const double &c2 = H(0,0);
	const double &bc = H(0,1);
	const double &ac = H(0,2);
	const double &b2a2 = H(1,0);
	const double &c2b2 = H(1,1);
	const double &bcab = H(1,2);
	const double &ab = H(2,0);
	const double &a2 = H(2,1);
	const double &c2b2a2 = H(2,2);
	
	G(i,0) = c2*z[0]+bc*z[1]+ac*z[2];
	G(i,1) = bc*z[0]+c2b2*z[1]+bcab*z[2]+ac*z[3];
	G(i,T-2) = ac*z[T-4]+bcab*z[T-3]+b2a2*z[T-2]+ab*z[T-1];
	G(i,T-1) = ac*z[T-3]+ab*z[T-2]+a2*z[T-1];
	for (int f = 2; f < T-2; ++f){
	  G(i,f) = ac*z[f-2]+bcab*z[f-1]+c2b2a2*z[f]+bcab*z[f+1]+ac*z[f+2];
	}
  }
  if (_keyframes.size() > 0)
	G = G*_insertKeyZ.transpose();
  assert_eq(G.size(),dim());
  memcpy(g,&G(0,0),sizeof(double)*G.size());
}

void MtlOptEnergy::setZ(const MatrixXd &Z){

  const int T = getT();
  const int r = reducedDim();
  assert_eq(Z.rows(),r);
  assert_eq(Z.cols(),T);
  const MatrixXd &Z0 = Z.block(0,0,r,T-2);
  const MatrixXd &Z2 = Z.block(0,2,r,T-2);

  _Z1 = Z.block(0,1,r,T-2);
  _Z2v = (Z2-_Z1)/_h;
  _Z2a = (Z2-2*_Z1+Z0)/(_h*_h);
}

double MtlOptEnergy::fun(const double *x){

  computeKD(x,_K,_D);
  const double objValue = (_Z2a+_D*_Z2v+_K*_Z1).norm();
  return objValue*objValue*0.5f;
}

void MtlOptEnergy::grad(const double *x,double *g){
  
  const int r = reducedDim();
  const double &ak = x[0];
  const double &am = x[1];
  const MatrixXd &A = Map<MatrixXd>(const_cast<double*>(&x[2]),r,r);
  _K = A.transpose()*A;

  {// grad ak, am
	const MatrixXd F1 = _K*_Z1+_Z2a;
	MatrixXd Fm = ak*_K*_Z2v+F1;
	MatrixXd Gm = _K*_Z2v;
	MatrixXd Hm = am*_Z2v+F1;
	const VectorXd &E = Map<VectorXd>(&_Z2v(0,0),_Z2v.size());
	const VectorXd &F = Map<VectorXd>(&Fm(0,0),Fm.size());
	const VectorXd &G = Map<VectorXd>(&Gm(0,0),Gm.size());
	const VectorXd &H = Map<VectorXd>(&Hm(0,0),Hm.size());
	g[1] = (E.transpose()*(am*E+F))(0,0); // grad am
	g[0] = (G.transpose()*(ak*G+H))(0,0); // grad ak
  }

  {// grad A.
	const MatrixXd &E = _Z2v;
	const MatrixXd &At = A.transpose();
	const double &J = ak;
	const MatrixXd &L = _Z1;
	const MatrixXd M = am*_Z2v+_Z2a;
	const MatrixXd &Et = E.transpose();
	const double &Jt = J;
	const MatrixXd &Lt = L.transpose();
	const MatrixXd &Mt = M.transpose();
	const MatrixXd gA = (At*A*L*Lt*At).transpose()+(A*L*Et*At*A*Jt)+(A*L*Mt)+
	  (A*At*A*L*Lt)+(L*Et*At*A*Jt*At).transpose()+(A*M*Lt)+
	  (A*Jt*At*A*L*Et)+(E*Et*At*A*Jt*J*At).transpose()+(A*Jt*M*Et)+
	  (Jt*At*A*L*Et*At).transpose()+(A*E*Et*At*A*Jt*J)+(A*E*Mt*J);
	memcpy(&g[2],&gA(0,0),sizeof(double)*r*r);
  }
}
