#include <BaseMtlOptEnergy.h>
using namespace LSW_ANI_EDITOR;

void BaseCtrlForceEnergy::setTotalFrames(const int T){
  assert_ge(T,3);
  _T = T;
  _keyfIndex.resize(T);
  for (int i = 0; i < T; ++i)
	_keyfIndex[i] = -1;
  clearPartialCon();
}

void BaseCtrlForceEnergy::setKD(const MatrixXd &K,const MatrixXd &D){

  assert_eq(K.rows(),reducedDim());
  assert_eq(K.cols(),reducedDim());
  assert_eq(D.rows(),reducedDim());
  assert_eq(D.cols(),reducedDim());

  // diag K and D
  SelfAdjointEigenSolver<MatrixXd> eigenK(K);
  _Lambda = eigenK.eigenvalues();
  const MatrixXd diagD = _U.transpose()*D*_U;
  _diagD.resize(reducedDim());
  for (int i = 0; i < reducedDim(); ++i)
	_diagD[i] = diagD(i,i);
  if (diagD.norm() > 0.0){
	// assert_lt( (diagD.norm()-_diagD.norm())/diagD.norm(), 1e-4 ); ///@todo
  }else{
	assert_lt( D.norm()-_diagD.norm(), 1e-8 );
  }

  const MatrixXd &U = eigenK.eigenvectors();
  _U = _U*U;

  this->precompute();
}

void BaseCtrlForceEnergy::addKeyframe(const VectorXd &unRotZk, const int f){

  assert_in(f,0,getT()-1);
  assert_eq(unRotZk.size(),reducedDim());
  const VectorXd zk = _U.transpose()*unRotZk;
  if (_keyfIndex[f] >= 0){
	_keyZ[ _keyfIndex[f] ] = zk;
  }else{
	_keyZ.push_back(zk);
	_keyframes.push_back(f);
	_keyfIndex[f] = _keyframes.size()-1;
  }
}

void BaseCtrlForceEnergy::setKeyframes(const vector<VectorXd> &keyZ, const vector<int> &keyframes){

  clearKeyframes();
  assert_eq(keyZ.size(),keyframes.size());
  for (size_t i = 0; i < keyframes.size(); ++i)
	addKeyframe(keyZ[i],keyframes[i]);
}

void BaseCtrlForceEnergy::setKeyframes(const MatrixXd &keyZ, const vector<int> &keyframes){
  
  clearKeyframes();
  assert_eq(keyZ.cols(),keyframes.size());
  for (size_t i = 0; i < keyframes.size(); ++i)
	addKeyframe(keyZ.col(i),keyframes[i]);
}

void BaseCtrlForceEnergy::setPartialCon(vector<int>&conF,vector<vector<int> >&conN,vector<VectorXd>&uc){

  clearPartialCon();
  _conFrames = conF;
  _conNodes = conN;
  _uc = uc;
  for (int i = 0; i < _conFrames.size(); ++i){
	const int f = _conFrames[i];
	assert_in(f,0,_confIndex.size()-1);
	_confIndex[f] = i;
  }
}

void BaseCtrlForceEnergy::clearPartialCon(){

  _conFrames.clear();
  _conNodes.clear();
  _uc.clear();
  _confIndex.resize(_T);
  for (int i = 0; i < _T; ++i)
	_confIndex[i] = -1;
}

void BaseCtrlForceEnergy::clearKeyframes(){
  _keyframes.clear();
  _keyZ.clear();
  _keyfIndex.resize(getT());
  for (size_t i = 0; i < _keyfIndex.size(); ++i){
	_keyfIndex[i] = -1;
  }
}

void BaseMtlOptEnergy::setMtl(const VectorXd &Lambda,const double ak,const double am){

  assert_ge(ak,0.0f);
  assert_ge(am,0.0f);
  const int r = Lambda.size();
  _AkAmA.resize(r*r+2);
  _AkAmA[0] = ak;
  _AkAmA[1] = am;
  _AkAmA.tail(r*r).setZero();
  int c = 0;
  for (int i = 0; i < r; ++i){
	assert_ge(Lambda[i],0.0f);
	_AkAmA[2+c] = sqrt(Lambda[i]);
	c += (r+1);
  }
  _K.resize(r,r);
  _K.setZero();
}

void BaseMtlOptEnergy::computeKD(const double *x,MatrixXd &K,MatrixXd &D)const{
  const int r = reducedDim();
  const double &ak = x[0];
  const VectorXd &am = VectorXd::Ones(r)*x[1];
  const MatrixXd &A = Map<MatrixXd>(const_cast<double*>(&x[2]),r,r);
  K = A.transpose()*A;
  D = am.asDiagonal();
  D += ak*K;
}
