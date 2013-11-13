#ifndef _MTLOPTENERGYAD_H_
#define _MTLOPTENERGYAD_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <CASADITools.h>
#include <RedRSWarperAD.h>
using namespace std;
using namespace Eigen;
using namespace CASADI;
using namespace LSW_WARPING;
using CasADi::SX;
using CasADi::SXMatrix;

class RedSpaceTimeEnergyAD{
  
public:
  RedSpaceTimeEnergyAD(){
	_penaltyCon = 1.0f;
	_usePartialCon = true;
  }
  RedSpaceTimeEnergyAD(pRedRSWarperAD warper):_warper(warper){
	_penaltyCon = 1.0f;
	_usePartialCon = true;
  }
  void setT(const int T){_T = T;}
  template<class SCALAR> 
  void setTimestep(const SCALAR h){_h = h;}
  void setDamping(const double ak, const double am, const VectorXd &lambda){
	setDamping((VectorXd)(am*VectorXd::Ones(lambda.size())+ak*lambda));
  }
  void setDamping(const VectorXd &d){
	const MatrixXd D = d.asDiagonal();
	convert(D,_D);
  }
  void setDamping(const VSX &d){
	_D = makeEyeMatrix(d);
  }
  void setDamping(const SXMatrix &D){
	assert_eq(D.size1(),D.size2());
	_D = D;
  }
  void setDamping(const MatrixXd &D){
	assert_eq(D.rows(),D.cols());
	convert(D,_D);
  }
  template<class VECTOR>
  void setDamping(const VSX &ak,const VSX &am, const VECTOR &lambda){
	assert_eq(ak.size(),am.size());
	assert_eq(ak.size(),lambda.size());
	VSX d = am;
	for (int i = 0; i < ak.size(); ++i)
	  d[i] += ak[i]*lambda[i];
	setDamping(d);
  }
  void setK(const SXMatrix &K){ assert_eq(K.size1(),K.size2());  _K = K;}
  void setK(const MatrixXd &K){ assert_eq(K.rows(),K.cols()); convert(K,_K); }
  template<class VECTOR>
  void setK(const VECTOR &v){_K = makeEyeMatrix(v);}
  template<class VECTOR>
  void setKeyframes(const MatrixXd &Kz, const VECTOR &fid){
  	assert_eq(Kz.cols(),(int)fid.size());
  	_keyZ = Kz;
  	_keyId.resize(fid.size());
  	for (int i = 0; i < fid.size(); ++i)
  	  _keyId[i] = fid[i];
  }

  void setWarper(pRedRSWarperAD warper){
	_warper = warper;
  }
  void usePartialCon(const bool use){
	_usePartialCon = use;
  }
  void setPenaltyCon(const double penalty){
	assert_gt(penalty,0.0f);
	_penaltyCon = penalty;
  }
  void setPartialCon(const vector<int>&conF,
					 const vector<vector<int> >&conN,
					 const vector<VectorXd>&uc){
	_conFrames = conF;
	_conNodes = conN;
	_uc = uc;
  }

  void assembleEnergy();

  const SXMatrix &getEnergy()const{return _energy;}
  const VSX &getVarZ()const{return _varZ;}
  int reducedDim()const{
  	return _K.size1();
  }
  int getT()const{
	return _T;
  }
  int isKeyframe(const int f)const{
	return isKeyframe(_keyId,f);
  }
  void insertKeyframes(MatrixXd &Z)const;
  template<class VECTOR>
  static int isKeyframe(const VECTOR& kid,const int f){
  	int k = -1;
  	for (int i = 0; i < kid.size(); ++i){
  	  if(f == kid[i]){
  		k = i;
  		break;
  	  }
  	}
  	return k;
  }
  
private:
  double _penaltyCon;
  bool _usePartialCon;
  int _T;
  SX _h;
  SXMatrix _D;
  SXMatrix _K;
  SXMatrix _energy;
  VSX _varZ;

  pRedRSWarperAD _warper;
  MatrixXd _keyZ;
  VectorXi _keyId;
  vector<int> _conFrames;
  vector<vector<int> > _conNodes;
  vector<VectorXd> _uc;
};

class MtlDataModel{
  
public:
  void setT(const int T){
	this->T = T;
  }
  void setTimestep(const double h){
	this->h = h;
  }
  void setDamping(double ak, double am, const VectorXd &lambda){
	Am = VectorXd::Ones(lambda.size())*am;
	Ak = VectorXd::Ones(lambda.size())*ak;
	D = Am.asDiagonal();
	D += lambda.asDiagonal()*ak;
  }
  void setLambda(const VectorXd &lambda){
	K = lambda.asDiagonal();
  }
  template<class VECTOR>
  void setKeyframes(const MatrixXd &Z, const VECTOR &kid){
	assert_eq(kid.size(),Z.cols());
	this->keyZ = Z;
	this->keyId.resize(kid.size());
	for (int i = 0; i < kid.size(); ++i)
	  this->keyId[i] = kid[i];
  }
  void setSubZ(const MatrixXd &Z){
	subZ = Z;
  }
  void print()const{

	cout<< "K \n" << K << endl;
	SelfAdjointEigenSolver<MatrixXd> eigenK(K);
	cout<< "eigen(K): " << eigenK.eigenvalues().transpose() << endl;

	cout<< "D \n" << D << endl;
	SelfAdjointEigenSolver<MatrixXd> eigenD(D);
	cout<< "eigen(D): " << eigenD.eigenvalues().transpose() << endl;
  }
  MatrixXd getRotZ()const{
	return getUt()*Z;
  }
  MatrixXd getUt()const{
	SelfAdjointEigenSolver<MatrixXd> eigenK(K);
	return eigenK.eigenvectors().transpose();
  }

  void setPenaltyCon(const double penalty){
	assert_gt(penalty,0.0f);
	penaltyCon = penalty;
  }
  void setPartialCon(const vector<int>&conF,
					 const vector<vector<int> >&conN,
					 const vector<VectorXd>&uc){
	conFrames = conF;
	conNodes = conN;
	this->uc = uc;
  }
  void setWarper(pRedRSWarperAD w){
	warper = w;
  }

public:
  int T;
  double h;
  MatrixXd keyZ;
  VectorXi keyId;

  MatrixXd subZ;
  MatrixXd Z;
  MatrixXd K;
  MatrixXd D;
  VectorXd Ak;
  VectorXd Am;

  pRedRSWarperAD warper;
  double penaltyCon;
  vector<int> conFrames;
  vector<vector<int> > conNodes;
  vector<VectorXd> uc;
};
typedef boost::shared_ptr<MtlDataModel> pMtlDataModel;

class MtlOptimizer{
  
 public:
  MtlOptimizer(pMtlDataModel m):_model(m){
	useHessian = true;
	tol = 1e-3;
	maxIt = 500;
  }
  void setInitialValue(const VectorXd &x0){
	_rlst = x0;
  }
  void setUseHessian(const bool useH){
	useHessian = useH;
  }
  void setTol(const double t){
	assert_gt(t,0.0f);
	tol = t;
  }
  void setMaxIt(const int mit){
	maxIt = mit;
  }
  virtual void optimize();
  
 protected:
  void resetEnergy(const bool useAllZ = true);
  virtual void optimizationBegin(){}
  virtual void optimizationEnd(){}
  virtual const VSX &getVariable()const = 0;
  virtual bool getInitValue(VectorXd &x0, const int size)const; 
  
  static void produceSymetricMat(const string name,const int dim,VSX&s,SXMatrix &SM);
  static int symIndex(int r,int c);
  virtual bool getLowBound(vector<double> &lower)const{
	return false;
  }

 protected:
  pMtlDataModel _model;
  RedSpaceTimeEnergyAD _energy;
  CasADi::SXFunction _fun;
  CasADi::IpoptSolver _solver;
  VectorXd _rlst;
  
  bool useHessian;
  double tol;
  int maxIt;
};
typedef boost::shared_ptr<MtlOptimizer> pMtlOptimizer;

// optimize reduced coordinates
class ZOptimizer:public MtlOptimizer{
  
public:
  ZOptimizer(pMtlDataModel m):MtlOptimizer(m){

	const MatrixXd &Z = m->subZ;
	assert_gt(Z.size(),0);
	MatrixXd zz = Z;
	_rlst = Map<VectorXd>( &(zz(0,0)),zz.size() );
  }

protected:
  void optimizationBegin(){
	_energy.usePartialCon(true);
	resetEnergy(false);
  }
  void optimizationEnd(){
	const int r = _energy.reducedDim();
	const int subT = _model->T - _model->keyId.size();
	assert_eq(r*subT, _rlst.size());
	_model->subZ = Map<MatrixXd>(&_rlst[0],r,subT);
	_model->Z = _model->subZ;
	_energy.insertKeyframes(_model->Z);
  }
  const VSX &getVariable()const{
	return _energy.getVarZ();
  }

};

// optimize K
class LambdaOptimizer:public MtlOptimizer{

public: 
  LambdaOptimizer(pMtlDataModel m,const bool useAkAm = true):
	MtlOptimizer(m),_useAkAm(useAkAm){
	const int r = m->K.rows();
	setInitK(m->K);
	_lambda = makeSymbolic(r,"la");
	_K = makeEyeMatrix(_lambda);
  }
  
protected:
  void setInitK(const MatrixXd &K){
	assert_eq(K.rows(),K.cols());
	const int r = K.rows();
	_rlst.resize(r);
	for (int i = 0; i < r; ++i)
	  _rlst[i] = K(i,i);
  }
  void optimizationBegin(){
	_energy.usePartialCon(false);
	resetEnergy();
	const SXMatrix &K = _K;
	_energy.setK(K);
	if(_useAkAm){
	  const SXMatrix Ak = makeEyeMatrix(convert(_model->Ak));
	  const SXMatrix Am = makeEyeMatrix(convert(_model->Am));
	  const SXMatrix D = Am + Ak.mul(K);
	  _energy.setDamping(D);
	}else{
	  _energy.setDamping(_model->D);
	}
  }
  void optimizationEnd(){
	const int r = _energy.reducedDim();
	MatrixXd &K = _model->K;
	K.resize(r,r);
	K.setZero();
	for (int i = 0; i < r; ++i)
	  K(i,i) = _rlst[i];
  }
  const VSX &getVariable()const{
	return _lambda;
  }
  bool getLowBound(vector<double> &lower)const{
	lower = vector<double>(this->getVariable().size(),0.0f);
	return true;
  }

private:
  const bool _useAkAm;
  VSX _lambda;
  SXMatrix _K;
};

class AtAOptimizer:public MtlOptimizer{

public:
  AtAOptimizer(pMtlDataModel m,const bool useAkAm = true):
	MtlOptimizer(m),_useAkAm(useAkAm){

	const int r = _model->K.rows();
	_A = makeSymbolic(r*r,"a");
	const SXMatrix A = convert(_A,r);
	_K = CasADi::trans(A).mul(A);
	setInitVal(m->K);
  }

protected:
  void setInitVal(const MatrixXd &K){
	const int r = K.rows();
	assert_gt(r,0);
	assert_eq(K.cols(),r);
	MatrixXd KK(r,r);
	KK.setZero();
	for (int i = 0; i < r; ++i)
	  KK(i,i) = sqrt( K(i,i) );
	_rlst = Map<VectorXd>(&KK(0,0),r*r);
  }
  void optimizationBegin(){
	_energy.usePartialCon(false);
	resetEnergy();
	const SXMatrix &K = _K;
	_energy.setK(K);
	if(_useAkAm){
	  const SXMatrix Ak = makeEyeMatrix(convert(_model->Ak));
	  const SXMatrix Am = makeEyeMatrix(convert(_model->Am));
	  const SXMatrix D = Am + Ak.mul(K);
	  _energy.setDamping(D);
	}else{
	  _energy.setDamping(_model->D);
	}
  }
  void optimizationEnd(){

	const int r = _energy.reducedDim();
	assert_eq(_rlst.size(),r*r);
	const MatrixXd A = Map<MatrixXd>(&_rlst[0],r,r);
	_model->K = A.transpose()*A;
  }
  const VSX &getVariable()const{
	return _A;
  }

private:
  const bool _useAkAm;
  VSX _A;
  SXMatrix _K;
};

class KOptimizer:public MtlOptimizer{

public: 
  KOptimizer(pMtlDataModel m,const bool useAkAm = true):MtlOptimizer(m),_useAkAm(useAkAm){

	setInitK(m->K);
	produceSymetricMat("k",_model->K.rows(),_kx,_K);
  }
  
protected:
  void setInitK(const MatrixXd &K){

	const int r = K.rows();
	const int n = r*(r+1)/2;
	assert_eq(K.rows(),K.cols());

	_rlst.resize(n);
	for (int i = 0; i < r; ++i)
	  for (int j = 0; j < r; ++j)
		_rlst[symIndex(i,j)] = K(i,j);
  }
  void optimizationBegin(){

	_energy.usePartialCon(false);
	resetEnergy();
	const SXMatrix &K = _K;
	_energy.setK(K);
	if(_useAkAm){
	  const SXMatrix Ak = makeEyeMatrix(convert(_model->Ak));
	  const SXMatrix Am = makeEyeMatrix(convert(_model->Am));
	  const SXMatrix D = Am + Ak.mul(K);
	  _energy.setDamping(D);
	}else{
	  _energy.setDamping(_model->D);
	}
  }
  void optimizationEnd(){

	const int r = _energy.reducedDim();
	const int n = r*(r+1)/2;
	assert_eq(_rlst.size(),n);

	MatrixXd &K = _model->K;
	K.resize(r,r);
	for (int i = 0; i < r; ++i)
	  for (int j = 0; j < r; ++j)
		K(i,j) = _rlst[symIndex(i,j)];
  }
  const VSX &getVariable()const{
	return _kx;
  }

private:
  const bool _useAkAm;
  VSX _kx;
  SXMatrix _K;
};

// optimize damping
class akamOptimizer:public MtlOptimizer{

public:
  akamOptimizer(pMtlDataModel m):MtlOptimizer(m){

	setInitAkAm(m->Ak[0], m->Am[0]);
	const int r = _model->K.rows();
	assert_gt(r,0);
	const SX ak("ak"), am("am");
	for (int i = 0; i < r; ++i){
	  _Ak.push_back(ak);
	  _Am.push_back(am);
	}
	_akam.push_back(ak);
	_akam.push_back(am);
  }
  
protected:
  virtual void setInitAkAm(const double &Ak,const double &Am){
	_rlst.resize(2);
	_rlst[0] = Ak;
	_rlst[1] = Am;
  }
  virtual void optimizationBegin(){

	_energy.usePartialCon(false);
	resetEnergy();
	const int r = _energy.reducedDim();
	const SXMatrix K = convert(_model->K);
	const SXMatrix Ak = makeEyeMatrix(_Ak);
	const SXMatrix Am = makeEyeMatrix(_Am);
	const SXMatrix D = Am + Ak.mul(K);
	_energy.setDamping(D);
  }
  virtual void optimizationEnd(){

	const int r = _energy.reducedDim();
	const MatrixXd &K = _model->K;
	assert_eq(_rlst.size(), 2);
	const double ak = _rlst[0];
	const double am = _rlst[1];
	_model->D = am*MatrixXd::Identity(r,r)+ak*K;
	_model->Ak = VectorXd::Ones(r)*ak;
	_model->Am = VectorXd::Ones(r)*am;
  }
  const VSX &getVariable()const{
	return _akam;
  }
  bool getLowBound(vector<double> &lower)const{
	lower.resize(2);
	lower[0] = 0.0f;
	lower[1] = 0.0f;
	return true;
  }

protected:
  VSX _Ak;
  VSX _Am;
  VSX _akam;
};

class AkAmOptimizer:public MtlOptimizer{

public:
  AkAmOptimizer(pMtlDataModel m):MtlOptimizer(m){
	
	setInitAkAm(m->Ak, m->Am);
	const int r = _model->K.rows();
	assert_gt(r,0);
	_Ak = makeSymbolic(r,"ak");
	_Am = makeSymbolic(r,"am");
	_AkAm.clear();
	_AkAm.insert(_AkAm.end(), _Ak.begin(), _Ak.end());
	_AkAm.insert(_AkAm.end(), _Am.begin(), _Am.end());
  }
  
protected:
  virtual void setInitAkAm( const VectorXd &Ak, const VectorXd &Am){

	assert_eq(Ak.size(), Am.size());
	const int r = Ak.size();
	_rlst.resize(Ak.size()*2);
	_rlst.head(r) = Ak;
	_rlst.tail(r) = Am;
  }
  virtual void optimizationBegin(){

	_energy.usePartialCon(false);
	resetEnergy();
	const int r = _energy.reducedDim();
	const SXMatrix K = convert(_model->K);
	const SXMatrix Ak = makeEyeMatrix(_Ak);
	const SXMatrix Am = makeEyeMatrix(_Am);
	const SXMatrix D = Am + Ak.mul(K);
	_energy.setDamping(D);
  }
  virtual void optimizationEnd(){

	const int r = _energy.reducedDim();
	const MatrixXd &K = _model->K;
	assert_ge(_rlst.size(), r*2);
	const VectorXd &Ak = _rlst.head(r);
	const VectorXd &Am = _rlst.segment(r,r);
	MatrixXd D = Am.asDiagonal();
	D = D + Ak.asDiagonal()*K;
	_model->D = D;
	_model->Ak = Ak;
	_model->Am = Am;
  }
  const VSX &getVariable()const{
	return _AkAm;
  }

protected:
  VSX _Ak;
  VSX _Am;
  VSX _AkAm;
};

// optimize both K+damping
class KAtAmOptimizer:public AkAmOptimizer{

public: 
  KAtAmOptimizer(pMtlDataModel m):AkAmOptimizer(m){

	setInitK(m->K);
	produceSymetricMat("k",_model->K.rows(),_kx,_K);
	_AkAmKx = _AkAm;
	_AkAmKx.insert(_AkAmKx.end(), _kx.begin(), _kx.end());
  }
  
protected:
  void setInitK(const MatrixXd &K){

	const int r = K.rows();
	const int n = r*(r+1)/2;
	assert_eq(K.rows(),K.cols());
	VectorXd temptrlst(r*2+n);
	assert_eq(_rlst.size(),r*2);
	temptrlst.head(r*2) = _rlst;
	for (int i = 0; i < r; ++i)
	  for (int j = 0; j < r; ++j)
		temptrlst[r*2+symIndex(i,j)] = K(i,j);
	_rlst = temptrlst;
  }
  void optimizationBegin(){

	_energy.usePartialCon(false);
	resetEnergy();
	const SXMatrix &K = _K;
	_energy.setK(K);
	const SXMatrix Ak = makeEyeMatrix(_Ak);
	const SXMatrix Am = makeEyeMatrix(_Am);
	const SXMatrix D = Am + Ak.mul(_K);
	_energy.setDamping(D);
  }
  void optimizationEnd(){

	const int r = _energy.reducedDim();
	const int n = r*(r+1)/2;
	assert_eq(_rlst.size(),n+2*r);

	MatrixXd &K = _model->K;
	K.resize(r,r);
	for (int i = 0; i < r; ++i)
	  for (int j = 0; j < r; ++j)
		K(i,j) = _rlst[r*2+symIndex(i,j)];

	AkAmOptimizer::optimizationEnd();
  }
  const VSX &getVariable()const{
	return _AkAmKx;
  }
  bool getLowBound(vector<double> &lower)const{

	lower = vector<double>(getVariable().size(),-std::numeric_limits<double>::infinity());
	for (int i = 0; i < _AkAm.size(); ++i){
	  lower[i] = 0.0f;
	}
	return true;
  }

private:
  VSX _kx;
  SXMatrix _K;
  VSX _AkAmKx;
};

class AtAAkAmOptimizer:public MtlOptimizer{

public:
  AtAAkAmOptimizer(pMtlDataModel m):MtlOptimizer(m){
	
	const int r = _model->K.rows();
	const VSX a = makeSymbolic(r*r,"a");
	const VSX ak = makeSymbolic(r,"ak");
	const VSX am = makeSymbolic(r,"am");
	_a_ak_am.clear();
	_a_ak_am.insert( _a_ak_am.end(),a.begin(),a.end() );
	_a_ak_am.insert( _a_ak_am.end(),ak.begin(),ak.end() );
	_a_ak_am.insert( _a_ak_am.end(),am.begin(),am.end() );

	_Ak = makeEyeMatrix( ak );
	_Am = makeEyeMatrix( am );
	const SXMatrix A = convert( a,r );
	_K = CasADi::trans(A).mul(A);
	setInitVal(m->K,m->Ak,m->Am);
  }

protected:
  void setInitVal(const MatrixXd &K,const VectorXd &Ak,const VectorXd &Am){

	const int r = Ak.size();
	assert_gt(r,0);
	assert_eq(Am.size(), r);
	assert_eq(K.rows(),r);
	assert_eq(K.cols(),r);

	_rlst.resize(r*r+r*2);
	MatrixXd KK(r,r);
	KK.setZero();
	for (int i = 0; i < r; ++i)
	  KK(i,i) = sqrt( K(i,i) );
	_rlst.head(r*r) = Map<VectorXd>(&KK(0,0),r*r);
	_rlst.segment(r*r,r) = Ak;
	_rlst.tail(r) = Am;
  }
  void optimizationBegin(){

	_energy.usePartialCon(false);
	resetEnergy();
	const SXMatrix &K = _K;
	_energy.setK(K);
	const SXMatrix D = _Am + _Ak.mul(K);
	_energy.setDamping(D);
  }
  void optimizationEnd(){

	const int r = _energy.reducedDim();
	assert_eq(_rlst.size(),r*r+r*2);
	const MatrixXd A = Map<MatrixXd>(&_rlst[0],r,r);
	_model->K = A.transpose()*A;
	_model->Ak = Map<VectorXd>(&_rlst[r*r],r);
	_model->Am = Map<VectorXd>(&_rlst[r*r+r],r);
	_model->D = _model->Am.asDiagonal();
	_model->D += _model->Ak.asDiagonal()*_model->K;
  }
  const VSX &getVariable()const{
	return _a_ak_am;
  }
  bool getLowBound(vector<double> &lower)const{
	const int r = _energy.reducedDim();
	lower = vector<double>(r*r+2*r,-std::numeric_limits<double>::infinity());
	for (int i = r*r; i < lower.size(); ++i)
	  lower[i] = 0.0f;
	return true;
  }

private:
  VSX _a_ak_am;
  SXMatrix _K;
  SXMatrix _Ak;
  SXMatrix _Am;
};

class AtAakamOptimizer:public MtlOptimizer{

public:
  AtAakamOptimizer(pMtlDataModel m):MtlOptimizer(m){
	
	const int r = _model->K.rows();
	const VSX a = makeSymbolic(r*r,"a");
	const SX ak("ak");
	const SX am("am");
	_a_ak_am.clear();
	_a_ak_am.insert( _a_ak_am.end(),a.begin(),a.end() );
	_a_ak_am.push_back(ak);
	_a_ak_am.push_back(am);

	VSX Ak,Am;
	for (int i = 0; i < r; ++i){
	  Ak.push_back(ak);
	  Am.push_back(am);
	}
	
	_Ak = CASADI::makeEyeMatrix(Ak);
	_Am = CASADI::makeEyeMatrix(Am);
	const SXMatrix A = convert( a,r );
	_K = CasADi::trans(A).mul(A);

	setInitVal(m->K,m->Ak[0],m->Am[0]);
  }

protected:
  void setInitVal(const MatrixXd &K,const double ak,const double &am){

	const int r = K.rows();
	assert_gt(r,0);
	assert_eq(K.rows(),r);
	assert_eq(K.cols(),r);

	_rlst.resize(r*r+2);
	MatrixXd KK(r,r);
	KK.setZero();
	for (int i = 0; i < r; ++i)
	  KK(i,i) = sqrt( K(i,i) );
	_rlst.head(r*r) = Map<VectorXd>(&KK(0,0),r*r);
	_rlst[r*r] = ak;
	_rlst[r*r+1] = am;
  }
  void optimizationBegin(){

	_energy.usePartialCon(true);
	resetEnergy();
	const SXMatrix &K = _K;
	_energy.setK(K);
	const SXMatrix D = _Am + _Ak.mul(K);
	_energy.setDamping(D);
  }
  void optimizationEnd(){

	const int r = _energy.reducedDim();
	assert_eq(_rlst.size(),r*r+2);
	const MatrixXd A = Map<MatrixXd>(&_rlst[0],r,r);
	_model->K = A.transpose()*A;
	_model->Ak.resize(r);
	_model->Am.resize(r);
	for (int i = 0; i < r; ++i){
	  _model->Ak[i] = _rlst[r*r];
	  _model->Am[i] = _rlst[r*r+1];
	}
	_model->D = _model->Am.asDiagonal();
	_model->D += _model->Ak.asDiagonal()*_model->K;
  }
  const VSX &getVariable()const{
	return _a_ak_am;
  }
  bool getLowBound(vector<double> &lower)const{
	const int r = _energy.reducedDim();
	lower = vector<double>(r*r+2,-std::numeric_limits<double>::infinity());
	for (int i = r*r; i < lower.size(); ++i)
	  lower[i] = 0.0f;
	return true;
  }

private:
  VSX _a_ak_am;
  SXMatrix _K;
  SXMatrix _Ak;
  SXMatrix _Am;
};

#endif /* _MTLOPTENERGYAD_H_ */
