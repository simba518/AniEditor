#ifndef _MTLOPT_H_
#define _MTLOPT_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <CASADITools.h>
using namespace std;
using namespace Eigen;
using namespace CASADI;
using CasADi::SX;
using CasADi::SXMatrix;

class RedSpaceTimeEnergyAD{
  
public:
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

  void assembleEnergy(){

  	const int T = _T;
  	const int r = reducedDim();
  	const VSX vz = makeSymbolic(T*r,"z");
  	VMatSX z(T);
  	_varZ.clear();
  	for (int i = 0; i < T; ++i){
  	  const int k = isKeyframe(_keyId,i);
  	  if(k >= 0){
  		CASADI::convert((VectorXd)(_keyZ.col(k)),z[i]);
  	  }else{
  		const VSX zi(vz.begin()+r*i,vz.begin()+r*(i+1));
  		assert_eq(zi.size(),r);
  		z[i] = zi;
  		_varZ.insert(_varZ.end(),zi.begin(),zi.end());
  	  }
  	}

  	_energy = 0;
  	for (int i = 1; i < T-1; ++i){
  	  const SXMatrix za = (z[i+1]-z[i]*2.0f+z[i-1])/(_h*_h);
  	  const SXMatrix zv = _D.mul(z[i+1]-z[i])/(_h);
  	  const SXMatrix diff = za+zv+_K.mul(z[i]);
  	  for (int j = 0; j < r; ++j)
  		_energy += diff.elem(j,0)*diff.elem(j,0);
  	}
  	_energy = _energy/2;
  }

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
  int _T;
  SX _h;
  SXMatrix _D;
  SXMatrix _K;
  MatrixXd _keyZ;
  VectorXi _keyId;
  SXMatrix _energy;
  VSX _varZ;
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
  void setKeyframes(const MatrixXd &Z, const VectorXi &kid){
	this->keyZ = Z;
	this->keyId = kid;
  }
  void setSubZ(const MatrixXd &Z){
	subZ = Z;
  }
  void print()const{

	EigenSolver<MatrixXd> eigenK(K);
	cout<< "eigen(K): " << eigenK.eigenvalues().transpose() << endl;

	EigenSolver<MatrixXd> eigenD(D);
	cout<< "eigen(D): " << eigenD.eigenvalues().transpose() << endl;
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
};

class MtlOptimizer{
  
 public:
  MtlOptimizer(MtlDataModel &m):_model(m){}
  void setInitialValue(const VectorXd &x0){
	_rlst = x0;
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
  MtlDataModel &_model;
  RedSpaceTimeEnergyAD _energy;
  CasADi::SXFunction _fun;
  CasADi::IpoptSolver _solver;
  VectorXd _rlst;
};

class ZOptimizer:public MtlOptimizer{
  
public:
  ZOptimizer(MtlDataModel &m):MtlOptimizer(m){
	setInitZ(m.subZ);
  }

protected:
  void setInitZ(const MatrixXd &Z){
	assert_gt(Z.size(),0);
	MatrixXd zz = Z;
	_rlst = Map<VectorXd>( &(zz(0,0)),zz.size() );
  }
  void optimizationBegin(){
	resetEnergy(false);
  }
  void optimizationEnd(){
	const int r = _energy.reducedDim();
	const int subT = _model.T - _model.keyId.size();
	assert_eq(r*subT, _rlst.size());
	_model.subZ = Map<MatrixXd>(&_rlst[0],r,subT);
	_model.Z = _model.subZ;
	_energy.insertKeyframes(_model.Z);
  }
  const VSX &getVariable()const{
	return _energy.getVarZ();
  }

};

class KOptimizer:public MtlOptimizer{

public: 
  KOptimizer(MtlDataModel &m,const bool useAkAm = true):MtlOptimizer(m),_useAkAm(useAkAm){

	setInitK(m.K);
	produceSymetricMat("k",_model.K.rows(),_kx,_K);
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

	resetEnergy();
	const SXMatrix &K = _K;
	_energy.setK(K);
	if(_useAkAm){
	  const SXMatrix Ak = makeEyeMatrix(convert(_model.Ak));
	  const SXMatrix Am = makeEyeMatrix(convert(_model.Am));
	  const SXMatrix D = Am + Ak.mul(K);
	  _energy.setDamping(D);
	}else{
	  _energy.setDamping(_model.D);
	}
  }
  void optimizationEnd(){

	const int r = _energy.reducedDim();
	const int n = r*(r+1)/2;
	assert_eq(_rlst.size(),n);

	MatrixXd &K = _model.K;
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

class AkAmOptimizer:public MtlOptimizer{

public:
  AkAmOptimizer(MtlDataModel &m):MtlOptimizer(m){
	
	setInitAkAm(m.Ak, m.Am);
	const int r = _model.K.rows();
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

	resetEnergy();
	const int r = _energy.reducedDim();
	const SXMatrix K = convert(_model.K);
	const SXMatrix Ak = makeEyeMatrix(_Ak);
	const SXMatrix Am = makeEyeMatrix(_Am);
	const SXMatrix D = Am + Ak.mul(K);
	_energy.setDamping(D);
  }
  virtual void optimizationEnd(){

	const int r = _energy.reducedDim();
	const MatrixXd &K = _model.K;
	assert_ge(_rlst.size(), r*2);
	const VectorXd &Ak = _rlst.head(r);
	const VectorXd &Am = _rlst.segment(r,r);
	MatrixXd D = Am.asDiagonal();
	D = D + Ak.asDiagonal()*K;
	_model.D = D;
	_model.Ak = Ak;
	_model.Am = Am;
  }
  const VSX &getVariable()const{
	return _AkAm;
  }

protected:
  VSX _Ak;
  VSX _Am;
  VSX _AkAm;
};

class KAtAmOptimizer:public AkAmOptimizer{

public: 
  KAtAmOptimizer(MtlDataModel &m):AkAmOptimizer(m){

	setInitK(m.K);
	produceSymetricMat("k",_model.K.rows(),_kx,_K);
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

	MatrixXd &K = _model.K;
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

#endif /* _MTLOPT_H_ */
