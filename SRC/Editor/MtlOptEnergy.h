#ifndef _MTLOPTENERGY_H_
#define _MTLOPTENERGY_H_

#include <boost/shared_ptr.hpp>
#include <CASADITools.h>
#include <RedRSWarperExt.h>
#include <assertext.h>
using namespace LSW_WARPING;

namespace LSW_ANI_EDITOR{

  class BaseFunGrad{
	
  public:
	virtual int dim()const = 0;
	virtual void init(double *x,const int n) = 0;
	virtual double fun(const double *x) = 0;
	virtual void grad(const double *x,double *g) = 0;
	virtual void setRlst(const double *x, const double objValue) = 0;
  };
  typedef boost::shared_ptr<BaseFunGrad> pBaseFunGrad;
  
  template<typename WARPER_POINTER>
  class CtrlForceEnergyT: public BaseFunGrad{
	
  public:
	void setRedWarper(WARPER_POINTER warper){
	  _warper = warper;
	}
	void setTimestep(const double h){
	  assert_gt(h,0.0f);
	  _h = h;
	}
	void setTotalFrames(const int T){

	  assert_ge(T,3);
	  _T = T;
	  _keyfIndex.resize(T);
	  for (int i = 0; i < T; ++i)
		_keyfIndex[i] = -1;
	  clearPartialCon();
	}
	void setPenaltyCon(const double pc){
	  assert_ge(pc,0.0f);
	  _penaltyCon = pc;
	}
	void setMtl(const VectorXd &Lambda,const VectorXd &diagD){
	  _Lambda = Lambda;
	  _diagD = diagD;
	  _U = MatrixXd::Identity(Lambda.size(),Lambda.size());
	}
	void setKD(const MatrixXd &K,const MatrixXd &D){

	  assert_eq(K.rows(),reducedDim());
	  assert_eq(K.cols(),reducedDim());
	  assert_eq(D.rows(),reducedDim());
	  assert_eq(D.cols(),reducedDim());

	  SelfAdjointEigenSolver<MatrixXd> eigenK(K);
	  _Lambda = eigenK.eigenvalues();
	  _U = eigenK.eigenvectors();
	  const MatrixXd diagD = _U.transpose()*D*_U;
	  _diagD.resize(reducedDim());
	  for (int i = 0; i < reducedDim(); ++i)
		_diagD[i] = D(i,i);
	  assert_lt( D.norm()-_diagD.norm(), 1e-8 );
	}
	void setZ0(const VectorXd &z0){
	  _z0 = z0;
	}
	void setV0(const VectorXd &v0){
	  _v0 = v0;
	}
	void precompute(){
  
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

	void setKeyframes(const VectorXd &zk, const int f){
	  if (_keyfIndex[f] >= 0){
	  	_keyZ[ _keyfIndex[f] ] = zk;
	  }else{
		_keyZ.push_back(zk);
	  	_keyframes.push_back(f);
	  	_keyfIndex[f] = _keyframes.size()-1;
	  }
	}
	void setPartialCon(vector<int>&conF,vector<vector<int> >&conN,vector<VectorXd>&uc){

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
	void clearPartialCon(){

	  _conFrames.clear();
	  _conNodes.clear();
	  _uc.clear();
	  _confIndex.resize(_T);
	  for (int i = 0; i < _T; ++i)
		_confIndex[i] = -1;
	}

	int dim()const{
	  const int n = (getT()-1)*reducedDim();
	  assert_gt(n,0);
	  return n;
	}
	void init(double *x,const int n){
	  assert_gt(n,0);
	  assert_eq(n,dim());
	  if (_w.size() != dim()){
		_w.resize(dim());
		_w.setZero();
	  }
	  memcpy(x,&_w[0],n*sizeof(double));
	}
	double fun(const double *x){
	  const VectorXd &w = Map<VectorXd>(const_cast<double*>(x),dim());
	  double objValue = 0.0f;
	  objValue = w.norm();
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
	void grad(const double *x,double *g){
  
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
	void setRlst(const double *x, const double objValue){
	  assert_gt(dim(),0);
	  assert_eq(_w.size(),dim());
	  memcpy(&_w[0],x,dim()*sizeof(double));
	  _objValue = objValue;
	}

	int reducedDim()const{
	  return _Lambda.size();
	}
	int getT()const{
	  return _T;
	}
	bool isKeyframe(const int f)const{
	  return f>=0 && f < _T && _keyfIndex[f] > 0;
	}
	const VectorXd &getLambda()const{
	  return _Lambda;
	}
	const MatrixXd &getU()const{
	  return _U;
	}
	const double getObjValue()const{
	  return _objValue;
	}
	void forward(const VectorXd &w,MatrixXd &V,MatrixXd &Z){

	  const int T = getT();
	  const int r = reducedDim();
	  assert_eq(_v0.size(),r);
	  assert_eq(_z0.size(),r);
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
	void forward(MatrixXd &V,MatrixXd &Z){
	  assert_gt(_w.size(),0);
	  forward(_w,V,Z);
	}

	const Matrix<double,-1,2> &getG()const{
	  return _G;
	}
	const Matrix<double,-1,1> &getS()const{
	  return _S;
	}
	
  protected:
	void adjoint(MatrixXd &R)const{
  
	  const int r = reducedDim();
	  const int T = getT();
	  R.resize(r*2,T-1);
	  R.col(T-2).setZero();
	  if (_confIndex[T-1] >= 0){
		const int kk = _confIndex[T-1];
		for (int c = 0; c < r; ++c)
		  R(c*2+1,T-2) = _pEpz(c,kk);
	  }

	  for (int f = T-3; f >= 0; --f){
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
		if (_confIndex[f_1] >= 0){
		  const int kk = _confIndex[f_1];
		  for (int c = 0; c < r; ++c)
			R(c*2+1,f) += _pEpz(c,kk);
		}
	  }
	}
	void compute_pEpz(MatrixXd &pEpz)const{

	  VectorXd g;
	  if (_conFrames.size() > 0){
		const int r = reducedDim();
		pEpz.resize(r,_conFrames.size());
		for (size_t i = 0; i < _conFrames.size(); ++i){
		  const int f = _conFrames[i];
		  const VectorXd &z = _Z.col(f);
		  _warper->jacobian(z,f,_conNodes[i],_uc[i],g);
		  pEpz.col(i) = g;
		}
		pEpz *= _penaltyCon;
	  }
	}
	
  private:
	WARPER_POINTER _warper;
	double _h;
	int _T;
	double _penaltyCon;
	VectorXd _Lambda; 
	VectorXd _diagD;
	MatrixXd _U; // K=U^T*Lambda*U
	double _objValue;
	VectorXd _w;
	VectorXd _z0;
	VectorXd _v0;

	Matrix<double,-1,2> _G;
	Matrix<double,-1,1> _S;

	// constraints
	vector<int> _keyframes;
	vector<int> _keyfIndex;
	vector<VectorXd> _keyZ;
	vector<int> _confIndex;
	vector<int> _conFrames;
	vector<vector<int> > _conNodes;
	vector<VectorXd> _uc;

	// tempt data
	MatrixXd _V, _Z, _pEpz, _r;
  };
  typedef CtrlForceEnergyT<pRedRSWarperExt> CtrlForceEnergy;
  typedef boost::shared_ptr<CtrlForceEnergy> pCtrlForceEnergy;

  class MtlOptEnergy: public BaseFunGrad{
	
  public:
	void setTotalFrames(const int T){
	  assert_ge(T,3);
	  _T = T;
	}
	void setTimestep(const double h){
	  assert_gt(h,0.0f);
	  _h = h;
	}
	void setMtl(const VectorXd &Lambda,const VectorXd &diagD){
	  
	}
	void setVZ(const MatrixXd &V, const MatrixXd &Z){
	  assert_eq(V.rows(),Z.rows());
	  assert_eq(V.cols(),Z.cols());
	  _V = V;
	  _Z = Z;
	}

	int dim()const{
	  const int r = reducedDim();
	  assert_gt(r,0);
	  return r*r+r*2;
	}
	void init(double *x,const int n){
	  
	}
	double fun(const double *x){
	  
	}
	void grad(const double *x,double *g){
	  
	}
	void setRlst(const double *x, const double objValue){
	  
	}

	MatrixXd &getK()const{
	  
	}
	MatrixXd &getD()const{
	  
	}
	int reducedDim()const{
	  return _Z.rows();
	}

  private:
	int _T;
	double _h;
	VectorXd _AkAmA;
	MatrixXd _V;
	MatrixXd _Z;
  };
  typedef boost::shared_ptr<MtlOptEnergy> pMtlOptEnergy;


  class NoWarp{
	// u(z) = W*z
  public:
	NoWarp(const MatrixXd &W){
	  this->W = W;
	}
	void warp(const VectorXd &z,int frame_id,const vector<int> &nodes,VectorXd &ui){
	  assert_eq(z.size(),W.cols());
	  ui = block(W,nodes)*z;
	}
	void jacobian(const VectorXd &z,int frame_id,const vector<int> &nodes,
				  const VectorXd &uc,VectorXd &g){
	  const MatrixXd Wb = block(W,nodes);
	  g = Wb.transpose()*(Wb*z-uc);
	}
	CasADi::SXMatrix warp(const CasADi::SXMatrix&z,const vector<int> &nodes){
	  assert_eq(z.size1(),W.cols());
	  assert_eq(z.size2(),1);
	  CasADi::SXMatrix sWi;  
	  CASADI::convert(this->block(W,nodes),sWi);
	  return sWi.mul(z);
	}
  
  protected:
	MatrixXd block(const MatrixXd &W,const vector<int> &nodes)const{

	  MatrixXd Wi(nodes.size()*3, W.cols());
	  for (size_t i = 0; i < nodes.size(); ++i){
		const int j3 = nodes[i]*3;
		assert_le(j3,W.rows()-3);
		Wi.block( i*3,0,3, W.cols() ) = W.block( j3,0,3,W.cols() );
	  }
	  return Wi;
	}
  
  private:
	MatrixXd W;
  };
  typedef boost::shared_ptr<NoWarp> pNoWarp;
  typedef CtrlForceEnergyT<pNoWarp> CtrlForceEnergyNoWarp;
  typedef boost::shared_ptr<CtrlForceEnergyNoWarp> pCtrlForceEnergyNoWarp;

}//end of namespace

#endif /*_MTLOPTENERGY_H_*/
