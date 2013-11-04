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
	void setPenaltyCon(const double pc,const double pk){
	  assert_ge(pc,0.0f);
	  assert_ge(pk,0.0f);
	  _penaltyCon = pc;
	  _penaltyKey = pk;
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
		_diagD[i] = diagD(i,i);
	  if (diagD.norm() > 0.0){
		assert_lt( (diagD.norm()-_diagD.norm())/diagD.norm(), 1e-4 );
	  }else{
		assert_lt( D.norm()-_diagD.norm(), 1e-8 );
	  }
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

	void addKeyframe(const VectorXd &zk, const int f){
	  assert_in(f,0,getT()-1);
	  assert_eq(zk.size(),reducedDim());
	  if (_keyfIndex[f] >= 0){
	  	_keyZ[ _keyfIndex[f] ] = zk;
	  }else{
		_keyZ.push_back(zk);
	  	_keyframes.push_back(f);
	  	_keyfIndex[f] = _keyframes.size()-1;
	  }
	}
	void setKeyframes(const vector<VectorXd> &keyZ, const vector<int> &keyframes){

	  clearKeyframes();
	  assert_eq(keyZ.size(),keyframes.size());
	  for (size_t i = 0; i < keyframes.size(); ++i)
		addKeyframe(keyZ[i],keyframes[i]);
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
	void clearKeyframes(){
	  _keyframes.clear();
	  _keyZ.clear();
	  _keyfIndex.resize(getT());
	  for (size_t i = 0; i < _keyfIndex.size(); ++i){
		_keyfIndex[i] = -1;
	  }
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
	  forward(w,_V,_Z);

	  {	// control forces
		objValue = w.norm();
		objValue = objValue*objValue;
	  }

	  { // keyframes
		for (size_t i = 0; i < _keyframes.size(); ++i){
		  const VectorXd &zf = _Z.col(_keyframes[i]);
		  const VectorXd &zk = _keyZ[i];
		  const double znorm = (zf-zk).norm();
		  objValue += _penaltyKey*znorm*znorm;
		}
	  }

	  { // partial constraints
		static VectorXd ui;
		for (size_t i = 0; i < _conFrames.size(); ++i){
		  const int f = _conFrames[i];
		  const VectorXd &z = _Z.col(f);
		  _warper->warp(z,f,_conNodes[i],ui);
		  assert_eq(_uc[i].size(),ui.size());
		  const double n = (_uc[i]-ui).norm();
		  objValue += _penaltyCon*n*n;
		}
	  }
	  return objValue*0.5f;
	}
	void grad(const double *x,double *g){
  
	  const VectorXd &w = Map<VectorXd>(const_cast<double*>(x),dim());
	  forward(w,_V,_Z);
	  compute_pEpz(_pEpzCon,_pEpzKey);
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
		  R(c*2+1,T-2) = _pEpzCon(c,kk);
	  }
	  if (_keyfIndex[T-1] >= 0){
		const int kk = _keyfIndex[T-1];
		for (int c = 0; c < r; ++c)
		  R(c*2+1,T-2) += _pEpzKey(c,kk);
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
			R(c*2+1,f) += _pEpzCon(c,kk);
		}
		if(_keyfIndex[f_1] >= 0){
		  const int kk = _keyfIndex[f_1];
		  for (int c = 0; c < r; ++c)
			R(c*2+1,f) += _pEpzKey(c,kk);
		}
	  }
	}
	void compute_pEpz(MatrixXd &pEpzCon,MatrixXd &pEpzKey)const{

	  { // keyframes
		const int r = reducedDim();
		pEpzKey.resize(r,_keyframes.size());
		for (size_t i = 0; i < _keyframes.size(); ++i){
		  const int f = _keyframes[i];
		  const VectorXd &zf = _Z.col(f);
		  pEpzKey.col(i) = zf-_keyZ[i];
		}
		pEpzKey *= _penaltyKey;
	  }

	  { // partial constraints
		static VectorXd g;
		if (_conFrames.size() > 0){
		  const int r = reducedDim();
		  pEpzCon.resize(r,_conFrames.size());
		  for (size_t i = 0; i < _conFrames.size(); ++i){
			const int f = _conFrames[i];
			const VectorXd &z = _Z.col(f);
			_warper->jacobian(z,f,_conNodes[i],_uc[i],g);
			pEpzCon.col(i) = g;
		  }
		  pEpzCon *= _penaltyCon;
		}
	  }
	}
	
  private:
	WARPER_POINTER _warper;
	double _h;
	int _T;
	double _penaltyCon;
	double _penaltyKey;
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
	MatrixXd _V, _Z, _pEpzCon, _pEpzKey, _r;
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
	void setMtl(const VectorXd &Lambda,const double ak,const double am){
	  assert_ge(ak,0.0f);
	  assert_ge(am,0.0f);
	  const int r = Lambda.size();
	  _AkAmA.resize(r+r+r*r);
	  _AkAmA.head(r) = ak*VectorXd::Ones(r);
	  _AkAmA.segment(r,r) = am*VectorXd::Ones(r);
	  _AkAmA.tail(r*r).setZero();
	  int c = 0;
	  for (int i = 0; i < r; ++i){
		assert_ge(Lambda[i],0.0f);
		_AkAmA[r*2+c] = sqrt(Lambda[i]);
		c += (r+1);
	  }
	}
	void setVZ(const MatrixXd &V, const MatrixXd &Z){
	  assert_eq(V.rows(),Z.rows());
	  assert_eq(V.cols(),Z.cols());
	  assert_gt(V.cols(),1);
	  assert_gt(_h,0.0f);
	  const int T = V.cols();
	  _hV1 = V.rightCols(T-1)*_h;
	  _hZ1 = Z.rightCols(T-1)*_h;
	  _V1_V0 = V.rightCols(T-1)-V.leftCols(T-1);
	}

	int dim()const{
	  const int r = reducedDim();
	  assert_gt(r,0);
	  return r*r+r*2;
	}
	void init(double *x,const int n){
	  assert_eq(n,dim());
	  assert_gt(n,0);
	  memcpy(x,&_AkAmA[0],n*sizeof(double));
	}
	double fun(const double *x){
	  computeKD(x,_K,_D);
	  const double E = (_D*_hV1+_K*_hZ1+_V1_V0).norm();
	  return E*E*0.5f;
	}
	void grad(const double *x,double *g){

	  const int r = reducedDim();
	  const VectorXd &ak = Map<VectorXd>(const_cast<double*>(x),r);
	  const VectorXd &am = Map<VectorXd>(const_cast<double*>(&x[r]),r);
	  const MatrixXd &A = Map<MatrixXd>(const_cast<double*>(&x[r*2]),r,r);
	  _K = A.transpose()*A;

	  {// grad ak, am
		const MatrixXd &E = _hV1;
		const MatrixXd F1 = _K*_hZ1+_V1_V0;
		const MatrixXd F = ak.asDiagonal()*_K*_hV1+F1;
		const MatrixXd G = _K*_hV1;
		const MatrixXd H = am.asDiagonal()*_hV1+F1;
		for (int i = 0; i < r; ++i){
		  const VectorXd &rowE = E.row(i);
		  const VectorXd &rowF = F.row(i);
		  const VectorXd &rowG = G.row(i);
		  const VectorXd &rowH = H.row(i);
		  g[i+r] = (rowE.transpose()*rowF+am[i]*rowE.transpose()*rowE)(0,0); // grad(am)
		  g[i] = (rowG.transpose()*rowH+ak[i]*rowG.transpose()*rowG)(0,0); // grad(ak)
		}
	  }

	  {// grad A.
		const MatrixXd &E = _hV1;
		const MatrixXd &At = A.transpose();
		const MatrixXd &J = ak.asDiagonal();
		const MatrixXd &L = _hZ1;
		const MatrixXd M = am.asDiagonal()*_hV1+_V1_V0;
		const MatrixXd &Et = E.transpose();
		const MatrixXd &Jt = J.transpose();
		const MatrixXd &Lt = L.transpose();
		const MatrixXd &Mt = M.transpose();
		const MatrixXd gA = (At*A*L*Lt*At).transpose()+(A*L*Et*At*A*Jt)+(A*L*Mt)+
		  (A*At*A*L*Lt)+(L*Et*At*A*Jt*At).transpose()+(A*M*Lt)+
		  (A*Jt*At*A*L*Et)+(E*Et*At*A*Jt*J*At).transpose()+(A*Jt*M*Et)+
		  (Jt*At*A*L*Et*At).transpose()+(A*E*Et*At*A*Jt*J)+(A*E*Mt*J);
		memcpy(&g[2*r],&gA(0,0),sizeof(double)*r*r);
	  }
	}
	void setRlst(const double *x, const double objValue){
	  memcpy(&_AkAmA[0],x,dim()*sizeof(double));
	  computeKD(x,_K,_D);
	}

	const MatrixXd &getK()const{
	  return _K;
	}
	const MatrixXd &getD()const{
	  return _D;
	}
	int reducedDim()const{
	  return _hZ1.rows();
	}
	const VectorXd &getX()const{
	  return _AkAmA;
	}

  protected:
	void computeKD(const double *x,MatrixXd &K,MatrixXd &D)const{
	  const int r = reducedDim();
	  const VectorXd &ak = Map<VectorXd>(const_cast<double*>(x),r);
	  const VectorXd &am = Map<VectorXd>(const_cast<double*>(&x[r]),r);
	  const MatrixXd &A = Map<MatrixXd>(const_cast<double*>(&x[r*2]),r,r);
	  K = A.transpose()*A;
	  D = am.asDiagonal();
	  D += ak.asDiagonal()*K;
	}

  protected:
	int _T;
	double _h;
	VectorXd _AkAmA;
	MatrixXd _hV1;
	MatrixXd _hZ1;
	MatrixXd _V1_V0;
	MatrixXd _K;
	MatrixXd _D;
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
