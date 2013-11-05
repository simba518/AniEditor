#ifndef _MTLOPTENERGYADJ_H_
#define _MTLOPTENERGYADJ_H_

#include <BaseMtlOptEnergy.h>

namespace LSW_ANI_EDITOR{
  
  template<typename WARPER_POINTER>
  class CtrlForceEnergyAdjT:public BaseCtrlForceEnergy,public BaseOptFun{
	
  public:
	void setRedWarper(WARPER_POINTER warper){
	  _warper = warper;
	}
	void setKD(const MatrixXd &K,const MatrixXd &D){

	  BaseCtrlForceEnergy::setKD(K,D);

	  // rotate w,z0,v0,_U
	  SelfAdjointEigenSolver<MatrixXd> eigenK(K);
	  const MatrixXd &U = eigenK.eigenvectors();
	  assert_eq(_w.size(),dim());
	  assert_gt(_w.size(),0);
	  MatrixXd W = U.transpose()*Map<MatrixXd>(&_w[0],reducedDim(),getT()-1);
	  _w = Map<VectorXd>(&W(0,0),W.size());
	  _z0 = U.transpose()*_z0;
	  _v0 = U.transpose()*_v0;
	  for (int i = 0; i < _keyZ.size(); ++i){
		_keyZ[i] = U.transpose()*_keyZ[i];
	  }
	}
	void setZ0(const VectorXd &z0){
	  assert_eq(z0.size(),_U.rows());
	  _z0 = _U.transpose()*z0;
	}
	void setV0(const VectorXd &v0){
	  assert_eq(v0.size(),_U.rows());
	  _v0 = _U.transpose()*v0;
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
		  const VectorXd z = _U*_Z.col(f);
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
			const VectorXd z = _U*_Z.col(f);
			_warper->jacobian(z,f,_conNodes[i],_uc[i],g);
			pEpzCon.col(i) = _U.transpose()*g;
		  }
		  pEpzCon *= _penaltyCon;
		}
	  }
	}
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
	
  private:
	WARPER_POINTER _warper;
	VectorXd _w;
	VectorXd _z0;
	VectorXd _v0;

	Matrix<double,-1,2> _G;
	Matrix<double,-1,1> _S;

	// tempt data
	MatrixXd _V, _Z, _pEpzCon, _pEpzKey, _r;
  };

  typedef CtrlForceEnergyAdjT<pRedRSWarperExt> CtrlForceEnergyAdj;
  typedef boost::shared_ptr<CtrlForceEnergyAdj> pCtrlForceEnergyAdj;
  typedef CtrlForceEnergyAdjT<pNoWarp> CtrlForceEnergyAdjNoWarp;
  typedef boost::shared_ptr<CtrlForceEnergyAdjNoWarp> pCtrlForceEnergyAdjNoWarp;

  class MtlOptEnergyAdj: public BaseMtlOptEnergy{
	
  public:
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
	double fun(const double *x){
	  computeKD(x,_K,_D);
	  const double E = (_D*_hV1+_K*_hZ1+_V1_V0).norm();
	  return E*E*0.5f/(_h*_h);
	}
	void grad(const double *x,double *g){

	  const int r = reducedDim();
	  const double &ak = x[0];
	  const double &am = x[1];
	  const MatrixXd &A = Map<MatrixXd>(const_cast<double*>(&x[2]),r,r);
	  _K = A.transpose()*A;

	  {// grad ak, am
		const MatrixXd F1 = _K*_hZ1+_V1_V0;
		MatrixXd Fm = ak*_K*_hV1+F1;
		MatrixXd Gm = _K*_hV1;
		MatrixXd Hm = am*_hV1+F1;

		const VectorXd &E = Map<VectorXd>(&_hV1(0,0),_hV1.size());
		const VectorXd &F = Map<VectorXd>(&Fm(0,0),Fm.size());
		const VectorXd &G = Map<VectorXd>(&Gm(0,0),Gm.size());
		const VectorXd &H = Map<VectorXd>(&Hm(0,0),Hm.size());
		
		g[1] = (E.transpose()*(am*E+F))(0,0)/(_h*_h); // grad am
		g[0] = (G.transpose()*(ak*G+H))(0,0)/(_h*_h); // grad ak
	  }

	  {// grad A.
		const MatrixXd &E = _hV1;
		const MatrixXd &At = A.transpose();
		const double &J = ak;
		const MatrixXd &L = _hZ1;
		const MatrixXd M = am*_hV1+_V1_V0;
		const MatrixXd &Et = E.transpose();
		const double &Jt = J;
		const MatrixXd &Lt = L.transpose();
		const MatrixXd &Mt = M.transpose();
		MatrixXd gA = (At*A*L*Lt*At).transpose()+(A*L*Et*At*A*Jt)+(A*L*Mt)+
		  (A*At*A*L*Lt)+(L*Et*At*A*Jt*At).transpose()+(A*M*Lt)+
		  (A*Jt*At*A*L*Et)+(E*Et*At*A*Jt*J*At).transpose()+(A*Jt*M*Et)+
		  (Jt*At*A*L*Et*At).transpose()+(A*E*Et*At*A*Jt*J)+(A*E*Mt*J);
		gA /= (_h*_h);
		memcpy(&g[2],&gA(0,0),sizeof(double)*r*r);
	  }
	}

  protected:
	MatrixXd _hV1;
	MatrixXd _hZ1;
	MatrixXd _V1_V0;
  };
  typedef boost::shared_ptr<MtlOptEnergyAdj> pMtlOptEnergyAdj;

}//end of namespace

#endif /* _MTLOPTENERGYADJ_H_ */
