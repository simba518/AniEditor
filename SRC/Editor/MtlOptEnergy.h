#ifndef _MTLOPTENERGY_H_
#define _MTLOPTENERGY_H_

#include <boost/shared_ptr.hpp>
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
  
  class CtrlForceEnergy: public BaseFunGrad{
	
  public:
	void setRedWarper(pRedRSWarperExt warper){
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
	  for (int i = 0; i < T; ++i){
		_keyfIndex[i] = -1;
		_confIndex[i] = -1;
	  }
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
	void precompute();

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
	  _conFrames = conF;
	  _conNodes = conN;
	  _uc = uc;
	}
	void clearPartialCon(){
	  _conFrames.clear();
	  _conNodes.clear();
	  _uc.clear();
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
	double fun(const double *x);
	void grad(const double *x,double *g);
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
	void forward(const VectorXd &w,MatrixXd &V,MatrixXd &Z);
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
	void adjoint(MatrixXd &r)const;
	void compute_pEpz(MatrixXd &pEpz)const;
	
  private:
	pRedRSWarperExt _warper;
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

}//end of namespace

#endif /*_MTLOPTENERGY_H_*/
