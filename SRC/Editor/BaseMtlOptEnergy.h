#ifndef _BASEMTLOPTENERGY_H_
#define _BASEMTLOPTENERGY_H_

#include <boost/shared_ptr.hpp>
#include <CASADITools.h>
#include <RedRSWarperExt.h>
#include <assertext.h>
#include <BaseOptFun.h>
#include <NoWarp.h>
using namespace LSW_WARPING;

namespace LSW_ANI_EDITOR{

  class BaseCtrlForceEnergy{
	
  public:
	void setTimestep(const double h){
	  assert_gt(h,0.0f);
	  _h = h;
	}
	void setTotalFrames(const int T);
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
	  precompute();
	}
	void setKD(const MatrixXd &K,const MatrixXd &D);

	void addKeyframe(const VectorXd &unRotZk, const int f);
	void setKeyframes(const vector<VectorXd> &keyZ, const vector<int> &keyframes);
	void setPartialCon(vector<int>&conF,vector<vector<int> >&conN,vector<VectorXd>&uc);
	void clearPartialCon();
	void clearKeyframes();

	int reducedDim()const{
	  return _Lambda.size();
	}
	int getT()const{
	  return _T;
	}
	bool isKeyframe(const int f)const{
	  return f>=0 && f < _T && _keyfIndex[f] > 0;
	}
	const vector<int> &getKeyframe()const{
	  return _keyframes;
	}
	const vector<VectorXd> &getKeyZ()const{
	  return _keyZ;
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
	
  protected:
	virtual void precompute(){
	
	}
	
  protected:
	double _h;
	int _T;
	double _penaltyCon;
	double _penaltyKey;
	VectorXd _Lambda; 
	VectorXd _diagD;
	MatrixXd _U; // K=U^T*Lambda*U
	double _objValue;

	// constraints
	vector<int> _keyframes;
	vector<int> _keyfIndex;
	vector<VectorXd> _keyZ;
	vector<int> _confIndex;
	vector<int> _conFrames;
	vector<vector<int> > _conNodes;
	vector<VectorXd> _uc;
  };

  class BaseMtlOptEnergy:public BaseOptFun{

  public:
	void setTotalFrames(const int T){
	  assert_ge(T,3);
	  _T = T;
	}
	void setTimestep(const double h){
	  assert_gt(h,0.0f);
	  _h = h;
	}
	void setMtl(const VectorXd &Lambda,const double ak,const double am);

	int dim()const{
	  const int r = reducedDim();
	  assert_gt(r,0);
	  return r*r+2;
	}
	void init(double *x,const int n){
	  assert_eq(n,dim());
	  assert_gt(n,0);
	  memcpy(x,&_AkAmA[0],n*sizeof(double));
	}
	void bounds(double *x_l,double *x_u,const int n){
	  BaseOptFun::bounds(x_l,x_u,n);
	  x_l[0] = 0.0f;
	  x_l[1] = 0.0f;
	}
	void setRlst(const double *x, const double objValue){
	  _objValue = objValue;
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
	  return _K.rows();
	}
	const VectorXd &getX()const{
	  return _AkAmA;
	}
	const double getObjValue()const{
	  return _objValue;
	}

  protected:
	void computeKD(const double *x,MatrixXd &K,MatrixXd &D)const;

  protected:
	double _objValue;
	int _T;
	double _h;
	VectorXd _AkAmA;
	MatrixXd _K;
	MatrixXd _D;
  };
  
}

#endif /* _BASEMTLOPTENERGY_H_ */
