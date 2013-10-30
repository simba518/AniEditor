#ifndef _MASIMULATORAD_H_
#define _MASIMULATORAD_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <symbolic/sx/sx.hpp>
#include <assertext.h>
using CasADi::SXMatrix;
using CasADi::SX;
using namespace Eigen;
using namespace std;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class MASimulatorAD modal-analysis simulator using auto-diff library.
   * 
   */
  class MASimulatorAD{

  public:
	// set parameters
	void setTimeStep(const SX &h){
	  _h = h;
	}
 	template<typename VECTOR>
	void setEigenValues(const VECTOR &eval){
	  _lambda.resize(eval.size());
	  for (int i = 0; i < (int)eval.size(); ++i)
		_lambda[i] = eval[i];
	}
	void setEigenValue(const SX eval){
	  // only one mode.
	  _lambda.resize(1);
	  _lambda[0] = eval;
	}
 	template<typename VECTOR>
	void setStiffnessDamping(const VECTOR &s){
	  _alphaK.resize(s.size());
	  for (int i = 0; i < (int)s.size(); ++i)
		_alphaK[i] = s[i];
	}
 	template<typename VECTOR>
 	void setMassDamping(const VECTOR &s){
	  _alphaM.resize(s.size());
	  for (int i = 0; i < (int)s.size(); ++i)
		_alphaM[i] = s[i];
	}
	void setStiffnessDamping(const SX s){
	  assert_gt(_lambda.size(),0);
	  _alphaK.resize(_lambda.size());
	  for (size_t i = 0; i < _lambda.size(); ++i)
		_alphaK[i] = s;
	}
 	void setMassDamping(const SX s){
	  assert_gt(_lambda.size(),0);
	  _alphaM.resize(_lambda.size());
	  for (size_t i = 0; i < _lambda.size(); ++i)
		_alphaM[i] = s;
	}
	void setStiffnessDamping(const double s){
	  setStiffnessDamping(SX(s));
	}
 	void setMassDamping(const double s){
	  setMassDamping(SX(s));
	}
	void setStiffnessDamping(const float s){
	  setStiffnessDamping(SX(s));
	}
 	void setMassDamping(const float s){
	  setMassDamping(SX(s));
	}
 	template<typename VECTOR>
	void setIntialStatus(const VECTOR &v0,const VECTOR &z0){
	  assert_eq(v0.size(),z0.size());
	  assert_eq(v0.size(),reducedDim());
	  _v.resize(v0.size());
	  _z.resize(z0.size());
	  for (int i = 0; i < v0.size(); ++i){
		_v[i] = v0[i];
		_z[i] = z0[i];
	  }
	}
	void removeVelecity(){
	  for (size_t i = 0; i < _v.size(); ++i)
		_v[i] = 0;
	}

	const SX &getTimeStep(){
	  return _h;
	}
	const vector<SX> &getStiffnessDamping(){
	  return _alphaK;
	}
 	const vector<SX> &getMassDamping(){
	  return _alphaM;
	}
	const vector<SX> &getEigenValues(){
	  return _lambda;
	}

	// simulate
	template<typename VECTOR>
	void forward(const VECTOR &w){

	  const int r = reducedDim();
	  assert_eq((int)w.size(),r);
	  const vector<SX> v0 = _v;
	  const vector<SX> z0 = _z;
	  SX G[2][2],S[2];
	  for (int i = 0; i < r; ++i){
		getImpIntegMat(G,S,i);
		_v[i] = G[0][0]*v0[i] + G[0][1]*z0[i] + S[0]*w[i];
		_z[i] = G[1][0]*v0[i] + G[1][1]*z0[i] + S[1]*w[i];
	  }
	}

	template<typename VECTOR>
	void forward(const vector<VECTOR> &ws,vector<VECTOR>&V,vector<VECTOR>&Z){

	  V.resize(ws.size()+1);
	  Z.resize(ws.size()+1);
	  getV(V[0]);
	  getZ(Z[0]);
	  for (int i = 0; i < ws.size(); ++i){
		this->forward(ws[i]);
		getV(V[i+1]);
		getZ(Z[i+1]);
	  }
	}

	// get data
	int reducedDim()const{
	  return _lambda.size();
	}
	const vector<SX> &getV()const{
	  return _v;
	}
	const vector<SX> &getZ()const{
	  return _z;
	}
	template<typename VECTOR> 
	void getV(VECTOR &v)const{
	  v.resize(_v.size());
	  for (int i = 0; i < v.size(); ++i) 
		v[i] = _v[i].getValue();
	}
	template<typename VECTOR>
	void getZ(VECTOR &z)const{
	  z.resize(_z.size());
	  for (int i = 0; i < z.size(); ++i) 
		z[i] = _z[i].getValue();
	}
	const VectorXd getEigenV()const{
	  VectorXd v(_v.size());
	  for (int i = 0; i < v.size(); ++i)
		v[i] = _v[i].getValue();
	  return v;
	}
	const VectorXd getEigenZ()const{
	  VectorXd z(_z.size());
	  for (int i = 0; i < z.size(); ++i)
		z[i] = _z[i].getValue();
	  return z;
	}

	void getImpIntegMat(SX G[2][2],SX S[2],const int mode_id)const{

	  assert_in(mode_id,0,reducedDim()-1);
	  assert_eq((int)_lambda.size(),reducedDim());
	  assert_eq((int)_alphaM.size(),reducedDim());
	  assert_eq((int)_alphaK.size(),reducedDim());

	  const SX lambda = _lambda[mode_id];
	  const SX alpha_k = _alphaK[mode_id];
	  const SX alpha_m = _alphaM[mode_id];
	  const SX gamma = 1/( (_h*alpha_k+_h*_h)*lambda + _h*alpha_m + 1);
	  G[0][0] = gamma;
	  G[0][1] = -_h*lambda*gamma;
	  G[1][0] = _h*gamma;
	  G[1][1] = 1-_h*_h*lambda*gamma;
	  S[0] = _h*gamma;
	  S[1] = _h*_h*gamma;
	}
	
  private:	
	SX _h;  //time step.	
	vector<SX> _alphaK; // stiffness damping.
	vector<SX> _alphaM; // mass damping.
	vector<SX> _lambda; // eigen values.
	vector<SX> _v; // velocity.
	vector<SX> _z; // reduced coordinates.
  };
  
  typedef boost::shared_ptr<MASimulatorAD> pMASimulatorAD;
  
}//end of namespace

#endif /*_MASIMULATORAD_H_*/
