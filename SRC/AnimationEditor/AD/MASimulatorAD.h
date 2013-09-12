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
	void setStiffnessDamping(const vector<SX> &s){
	  _alphaK = s;
	}
 	void setMassDamping(const vector<SX> &s){
	  _alphaM = s;
	}
	void setEigenValues(const vector<SX> &eval){
	  _lambda = eval;
	}
	void setEigenValues(const VectorXd &eval){
	  _lambda.resize(eval.size());
	  for (int i = 0; i < eval.size(); ++i){
		_lambda[i] = eval[i];
	  }
	}
	void setIntialStatus(const VectorXd &v0,const VectorXd &z0){
	  assert_eq(v0.size(),z0.size());
	  assert_eq(v0.size(),reducedDim());
	  _v.resize(v0.size());
	  _z.resize(z0.size());
	  for (int i = 0; i < v0.size(); ++i){
		_v[i] = v0[i];
		_z[i] = z0[i];
	  }
	}
	void setIntialStatus(const vector<SX> &v0,const vector<SX> &z0){
	  assert_eq(v0.size(),z0.size());
	  assert_eq(v0.size(),reducedDim());
	  _v = v0;
	  _z = z0;
	}
	void removeVelecity(){
	  for (size_t i = 0; i < _v.size(); ++i){
		_v[i] = 0;
	  }
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
	void forward(const vector<SX> &w){

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

  protected:
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
