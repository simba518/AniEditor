#ifndef _MTLOPT_H_
#define _MTLOPT_H_

#include <eigen3/Eigen/Dense>
#include <CASADITools.h>
using namespace Eigen;
using namespace CASADI;
using CasADi::SX;
using CasADi::SXMatrix;

class RedSpaceTimeEnergyAD{
  
public:
  void setT(const int T){_T = T;}
  template<class SCALAR> 
  void setTimestep(const SCALAR h){_h = h;}
  template<class SCALAR>
  void setDamping(const SCALAR ak,const SCALAR am){ _alphaK=ak ; _alphaM=am;}
  void setK(const SXMatrix &K){ assert_eq(_K.size1(),_K.size2());  _K = K;}
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
	  const int k = isKeyframe(i);
	  if(k >= 0){
		CASADI::convert((VectorXd)(_keyZ.col(k)),z[i]);
	  }else{
		const VSX zi(vz.begin()+r*i,vz.begin()+r*(i+1));
		assert_eq(zi.size(),r);
		z[i] = zi;
		_varZ.insert(_varZ.end(),zi.begin(),zi.end());
	  }
	}

	const SXMatrix I = makeEyeMatrix(VectorXd::Ones(r));
	const SXMatrix D = I*_alphaM + _K*_alphaK;
	_energy = 0;
	for (int i = 1; i < T-1; ++i){
	  const SXMatrix za = (z[i+1]-z[i]*2.0f+z[i-1])/(_h*_h);
	  const SXMatrix zv = D.mul(z[i+1]-z[i])/(_h);
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

protected:
  int isKeyframe(const int f)const{
	int k = -1;
	for (int i = 0; i < _keyId.size(); ++i){
	  if(f == _keyId[i]){
		k = i;
		break;
	  }
	}
	return k;
  }
  
private:
  int _T;
  SX _h;
  SX _alphaK;
  SX _alphaM;
  SXMatrix _K;
  MatrixXd _keyZ;
  VectorXi _keyId;
  SXMatrix _energy;
  VSX _varZ;
};

#endif /* _MTLOPT_H_ */
