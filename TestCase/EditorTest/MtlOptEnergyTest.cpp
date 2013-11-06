#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MtlOptEnergy.h>
#include <MtlOptEnergyAD.h>
using namespace Eigen;
using namespace LSW_ANI_EDITOR;

class MtlOptDataModeTest{
  
public:
  MtlOptDataModeTest(){
	r = 2;
	lambda = VectorXd::Ones(r)*10.0f+VectorXd::Random(r);
	ak = 0.1f;
	am = 0.2f;
	penC = 220;
	penK = 1220;
	T = 5;
	h = 0.15f;
	Z = MatrixXd::Random(r,T);
	keyframe.push_back(0);
	keyframe.push_back(T/2);
	keyframe.push_back(T-1);
	Kz = MatrixXd::Random(r,keyframe.size())*1.5f;
  }
  void clearKeyframes(){
	keyframe.clear();
	Kz.resize(0,0);
  }
  void init(MtlOptEnergy &fun){
	fun.setTotalFrames(T);
	fun.setTimestep(h);
	fun.setMtl(lambda,ak,am);
	fun.setZ(Z);
  }
  void init(RedSpaceTimeEnergyAD &fun){

	fun.setT(T);
	fun.setTimestep(h);
	fun.setDamping(ak,am,lambda);
	fun.setK(lambda);
	fun.setKeyframes(Kz,keyframe);
  }
  void init(CtrlForceEnergy &fun){
	
	fun.setTotalFrames(T);
	fun.setTimestep(h);
	fun.setMtl(lambda,ak*lambda+am*VectorXd::Ones(lambda.size()));
	fun.setKeyframes(Kz,keyframe);
  }
  
public:
  VectorXd lambda;
  int r;
  double ak;
  double am;
  double penC;
  double penK;
  int T;
  double h;
  MatrixXd Z;
  MatrixXd Kz;
  vector<int> keyframe;
};

class MtlOptEnergyADTest{
  
public:
  MtlOptEnergyADTest(RedSpaceTimeEnergyAD &ad,const MatrixXd &Z){

	const int r = ad.reducedDim();
	const VSX a = makeSymbolic(r*r,"a");
	const SX ak("ak");
	const SX am("am");
	VSX _a_ak_am;
	_a_ak_am.push_back(ak);
	_a_ak_am.push_back(am);
	_a_ak_am.insert( _a_ak_am.end(), a.begin(), a.end() );
	VSX Ak,Am;
	for (int i = 0; i < r; ++i){
	  Ak.push_back(ak);
	  Am.push_back(am);
	}
	const SXMatrix _Ak = CASADI::makeEyeMatrix(Ak);
	const SXMatrix _Am = CASADI::makeEyeMatrix(Am);
	const SXMatrix A = convert( a,r );
	const SXMatrix _K = CasADi::trans(A).mul(A);
	const SXMatrix D = _Am + _Ak.mul(_K);
	ad.setDamping(D);
	ad.setK(_K);

	vector<int> keyframe;
	for (int i = 0; i < Z.cols(); ++i){
	  keyframe.push_back(i);
	}
	ad.setKeyframes(Z,keyframe);

	ad.assembleEnergy();
	
	_fun = CasADi::SXFunction(_a_ak_am,ad.getEnergy());
	_fun.init();
	const SXMatrix g = _fun.jac();
	_grad = CasADi::SXFunction(_a_ak_am,g);
	_grad.init();
  }
  double fun(const VectorXd &x){
	VectorXd rlst;
	CASADI::evaluate(_fun,x,rlst);
	return rlst[0];
  }
  VectorXd grad(const VectorXd &x){
	VectorXd rlst;
	CASADI::evaluate(_grad,x,rlst);
	return rlst;
  }

private:
  CasADi::SXFunction _fun;
  CasADi::SXFunction _grad;
};

class MtlOptForcesADTest{
  
public:
  MtlOptForcesADTest(RedSpaceTimeEnergyAD &ad){

	ad.assembleEnergy();
	_fun = CasADi::SXFunction(ad.getVarZ(),ad.getEnergy());
	_fun.init();
	const SXMatrix g = _fun.jac();
	_grad = CasADi::SXFunction(ad.getVarZ(),g);
	_grad.init();
  }
  double fun(const VectorXd &x){
	VectorXd rlst;
	CASADI::evaluate(_fun,x,rlst);
	return rlst[0];
  }
  VectorXd grad(const VectorXd &x){
	VectorXd rlst;
	CASADI::evaluate(_grad,x,rlst);
	return rlst;
  }

private:
  CasADi::SXFunction _fun;
  CasADi::SXFunction _grad;
};

BOOST_AUTO_TEST_SUITE(MtlOptEnergyTest)

BOOST_AUTO_TEST_CASE(testMtlOpt){

  MtlOptDataModeTest data;
  MtlOptEnergy optfun;
  RedSpaceTimeEnergyAD optfunAD;

  data.init(optfun);
  data.init(optfunAD);
  
  const VectorXd x = VectorXd::Random(optfun.dim());
  const double obj1 = optfun.fun(&x[0]);
  VectorXd g(optfun.dim());
  optfun.grad(&x[0],&g[0]);
  
  MtlOptEnergyADTest enAD(optfunAD,data.Z);
  const double obj2 = enAD.fun(x);
  const VectorXd g1 = enAD.grad(x);
  
  ASSERT_EQ_TOL(obj1,obj2,1e-10);
  ASSERT_EQ(g.size(),g1.size());
  ASSERT_EQ_SMALL_VEC_TOL(g,g1,g.size(),1e-12);
}

BOOST_AUTO_TEST_CASE(testZOpt){

  MtlOptDataModeTest data;
  CtrlForceEnergy optfun;
  RedSpaceTimeEnergyAD optfunAD;

  data.init(optfun);
  data.init(optfunAD);
  
  const VectorXd x = VectorXd::Random(optfun.dim());
  const double obj1 = optfun.fun(&x[0]);
  VectorXd g(optfun.dim());
  optfun.grad(&x[0],&g[0]);
  
  MtlOptForcesADTest enAD(optfunAD);
  const double obj2 = enAD.fun(x);
  const VectorXd g2 = enAD.grad(x);
  
  ASSERT_EQ_TOL(obj1,obj2,1e-10);
  ASSERT_EQ(g.size(),g2.size());
  ASSERT_EQ_SMALL_VEC_TOL(g,g2,g.size(),1e-10);
}

BOOST_AUTO_TEST_CASE(testZOptNoKey){

  MtlOptDataModeTest data;
  data.clearKeyframes();
  CtrlForceEnergy optfun;
  RedSpaceTimeEnergyAD optfunAD;

  data.init(optfun);
  data.init(optfunAD);
  
  const VectorXd x = VectorXd::Random(optfun.dim());
  const double obj1 = optfun.fun(&x[0]);
  VectorXd g(optfun.dim());
  optfun.grad(&x[0],&g[0]);
  
  MtlOptForcesADTest enAD(optfunAD);
  const double obj2 = enAD.fun(x);
  const VectorXd g2 = enAD.grad(x);
  
  ASSERT_EQ_TOL(obj1,obj2,1e-10);
  ASSERT_EQ(g.size(),g2.size());
  ASSERT_EQ_SMALL_VEC_TOL(g,g2,g.size(),1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
