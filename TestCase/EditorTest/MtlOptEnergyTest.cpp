#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MtlOptEnergy.h>
#include <MtlOptIpopt.h>
#include <MatrixTools.h>
#include <MASimulatorAD.h>
using namespace Eigen;
using namespace LSW_ANI_EDITOR;

class TestIpoptObjfun:public BaseFunGrad{
  // (x0-2)^2 + (x1-1)^2
public:
  int dim()const{
	return 2;
  }
  void init(double *x,const int n){
	x[0] = 0;
	x[1] = 0;
  }
  double fun(const double *x){
	return (x[0]-2)*(x[0]-2) + (x[1]-1)*(x[1]-1);
  }
  void grad(const double *x,double *g){
	g[0] = 2*(x[0]-2);
	g[1] = 2*(x[1]-1);
  }
  void setRlst(const double *x, const double objValue){
	_objval = objValue;
	_x[0] = x[0];
	_x[1] = x[1];
  }

  const Vector2d &getOptX()const{
	return _x;
  }
  double getObjVal()const{
	return _objval;
  }

private:
  double _objval;
  Vector2d _x;
};
typedef boost::shared_ptr<TestIpoptObjfun> pTestIpoptObjfun;

BOOST_AUTO_TEST_SUITE(MtlOptEnergyTest)

BOOST_AUTO_TEST_CASE(testSimpleFunOpt){

  pTestIpoptObjfun energy = pTestIpoptObjfun(new TestIpoptObjfun());
  NoConIpoptSolver solver(energy);
  solver.setTol(1e-8);
  solver.setMaxIt(100);
  TEST_ASSERT( solver.initialize() );
  TEST_ASSERT( solver.solve() );
  
  const double optObjVal = 0.0f;
  Vector2d optX;
  optX << 2.0f,1.0f;
  ASSERT_EQ_TOL( energy->getObjVal(),optObjVal,1e-8);
  ASSERT_EQ_SMALL_VEC_TOL( energy->getOptX(),optX,2,1e-8);
}

BOOST_AUTO_TEST_CASE(testCtrlForward){
 
  const int r = 2;
  const int T = 4;
  const double h = 3.0;
  double alpha_k = 0.1f, alpha_m = 0.2f;
  Vector2d lambda,v0,z0;
  lambda << 0.25f, 2.5f;
  v0 << 1,2;
  z0 << 3,4;
  VectorXd w = VectorXd::Random(r*T-r);

  CtrlForceEnergy ctrlF;
  ctrlF.setTimestep(h);
  ctrlF.setTotalFrames(T);
  const VectorXd diagD = alpha_k*lambda+VectorXd::Ones(r)*alpha_m;
  ctrlF.setMtl(lambda,diagD);
  ctrlF.setV0(v0);
  ctrlF.setZ0(z0);
  ctrlF.precompute();

  MatrixXd V,Z;
  ctrlF.forward(w,V,Z);
  
  MASimulatorAD simulator;
  simulator.setTimeStep(h);
  simulator.setEigenValues<VectorXd>(lambda);
  simulator.setStiffnessDamping(SX(alpha_k));
  simulator.setMassDamping(SX(alpha_m));
  simulator.setIntialStatus(v0,z0);

  for (int i = 0; i < r; ++i){

	SX sG[2][2], sS[2];
	simulator.getImpIntegMat(sG,sS,i);
	Eigen::Matrix<double,2,2> G;
	Eigen::Matrix<double,2,1> S;
	G(0,0) = sG[0][0].getValue();
	G(0,1) = sG[0][1].getValue();
	G(1,0) = sG[1][0].getValue();
	G(1,1) = sG[1][1].getValue();
	S(0,0) = sS[0].getValue();
	S(1,0) = sS[1].getValue();

	ASSERT_EQ(G,(ctrlF.getG().block<2,2>(i*2,0)));
	ASSERT_EQ(S,(ctrlF.getS().block<2,1>(i*2,0)));
  }

  vector<VectorXd> w2,v2,z2;
  EIGEN3EXT::convert(Map<MatrixXd>(&w[0],r,T-1),w2);
  simulator.forward(w2,v2,z2);

  MatrixXd V2,Z2;
  EIGEN3EXT::convert(v2,V2);
  EIGEN3EXT::convert(z2,Z2);
  
  ASSERT_EQ(V.size(),V2.size());
  ASSERT_EQ(Z.size(),Z2.size());
  ASSERT_EQ(V,V2);
  ASSERT_EQ(Z,Z2);
    
}

BOOST_AUTO_TEST_CASE(testCtrlObjValue){
 
  CtrlForceEnergy ctrlF;
  
}

BOOST_AUTO_TEST_CASE(testCtrlGrad){
 
  CtrlForceEnergy ctrlF;
  
}

BOOST_AUTO_TEST_CASE(testCtrlOpt){
 
  CtrlForceEnergy ctrlF;
  
}

BOOST_AUTO_TEST_SUITE_END()
