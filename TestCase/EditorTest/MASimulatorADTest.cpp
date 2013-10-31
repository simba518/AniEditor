#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <CASADITools.h>
#include <MASimulatorAD.h>
using namespace Eigen;
using namespace LSW_ANI_EDITOR;

BOOST_AUTO_TEST_SUITE(MASimulatorADTest)

BOOST_AUTO_TEST_CASE(testAll){
  
  const int r = 2;
  const double h = 3.0;
  vector<SX> alpha_k(r), alpha_m(r), lambda(r);
  alpha_k[0] = 0.1;
  alpha_k[1] = 0.1;

  alpha_m[0] = 0.2;
  alpha_m[1] = 0.2;

  lambda[0] = 0.25;
  lambda[1] = 2.5;

  VectorXd v0(r);
  VectorXd z0(r);
  v0 << 1,2;
  z0 << 3,4;

  vector<SX> w(r);
  w[0] = 100;
  w[1] = 20;

  MASimulatorAD simulator;
  simulator.setTimeStep(h);
  simulator.setStiffnessDamping<vector<SX> >(alpha_k);
  simulator.setMassDamping(alpha_m);
  simulator.setEigenValues(lambda);
  simulator.setIntialStatus(v0,z0);
  ASSERT_EQ( simulator.reducedDim(),r );

  simulator.forward(w);
  VectorXd z(2),v(2);
  z[0] = simulator.getZ()[0].getValue();
  z[1] = simulator.getZ()[1].getValue();
  v[0] = simulator.getV()[0].getValue();
  v[1] = simulator.getV()[1].getValue();

  VectorXd correctZ(2),correctV(2);  
  correctZ << 231.343949044586, 7.863179074446681;
  correctV << 76.11464968152866, 1.287726358148893;

  ASSERT_EQ_SMALL_VEC_TOL( z, correctZ,2, 1e-10);
  ASSERT_EQ_SMALL_VEC_TOL( v, correctV,2, 1e-10);
}

BOOST_AUTO_TEST_CASE(testADMethod){
 
  vector<SXMatrix> sW,sV,sZ;
  vector<VectorXd> W,V,Z;
  const int T = 10;
  const int r = 3;
  W.resize(T-1);
  sW.resize(T-1);
  for (int i = 0; i < W.size(); ++i){
	W[i] = VectorXd::Random(r);
    sW[i] = CASADI::convert(W[i]);
  }

  const VectorXd z0 = VectorXd::Random(r);
  const VectorXd v0 = VectorXd::Random(r);
  
  const double h = 3.0;
  const double alpha_k = 0.1;
  const double alpha_m = 0.2;
  const VectorXd lambda = VectorXd::Random(r)+10.0f*VectorXd::Ones(r);

  MASimulatorAD simulator;
  simulator.setTimeStep(h);
  simulator.setEigenValues(lambda);
  simulator.setStiffnessDamping(alpha_k);
  simulator.setMassDamping(alpha_m);
  simulator.setIntialStatus(v0,z0);
  simulator.forward(W,V,Z);

  simulator.setIntialStatus(v0,z0);
  simulator.forward(sW,sV,sZ);

  ASSERT_EQ(Z.size(),sZ.size());
  ASSERT_EQ(V.size(),sV.size());
  ASSERT_EQ(Z.size(),T);
  ASSERT_EQ(V.size(),T);

  for (int j = 0; j < T-1; ++j){
	const VectorXd wi = CASADI::convert2Vec<double>(sW[j]);
	ASSERT_EQ(wi.size(),W[j].size());
	ASSERT_EQ_SMALL_VEC_TOL(wi,W[j],wi.size(),1e-12);	
  }

  for (int j = 0; j < T; ++j){

	const VectorXd zi = CASADI::convert2Vec<double>(sZ[j]);
	const VectorXd vi = CASADI::convert2Vec<double>(sV[j]);

	ASSERT_EQ(zi.size(),Z[j].size());
	ASSERT_EQ(vi.size(),V[j].size());
	ASSERT_EQ_SMALL_VEC_TOL(zi,Z[j],zi.size(),1e-12);
	ASSERT_EQ_SMALL_VEC_TOL(vi,V[j],vi.size(),1e-12);
  }
}

BOOST_AUTO_TEST_SUITE_END()
