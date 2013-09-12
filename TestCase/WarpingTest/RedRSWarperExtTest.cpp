#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <RedRSWarperExt.h>
#include <RedRSWarperAD.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(RedRSWarperExtTest)

BOOST_AUTO_TEST_CASE(jacobianTest){

  const double Sqrt_V = 1.5f;
  const Eigen::Matrix<double,9,1> yj = Eigen::MatrixXd::Random(9,1);
  Eigen::Matrix<double,9,9> rlstAD;
  Eigen::Matrix<double,9,9> rlst;

  LSW_WARPING::RedRSWarperGradAD gradAD;
  LSW_WARPING::RedRSWarperGrad grad;
  for (int i = 0; i < 144*5*20; ++i){
	gradAD.evaluateGrad(Sqrt_V, yj, rlstAD);
  }

  for (int i = 0; i < 144*5*20; ++i){
	grad.grad(Sqrt_V, &yj(0,0), &rlst(0,0));
  }

  ASSERT_EQ_SMALL_MAT_TOL (rlst,rlstAD,1e-8);
}

BOOST_AUTO_TEST_CASE(jacobianZTest){

  const VectorXd y = VectorXd::Random(9*100);
  vector<double> Sqrt_V(y.size()/9);
  for (size_t i = 0; i < Sqrt_V.size(); ++i){
    Sqrt_V[i] = i+1.0f;
  }
  const int nodeId = 2;
  const MatrixXd hatW = MatrixXd::Random(y.size(),80);
  const MatrixXd BP = MatrixXd::Random(2000*3,y.size());
  Eigen::Matrix<double,3,-1> Jz, autoJz;
  
  LSW_WARPING::RedRSWarperGradAD gradAD;
  LSW_WARPING::RedRSWarperGrad grad;
  for (int i = 0; i < 5*100; ++i){
	gradAD.jacobian(y,Sqrt_V,nodeId,hatW,BP,autoJz);
  }

  for (int i = 0; i < 5*100; ++i){
	grad.jacobian(y,Sqrt_V,nodeId,hatW,BP,Jz);
  }
  ASSERT_EQ_SMALL_MAT_TOL (Jz,autoJz,1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
