#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <CASADITools.h>
#include <ComputeBjAD.h>
#include <ComputeBj.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(ComputeBjADTest)

BOOST_AUTO_TEST_CASE(computeTest){

  const VectorXd eyj = VectorXd::Random(9);
  const double Sqrt_Vj = 1.5f;

  const CasADi::SXMatrix yj=CASADI::convert(eyj);
  CasADi::SXMatrix bj;
  ComputeBjAD::compute(yj,Sqrt_Vj,bj);
  const VectorXd ebj = CASADI::convert2Vec<double>(bj);

  VectorXd b2;
  ComputeBj::compute(eyj,Sqrt_Vj,b2);

  ASSERT_EQ(b2.size(),9);
  ASSERT_EQ(ebj.size(),9);
  ASSERT_EQ_SMALL_VEC(ebj,b2,9);
}

BOOST_AUTO_TEST_SUITE_END()
